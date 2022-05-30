# Header -------------------------------------------------------------
#
# Replicate simulation results from Lin et al. (1997):
# <https://dlin.web.unc.edu/wp-content/uploads/sites/1568/2013/04/LinEA97.pdf>.
#
# Reference:
#
# Lin DY, Feuer EJ, Etzioni R, Wax Y. Estimating medical costs from
# incomplete follow-up data. Biometrics. 1997 Jun;53(2):419-34. PMID:
# 9192444.
#
# John Sahrmann
# 20220530


# Setup --------------------------------------------------------------

library(dplyr)
library(purrr)
library(survival)


# Constant definitions -----------------------------------------------

n_pt <- 1000
study_length <- 10
exp_surv_mean <- 6
base_costs <- c(min = 1000, max = 3000)
diagn_costs <- c(min = 5000, max = 15000)
death_costs <- c(min = 10000, max = 30000)
probab_cens_low <- c(rep(.05, times = study_length-1), .55)
probab_cens_moderat <- c(rep(.08, times = study_length-1), .28)


# Function definitions -----------------------------------------------

total_costs_ea <- function() {
}

# Survival time generating functions.
unif_surv <- function() {
  runif(n_pt, 0, study_length)
}
unif_exp <- function() {
  rexp(n_pt, 1 / exp_surv_mean)
}

# Censoring time generating functions.
low_cens_start <- function() {
  sample(1:study_length, n_pt, replace = TRUE, probab_cens_low)
}
moderat_cens_start <- function() {
  sample(1:study_length, n_pt, replace = TRUE, probab_cens_moderat)
}


# Simulations --------------------------------------------------------

set.seed(577522)

simul_cost_hist <- function(f_surv, f_cens) {
  pt_ds <- dplyr::tibble(
    id = 1:n_pt,
    surv_time = f_surv(),
    cens_time = f_cens(),
    fu_time = pmin(surv_time, cens_time),
    fu_status = as.integer(surv_time == fu_time)
  )

  interval_ds <- dplyr::tibble(
    id = rep(pt_ds$id, each = study_length),
    i = rep(0:(study_length-1), times = n_pt),
    j = i + 1,
    surv_time = rep(pt_ds$surv_time, each = study_length),
    cens_time = rep(pt_ds$cens_time, each = study_length),
    fu_time = rep(pt_ds$fu_time, each = study_length),
    fu_status = rep(pt_ds$fu_status, each = study_length),
    status = dplyr::case_when(
      fu_status == 1 & surv_time < i                     ~ "y",
      fu_status == 0 & cens_time <= i                    ~ "d",
      fu_status == 1 & i < surv_time & surv_time < j     ~ "x",
      fu_status == 1 & j < surv_time & surv_time < j + 1 ~ "w",
      fu_status == 0 & i < cens_time & cens_time < j
      & surv_time - cens_time < 1                        ~ "cw",
      fu_status == 0 & i <= cens_time & cens_time <= j   ~ "c",
      TRUE                                               ~ "o"
    ),
    under_obs_at_interval_start = fu_time > i,
    under_obs_at_interval_end = fu_time > j,
    propor_base =
      ifelse(status == "y", 0,
        ifelse(status == "d", 0,
          ifelse(status %in% c("c", "cw", "x"), fu_time - i,
            1))),
    propor_diagn = ifelse(i == 0, propor_base, 0),
    propor_death =
      ifelse(status == "x", fu_time - i,
        ifelse(status == "w", j - (fu_time - 1),
          ifelse(status == "cw", cens_time - (surv_time - 1),
            0))),
    base_costs_raw = runif(
      n_pt * study_length, base_costs["min"], base_costs["max"]),
    diagn_costs_raw = rep(
      runif(n_pt, diagn_costs["min"], diagn_costs["max"]),
      each = study_length),
    death_costs_raw = rep(
      runif(n_pt, death_costs["min"], death_costs["max"]),
      each = study_length),
    base_costs = ifelse(
      fu_status == 0 & cens_time <= i,
      NA, base_costs_raw * propor_base),
    diagn_costs = ifelse(
      fu_status == 0 & cens_time <= i,
      NA, diagn_costs_raw * propor_diagn),
    death_costs = ifelse(
      fu_status == 0 & cens_time <= i,
      NA, death_costs_raw * propor_death),
    total_costs = ifelse(
      fu_status == 0 & cens_time <= i,
      NA, base_costs + diagn_costs + death_costs),
    propor_base_no_cens =
      ifelse(surv_time < i, 0,
        ifelse(i < surv_time & surv_time < j, surv_time - i,
          1)),
    propor_diagn_no_cens = ifelse(i == 0, propor_base_no_cens, 0),
    propor_death_no_cens =
      ifelse(i < surv_time & surv_time < j, surv_time - i,
        ifelse(j < surv_time & surv_time < j + 1, j - (surv_time - 1),
          0)),
    base_costs_no_cens = base_costs_raw * propor_base_no_cens,
    diagn_costs_no_cens = diagn_costs_raw * propor_diagn_no_cens,
    death_costs_no_cens = death_costs_raw * propor_death_no_cens,
    total_costs_no_cens = (
      base_costs_no_cens + diagn_costs_no_cens + death_costs_no_cens)
  )

  list(pt_ds = pt_ds, interval_ds = interval_ds)
}

res <- simul_cost_hist(unif_surv, low_cens_start)


readr::write_csv(res$interval_ds, "../output/interval_ds.csv")

pt_total_costs <- res$interval_ds %>%
  dplyr::group_by(id) %>%
  dplyr::summarise(fu_total_costs_no_cens = sum(total_costs_no_cens))
mean(pt_total_costs$fu_total_costs_no_cens)


## readr::write_csv(dat, "../output/dat.csv")

crude_total_costs <- interval_ds %>%
  dplyr::group_by(id) %>%
  summarise(total_costs = sum(cc, na.rm = TRUE))
mean(crude_total_costs$total_costs)


# $\hat{S}_k$

km <- summary(survival::survfit(survival::Surv(fu_time, fu_status) ~ 1, data = person_ds))
km_estim <- dplyr::tibble(
  time = km$time, surv = km$surv, interval_end = ceiling(time), interval = interval_end + 1
)

probab_surv_to_interval_start <- dplyr::tibble(
  interval = 1, surv_interval_start = 1
) %>%
  dplyr::bind_rows(
    km_estim %>%
      dplyr::filter(interval <= study_length) %>%
      dplyr::group_by(interval) %>%
      dplyr::summarise(surv_interval_start = dplyr::last(surv))
  )

# $\hat{E}_k$

estim_total_costs_interval <- interval_ds2 %>%
  dplyr::filter(!is.na(cc)) %>%
  dplyr::group_by(j) %>%
  dplyr::summarise(
    total_costs_interval = (
      sum(under_obs_at_interval_start * cc)
      / sum(under_obs_at_interval_start))
  )

# $\hat{E}_A$

estim_total_costs <- probab_surv_to_interval_start %>%
  dplyr::inner_join(estim_total_costs_interval, by = c(interval = "j")) %>%
  dplyr::mutate(
    surv_weighted_total_costs_interval = surv_interval_start * total_costs_interval
  ) %>%
  dplyr::pull(surv_weighted_total_costs_interval) %>%
  sum()
