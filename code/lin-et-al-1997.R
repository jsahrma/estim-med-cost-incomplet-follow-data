## Header ------------------------------------------------------------
##
## Replicate simulation results from Lin et al. (1997):
## <https://dlin.web.unc.edu/wp-content/uploads/sites/1568/2013/04/LinEA97.pdf>.
##
## Reference:
##
## Lin DY, Feuer EJ, Etzioni R, Wax Y. Estimating medical costs from
## incomplete follow-up data. Biometrics. 1997 Jun;53(2):419-34. PMID:
## 9192444.
##
## John Sahrmann
## 20220530


## Setup -------------------------------------------------------------

library(dplyr)
library(purrr)
library(survival)


## Constant definitions ----------------------------------------------

n_simul <- 50000
n_pt <- 1000
study_length <- 10
exp_surv_mean <- 6
base_costs <- c(min = 1000, max = 3000)
diagn_costs <- c(min = 5000, max = 15000)
death_costs <- c(min = 10000, max = 30000)
probab_cens_low <- c(rep(.05, times = study_length-1), .55)
probab_cens_moderat <- c(rep(.08, times = study_length-1), .28)


# Function definitions -----------------------------------------------

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

total_costs_true <- function(data) {
  data %>%
    dplyr::group_by(id) %>%
    dplyr::summarise(pt_total_costs = sum(total_costs_no_cens)) %>%
    dplyr::pull(pt_total_costs) %>%
    mean()
}

total_costs_full_sampl <- function(data) {
  data %>%
    dplyr::group_by(id) %>%
    dplyr::summarise(
      pt_total_costs = sum(total_costs, na.rm = TRUE)
    ) %>%
    dplyr::pull(pt_total_costs) %>%
    mean()
}

total_costs_uncens_sampl <- function(data) {
  data %>%
    dplyr::filter(fu_status == 1) %>%
    dplyr::group_by(id) %>%
    dplyr::summarise(pt_total_costs = sum(total_costs)) %>%
    dplyr::pull(pt_total_costs) %>%
    mean()
}

total_costs_ea <- function(pt_ds, interval_ds) {
  # Get Kaplan-Meier estimates for survival.
  km <- summary(
    survival::survfit(
      survival::Surv(fu_time, fu_status) ~ 1, data = pt_ds))
  km_estim <- dplyr::tibble(
    time = km$time, surv = km$surv, interval_end = ceiling(time),
    interval = interval_end + 1
  )
  # Produce a data set containing the survival estimates at the start
  # of each interval. For the first interval, this is simply one. For
  # subsequent intervals, use the survival estimate at the time point
  # closest to, but before, the start of the interval.
  probab_surv_to_interval_start <- dplyr::tibble(
    interval = 1, surv_interval_start = 1
  ) %>%
    dplyr::bind_rows(
      km_estim %>%
        dplyr::filter(interval <= study_length) %>%
        dplyr::group_by(interval) %>%
        dplyr::summarise(surv_interval_start = dplyr::last(surv))
    )

  # Produce a data set containing interval-specific total costs
  # conditional on patients surviving to the start of each interval.
  total_costs_per_interval <- interval_ds %>%
    dplyr::filter(!is.na(total_costs)) %>%
    dplyr::group_by(j) %>%
    dplyr::summarise(
      total_costs_in_interval = (
        sum(under_obs_at_interval_start * total_costs)
        / sum(under_obs_at_interval_start))
  )

  # Compute the estimator.
  probab_surv_to_interval_start %>%
    dplyr::inner_join(
      total_costs_per_interval, by = c(interval = "j")
    ) %>%
    dplyr::mutate(
      surv_weighted_total_costs_in_interval = (
        surv_interval_start * total_costs_in_interval)
  ) %>%
  dplyr::pull(surv_weighted_total_costs_in_interval) %>%
  sum()
}

total_costs_eb <- function(pt_ds, interval_ds) {
  # Get Kaplan-Meier estimates for survival.
  km <- summary(
    survival::survfit(
      survival::Surv(fu_time, fu_status) ~ 1, data = pt_ds))
  km_estim <- dplyr::tibble(
    time = km$time, surv = km$surv, interval_end = ceiling(time),
    interval = interval_end + 1
  )
  # Produce a data set containing the survival estimates at the start
  # of each interval. For the first interval, this is simply one. For
  # subsequent intervals, use the survival estimate at the time point
  # closest to, but before, the start of the interval.
  probab_surv_to_interval_start <- dplyr::tibble(
    interval = 1, surv_interval_start = 1
  ) %>%
    dplyr::bind_rows(
      km_estim %>%
        dplyr::filter(interval <= study_length) %>%
        dplyr::group_by(interval) %>%
        dplyr::summarise(surv_interval_start = dplyr::last(surv))
    )

  # Produce a data set containing interval-specific total costs
  # conditional on patients surviving to just after the start of each
  # interval.
  n_pt_under_obs_at_interval_end <- tapply(
    interval_ds$under_obs_at_interval_end,
    interval_ds$j,
    sum
  )
  any_pt_in_last_interval <- (
    n_pt_under_obs_at_interval_end[[study_length]] > 0)
  if (any_pt_in_last_interval) {
    total_costs_per_interval <- interval_ds %>%
      dplyr::filter(!is.na(total_costs)) %>%
      dplyr::group_by(j) %>%
      dplyr::summarise(
        total_costs_in_interval = (
          sum(under_obs_at_interval_end * total_costs)
          / sum(under_obs_at_interval_end))
      )
  } else {
    total_costs_per_interval <- interval_ds %>%
      dplyr::filter(j < study_length) %>%
      dplyr::filter(!is.na(total_costs)) %>%
      dplyr::group_by(j) %>%
      dplyr::summarise(
        total_costs_in_interval = (
          sum(under_obs_at_interval_end * total_costs)
          / sum(under_obs_at_interval_end))
      )
  }

  # Compute the estimator.
  probab_surv_to_interval_start %>%
    dplyr::inner_join(
      total_costs_per_interval, by = c(interval = "j")
    ) %>%
    dplyr::mutate(
      surv_weighted_total_costs_in_interval = (
        surv_interval_start * total_costs_in_interval)
  ) %>%
  dplyr::pull(surv_weighted_total_costs_in_interval) %>%
  sum()
}


# Simulations --------------------------------------------------------

simul_pt <- function(f_surv, f_cens) {
  dplyr::tibble(
    id = 1:n_pt,
    surv_time = f_surv(),
    cens_time = f_cens(),
    fu_time = pmin(surv_time, cens_time),
    fu_status = as.integer(surv_time == fu_time)
  )
}


simul_cost_hist <- function(data) {
  dplyr::tibble(
    id = rep(data$id, each = study_length),
    i = rep(0:(study_length-1), times = n_pt),
    j = i + 1,
    surv_time = rep(data$surv_time, each = study_length),
    cens_time = rep(data$cens_time, each = study_length),
    fu_time = rep(data$fu_time, each = study_length),
    fu_status = rep(data$fu_status, each = study_length),
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
    # under_obs_at_interval_start = !(status %in% c("d", "y")),
    # under_obs_at_interval_end = !(status %in% c("c", "d", "x", "y")),
    under_obs_at_interval_start = !(status %in% c("d")),
    under_obs_at_interval_end = !(status %in% c("c", "d")),
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
}


set.seed(577522)

pt_ds <- simul_pt(unif_surv, low_cens_start)
cost_hist_ds <- simul_cost_hist(pt_ds)

E <- total_costs_true(cost_hist_ds)
E_F <- total_costs_full_sampl(cost_hist_ds)
E_U <- total_costs_uncens_sampl(cost_hist_ds)
E_A <- total_costs_ea(pt_ds, cost_hist_ds)
E_B <- total_costs_eb(pt_ds, cost_hist_ds)

E; E_F; E_U; E_A; E_B

res <- simul_cost_hist(unif_surv, low_cens_start)


readr::write_csv(res$interval_ds, "../output/interval_ds.csv")

total_costs_true(res$interval_ds)
total_costs_ea(res$pt_ds, res$interval_ds)
total_costs_eb(res$pt_ds, res$interval_ds)
