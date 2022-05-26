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
# 20220526


# Setup --------------------------------------------------------------

library(dplyr)
library(purrr)
library(survival)


# Function definitions -----------------------------------------------

total_costs_ea <- function() {
}


# Constant definitions -----------------------------------------------

npt <- 100
study_length <- 10
exp_surv_mean <- 6
base_costs <- c(min = 1000, max = 3000)
diagn_costs <- c(min = 5000, max = 15000)
death_costs <- c(min = 10000, max = 30000)
probab_cens_low <- c(rep(.05, times = study_length-1), .55)
probab_cens_mod <- c(rep(.08, times = study_length-1), .28)


# Simulations --------------------------------------------------------

set.seed(577522)

person_ds <- dplyr::tibble(
  id = 1:npt,
  surv_time = runif(npt, 0, study_length),
  ## surv_time <- rexp(npt, 1 / exp_surv_mean),
  cens_time = sample(1:10, npt, replace = TRUE, probab_cens_low),
  fu_time = pmin(surv_time, cens_time),
  fu_status = as.integer(surv_time == fu_time)
)

person_ds <- dplyr::bind_rows(
  dplyr::tibble(
    id = 1, surv_time = 5.5, cens_time = 4.75,
    fu_time = pmin(surv_time, cens_time),
    fu_status = as.integer(surv_time == fu_time)
    ),
  dplyr::filter(person_ds, id > 1)
)

interval_ds <- dplyr::tibble(
  id = rep(person_ds$id, each = study_length),
  i = rep(0:(study_length-1), times = npt),
  j = i + 1,
  surv_time = rep(person_ds$surv_time, each = study_length),
  cens_time = rep(person_ds$cens_time, each = study_length),
  fu_time = rep(person_ds$fu_time, each = study_length),
  fu_status = rep(person_ds$fu_status, each = study_length),
  status = dplyr::case_when(
    fu_status == 1 & surv_time < i ~ "y",
    fu_status == 0 & cens_time < i ~ "d",
    fu_status == 1 & i < surv_time & surv_time < j ~ "x",
    fu_status == 1 & j < surv_time & surv_time < j + 1 ~ "w",
    fu_status == 0 & i < cens_time & cens_time < j & surv_time - cens_time < 1 ~ "cw",
    fu_status == 0 & i < cens_time & cens_time <= j ~ "c",
    TRUE ~ "o"
  ),
  pb =
    ifelse(status == "y", 0,
      ifelse(status == "d", 0,
        ifelse(status %in% c("c", "cw", "x"), fu_time - i,
          1))),
  pd = ifelse(i == 0, pb, 0),
  px =
    ifelse(status == "x", fu_time - i,
      ifelse(status == "w", j - (fu_time - 1),
        ifelse(status == "cw", cens_time - (surv_time - 1),
          0))),
  cb_ = runif(npt * study_length, base_costs["min"], base_costs["max"]),
  cd_ = rep(runif(npt, diagn_costs["min"], diagn_costs["max"]), each = study_length),
  cx_ = rep(runif(npt, death_costs["min"], death_costs["max"]), each = study_length),
  cb = ifelse(fu_status == 0 & cens_time < i, NA, cb_ * pb),
  cd = ifelse(fu_status == 0 & cens_time < i, NA, cd_ * pd),
  cx = ifelse(fu_status == 0 & cens_time < i, NA, cx_ * px),
  cc = ifelse(fu_status == 0 & cens_time < i, NA, cb + cd + cx)
)

readr::write_csv(interval_ds, "../output/interval_ds.csv")

## readr::write_csv(dat, "../output/dat.csv")


# $\hat{S}_k

km <- survival::survfit(
  survival::Surv(fu_time, fu_status) ~ 1, data = person_ds) %>%
  summary()

Sk <- dplyr::tibble(ak = floor(km$time), survival = km$surv) %>%
  dplyr::group_by(ak) %>%
  dplyr::summarise(survival = dplyr::last(survival)) %>%
  pull(survival)
Sk <- c(1, Sk[-study_length])

kmdat <- data.frame(ak = floor(km$time), survival = km$surv)
Sk <- aggregate(kmdat, by = list(kmdat$ak), min)

Sk <- aggregate(survival ~ ak, data = kmdat, FUN = min)

km <- summary(survival::survfit(survival::Surv(fu_time, fu_status) ~ 1, data = person_ds))
km_estim <- data.frame(time = )
probab_surv_to_interval_start <- 

probab_surv_to_interval_start <- cbind(
  data.frame(interval = 1, surv = 1),
  data.frame(interval = km_estim$time))

km <- summary(survival::survfit(survival::Surv(fu_time, fu_status) ~ 1, data = person_ds))
km_estim <- dplyr::tibble(
  time = km$time, surv = km$surv, interval_end = ceiling(time), interval = interval_end + 1
)

probab_surv_to_interval_start <- km_estim %>%
  dplyr::group_by(interval) %>%
  dplyr::summarise(surv_interval_start = dplyr::last(surv))

probab_surv_to_interval_start <- dplyr::tibble(
  interval = 1, surv_interval_start = 1
) %>%
  dplyr::bind_rows(
    km_estim %>%
      dplyr::filter(interval <= study_length) %>%
      dplyr::group_by(interval) %>%
      dplyr::summarise(surv_interval_start = dplyr::last(surv))
  )

