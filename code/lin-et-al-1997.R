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
# 20220524


# Setup --------------------------------------------------------------

library(dplyr)
library(purrr)


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

surv_time <- runif(npt, 0, study_length)
## surv_time <- rexp(npt, 1 / exp_surv_mean)
cens_time <- sample(1:10, npt, TRUE, probab_cens_low)

dat <- tibble::tibble(
  i = rep(1:study_length, times = npt),
  j = i + 1,
  T = rep(surv_time, each = study_length),
  U = rep(cens_time, each = study_length),
  X = pmin(T, U),
  dt = as.integer(T <= U),
  stat = dplyr::case_when(
    dt == 1 & T < i ~ "z",
    dt == 0 & U < i ~ "z",
    dt == 1 & j < T & T < (j+1) ~ "w",
    dt == 1 & i < T & T < j ~ "x",
    dt == 0 & i < U & U < j ~ "c",
    TRUE ~ "o"
  ),
  cb_ = runif(npt * study_length, base_costs["min"], base_costs["max"]),
  cd_ = rep(runif(npt, diagn_costs["min"], diagn_costs["max"]), each = study_length),
  cx_ = rep(runif(npt, death_costs["min"], death_costs["max"]), each = study_length),
) %>%
  dplyr::mutate(
    pb =
      ifelse(stat %in% c("o", "w"), 1,
        ifelse(stat == "c", U - i,
          ifelse(stat == "x", T - i,
            0))),
    pd =
      ifelse(i == 1 & stat %in% c("o", "w"), 1,
        ifelse(i == 1 & stat == "c", U - i,
          ifelse(i == 1 & stat == "x", T - i,
            0))),
    px =
      ifelse(stat == "w", j - (T - 1),
        ifelse(stat == "x", T - i,
          0)),
    cb = cb_ * pb,
    cd = cd_ * pd,
    cx = cx_ * px,
    cc = cb + cd + cx
)

readr::write_csv(dat, "dat.csv")
