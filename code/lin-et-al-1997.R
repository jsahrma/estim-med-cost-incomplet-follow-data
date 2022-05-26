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


# Constant definitions -----------------------------------------------

npt <- 1e4
study_length <- 10
exp_surv_mean <- 6
base_costs <- c(min = 1000, max = 3000)
diagn_costs <- c(min = 5000, max = 15000)
death_costs <- c(min = 10000, max = 30000)
probab_cens_low <- c(rep(.05, times = study_length-1), .55)
probab_cens_mod <- c(rep(.08, times = study_length-1), .28)


# Simulations --------------------------------------------------------

set.seed(577522)

## surv_time <- runif(npt, 0, study_length)
surv_time <- rexp(npt, 1 / exp_surv_mean)
## cens_time <- sample(1:10, npt, TRUE, probab_cens_low)
cens_time <- rep(10, times = npt)

dat <- tibble::tibble(
  id = rep(1:npt, each = study_length),
  i = rep(0:(study_length-1), times = npt),
  j = i + 1,
  T = rep(surv_time, each = study_length),
  U = rep(cens_time, each = study_length),
  X = pmin(T, U),
  dt = as.integer(T <= U),
  stat = dplyr::case_when(
    X < i ~ "z",
    j < T & T < (j+1) ~ "w",
    i < T & T < j ~ "x",
    i < U & U < j ~ "c",
    TRUE ~ "o"
  ),
  cb_ = runif(npt * study_length, base_costs["min"], base_costs["max"]),
  cd_ = rep(runif(npt, diagn_costs["min"], diagn_costs["max"]), each = study_length),
  cx_ = rep(runif(npt, death_costs["min"], death_costs["max"]), each = study_length),
) %>%
  dplyr::mutate(
    pb =
      ifelse(stat %in% c("o", "w"), 1,
        ifelse(stat %in% c("c", "x"), X - i,
          0)),
    pd =
      ifelse(i == 0 & stat %in% c("o", "w"), 1,
        ifelse(i == 0 & stat %in% c("c", "x"), X - i,
          0)),
    px =
      ifelse(stat == "w", j - (T - 1),
        ifelse(stat == "x", T - i,
          0)),
    cb = cb_ * pb,
    cd = cd_ * pd,
    cx = cx_ * px,
    cc = cb + cd + cx
)

readr::write_csv(dat, "../output/dat.csv")

i_costs <- dat %>%
  dplyr::group_by(j) %>%
  dplyr::summarise(i_costs = mean(cc))
i_costs

# uniform distribution
##        j i_costs
##    <dbl>   <dbl>
##  1     1  13392.  13400
##  2     2   3699.   3700
##  3     3   3503.   3500
##  4     4   3306.   3300
##  5     5   3089.   3100
##  6     6   2902.   2900
##  7     7   2698.   2700
##  8     8   2502.   2500
##  9     9   2293.   2300
## 10    10   1103.   1100
# exponential distribution
##        j i_costs
##    <dbl>   <dbl>
##  1     1  13877.
##  2     2   3948.
##  3     3   3349.
##  4     4   2835.
##  5     5   2402.
##  6     6   2025.
##  7     7   1723.
##  8     8   1453.
##  9     9   1237.
## 10    10   1045.

total_costs <- dat %>%
  dplyr::group_by(id) %>%
  dplyr::summarise(total_cc = sum(cc)) %>%
  dplyr::pull(total_cc)
mean(total_costs)


(0.1 * 20000) + 2000 + 10000 +
.9*0.1*20000 + .9*2000 +
.8*0.1*20000 + .8*2000 +
.7*0.1*20000 + .7*2000 +
.6*0.1*20000 + .6*2000 +
.5*0.1*20000 + .5*2000 +
.4*0.1*20000 + .4*2000 +
.3*0.1*20000 + .3*2000 +
.2*0.1*20000 + .2*2000 +
.1*0.1*20000 + .1*2000

sum(12000, seq(3600, 2000, by = -200))

(0.1 * 20000) + 2000 + 10000 +
0.1*20000 + .9*2000 +
0.1*20000 + .8*2000 +
0.1*20000 + .7*2000 +
0.1*20000 + .6*2000 +
0.1*20000 + .5*2000 +
0.1*20000 + .4*2000 +
0.1*20000 + .3*2000 +
0.1*20000 + .2*2000 +
0.1*20000 + .1*2000

.9 * 12000 + .1 * 30000 +
  .8 * 2000 + .1 * 20000 +
  .7 * 2000 + .1 * 20000 +
  .6 * 2000 + .1 * 20000 +
  .5 * 2000 + .1 * 20000 +
  .4 * 2000 + .1 * 20000 +
  .3 * 2000 + .1 * 20000 +
  .2 * 2000 + .1 * 20000 +
  .1 * 2000 + .1 * 20000 +
  .0 * 2000 + .1 * 20000
# 39000

.9 * 12000 + .1 * 21000 +
  .8 * 2000 + .1 * 21000 +
  .7 * 2000 + .1 * 21000 +
  .6 * 2000 + .1 * 21000 +
  .5 * 2000 + .1 * 21000 +
  .4 * 2000 + .1 * 21000 +
  .3 * 2000 + .1 * 21000 +
  .2 * 2000 + .1 * 21000 +
  .1 * 2000 + .1 * 21000 +
  .0 * 2000 + .1 * 21000
# 39000

.9 * 12000 + .1 * 32000 +
  .8 * 2000 + .1 * 20000 +
  .7 * 2000 + .1 * 20000 +
  .6 * 2000 + .1 * 20000 +
  .5 * 2000 + .1 * 20000 +
  .4 * 2000 + .1 * 20000 +
  .3 * 2000 + .1 * 20000 +
  .2 * 2000 + .1 * 20000 +
  .1 * 2000 + .1 * 20000 +
  .0 * 2000 + .1 * 20000
# 39200

.9 * 12000 + .1 * 32000 +
  .8 * 2000 + .1 * 22000 +
  .7 * 2000 + .1 * 22000 +
  .6 * 2000 + .1 * 22000 +
  .5 * 2000 + .1 * 22000 +
  .4 * 2000 + .1 * 22000 +
  .3 * 2000 + .1 * 22000 +
  .2 * 2000 + .1 * 22000 +
  .1 * 2000 + .1 * 22000 +
  .0 * 2000 + .1 * 22000
# 41000

.9 * 12000 + .1 * 16000 +
  .8 * 2000 + .1 * 21000 +
  .7 * 2000 + .1 * 21000 +
  .6 * 2000 + .1 * 21000 +
  .5 * 2000 + .1 * 21000 +
  .4 * 2000 + .1 * 21000 +
  .3 * 2000 + .1 * 21000 +
  .2 * 2000 + .1 * 21000 +
  .1 * 2000 + .1 * 21000 +
  .0 * 2000 + .1 * 21000
# 38500

seq(.9, 0, -.1) * c(12000, rep(2000, 9)) + rep(.1, 10) * c(16000, rep(21000, 9))

.9 * 12000 + .1 * 26000 +
  .8 * 2000 + .1 * 21000 +
  .7 * 2000 + .1 * 21000 +
  .6 * 2000 + .1 * 21000 +
  .5 * 2000 + .1 * 21000 +
  .4 * 2000 + .1 * 21000 +
  .3 * 2000 + .1 * 21000 +
  .2 * 2000 + .1 * 21000 +
  .1 * 2000 + .1 * 21000 +
  .0 * 2000 + .1 * 11000
# 38500

seq(.9, 0, -.1) * c(12000, rep(2000, 9)) + rep(.1, 10) * c(26000, rep(21000, 8), 11000)

5*2000 + 20000 + .9*10000
# 39000
5*2000 + .95*20000 + .95*10000
# 38500

sum((1 - punif(1:10, min = 0, max = 10)) * c(12000, rep(2000, 9))) +
sum(dunif(1:10, min = 0, max = 10) * c(30000, rep(20000, 9)))

sum((1 - pexp(1:10, 1 / exp_surv_mean)) * c(12000, rep(2000, 9))) +
sum(dexp(1:10, 1 / exp_surv_mean) * c(30000, rep(20000, 9)))


6*2000 + pexp(10, 1/exp_surv_mean)*20000 + (1-pexp(1, 1/exp_surv_mean))*10000

sum((1 - pexp(1:10, 1 / exp_surv_mean)) * c(12000, rep(2000, 9))) +
sum(dexp(1:10, 1 / exp_surv_mean) * rep(21000, 10))


# proportion of sample surviving to start of each interval
1 - pexp(1:10, rate = 1 / exp_surv_mean)

# proportion of sample dying in each interval
pexp(1:10, rate = 1 / exp_surv_mean) - pexp(0:9, rate = 1 / exp_surv_mean)



sum((1 - pexp(1:10, rate = 1 / exp_surv_mean)) * c(12000, rep(2000, 9)) +
(pexp(1:10, rate = 1 / exp_surv_mean) - pexp(0:9, rate = 1 / exp_surv_mean)) * 20000)

sum((1 - pexp(1:10, rate = 1 / exp_surv_mean)) * c(12000, rep(2000, 9)) +
(pexp(1:10, rate = 1 / exp_surv_mean) - pexp(0:9, rate = 1 / exp_surv_mean)) * 21000)


(1 - pexp(1:10, rate = 1 / exp_surv_mean)) * c(12000, rep(2000, 9)) +
(pexp(1:10, rate = 1 / exp_surv_mean) - pexp(0:9, rate = 1 / exp_surv_mean)) * 20000 +
(pexp(1:10, rate = 1 / exp_surv_mean) - pexp(0:9, rate = 1 / exp_surv_mean)) * 2000

sum((1 - pexp(1:10, rate = 1 / exp_surv_mean)) * c(12000, rep(2000, 9)) + (pexp(1:10, rate = 1 / exp_surv_mean) - pexp(0:9, rate = 1 / exp_surv_mean)) * 20000)


sum(
  (1 - pexp(1:10, rate = 1 / exp_surv_mean)) * c(12000, rep(2000, 9)) +
(pexp(1:10, rate = 1 / exp_surv_mean) - pexp(0:9, rate = 1 / exp_surv_mean)) * 20000 +
(pexp(1:10, rate = 1 / exp_surv_mean) - pexp(0:9, rate = 1 / exp_surv_mean)) * 1000
)

sum(
  (1 - pexp(1:10, rate = 1 / exp_surv_mean)) * c(12000, rep(2000, 9)) +
(pexp(1:10, rate = 1 / exp_surv_mean) - pexp(0:9, rate = 1 / exp_surv_mean)) * 20000 +
(1 - (pexp(1:10, rate = 1 / exp_surv_mean) - pexp(0:9, rate = 1 / exp_surv_mean))) * 2000
)
