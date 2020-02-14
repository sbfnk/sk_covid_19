library("tidyverse")
library("readxl")
library("EpiEstim")
library("bayestestR")
library("HDInterval")

## Read data, determine import stats and reporting delay
df <- read_excel("South Korean case summary.xlsx",
                 col_types = c("numeric", "text", "numeric", "text", "text",
                               "text", "date", "date", "date", "text")) %>%
  mutate(import_status = if_else(is.na(Imported), "local", "imported"), 
         report_delay = as.integer(as.Date(`Dx date`) - as.Date(`Sx onset`)))

## A few cases don't have a symptom date (either unknown or no symptoms). For
## these we estimate the onset of infectiousness by assuming the delay to
## reporting is the same as the mean delay. This could be done better to fully
## account for uncertainty, but at the moment it is only a few cases so this
## will not matter much. 
mean_delay <- round(mean(df$report_delay, na.rm = TRUE))
df_with_estimated_onsets <- df %>%
  mutate(`Sx onset` =
           if_else(is.na(report_delay), `Dx date` - mean_delay, `Sx onset`))

## Create time series of local and imported cases
ts <- df_with_estimated_onsets %>%
  mutate(dates = as.Date(`Sx onset`)) %>%
  group_by(dates) %>%
  summarise(local = sum(is.na(Imported)), imported = sum(!is.na(Imported))) %>%
  ungroup() %>%
  complete(dates = seq(min(dates), max(dates), by = "day"),
           fill = list(local = 0, imported = 0))

## Estimate time-varying reproduction number using serial interval distribution
## from Li et al., 2020, NEJM, https://doi.org/10.1056/NEJMoa2001316
x <- estimate_R(ts, method = "parametric_si",
                config = make_config(list(mean_si = 7.5, std_si = 3.4)))

## Plot - note that early on (left part of the plot) there is not enough data to
## inform R so the confidence intervals are very wide; as more data are used
## they become narrower. 
plot(x)

## To get the current estimate of R, extract the last data point
last_datapoint <- unlist(x$R[nrow(x$R), ])

## The function estimate_R assumes that the reproduction number is gamma
## distributed - to get estimates of most likely R, confidence interval and
## probability of R>1 we retrieve the parameters of the gamma distribution
## (shape and rate).
shape <- (last_datapoint[["Mean(R)"]] / last_datapoint[["Std(R)"]])^2
rate <- last_datapoint[["Mean(R)"]] / last_datapoint[["Std(R)"]]^2

## Generate 10000 samples from the gamma distribution to use for estimation
samples <- rgamma(n = 10000, shape = shape, rate = rate)

## Best guess, i.e. the maximum a posteriori estimate
best <- round(as.numeric(map_estimate(samples)), 1)
best
## [1] 0.5

## Probability of R > 1
prob_gt_1 <- round(sum(samples > 1) / length(samples), 2)
round(prob_gt_1, 2)
## [1] 0.08

## 95 percent high-density interval
ci_95 <- as.numeric(hdi(samples))
round(ci_95, 1)
## [1] 0.2 1.1
