# sk_covid_19

> **Archived.** No longer maintained; last updated February 2020. The case data file it reads is not included.

A single script from February 2020 estimating the time-varying reproduction number of COVID-19 in South Korea from a case line list.

- `sk_ncov.r`: reads `South Korean case summary.xlsx`, imputes missing symptom onset dates from the mean reporting delay, separates local and imported cases, and estimates R with `EpiEstim` using the serial interval from Li et al. (2020, NEJM, https://doi.org/10.1056/NEJMoa2001316). It then summarises the latest estimate, with its credible interval and the probability that R is above 1.
