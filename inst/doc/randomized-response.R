## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  message = FALSE,
  warning = FALSE,
  fig.width = 6,
  fig.height = 3.5,
  fig.align = "center"
)

## ----confidential-data--------------------------------------------------------
library(dapper)
library(ggplot2)

cells <- data.frame(
  sex = c(1, 1, 0, 0),
  status = c(1, 0, 1, 0)
)
counts <- c(1198, 1493, 557, 1278)
cnf_df <- cells[rep(seq_len(nrow(cells)), times = counts), ]

set.seed(1)
n <- 400
ix <- sample(seq_len(nrow(cnf_df)), n, replace = FALSE)
cnf_df <- cnf_df[ix, ]

## ----randomized-release-------------------------------------------------------
ri <- as.logical(rbinom(2 * n, 1, 1/2))
ra <- rbinom(sum(ri), 1, 1/2)

sdp <- c(as.matrix(cnf_df))
sdp[ri] <- ra
prv_df <- data.frame(sex = sdp[seq_len(n)], status = sdp[n + seq_len(n)])

## ----admission-tables---------------------------------------------------------
admission_table <- function(x) {
  table(
    Sex = factor(x$sex, levels = c(0, 1), labels = c("Female", "Male")),
    Status = factor(x$status, levels = c(1, 0),
                    labels = c("Admitted", "Rejected"))
  )
}

knitr::kable(admission_table(cnf_df), caption = "Confidential counts")
knitr::kable(admission_table(prv_df), caption = "Privatized counts")

## ----latent-component---------------------------------------------------------
latent_f <- function(theta) {
  tl <- list(c(1, 1), c(1, 0), c(0, 1), c(0, 0))
  rs <- sample(tl, n, replace = TRUE, prob = theta)
  do.call(rbind, rs)
}

## ----posterior-component------------------------------------------------------
posterior_f <- function(dmat, theta) {
  sex <- dmat[, 1]
  status <- dmat[, 2]
  x <- c(
    sum(sex & status),
    sum(sex & !status),
    sum(!sex & status),
    sum(!sex & !status)
  )
  t1 <- rgamma(4, shape = x + 1, rate = 1)
  t1 / sum(t1)
}

## ----statistic-component------------------------------------------------------
statistic_f <- function(xi, sdp, i) {
  n <- length(sdp) %/% 2
  sum(xi == sdp[c(i, n + i)])
}

## ----mechanism-component------------------------------------------------------
mechanism_f <- function(sdp, sx) {
  sx * log(3/4) + (length(sdp) - sx) * log(1/4)
}

## ----sample-posterior---------------------------------------------------------
dmod <- new_privacy(
  posterior_f = posterior_f,
  latent_f = latent_f,
  mechanism_f = mechanism_f,
  statistic_f = statistic_f,
  npar = 4,
  varnames = c("pi_11", "pi_10", "pi_01", "pi_00")
)

dp_out <- dapper_sample(
  dmod,
  sdp = sdp,
  seed = 123,
  niter = 6000,
  warmup = 1000,
  chains = 4,
  init_par = rep(0.25, 4)
)

## ----parameter-summary--------------------------------------------------------
summary(dp_out)

## ----trace-plot, fig.height=4, fig.cap="Trace plots for the four cell probabilities. Colors distinguish the four chains.", fig.alt="Four panels show the sampled cell probabilities over iterations, with a different color for each chain."----
plot(dp_out) +
  scale_color_manual(values = c("#00468B", "#ED0000", "#42B540", "#925E9F")) +
  theme_bw(base_size = 11)

## ----odds-ratio---------------------------------------------------------------
odds_ratio_draws <- posterior::mutate_variables(
  dp_out$chain,
  odds_ratio = (pi_11 * pi_00) / (pi_10 * pi_01)
)
odds_ratio_draws <- posterior::subset_draws(
  odds_ratio_draws, variable = "odds_ratio"
)
posterior::summarise_draws(odds_ratio_draws)

## ----odds-ratio-density, fig.height=2.8, fig.cap="Privacy-aware posterior distribution of the odds ratio, displayed from 0 to 10.", fig.alt="Posterior density of the odds ratio, with a shaded central 80 percent interval and a dashed vertical reference line at one."----
bayesplot::mcmc_areas(
  odds_ratio_draws,
  prob = 0.8,
  prob_outer = 0.95,
  point_est = "none"
) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.5) +
  coord_cartesian(xlim = c(0, 10)) +
  labs(x = "Odds ratio", y = NULL) +
  theme_bw(base_size = 11) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

## ----comparison-draws---------------------------------------------------------
set.seed(1)
confidential_data <- as.matrix(cnf_df)
cps <- t(replicate(20000, posterior_f(confidential_data, NULL)))
odds_ratio_conf <- (cps[, 1] * cps[, 4]) / (cps[, 2] * cps[, 3])

set.seed(1)
noisy_data <- matrix(sdp, ncol = 2)
cps <- t(replicate(20000, posterior_f(noisy_data, NULL)))
odds_ratio_noisy <- (cps[, 1] * cps[, 4]) / (cps[, 2] * cps[, 3])

odds_ratio_private <- posterior::extract_variable(
  odds_ratio_draws, "odds_ratio"
)

## ----comparison-densities-----------------------------------------------------
posterior_density <- function(draws, data, analysis) {
  estimate <- density(draws, from = 0, to = 8)
  data.frame(odds_ratio = estimate$x, density = estimate$y,
             data = data, analysis = analysis)
}

comparison <- rbind(
  posterior_density(odds_ratio_conf, "Confidential", "Privacy-aware"),
  posterior_density(odds_ratio_private, "Privatized", "Privacy-aware"),
  posterior_density(odds_ratio_conf, "Confidential", "Naïve"),
  posterior_density(odds_ratio_noisy, "Privatized", "Naïve")
)
comparison$analysis <- factor(comparison$analysis,
                              levels = c("Privacy-aware", "Naïve"))

## ----posterior-comparison, fig.cap="Odds-ratio posteriors for confidential data (blue) and privatized data (red). The left panel accounts for randomized response; the right panel ignores it.", fig.alt="Two density panels compare confidential and privatized analyses. The privacy-aware posterior is wider, while the naïve posterior is concentrated closer to one."----
data_colors <- c(Confidential = "#00468B", Privatized = "#ED0000")

ggplot(comparison, aes(x = odds_ratio, y = density, fill = data, color = data)) +
  geom_area(position = "identity", alpha = 0.2, color = NA) +
  geom_line(linewidth = 0.6) +
  facet_wrap(~ analysis) +
  coord_cartesian(xlim = c(0, 8)) +
  scale_fill_manual(name = "Data", values = data_colors) +
  scale_color_manual(name = "Data", values = data_colors) +
  labs(x = "Odds ratio", y = "Density") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

