m_extent_prior <- ulam(
  alist(
    E ~ dexp(lambda),
    log(lambda) <- aE + bEU*U + bEA*A + bES*S + bEC*C,
    aE  ~ dnorm(0, 1),
    c(bEU, bEA, bES, bEC) ~ dnorm(0, 0.5)
  ),
  data = dat_extent,
  chains = 1, cores = 1,
  sample_prior = TRUE
)

post_prior_E <- extract.samples(m_extent_prior)
lambda_prior <- exp(
  post_prior_E$aE +
  post_prior_E$bEU * dat_extent$U +
  post_prior_E$bEA * dat_extent$A +
  post_prior_E$bES * dat_extent$S +
  post_prior_E$bEC * dat_extent$C
)

# Simulate prior predictive extent
extent_prior <- rexp(length(lambda_prior), rate = lambda_prior)

# Compare to observed
par(mfrow = c(1, 2))
hist(d$extent, main = "Observed extent", xlab = "Extent")
hist(extent_prior, main = "Prior predictive extent", xlab = "Extent")
par(mfrow = c(1, 1))

m_quality_prior <- ulam(
  alist(
    Q ~ dordlogit(phi, cutpoints),
    phi <- aQ + bQU*U + bQA*A + bQS*S + bQC*C,
    aQ  ~ dnorm(0, 1),
    c(bQU, bQA, bQS, bQC) ~ dnorm(0, 0.5),
    cutpoints ~ dnorm(0, 1.5)
  ),
  data = dat_quality,
  chains = 1, cores = 1,
  sample_prior = TRUE
)

# Use link() to simulate prior predictive probabilities
Q_prior <- sim(m_quality_prior, data = dat_quality, n = 1000)

# Compare marginal distributions (e.g., barplot of rowMeans)
obs_Q <- table(d$quality) / nrow(d)
prior_Q <- table(Q_prior) / length(Q_prior)

barplot(rbind(obs_Q, prior_Q),
        beside = TRUE,
        legend.text = c("Observed", "Prior predictive"),
        main = "Quality: observed vs prior predictive")