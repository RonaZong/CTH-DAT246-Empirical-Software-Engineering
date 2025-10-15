set.seed(2025)
N <- 500

# Simulate basic covariates
Engineer_age <- rnorm(N, mean = 35, sd = 8)
Engineer_experience <- pmax(0, Engineer_age - rnorm(N, 22, 4))
Engineer_skill <- pmin(5, pmax(1, round(3 + 0.02 * Engineer_experience + rnorm(N, 0, 0.6), 2)))
Issue_complexity <- sample(1:5, N, replace = TRUE, prob = c(0.15,0.25,0.3,0.2,0.1))
AI_familiarity <- pmin(5, pmax(1, round(2 + (40 - Engineer_age)/20 + rnorm(N,0,0.8))))

# Probability of using AI increases with familiarity and decreases with age
logit_p_ai <- -1.2 + 0.6 * (AI_familiarity - 3) - 0.03 * (Engineer_age - 35)
AI_prob <- plogis(logit_p_ai)
AI_used <- rbinom(N, 1, AI_prob)

# Outcomes: CR_length and CR_quality (encode DAG A effects)
# CR_length (characters) - baseline around 1500, increase with complexity and decrease slightly with skill; AI may shorten length
CR_length_mean <- 1500 + 250 * (Issue_complexity - 3) - 80 * (Engineer_skill - 3) - 120 * AI_used
CR_length <- rpois(N, pmax(50, CR_length_mean))

# CR_quality (1-5) - baseline 3, increases with skill and experience, decreases with complexity; small positive or negative AI effect
CR_quality_mean <- 3 + 0.3 * (Engineer_skill - 3) + 0.01 * (Engineer_experience - 5) - 0.15 * (Issue_complexity - 3) + 0.1 * AI_used
CR_quality <- pmin(5, pmax(1, round(rnorm(N, CR_quality_mean, 0.6),2)))

# CR_relevance (perceived) - depends on length, quality and aversion to AI
logit_relevance <- -0.5 + 0.002 * (CR_length - 1500) + 0.8 * (CR_quality - 3) - 1.0 * AI_used
CR_relevance <- rbinom(N, 1, plogis(logit_relevance))

# Put into a data frame
df <- data.frame(
  CR_length, CR_quality, AI_used, Issue_complexity, Engineer_experience,
  Engineer_age, Engineer_skill, AI_familiarity, CR_relevance
)

print(summary(df))
write.csv(head(df,20), file = "Assignments/sim_head.csv", row.names = FALSE)
cat("Wrote Assignments/sim_head.csv (first 20 rows).\n")
