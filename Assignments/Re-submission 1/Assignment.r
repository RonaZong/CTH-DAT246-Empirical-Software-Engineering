# ============================================================
# DAT246 — DAG-driven Bayesian models (rethinking / ulam)
# Outcomes: extent (continuous), quality (ordinal), result (binary)
# ============================================================

# --- packages
library(rethinking)   # ulam, standardize, precis, PSIS/WAIC helpers
library(tidyverse)
library(ggplot2)

# --- load data
reviewers <- read.csv("Assignments/reviewers.csv", stringsAsFactors = FALSE)
reviews <- read.csv("Assignments/reviews.csv", stringsAsFactors = FALSE)
str(reviewers)
str(reviews)
head(reviewers)
head(reviews)

# --- try a sensible merge (adjust if needed)
df <- reviews %>% left_join(reviewers, by = "reviewer.id")

# Alternative: inner join to keep only matched rows
# dat <- reviews %>%
#   inner_join(reviewers, by = "reviewer.id") %>%
#   mutate(
#     used = as.integer(used.cr.technology),      # 0/1
#     extent_log = log1p(extent),                 # positive -> log1p
#     # z-scales for continuous predictors
#     complexity_z = as.numeric(scale(complexity)),
#     skill_z      = as.numeric(scale(skill)),
#     age_z        = as.numeric(scale(age)),
#     extent_z     = as.numeric(scale(extent_log)),
#     quality_int  = as.integer(quality),         # 1..5
#     quality_z    = as.numeric(scale(quality))   # for result model convenience
#   )

# Alternative: try common id-like columns
# common <- intersect(names(reviews), names(reviewers))
# key <- dplyr::case_when(
#   "reviewer_id" %in% common ~ "reviewer_id",
#   "reviewer"    %in% common ~ "reviewer",
#   TRUE ~ common[grepl("id|reviewer", tolower(common))][1]
# )
# if (!is.na(key) && !is.null(key)) {
#   df <- left_join(reviews, reviewers, by = key)
# } else {
#   df <- reviews
# }

# Rename to match DAG (adapt to your real column names!)
df <- df %>%
  rename(
    age        = age,                 # reviewer age
    skill      = skill,               # reviewer skill/expertise
    complexity = complexity,          # change complexity
    used       = used.cr.technology,  # AI usage (0/1)
    extent     = extent,              # review length / duration
    quality    = quality,             # quality score
    result     = result               # approval (0/1)
  )

# Coerce binaries and standardize predictors (z = (x-mean)/sd)
df <- df %>%
  mutate(
    used   = as.integer(used),
    result = as.integer(result),
    z_age        = as.numeric(scale(age)),
    z_skill      = as.numeric(scale(skill)),
    z_complexity = as.numeric(scale(complexity)),
    z_extent     = as.numeric(scale(extent)),
    z_quality    = as.numeric(scale(quality))
  )


summary(df[,c("age","skill","complexity","extent","quality")])
table(used = df$used, result = df$result)

ggplot(df, aes(x=factor(used))) + geom_bar() + labs(x="AI used", y="Count")
ggplot(df, aes(x=factor(result))) + geom_bar() + labs(x="Approval", y="Count")
ggplot(df, aes(z_complexity, z_extent)) + geom_point(alpha=.3) + geom_smooth(method="lm", se=FALSE)
ggplot(df, aes(z_skill, z_quality)) + geom_point(alpha=.3) + geom_smooth(method="lm", se=FALSE)


library(rethinking)

# EXTENT
dl_extent <- list(
  extent = as.numeric(scale(d$extent)),
  age = d$z_age, skill = d$z_skill,
  complexity = d$z_complexity, used = d$used
)

m_extent_prior <- ulam(
  alist(
    extent ~ dnorm(mu, sigma),
    mu <- a + b_age*age + b_skill*skill + b_complexity*complexity + b_used*used,
    a ~ dnorm(0, 1),
    c(b_age, b_skill, b_complexity, b_used) ~ dnorm(0, 0.5),
    sigma ~ dexp(1)
  ),
  data = dl_extent, chains=1, iter=500, sample=TRUE, dofit=FALSE
)
# Prior predictive: simulate with extract.prior() then link()

# --- keep a modeling frame
keep <- c(age_col, skill_col, complexity_col, ai_col, extent_col, quality_col, result_col)
keep <- keep[!is.na(keep)]
df <- d |> dplyr::select(all_of(keep)) |> dplyr::rename(
  age    = !!age_col,
  skill  = !!skill_col,
  complexity = !!complexity_col,
  AI     = !!ai_col,
  extent = !!extent_col,
  quality= !!quality_col,
  result = !!result_col
)


# --- coerce numerics + handle ordinal quality
numify <- function(x) suppressWarnings(as.numeric(x))
for (v in c("age","skill","complexity","extent")) if (v %in% names(df)) df[[v]] <- numify(df[[v]])

# Quality as ordered factor 1<2<...<K (K typically 5)
if ("quality" %in% names(df)) {
  # try to coerce to {1,...,K}
  qnum <- suppressWarnings(as.numeric(df$quality))
  if (all(!is.na(qnum))) {
    qK <- sort(na.omit(unique(qnum)))
    df$quality <- factor(qnum, levels = sort(qK), ordered = TRUE)
  } else {
    # if text, map common labels
    lab <- tolower(as.character(df$quality))
    map <- c("very low"=1,"low"=2,"ok"=3,"good"=4,"excellent"=5)
    q <- ifelse(lab %in% names(map), unname(map[lab]), NA)
    df$quality <- factor(q, levels = sort(na.omit(unique(q))), ordered = TRUE)
  }
}

# --- drop rows with missing essentials per model later
summary(df)

# --- prepare data list for modeling
data_list <- list(
  extent = d$extent,
  skill = standardize(d$skill),
  complexity = standardize(d$complexity),
  used_cr_tech = d$used.cr.technology,
  age = standardize(d$age)
)

# ============================================================
# 1) EXTENT model: Normal (identity link)
#    mu = a + b_skill*skill + b_complexity*complexity + b_ai*AI + b_age*age
# ============================================================


d_extent <- df |> dplyr::select(extent, skill, complexity, AI, age) |> na.omit()
dat1 <- list(
  extent = d_extent$extent,
  skill  = standardize(d_extent$skill),
  complexity = standardize(d_extent$complexity),
  AI     = as.numeric(d_extent$AI),
  age    = standardize(d_extent$age),
  N = nrow(d_extent)
)

m_extent <- ulam(
  alist(
    extent ~ dnorm(mu, sigma),
    mu <- a + bS*skill + bC*complexity + bAI*AI + bA*age,
    a ~ dnorm(0, 1),
    c(bS, bC, bAI, bA) ~ dnorm(0, 1),
    sigma ~ dexp(1)
  ),
  data = dat1, chains = 4, cores = 4, log_lik = TRUE
)

# ============================================================
# 2) QUALITY model: Ordered logistic
#    eta = a + b_skill*skill + b_complexity*complexity + b_ai*AI + b_age*age
#    quality ~ ordered_logistic(eta, kappa)  # kappa are ordered cutpoints
# ============================================================
d_quality <- df |> dplyr::select(quality, skill, complexity, AI, age) |> na.omit()
# Convert ordered factor to integer categories 1..K
q_y <- as.integer(d_quality$quality)
K <- length(unique(q_y))

dat2 <- list(
  quality = q_y,
  skill   = standardize(as.numeric(d_quality$skill)),
  complexity = standardize(as.numeric(d_quality$complexity)),
  AI      = as.numeric(d_quality$AI),
  age     = standardize(as.numeric(d_quality$age)),
  N = length(q_y),
  K = K
)

m_quality <- ulam(
  alist(
    quality ~ dordlogit(eta, kappa),
    eta <- a + bS*skill + bC*complexity + bAI*AI + bA*age,
    a ~ dnorm(0, 1),
    c(bS, bC, bAI, bA) ~ dnorm(0, 1),
    kappa ~ dnorm(0, 1.5)
  ),
  data = dat2,
  chains = 4, cores = 4, log_lik = TRUE,
  constraints = list(kappa = "ordered")  # enforce ordered cutpoints
)

# ============================================================
# 3) RESULT models (Binomial/Logit)
#    A) TOTAL effect of AI on approval: omit mediators extent/quality
#    B) DIRECT effect: condition on extent & quality
# ============================================================
# --- TOTAL
d_res_total <- df |> dplyr::select(result, AI, age, skill, complexity) |> na.omit()
dat3t <- list(
  result = as.integer(d_res_total$result),
  AI     = as.numeric(d_res_total$AI),
  age    = standardize(as.numeric(d_res_total$age)),
  skill  = standardize(as.numeric(d_res_total$skill)),
  complexity = standardize(as.numeric(d_res_total$complexity)),
  N = nrow(d_res_total)
)

m_result_total <- ulam(
  alist(
    result ~ dbinom(1, p),
    logit(p) <- a + bAI*AI + bA*age + bS*skill + bC*complexity,
    a ~ dnorm(0, 1),
    c(bAI, bA, bS, bC) ~ dnorm(0, 1)
  ),
  data = dat3t, chains = 4, cores = 4, log_lik = TRUE
)

# --- DIRECT (condition on mediators extent, quality)
# Join extent & quality; drop NAs row-wise
d_res_direct <- df |> dplyr::select(result, AI, extent, quality, age, skill, complexity) |> na.omit()
# For ordered quality, convert to numeric score for a simple direct model;
# (more purist: use a latent expectation from m_quality via link(), but numeric works fine in practice)
q_num <- if(is.ordered(d_res_direct$quality)) as.numeric(d_res_direct$quality) else as.numeric(d_res_direct$quality)

dat3d <- list(
  result = as.integer(d_res_direct$result),
  AI     = as.numeric(d_res_direct$AI),
  extent = standardize(as.numeric(d_res_direct$extent)),
  quality= standardize(q_num),
  age    = standardize(as.numeric(d_res_direct$age)),
  skill  = standardize(as.numeric(d_res_direct$skill)),
  complexity = standardize(as.numeric(d_res_direct$complexity)),
  N = nrow(d_res_direct)
)

m_result_direct <- ulam(
  alist(
    result ~ dbinom(1, p),
    logit(p) <- a + bAI*AI + bE*extent + bQ*quality + bA*age + bS*skill + bC*complexity,
    a ~ dnorm(0, 1),
    c(bAI, bE, bQ, bA, bS, bC) ~ dnorm(0, 1)
  ),
  data = dat3d, chains = 4, cores = 4, log_lik = TRUE
)

# ============================================================
# Summaries & comparisons
# ============================================================
cat("\n=== EXTENT model ===\n");  print(precis(m_extent, 2))
cat("\n=== QUALITY model (ordered logit) ===\n"); print(precis(m_quality, 2))
cat("\n=== RESULT TOTAL model ===\n"); print(precis(m_result_total, 2))
cat("\n=== RESULT DIRECT model ===\n"); print(precis(m_result_direct, 2))

# LOO/PSIS comparisons (predictive ability)
cat("\n=== Compare RESULT models (PSIS-LOO) ===\n")
print(compare(m_result_total, m_result_direct, func = PSIS))

# ============================================================
# Interpretation helpers
# ============================================================

# Posterior probability of direction (e.g., P(bAI < 0))
post_total  <- extract.samples(m_result_total)
post_direct <- extract.samples(m_result_direct)

cat("\nP_total( AI effect < 0 ) = ",
    mean(post_total$bAI < 0), "\n")
cat("P_direct( AI effect < 0 | extent,quality ) = ",
    mean(post_direct$bAI < 0), "\n")

# Marginal effect size (median odds ratio for AI in direct model)
OR_direct <- exp(post_direct$bAI)
cat("Median OR_direct(AI) = ", median(OR_direct),
    "  89% PI [", quantile(OR_direct,0.055), ", ", quantile(OR_direct,0.945), "]\n", sep="")

# Quality → Result (expected positive)
OR_quality <- exp(post_direct$bQ)
cat("Median OR(quality) = ", median(OR_quality),
    "  89% PI [", quantile(OR_quality,0.055), ", ", quantile(OR_quality,0.945), "]\n", sep="")

# Extent → Result (often negative or weak)
OR_extent <- exp(post_direct$bE)
cat("Median OR(extent) = ", median(OR_extent),
    "  89% PI [", quantile(OR_extent,0.055), ", ", quantile(OR_extent,0.945), "]\n", sep="")

# Quick checks for divergences / mixing
# (you can also use trankplot(), pairs() from rethinking)
max(rhat(m_extent)); min(neff(m_extent))
max(rhat(m_quality)); min(neff(m_quality))
max(rhat(m_result_total)); min(neff(m_result_total))
max(rhat(m_result_direct)); min(neff(m_result_direct))
