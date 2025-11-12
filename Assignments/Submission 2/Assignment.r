# ============================================================
# DAT246 — DAG-driven Bayesian models (rethinking / ulam)
# Outcomes: extent (continuous), quality (ordinal), result (binary)
# ============================================================

# --- packages
library(rethinking)   # ulam, standardize, precis, PSIS/WAIC helpers
library(tidyverse)

# --- load data
reviewers <- read.csv("Assignments/reviewers.csv", stringsAsFactors = FALSE)
reviews <- read.csv("Assignments/reviews.csv", stringsAsFactors = FALSE)

# --- try a sensible merge (adjust if needed)
common <- intersect(names(reviews), names(reviewers))
key <- dplyr::case_when(
  "reviewer_id" %in% common ~ "reviewer_id",
  "reviewer"    %in% common ~ "reviewer",
  TRUE ~ common[grepl("id|reviewer", tolower(common))][1]
)
if (!is.na(key) && !is.null(key)) {
  d <- left_join(reviews, reviewers, by = key)
} else {
  d <- reviews
}

# --- lowercase, snake_case to make life easier
names(d) <- names(d) |> tolower() |> gsub("\\s+", "_", x = _)

# ===============================
# 🧭 mapping your columns to the DAG
# ===============================
# If these guesses look wrong on your data, set them manually, e.g.:
# age_col      <- "age_years"
# skill_col    <- "skill"
# complexity_col <- "complexity"
# ai_col       <- "used.cr.technology"  # 0/1 (or factor we will coerce)
# extent_col   <- "extent"              # length/time tokens/words/etc.
# quality_col  <- "quality"             # 1..5 ordinal
# result_col   <- "result"              # approved vs discarded

guess <- function(patterns) {
  # find first column containing all substrings in 'patterns'
  ix <- which(vapply(names(d), function(nm){
    all(vapply(patterns, function(p) grepl(p, nm, ignore.case = TRUE), logical(1)))
  }, logical(1)))
  if (length(ix)) names(d)[ix[1]] else NA_character_
}

age_col        <- guess(c("age"))
skill_col      <- guess(c("skill"))
complexity_col <- guess(c("complex"))
ai_col         <- guess(c("used","cr")) %||% guess(c("ai","used")) %||% guess(c("tech","used"))
extent_col     <- guess(c("extent")) %||% guess(c("length")) %||% guess(c("time"))
quality_col    <- guess(c("quality")) %||% guess(c("clarity"))
result_col     <- guess(c("result")) %||% guess(c("approve")) %||% guess(c("decision","status"))

# If any of these are NA, set explicitly here:
# age_col        <- age_col        %||% "age"
# skill_col      <- skill_col      %||% "skill"
# complexity_col <- complexity_col %||% "complexity"
# ai_col         <- ai_col         %||% "used_cr_technology"
# extent_col     <- extent_col     %||% "extent"
# quality_col    <- quality_col    %||% "quality"
# result_col     <- result_col     %||% "result"

`%||%` <- function(a,b) if(!is.null(a) && !is.na(a) && nzchar(a)) a else b

message("Detected columns:")
print(list(age=age_col, skill=skill_col, complexity=complexity_col,
           ai=ai_col, extent=extent_col, quality=quality_col, result=result_col))

# --- basic coercions
if (!is.na(ai_col)) {
  d[[ai_col]] <- d[[ai_col]] |>
    tolower() |>
    dplyr::recode("yes"="1","true"="1","y"="1","ai"="1","used"="1",
                  "no"="0","false"="0","n"="0","not_used"="0",
                  .default=as.character(d[[ai_col]])) |>
    as.numeric()
}

if (!is.na(result_col)) {
  d[[result_col]] <- tolower(as.character(d[[result_col]]))
  d[[result_col]] <- dplyr::recode(d[[result_col]],
     "approved"="1","approve"="1","accepted"="1","accept"="1","merged"="1","kept"="1","true"="1","yes"="1",
     "rejected"="0","reject"="0","discarded"="0","discard"="0","removed"="0","false"="0","no"="0",
     .default=d[[result_col]])
  # if still not 0/1, fall back to numeric/coerce
  if (!all(unique(na.omit(d[[result_col]])) %in% c("0","1"))) {
    d[[result_col]] <- as.numeric(as.factor(d[[result_col]])) - 1 # map to {0,1}
  } else d[[result_col]] <- as.numeric(d[[result_col]])
}

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
