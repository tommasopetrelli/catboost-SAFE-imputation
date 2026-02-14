################################################################################
# Method B — Model I (Version 5)
# CatBoost MultiClass on Q8A bins + probabilistic imputation
#
#
# Goal of this script:
#   1) Train CatBoost MultiClass on observed Q8A (bins)
#   2) Run a time-aware wave holdout (diagnostic only)
#   3) Report comparison metrics (Table 10)
#   4) Fit the final model on all observed Q8A
#   5) Impute missing Q8A by stochastic draw from predicted probabilities
#   6) Merge imputed midpoints back to the “general” SAFE dataset

################################################################################

library(tidyverse)
library(haven)
library(catboost)
library(caret)
library(parallel)

set.seed(123)
t0 <- Sys.time()

# ==============================================================================
# BLOCK 0 — Paths + runtime knobs
# ==============================================================================
path_in       <- "C:/Università/Tesi 3/ecb.SAFE_microdata/safepanel_allrounds.dta"
path_out_data <- "C:/Università/Tesi 3/Dataset_Method_B_V5.dta"

path_out_dist <- "C:/Università/Tesi 3/Comparison Metrics/Method_B_Diag_dist_V5.dta"
path_out_tvd  <- "C:/Università/Tesi 3/Comparison Metrics/Method_B_Diag_tvd_V5.dta"

path_out_metrics_csv <- "C:/Università/Tesi 3/Comparison Metrics/Method_B_TestMetrics_V5.csv"
path_out_metrics_tex <- "C:/Università/Tesi 3/Comparison Metrics/Method_B_TestMetrics_V5.tex"

n_cores <- max(1L, detectCores())

# Fixed bin order used in the thesis (Table 10)
BIN_LEVELS <- c("1", "2", "5", "6", "4")

# Midpoints used for €-space evaluation and final merged output
MIDPOINTS_MAP <- c(
  "1" = 12500,
  "2" = 62500,
  "5" = 175000,
  "6" = 625000,
  "4" = 1500000
)

# ==============================================================================
# BLOCK 1 — Helpers (metrics + credibility tables)
# ==============================================================================

# Option C missingness-informativeness test (categorical target)
missingness_informative <- function(data, var, target, p_thresh = 0.20) {
  keep_rows <- !is.na(data[[target]])
  d <- data[keep_rows, , drop = FALSE]

  miss_flag <- is.na(d[[var]])
  if (all(!miss_flag)) return(TRUE)

  tab <- table(miss_flag, d[[target]])
  if (min(dim(tab)) < 2) return(TRUE)

  pval <- suppressWarnings(chisq.test(tab)$p.value)
  !is.na(pval) && pval < p_thresh
}

safe_logloss <- function(probs, y_true_cls, n_classes) {
  # probs: N x K (or vector length N*K), y_true_cls: 0..K-1
  eps <- 1e-15
  if (is.null(dim(probs))) probs <- matrix(probs, ncol = n_classes, byrow = TRUE)
  idx <- cbind(seq_len(nrow(probs)), y_true_cls + 1L)
  p_true <- pmax(probs[idx], eps)
  -mean(log(p_true))
}

quadratic_weighted_kappa <- function(truth_int, pred_int, k) {
  O <- table(factor(truth_int, levels = 1:k), factor(pred_int, levels = 1:k))
  O <- as.matrix(O)
  N <- sum(O)
  if (N == 0) return(NA_real_)

  row_m <- rowSums(O)
  col_m <- colSums(O)
  E <- outer(row_m, col_m) / N

  W <- outer(1:k, 1:k, function(i, j) ((i - j)^2) / ((k - 1)^2))
  1 - (sum(W * O) / sum(W * E))
}

rmse  <- function(p, y) sqrt(mean((p - y)^2, na.rm = TRUE))
mae   <- function(p, y) mean(abs(p - y), na.rm = TRUE)
medae <- function(p, y) median(abs(p - y), na.rm = TRUE)

bias_mean   <- function(p, y) mean(p - y, na.rm = TRUE)
bias_median <- function(p, y) median(p - y, na.rm = TRUE)

r2 <- function(p, y) {
  ok <- is.finite(p) & is.finite(y)
  y <- y[ok]; p <- p[ok]
  if (length(y) < 2) return(NA_real_)
  ss_res <- sum((y - p)^2)
  ss_tot <- sum((y - mean(y))^2)
  if (ss_tot == 0) return(NA_real_)
  1 - ss_res / ss_tot
}

spearman_fun <- function(p, y) {
  suppressWarnings(cor(p, y, method = "spearman", use = "complete.obs"))
}

snap_to_midpoint <- function(x, mps) {
  # Deterministic nearest-neighbour snap (ties resolved by first midpoint)
  mps[max.col(-abs(outer(x, mps, "-")), ties.method = "first")]
}

# Credibility tables: distributions and TVD (overall / by group)
make_dist <- function(df, y_levels, group_vars = NULL) {
  if (is.null(group_vars) || length(group_vars) == 0) {
    df <- df %>% mutate(grp = "ALL")
    group_vars <- "grp"
  }

  ct_true <- df %>%
    count(across(all_of(group_vars)), q8a_true_bin, name = "n_true") %>%
    rename(bin = q8a_true_bin)

  ct_pred <- df %>%
    count(across(all_of(group_vars)), q8a_pred_bin, name = "n_pred") %>%
    rename(bin = q8a_pred_bin)

  grid2 <- tidyr::crossing(
    df %>% distinct(across(all_of(group_vars))),
    bin = factor(y_levels, levels = y_levels, ordered = TRUE)
  )

  dist <- grid2 %>%
    left_join(ct_true, by = c(group_vars, "bin")) %>%
    left_join(ct_pred, by = c(group_vars, "bin")) %>%
    mutate(
      n_true = tidyr::replace_na(n_true, 0L),
      n_pred = tidyr::replace_na(n_pred, 0L)
    ) %>%
    group_by(across(all_of(group_vars))) %>%
    mutate(
      share_true = n_true / sum(n_true),
      share_pred = n_pred / sum(n_pred),
      diff_share = share_pred - share_true
    ) %>%
    ungroup()

  tvd_tbl <- dist %>%
    group_by(across(all_of(group_vars))) %>%
    summarise(
      n = sum(n_true),
      tvd = 0.5 * sum(abs(share_pred - share_true)),
      .groups = "drop"
    ) %>%
    left_join(
      df %>%
        group_by(across(all_of(group_vars))) %>%
        summarise(
          mean_true_int = mean(q8a_true_int, na.rm = TRUE),
          mean_pred_int = mean(q8a_pred_int, na.rm = TRUE),
          mean_shift    = mean_pred_int - mean_true_int,
          .groups = "drop"
        ),
      by = group_vars
    )

  list(dist = dist, tvd = tvd_tbl)
}

# Tail-risk diagnostics (Table 10)
tail_risk_metrics <- function(true_int, pred_int, true_mid, pred_mid_snap) {
  Kmax <- max(true_int, na.rm = TRUE)
  Kmin <- min(true_int, na.rm = TRUE)
  ae_eur <- abs(pred_mid_snap - true_mid)

  list(
    big_mistake_rate  = mean(abs(pred_int - true_int) >= 2, na.rm = TRUE),
    ae_p95_eur        = as.numeric(quantile(ae_eur, 0.95, na.rm = TRUE)),
    topbin_under_rate = {
      idx <- (true_int == Kmax)
      if (!any(idx, na.rm = TRUE)) NA_real_ else mean(pred_int[idx] < Kmax, na.rm = TRUE)
    },
    catastrophic_rate = mean(
      (true_int == Kmin & pred_int == Kmax) | (true_int == Kmax & pred_int == Kmin),
      na.rm = TRUE
    )
  )
}

# ==============================================================================
# BLOCK 2 — Load + clean data (q8a only)
#   * non-EU countries are dropped
#   * non-relevant waves are dropped
#   * ony observations representing the effective demand are kept
#   * vulnerable firms are dropped
#   * useless variables are dropped
#   * q8a == 9 (DK/NA) treated as missing
# ==============================================================================
dataset_A <- read_dta(path_in) %>%
  filter(
    d0 %in% c("BE","DK","DE","EE","IE","GR","ES","FR","HR","IT",
              "CY","LV","LT","LU","MT","NL","AT","PL","PT","RO",
              "SI","SK","FI","SE","BG","CZ","HU"),
    wave > 10,
    !wave %in% c(31, 33, 35),
    q7a_a %in% c(1, 2) | q32 %in% c(1, 2, 3, 4, 5, 6),
    (vulnerable != 1)
  ) %>%
  select(
    -c(wgtcommon, wgtentr, wgtoldentr, wgtoldcommon,
       ida, experiment_1, experiment_2, intdate),
    -any_of("q8a_rec")
  ) %>%
  arrange(permid, wave) %>%
  mutate(
    across(
      -c(permid, wave),
      ~ if (inherits(.x, "haven_labelled")) haven::as_factor(.x) else as.factor(.x)
    ),
    q8a = as.character(q8a),
    q8a = stringr::str_trim(q8a),
    q8a = dplyr::na_if(q8a, ""),
    q8a = stringr::str_extract(q8a, "\\d+"),
    q8a = ifelse(q8a == "9", NA_character_, q8a),
    q8a = factor(q8a)
  )

cat("Step 1 completed | Rows:", nrow(dataset_A),
    "| Columns:", ncol(dataset_A),
    "| Missing q8a:", sum(is.na(dataset_A$q8a)), "\n")

# ==============================================================================
# BLOCK 3 — Feature selection + NA encoding (hybrid rule)
# ==============================================================================
id_vars    <- c("permid", "wave", "d0")
target_var <- "q8a"
leak_vars  <- c("q8a")

# (1) Drop almost-empty variables
hard_drop <- names(which(colMeans(is.na(dataset_A)) > 0.90))
dataset_A <- dataset_A %>% select(-all_of(hard_drop))

# (2) Keep variables whose missingness is informative of q8a (Option C)
candidate_vars <- setdiff(names(dataset_A), c(id_vars, leak_vars))
vars_to_keep <- candidate_vars[
  sapply(candidate_vars, function(v) missingness_informative(dataset_A, v, target_var))
]

# Force keep q7a_a if present
vars_to_keep <- union(vars_to_keep, intersect("q7a_a", names(dataset_A)))

# Keep id + target + selected predictors
dataset_A <- dataset_A %>%
  select(all_of(c(id_vars, target_var, vars_to_keep)))

# (3) Explicit NA encoding for predictors only
pred_cols <- setdiff(names(dataset_A), c(id_vars, leak_vars))
dataset_A <- dataset_A %>%
  mutate(across(all_of(pred_cols), ~ if (is.factor(.x)) forcats::fct_explicit_na(.x, "MISSING") else .x))

cat("Step 2 completed | Columns after missingness control:", ncol(dataset_A), "\n")

# ==============================================================================
# BLOCK 4 — Model matrix prep (observed q8a only)
#   * Ordinal level order used in the thesis
#   * CatBoost MultiClass expects labels 0..K-1
# ==============================================================================
model_df <- dataset_A %>% filter(!is.na(.data[[target_var]]))

y_factor <- as.factor(model_df[[target_var]])
y_levels <- intersect(BIN_LEVELS, levels(y_factor))
n_classes <- length(y_levels)

y_ord_factor <- factor(y_factor, levels = y_levels, ordered = TRUE)
y_testable_int <- as.integer(y_ord_factor)           # 1..K (ordinal scale)
y_testable_cls <- y_testable_int - 1L                # 0..K-1 (CatBoost)

# Predictors (avoid leakage): drop permid + q8a
X_all <- model_df %>% select(-all_of(c("permid", leak_vars)))
feature_names <- names(X_all)

cat_features <- which(sapply(X_all, function(x) is.factor(x) || is.character(x))) - 1L

cat("Step 3 completed | Observed rows:", nrow(model_df),
    "| Classes:", n_classes,
    "| X columns:", ncol(X_all), "\n")

# ==============================================================================
# BLOCK 5 — Time-aware holdout split by wave (diagnostic only)
# ==============================================================================
uniq_waves <- sort(unique(model_df$wave))
test_n     <- max(1L, floor(0.20 * length(uniq_waves)))
test_waves <- tail(uniq_waves, test_n)

train_idx <- which(!model_df$wave %in% test_waves)
test_idx  <- which( model_df$wave %in% test_waves)

y_train_int <- y_testable_int[train_idx]
y_test_int  <- y_testable_int[test_idx]

y_train_cls <- y_testable_cls[train_idx]
y_test_cls  <- y_testable_cls[test_idx]

# Optional: inverse-frequency weights on training set (same as original V5)
freq_train <- table(y_train_cls)
w_map      <- median(as.numeric(freq_train)) / as.numeric(freq_train)
names(w_map) <- names(freq_train)

w_train <- as.numeric(w_map[as.character(y_train_cls)])
w_train[is.na(w_train)] <- 1

train_pool <- catboost.load_pool(
  X_all[train_idx, ], y_train_cls,
  cat_features = cat_features,
  weight       = w_train
)

test_pool <- catboost.load_pool(
  X_all[test_idx, ], y_test_cls,
  cat_features = cat_features
)

# ==============================================================================
# BLOCK 6 — One-shot grid search on holdout waves (diagnostic only)
#   * Same grid + early stopping as original V5
# ==============================================================================
grid <- expand.grid(
  depth         = c(5, 6, 8),
  learning_rate = c(0.03, 0.05),
  l2_leaf_reg   = c(3, 5)
)

results <- vector("list", nrow(grid))

for (i in seq_len(nrow(grid))) {
  params <- list(
    loss_function  = "MultiClass",
    eval_metric    = "MultiClass",
    iterations     = 3000,
    learning_rate  = grid$learning_rate[i],
    depth          = grid$depth[i],
    l2_leaf_reg    = grid$l2_leaf_reg[i],
    od_type        = "Iter",
    od_wait        = 50,
    random_seed    = 123,
    thread_count   = n_cores,
    use_best_model = TRUE,
    logging_level  = "Silent"
  )

  m <- catboost.train(train_pool, test_pool, params = params)
  probs_i <- catboost.predict(m, test_pool, prediction_type = "Probability")

  logloss_i <- safe_logloss(probs_i, y_test_cls, n_classes)

  best_iter_i <- tryCatch(catboost.get_best_iteration(m), error = function(e) NA_integer_)
  if (!is.na(best_iter_i)) best_iter_i <- best_iter_i + 1L
  if (is.na(best_iter_i))  best_iter_i <- params$iterations

  results[[i]] <- list(params = grid[i, ], best_iter = best_iter_i, logloss = logloss_i)

  cat(sprintf(
    "Grid %02d/%02d | Logloss %.5f | Iter %d | depth=%d lr=%.3f l2=%d\n",
    i, nrow(grid), logloss_i, best_iter_i,
    grid$depth[i], grid$learning_rate[i], grid$l2_leaf_reg[i]
  ))
}

best_id <- which.min(sapply(results, `[[`, "logloss"))
best    <- results[[best_id]]

cat("\nBest grid choice:\n")
print(best$params)
cat(sprintf("Best holdout logloss: %.5f | Best iteration: %d\n", best$logloss, best$best_iter))

# ==============================================================================
# BLOCK 7 — Final model (fit on ALL observed q8a)
# ==============================================================================
freq_all   <- table(y_testable_cls)
w_map_all  <- median(as.numeric(freq_all)) / as.numeric(freq_all)
names(w_map_all) <- names(freq_all)

w_all <- as.numeric(w_map_all[as.character(y_testable_cls)])
w_all[is.na(w_all)] <- 1

all_pool <- catboost.load_pool(
  X_all, y_testable_cls,
  cat_features = cat_features,
  weight       = w_all
)

final_params <- list(
  loss_function  = "MultiClass",
  eval_metric    = "MultiClass",
  iterations     = best$best_iter,
  learning_rate  = best$params$learning_rate,
  depth          = best$params$depth,
  l2_leaf_reg    = best$params$l2_leaf_reg,
  random_seed    = 123,
  thread_count   = n_cores,
  logging_level  = "Silent"
)

final_model <- catboost.train(all_pool, params = final_params)

# ==============================================================================
# BLOCK 8 — Holdout evaluation (ONLY Table 10 metrics)
#   A) Bin-space: exact / within-1 / QWK / macro-F1
#   B) Credibility: TVD overall, avg by wave, avg by country
#   C) € metrics on snapped midpoints + log1p metrics
#   D) Tail-risk diagnostics
# ==============================================================================

# Diagnostic refit (train -> early-stop on holdout), same as original V5
diag_params <- final_params
diag_params$iterations     <- 3000
diag_params$use_best_model <- TRUE
diag_params$od_type        <- "Iter"
diag_params$od_wait        <- 80

diag_model <- catboost.train(train_pool, test_pool, params = diag_params)

probs_test <- catboost.predict(diag_model, test_pool, prediction_type = "Probability")
if (is.null(dim(probs_test))) probs_test <- matrix(probs_test, ncol = n_classes, byrow = TRUE)

# Predicted bin = argmax probability (bin-space metrics)
pred_cls <- max.col(probs_test) - 1L
pred_int <- pred_cls + 1L

# ----- Bin-space metrics (Table 10) -----
acc_exact   <- mean(pred_int == y_test_int)
acc_within1 <- mean(abs(pred_int - y_test_int) <= 1L)
qwk         <- quadratic_weighted_kappa(y_test_int, pred_int, n_classes)

pred_factor <- factor(pred_int, levels = 1:n_classes, labels = y_levels, ordered = TRUE)
obs_factor  <- factor(y_test_int, levels = 1:n_classes, labels = y_levels, ordered = TRUE)
cm          <- confusionMatrix(pred_factor, obs_factor)

macro_f1 <- if (is.matrix(cm$byClass)) {
  mean(cm$byClass[, "F1"], na.rm = TRUE)
} else {
  unname(cm$byClass["F1"])
}

# ----- Credibility (TVD) metrics (Table 10) -----
eval_df <- model_df[test_idx, ] %>%
  mutate(
    q8a_true_int = y_test_int,
    q8a_pred_int = pred_int,
    q8a_true_bin = factor(y_test_int, levels = 1:n_classes, labels = y_levels, ordered = TRUE),
    q8a_pred_bin = factor(pred_int, levels = 1:n_classes, labels = y_levels, ordered = TRUE)
  )

overall    <- make_dist(eval_df, y_levels)
by_wave    <- make_dist(eval_df, y_levels, c("wave"))
by_country <- make_dist(eval_df, y_levels, c("d0"))

dist_all <- bind_rows(
  overall$dist    %>% mutate(scope = "overall", wave = NA_real_, d0 = NA_character_) %>% select(scope, wave, d0, everything()),
  by_wave$dist    %>% mutate(scope = "wave",   d0   = NA_character_)                %>% select(scope, wave, d0, everything()),
  by_country$dist %>% mutate(scope = "d0",     wave = NA_real_)                     %>% select(scope, wave, d0, everything())
)

tvd_all <- bind_rows(
  overall$tvd    %>% mutate(scope = "overall", wave = NA_real_, d0 = NA_character_) %>% select(scope, wave, d0, everything()),
  by_wave$tvd    %>% mutate(scope = "wave",   d0   = NA_character_)                 %>% select(scope, wave, d0, everything()),
  by_country$tvd %>% mutate(scope = "d0",     wave = NA_real_)                      %>% select(scope, wave, d0, everything())
)

write_dta(dist_all, path_out_dist)
write_dta(tvd_all,  path_out_tvd)

overall_tvd        <- overall$tvd$tvd[1]
avg_tvd_by_wave    <- mean(by_wave$tvd$tvd, na.rm = TRUE)
avg_tvd_by_country <- mean(by_country$tvd$tvd, na.rm = TRUE)

# ----- € metrics on snapped midpoints (Table 10) -----
true_code <- as.character(model_df$q8a[test_idx])
true_mid  <- as.numeric(MIDPOINTS_MAP[true_code])
mp_vec    <- as.numeric(MIDPOINTS_MAP[y_levels])

# Ensure probs_test columns align with y_levels (important for EV)
if (!is.null(colnames(probs_test)) && !all(colnames(probs_test) == y_levels)) {
  probs_test <- probs_test[, y_levels, drop = FALSE]
}

# Model I € prediction: expected midpoint from probabilities, then snap
pred_mid_ev   <- as.vector(probs_test %*% mp_vec)
pred_mid_snap <- snap_to_midpoint(pred_mid_ev, mp_vec)

metrics_eur <- list(
  rmse_eur        = rmse(pred_mid_snap, true_mid),
  mae_eur         = mae(pred_mid_snap, true_mid),
  medae_eur       = medae(pred_mid_snap, true_mid),
  spearman_eur    = spearman_fun(pred_mid_snap, true_mid),
  r2_eur          = r2(pred_mid_snap, true_mid),
  bias_mean_eur   = bias_mean(pred_mid_snap, true_mid),
  bias_median_eur = bias_median(pred_mid_snap, true_mid)
)

pred_log <- log1p(pred_mid_snap)
true_log <- log1p(true_mid)

metrics_log <- list(
  rmse_log1p        = rmse(pred_log, true_log),
  mae_log1p         = mae(pred_log, true_log),
  medae_log1p       = medae(pred_log, true_log),
  spearman_log1p    = spearman_fun(pred_log, true_log),
  r2_log1p          = r2(pred_log, true_log),
  bias_mean_log1p   = bias_mean(pred_log, true_log),
  bias_median_log1p = bias_median(pred_log, true_log)
)

# ----- Tail-risk diagnostics (Table 10) -----
tail <- tail_risk_metrics(
  true_int      = y_test_int,
  pred_int      = pred_int,
  true_mid      = true_mid,
  pred_mid_snap = pred_mid_snap
)

cat("\n================ TABLE 10 METRICS (HOLDOUT) — MODEL I (V5) ================\n")
cat(sprintf("Exact accuracy              : %.4f\n", acc_exact))
cat(sprintf("Within-1 accuracy           : %.4f\n", acc_within1))
cat(sprintf("Quadratic weighted kappa    : %.4f\n", qwk))
cat(sprintf("Macro F1                    : %.4f\n", macro_f1))
cat(sprintf("Overall TVD                 : %.4f\n", overall_tvd))
cat(sprintf("Avg TVD by wave             : %.4f\n", avg_tvd_by_wave))
cat(sprintf("Avg TVD by country          : %.4f\n", avg_tvd_by_country))
cat(sprintf("RMSE (€)                    : %.1f\n", metrics_eur$rmse_eur))
cat(sprintf("MAE (€)                     : %.1f\n", metrics_eur$mae_eur))
cat(sprintf("Median AE (€)               : %.1f\n", metrics_eur$medae_eur))
cat(sprintf("Spearman (€)                : %.4f\n", metrics_eur$spearman_eur))
cat(sprintf("R2 (€)                      : %.4f\n", metrics_eur$r2_eur))
cat(sprintf("Bias mean (€)               : %.1f\n", metrics_eur$bias_mean_eur))
cat(sprintf("Bias median (€)             : %.1f\n", metrics_eur$bias_median_eur))
cat(sprintf("RMSE (log1p)                : %.4f\n", metrics_log$rmse_log1p))
cat(sprintf("MAE (log1p)                 : %.4f\n", metrics_log$mae_log1p))
cat(sprintf("Median AE (log1p)           : %.4f\n", metrics_log$medae_log1p))
cat(sprintf("Spearman (log1p)            : %.4f\n", metrics_log$spearman_log1p))
cat(sprintf("R2 (log1p)                  : %.4f\n", metrics_log$r2_log1p))
cat(sprintf("Bias mean (log1p)           : %.4f\n", metrics_log$bias_mean_log1p))
cat(sprintf("Bias median (log1p)         : %.4f\n", metrics_log$bias_median_log1p))
cat(sprintf("Big mistake rate (|Δbin|>=2): %.6f\n", tail$big_mistake_rate))
cat(sprintf("AE p95 (€)                  : %.1f\n", tail$ae_p95_eur))
cat(sprintf("Top-bin underpred. rate     : %.6f\n", tail$topbin_under_rate))
cat(sprintf("Catastrophic flip rate      : %.6f\n", tail$catastrophic_rate))

# ==============================================================================
# BLOCK 9 — Probabilistic imputation on FULL dataset_A
#   * Missing q8a -> stochastic draw from predicted probabilities
# ==============================================================================
q8a_was_missing <- is.na(dataset_A$q8a)
miss_idx <- which(q8a_was_missing)

q8a_imputed_chr <- as.character(dataset_A$q8a)

if (length(miss_idx) > 0) {
  X_miss <- dataset_A[miss_idx, feature_names, drop = FALSE]
  miss_pool <- catboost.load_pool(X_miss, cat_features = cat_features)

  probs_miss <- catboost.predict(final_model, miss_pool, prediction_type = "Probability")
  if (is.null(dim(probs_miss))) probs_miss <- matrix(probs_miss, ncol = n_classes, byrow = TRUE)

  set.seed(123)  # keep deterministic stochastic imputations
  draw_one <- function(p) sample(y_levels, size = 1, prob = p)
  q8a_draws <- apply(probs_miss, 1, draw_one)

  q8a_imputed_chr[miss_idx] <- q8a_draws
}

dataset_A <- dataset_A %>%
  mutate(
    q8a          = factor(as.character(q8a), levels = y_levels, ordered = TRUE),
    q8a_imputed  = factor(q8a_imputed_chr,   levels = y_levels, ordered = TRUE),
    q8a_imp_flag = ifelse(q8a_was_missing, 1L, 0L)
  )

# ==============================================================================
# BLOCK 10 — Merge imputed midpoint back to the general SAFE dataset
# ==============================================================================
q8a_mid_tbl <- dataset_A %>%
  transmute(
    permid,
    wave,
    q8a_midpoints = as.numeric(MIDPOINTS_MAP[as.character(q8a_imputed)])
  )

dataset_general <- read_dta(path_in) %>%
  filter(
    d0 %in% c("BE","DK","DE","EE","IE","GR","ES","FR","HR","IT",
              "CY","LV","LT","LU","MT","NL","AT","PL","PT","RO",
              "SI","SK","FI","SE","BG","CZ","HU"),
    wave > 10,
    !wave %in% c(31, 33, 35)
  ) %>%
  arrange(permid, wave) %>%
  left_join(q8a_mid_tbl, by = c("permid", "wave"))

write_dta(dataset_general, path_out_data)
cat("\nSaved merged dataset:", path_out_data, "\n")

# ==============================================================================
# BLOCK 11 — Export Table 10 metrics table (CSV + LaTeX)
# ==============================================================================
test_waves_label <- paste(test_waves, collapse = ", ")

metrics_tbl <- tribble(
  ~section, ~metric, ~value,
  "Setup", "Holdout waves (by wave)", test_waves_label,
  "Setup", "N holdout observations", as.character(length(test_idx)),

  "Bin-space", "Exact accuracy", sprintf("%.4f", acc_exact),
  "Bin-space", "Within-1 accuracy", sprintf("%.4f", acc_within1),
  "Bin-space", "Quadratic weighted kappa", sprintf("%.4f", qwk),
  "Bin-space", "Macro F1", sprintf("%.4f", macro_f1),

  "Distribution (TVD)", "Overall TVD (bins)", sprintf("%.4f", overall_tvd),
  "Distribution (TVD)", "Avg. TVD by wave", sprintf("%.4f", avg_tvd_by_wave),
  "Distribution (TVD)", "Avg. TVD by country", sprintf("%.4f", avg_tvd_by_country),

  "€-space (snapped midpoints)", "RMSE (€)", sprintf("%.1f", metrics_eur$rmse_eur),
  "€-space (snapped midpoints)", "MAE (€)", sprintf("%.1f", metrics_eur$mae_eur),
  "€-space (snapped midpoints)", "Median AE (€)", sprintf("%.1f", metrics_eur$medae_eur),
  "€-space (snapped midpoints)", "Spearman", sprintf("%.4f", metrics_eur$spearman_eur),
  "€-space (snapped midpoints)", "R2", sprintf("%.4f", metrics_eur$r2_eur),
  "€-space (snapped midpoints)", "Bias mean", sprintf("%.1f", metrics_eur$bias_mean_eur),
  "€-space (snapped midpoints)", "Bias median", sprintf("%.1f", metrics_eur$bias_median_eur),

  "log1p(€) (snapped midpoints)", "RMSE (log)", sprintf("%.4f", metrics_log$rmse_log1p),
  "log1p(€) (snapped midpoints)", "MAE (log)", sprintf("%.4f", metrics_log$mae_log1p),
  "log1p(€) (snapped midpoints)", "Median AE (log)", sprintf("%.4f", metrics_log$medae_log1p),
  "log1p(€) (snapped midpoints)", "Spearman (log)", sprintf("%.4f", metrics_log$spearman_log1p),
  "log1p(€) (snapped midpoints)", "R2 (log)", sprintf("%.4f", metrics_log$r2_log1p),
  "log1p(€) (snapped midpoints)", "Bias mean (log)", sprintf("%.4f", metrics_log$bias_mean_log1p),
  "log1p(€) (snapped midpoints)", "Bias median (log)", sprintf("%.4f", metrics_log$bias_median_log1p),

  "Tail-risk (snapped midpoints)", "Big mistake rate (|Δbin| ≥ 2)", sprintf("%.8f", tail$big_mistake_rate),
  "Tail-risk (snapped midpoints)", "AE p95 (€)", sprintf("%.1f", tail$ae_p95_eur),
  "Tail-risk (snapped midpoints)", "Top-bin underprediction rate", sprintf("%.8f", tail$topbin_under_rate),
  "Tail-risk (snapped midpoints)", "Catastrophic flip rate (bottom↔top)", sprintf("%.8f", tail$catastrophic_rate)
)

readr::write_csv(metrics_tbl, path_out_metrics_csv)

tex_lines <- c(
  "\\begin{table}[!htbp]\\centering",
  "\\caption{Holdout comparison metrics for Method B (Model I: CatBoost MultiClass)}",
  "\\label{tab:methodB_table10_modelI}",
  "\\begin{tabular}{lll}",
  "\\hline",
  "Section & Metric & Value \\\\",
  "\\hline",
  paste0(metrics_tbl$section, " & ", metrics_tbl$metric, " & ", metrics_tbl$value, " \\\\") ,
  "\\hline",
  "\\end{tabular}",
  "\\end{table}"
)

writeLines(tex_lines, path_out_metrics_tex)

cat("\nSaved Table 10 metrics table:\n")
cat(" - CSV:", path_out_metrics_csv, "\n")
cat(" - TEX:", path_out_metrics_tex, "\n")

t1 <- Sys.time()
cat("\nTOTAL runtime:\n")
print(t1 - t0)
