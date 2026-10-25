################################################################################
# Method B — Model I
# CatBoost MultiClass on Q8A bins + probabilistic imputation
#
# Pipeline:
#   0)  Paths + runtime knobs
#   1)  Helpers (metrics + credibility tables + coerce_wave_to_int)
#   2)  Load + clean data
#       2.1  Load raw data, aggregate quarterly waves, deduplicate
#       2.2  Apply sample filters + q8a cleaning
#   3)  Feature selection + NA encoding
#   4)  Model matrix prep
#   5)  Time-aware holdout split
#   6)  Grid search (diagnostic)
#   7)  Final model (fit on all observed q8a)
#   8)  Holdout evaluation
#   9)  Probabilistic imputation
#   10) Merge imputed midpoints back to general SAFE dataset
#   11) Export Table 10 metrics (CSV + LaTeX)
#   12) Heatmap: log(1+gap) error by country × sector/size
#
# Revision log (details: docs/Changes_vs_thesis.docx)
#   C1  Class weights: switch CLASS_WEIGHT_MODE ("none" default / "corrected" /
#       "legacy"); prior-shift correction of probabilities before evaluation and
#       before the imputation draw.
#   C2  Distributional credibility (TVD, mean shift) computed on expected shares
#       = mean of predicted probabilities; argmax-based TVD (thesis) kept as legacy.
#   C3  Sector/size grouping variables corrected (d3_rec / d1_rec).
################################################################################

library(ggplot2)
library(patchwork)
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

path_raw        <- ""
path_out_data   <- ""

path_out_dist <- ""
path_out_tvd  <- ""

path_out_metrics_csv <- ""
path_out_metrics_tex <- ""

path_out_heatmap <- ""

n_cores <- max(1L, detectCores())

# Class-weighting mode (change C1 — see docs/Changes_vs_thesis.docx)
#   "none"      : no class weights. Predicted probabilities stay calibrated, so the
#                 probabilistic draw in Block 9 reproduces the true bracket shares.
#                 (default)
#   "corrected" : class weights w_k = median(n)/n_k in training (as in the thesis),
#                 then prior-shift correction p_k <- (p_k / w_k) / sum_j (p_j / w_j)
#                 before any evaluation or draw.
#   "legacy"    : thesis setting — class weights, no correction (reproduces thesis).
CLASS_WEIGHT_MODE <- "none"
stopifnot(CLASS_WEIGHT_MODE %in% c("none", "corrected", "legacy"))

BIN_LEVELS <- c("1", "2", "5", "6", "4")

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

# Coerce wave variable to integer regardless of storage type
# (numeric → direct cast; haven_labelled → via as_factor; other → character parse)
coerce_wave_to_int <- function(x) {
  if (is.numeric(x)) return(as.integer(x))
  if (!is.null(attr(x, "labels"))) {
    s <- as.character(haven::as_factor(x))
    if (all(grepl("^\\s*-?\\d+\\s*$", na.omit(s)))) return(as.integer(s))
  }
  s <- as.character(x)
  if (all(grepl("^\\s*-?\\d+\\s*$", na.omit(s)))) return(as.integer(s))
  as.integer(as.numeric(s))
}

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
  eps <- 1e-15
  if (is.null(dim(probs))) probs <- matrix(probs, ncol = n_classes, byrow = TRUE)
  idx    <- cbind(seq_len(nrow(probs)), y_true_cls + 1L)
  p_true <- pmax(probs[idx], eps)
  -mean(log(p_true))
}

quadratic_weighted_kappa <- function(truth_int, pred_int, k) {
  O <- as.matrix(table(factor(truth_int, levels = 1:k),
                       factor(pred_int,  levels = 1:k)))
  N <- sum(O)
  if (N == 0) return(NA_real_)
  
  E <- outer(rowSums(O), colSums(O)) / N
  W <- outer(1:k, 1:k, function(i, j) ((i - j)^2) / ((k - 1)^2))
  1 - (sum(W * O) / sum(W * E))
}

rmse        <- function(p, y) sqrt(mean((p - y)^2, na.rm = TRUE))
mae         <- function(p, y) mean(abs(p - y), na.rm = TRUE)
medae       <- function(p, y) median(abs(p - y), na.rm = TRUE)
bias_mean   <- function(p, y) mean(p - y, na.rm = TRUE)
bias_median <- function(p, y) median(p - y, na.rm = TRUE)

r2 <- function(p, y) {
  ok <- is.finite(p) & is.finite(y)
  y  <- y[ok]; p <- p[ok]
  if (length(y) < 2) return(NA_real_)
  ss_res <- sum((y - p)^2)
  ss_tot <- sum((y - mean(y))^2)
  if (ss_tot == 0) return(NA_real_)
  1 - ss_res / ss_tot
}

spearman_fun <- function(p, y) {
  suppressWarnings(cor(p, y, method = "spearman", use = "complete.obs"))
}

to_numeric_safe <- function(x) {
  if (inherits(x, "haven_labelled")) return(as.numeric(x))
  if (is.factor(x))    return(suppressWarnings(as.numeric(as.character(x))))
  if (is.character(x)) return(suppressWarnings(as.numeric(x)))
  suppressWarnings(as.numeric(x))
}

snap_to_midpoint <- function(x, mps) {
  mps[max.col(-abs(outer(x, mps, "-")), ties.method = "first")]
}

make_dist <- function(df, y_levels, group_vars = NULL) {
  if (is.null(group_vars) || length(group_vars) == 0) {
    df         <- df %>% mutate(grp = "ALL")
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
      n   = sum(n_true),
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

# --- Change C1: class weights and prior-shift correction ------------------------
# Returns the named weight map w_k = median(n)/n_k (names = class codes 0..K-1),
# or NULL when CLASS_WEIGHT_MODE == "none".
class_weight_map <- function(y_cls, mode = CLASS_WEIGHT_MODE) {
  if (mode == "none") return(NULL)
  freq  <- table(y_cls)
  w_map <- median(as.numeric(freq)) / as.numeric(freq)
  names(w_map) <- names(freq)
  w_map
}

# Row weights for catboost.load_pool (NULL = unweighted)
row_weights <- function(y_cls, w_map) {
  if (is.null(w_map)) return(NULL)
  w <- as.numeric(w_map[as.character(y_cls)])
  w[is.na(w)] <- 1
  w
}

# Undo the class weighting on predicted probabilities (only in "corrected" mode).
# Weighting class k by w_k multiplies its odds by w_k; dividing by w_k and
# renormalising restores probabilities on the original class frequencies.
correct_probs <- function(probs, w_map, n_classes, mode = CLASS_WEIGHT_MODE) {
  if (is.null(dim(probs))) probs <- matrix(probs, ncol = n_classes, byrow = TRUE)
  if (mode != "corrected" || is.null(w_map)) return(probs)
  w   <- as.numeric(w_map[as.character(0:(n_classes - 1L))])
  w[is.na(w)] <- 1
  adj <- sweep(probs, 2, w, "/")
  adj / rowSums(adj)
}

# --- Change C2: TVD on expected shares (mean of predicted probabilities) -----
# For each group g and bracket k:
#   share_pred_gk = (1/n_g) * sum_{i in g} p_hat_ik      (expected share under the draw)
#   share_true_gk = (1/n_g) * sum_{i in g} 1{k_i = k}
#   TVD_g         = 0.5 * sum_k |share_pred_gk - share_true_gk|
# This evaluates exactly the object produced by the probabilistic imputation
# (Block 9) without Monte Carlo noise; the thesis computed TVD on argmax labels.
make_dist_prob <- function(df, probs, y_levels, group_vars = NULL) {
  K     <- length(y_levels)
  pcols <- paste0("p_", seq_len(K))
  colnames(probs) <- pcols

  if (is.null(group_vars) || length(group_vars) == 0) {
    df         <- df %>% mutate(grp = "ALL")
    group_vars <- "grp"
  }
  d <- bind_cols(
    df %>% select(all_of(group_vars), q8a_true_int),
    as.data.frame(probs)
  )

  pred <- d %>%
    group_by(across(all_of(group_vars))) %>%
    summarise(across(all_of(pcols), mean), .groups = "drop") %>%
    tidyr::pivot_longer(all_of(pcols), names_to = "k", values_to = "share_pred") %>%
    mutate(k = as.integer(sub("p_", "", k)))

  true <- d %>%
    count(across(all_of(group_vars)), k = q8a_true_int, name = "n_true")

  dist <- tidyr::crossing(
    d %>% distinct(across(all_of(group_vars))),
    k = seq_len(K)
  ) %>%
    left_join(pred, by = c(group_vars, "k")) %>%
    left_join(true, by = c(group_vars, "k")) %>%
    mutate(n_true = tidyr::replace_na(n_true, 0L)) %>%
    group_by(across(all_of(group_vars))) %>%
    mutate(
      share_true = n_true / sum(n_true),
      diff_share = share_pred - share_true
    ) %>%
    ungroup() %>%
    mutate(bin = factor(y_levels[k], levels = y_levels, ordered = TRUE))

  tvd_tbl <- dist %>%
    group_by(across(all_of(group_vars))) %>%
    summarise(
      n             = sum(n_true),
      tvd           = 0.5 * sum(abs(share_pred - share_true)),
      mean_true_int = sum(k * share_true),
      mean_pred_int = sum(k * share_pred),
      mean_shift    = mean_pred_int - mean_true_int,
      .groups = "drop"
    )

  list(dist = dist %>% select(-k) %>% relocate(bin, .after = all_of(group_vars)),
       tvd  = tvd_tbl)
}

tail_risk_metrics <- function(true_int, pred_int, true_mid, pred_mid_snap) {
  Kmax   <- max(true_int, na.rm = TRUE)
  Kmin   <- min(true_int, na.rm = TRUE)
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
# BLOCK 2 — Load + clean data
# ==============================================================================

# --- 2.1  Load raw data, aggregate quarterly waves, deduplicate ---------------
#
#   Semiannual regime (waves 1–29): wave_agg = wave (unchanged)
#   Quarterly regime  (waves 30+):  consecutive quarter-pairs collapsed into
#                                   sequential semiannual indices:
#     30–31 (2024 Q1–Q2) → 30 (2024 H1)
#     32–33 (2024 Q3–Q4) → 31 (2024 H2)
#     34–35 (2025 Q1–Q2) → 32 (2025 H1)
#     36–37 (2025 Q3–Q4) → 33 (2025 H2)
#
#   Deduplication: one row per (permid, wave_agg); for the ~980 firms appearing
#   in both quarters of the same half-year, the earlier quarter's record is kept.

dat_agg <- read_dta(path_raw) %>%
  mutate(
    permid    = as.character(permid),
    orig_wave = coerce_wave_to_int(wave),
    wave_agg  = case_when(
      orig_wave %in% c(30L, 31L) ~ 30L,
      orig_wave %in% c(32L, 33L) ~ 31L,
      orig_wave %in% c(34L, 35L) ~ 32L,
      orig_wave %in% c(36L, 37L) ~ 33L,
      TRUE                       ~ orig_wave
    )
  ) %>%
  arrange(permid, wave_agg, orig_wave) %>%
  distinct(permid, wave_agg, .keep_all = TRUE)

attributes(dat_agg)[c("label", "datalabel", "notes")] <- NULL
# drop any leftover Stata file-level label if present

# --- 2.2  Apply sample filters + q8a cleaning  --------------------------------
#   Uses dat_agg directly — no additional disk read.
#   Same country/wave/financing filters and q8a recoding as original V5.

dataset_A_raw <- dat_agg %>%
  filter(
    d0 %in% c("BE","DK","DE","EE","IE","GR","ES","FR","HR","IT",
              "CY","LV","LT","LU","MT","NL","AT","PL","PT","RO",
              "SI","SK","FI","SE","BG","CZ","HU"),
    wave_agg > 10,
    q7a_a %in% c(1, 2) | q32 %in% c(1, 2, 3, 4, 5, 6)
  ) %>%
  select(
    -c(wgtcommon, wgtentr, wgtoldentr, wgtoldcommon,
       ida, experiment_1, experiment_2, intdate),
    -any_of("q8a_rec")
  ) %>%
  arrange(permid, wave_agg)

dataset_A <- dataset_A_raw %>%
  mutate(
    across(
      -c(permid, wave_agg),
      ~ if (inherits(.x, "haven_labelled")) haven::as_factor(.x) else as.factor(.x)
    ),
    q8a = as.character(q8a),
    q8a = stringr::str_trim(q8a),
    q8a = dplyr::na_if(q8a, ""),
    q8a = stringr::str_extract(q8a, "\\d+"),
    q8a = ifelse(q8a == "9", NA_character_, q8a),
    q8a = factor(q8a)
  )

cat("2.2 done | Rows:", nrow(dataset_A),
    "| Columns:", ncol(dataset_A),
    "| Missing q8a:", sum(is.na(dataset_A$q8a)), "\n")

# ==============================================================================
# BLOCK 3 — Feature selection + NA encoding (hybrid rule)
# ==============================================================================
id_vars    <- c("permid", "wave_agg", "d0")
target_var <- "q8a"
leak_vars  <- c("q8a")

hard_drop <- names(which(colMeans(is.na(dataset_A)) > 0.90))
dataset_A <- dataset_A %>% select(-all_of(hard_drop))

candidate_vars <- setdiff(names(dataset_A), c(id_vars, leak_vars))
vars_to_keep   <- candidate_vars[
  sapply(candidate_vars, function(v) missingness_informative(dataset_A, v, target_var))
]
vars_to_keep <- union(vars_to_keep, intersect("q7a_a", names(dataset_A)))

dataset_A <- dataset_A %>%
  select(all_of(c(id_vars, target_var, vars_to_keep)))

pred_cols <- setdiff(names(dataset_A), c(id_vars, leak_vars))
dataset_A <- dataset_A %>%
  mutate(across(all_of(pred_cols),
                ~ if (is.factor(.x)) forcats::fct_explicit_na(.x, "MISSING") else .x))

cat("Step 3 completed | Columns after missingness control:", ncol(dataset_A), "\n")

# ==============================================================================
# BLOCK 4 — Model matrix prep (observed q8a only)
# ==============================================================================
model_df <- dataset_A %>% filter(!is.na(.data[[target_var]]))

y_factor      <- as.factor(model_df[[target_var]])
y_levels      <- intersect(BIN_LEVELS, levels(y_factor))
n_classes     <- length(y_levels)

y_ord_factor   <- factor(y_factor, levels = y_levels, ordered = TRUE)
y_testable_int <- as.integer(y_ord_factor)
y_testable_cls <- y_testable_int - 1L

X_all         <- model_df %>% select(-all_of(c("permid", leak_vars)))
feature_names <- names(X_all)
cat_features  <- which(sapply(X_all, function(x) is.factor(x) || is.character(x))) - 1L

cat("Step 4 completed | Observed rows:", nrow(model_df),
    "| Classes:", n_classes,
    "| X columns:", ncol(X_all), "\n")

# ==============================================================================
# BLOCK 5 — Time-aware holdout split by wave_agg (diagnostic only)
# ==============================================================================
uniq_waves <- sort(unique(model_df$wave_agg))
test_n     <- max(1L, floor(0.20 * length(uniq_waves)))
test_waves <- tail(uniq_waves, test_n)

train_idx <- which(!model_df$wave_agg %in% test_waves)
test_idx  <- which( model_df$wave_agg %in% test_waves)

y_train_int <- y_testable_int[train_idx]
y_test_int  <- y_testable_int[test_idx]
y_train_cls <- y_testable_cls[train_idx]
y_test_cls  <- y_testable_cls[test_idx]

# Change C1: class weights only if CLASS_WEIGHT_MODE != "none"
w_map   <- class_weight_map(y_train_cls)
w_train <- row_weights(y_train_cls, w_map)

train_pool <- catboost.load_pool(
  X_all[train_idx, ], y_train_cls,
  cat_features = cat_features, weight = w_train
)
test_pool <- catboost.load_pool(
  X_all[test_idx, ], y_test_cls,
  cat_features = cat_features
)

# ==============================================================================
# BLOCK 6 — One-shot grid search on holdout waves (diagnostic only)
# ==============================================================================
grid <- expand.grid(
  depth         = c(5, 6, 8),
  learning_rate = c(0.03, 0.05),
  l2_leaf_reg   = c(3, 5)
)

results <- vector("list", nrow(grid))

for (i in seq_len(nrow(grid))) {
  params <- list(
    loss_function       = "MultiClass",
    eval_metric         = "MultiClass",
    iterations          = 3000L,
    learning_rate       = grid$learning_rate[i],
    depth               = grid$depth[i],
    l2_leaf_reg         = grid$l2_leaf_reg[i],
    od_type             = "Iter",
    od_wait             = 100,
    random_seed         = 123,
    thread_count        = n_cores,
    use_best_model      = TRUE,
    logging_level       = "Silent",
    allow_writing_files = FALSE
  )
  
  m         <- catboost.train(train_pool, test_pool, params = params)
  probs_i   <- catboost.predict(m, test_pool, prediction_type = "Probability")
  probs_i   <- correct_probs(probs_i, w_map, n_classes)   # change C1
  logloss_i <- safe_logloss(probs_i, y_test_cls, n_classes)
  
  best_iter_i <- tryCatch(catboost.get_best_iteration(m) + 1L,
                          error = function(e) params$iterations)
  
  results[[i]] <- list(params = grid[i, ], best_iter = best_iter_i, logloss = logloss_i)
  
  cat(sprintf("Grid %02d/%02d | Logloss %.5f | Iter %d | depth=%d lr=%.3f l2=%d\n",
              i, nrow(grid), logloss_i, best_iter_i,
              grid$depth[i], grid$learning_rate[i], grid$l2_leaf_reg[i]))
  
  rm(m, probs_i); gc()
}

best_id <- which.min(sapply(results, `[[`, "logloss"))
best    <- results[[best_id]]

cat("\nBest grid choice:\n"); print(best$params)
cat(sprintf("Best holdout logloss: %.5f | Best iteration: %d\n", best$logloss, best$best_iter))

# ==============================================================================
# BLOCK 7 — Final model (fit on ALL observed q8a)
# ==============================================================================
# Change C1: class weights only if CLASS_WEIGHT_MODE != "none"
w_map_all <- class_weight_map(y_testable_cls)
w_all     <- row_weights(y_testable_cls, w_map_all)

all_pool <- catboost.load_pool(
  X_all, y_testable_cls,
  cat_features = cat_features, weight = w_all
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
# BLOCK 8 — Holdout evaluation
# ==============================================================================
diag_params <- modifyList(final_params, list(
  iterations     = 3000L,
  use_best_model = TRUE,
  od_type        = "Iter",
  od_wait        = 100
))

diag_model <- catboost.train(train_pool, test_pool, params = diag_params)

probs_test <- catboost.predict(diag_model, test_pool, prediction_type = "Probability")
probs_test <- correct_probs(probs_test, w_map, n_classes)   # change C1 (no-op unless "corrected")

# Argmax labels: used ONLY for individual-level accuracy metrics (accuracy, QWK, F1,
# tail risk). Aggregate/distributional metrics use the probabilities (change C2).
pred_cls <- max.col(probs_test) - 1L
pred_int <- pred_cls + 1L

# Bin-space metrics
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

# Credibility (TVD) metrics
eval_df <- model_df[test_idx, ] %>%
  mutate(
    q8a_true_int = y_test_int,
    q8a_pred_int = pred_int,
    q8a_true_bin = factor(y_test_int, levels = 1:n_classes, labels = y_levels, ordered = TRUE),
    q8a_pred_bin = factor(pred_int,   levels = 1:n_classes, labels = y_levels, ordered = TRUE)
  )

# Change C3: grouping variables. In SAFE, D1 = size (employees) and D3 = sector
# (thesis, Ch. 2). The thesis script grouped "sector" by d1 and "size" by d2,
# so its "TVD by sector" was in fact TVD by size class. The recoded versions
# (d1_rec, d3_rec — also used by the heatmap in Block 12) are preferred.
SECTOR_VAR <- intersect(c("d3_rec", "d3"), names(dataset_A_raw))[1]
SIZE_VAR   <- intersect(c("d1_rec", "d1"), names(dataset_A_raw))[1]
extra_grp   <- na.omit(c(SECTOR_VAR, SIZE_VAR))
missing_grp <- setdiff(extra_grp, names(eval_df))
if (length(missing_grp) > 0) {
  eval_df <- eval_df %>%
    left_join(dataset_A_raw %>% select(permid, wave_agg, all_of(missing_grp)),
              by = c("permid", "wave_agg"))
}

group_specs <- list(
  overall  = NULL,
  wave_agg = "wave_agg",
  d0       = "d0",
  sector   = if (!is.na(SECTOR_VAR)) SECTOR_VAR else NA,
  size     = if (!is.na(SIZE_VAR))   SIZE_VAR   else NA
)
group_specs <- group_specs[!sapply(group_specs, function(g) length(g) == 1 && is.na(g))]

# Change C2: primary = expected shares (mean of probabilities);
#            legacy  = argmax labels (thesis definition), kept for comparison.
run_dist <- function(fun_type) {
  res <- lapply(group_specs, function(g) {
    if (fun_type == "expected") make_dist_prob(eval_df, probs_test, y_levels, g)
    else                        make_dist(eval_df, y_levels, g)
  })
  names(res) <- names(group_specs)
  res
}
dist_expected <- run_dist("expected")
dist_argmax   <- run_dist("argmax")

stack_scope <- function(res, what, method) {
  bind_rows(lapply(names(res), function(s) {
    x <- res[[s]][[what]] %>% mutate(scope = s, tvd_method = method)
    # harmonise group columns into one 'group' column for export
    gcols <- intersect(names(x), c("grp", "wave_agg", "d0", SECTOR_VAR, SIZE_VAR))
    x %>% mutate(group = if (length(gcols)) as.character(.data[[gcols[1]]]) else "ALL") %>%
      select(-any_of(gcols))
  })) %>% relocate(tvd_method, scope, group)
}

dist_all <- bind_rows(stack_scope(dist_expected, "dist", "expected_share"),
                      stack_scope(dist_argmax,   "dist", "argmax"))
tvd_all  <- bind_rows(stack_scope(dist_expected, "tvd",  "expected_share"),
                      stack_scope(dist_argmax,   "tvd",  "argmax"))

write_dta(dist_all %>% mutate(bin = as.character(bin)) %>%
            select(-any_of(c("q8a_true_bin", "q8a_pred_bin"))), path_out_dist)
write_dta(tvd_all, path_out_tvd)

summ_tvd <- function(res) {
  avg <- function(s, col) if (!is.null(res[[s]])) mean(res[[s]]$tvd[[col]], na.rm = TRUE) else NA_real_
  list(
    overall_tvd        = res$overall$tvd$tvd[1],
    avg_tvd_wave       = avg("wave_agg", "tvd"),
    avg_tvd_country    = avg("d0", "tvd"),
    avg_tvd_sector     = avg("sector", "tvd"),
    avg_tvd_size       = avg("size", "tvd"),
    overall_shift      = res$overall$tvd$mean_shift[1],
    avg_shift_wave     = avg("wave_agg", "mean_shift"),
    avg_shift_country  = avg("d0", "mean_shift"),
    avg_shift_sector   = avg("sector", "mean_shift"),
    avg_shift_size     = avg("size", "mean_shift")
  )
}
tvd_E <- summ_tvd(dist_expected)   # primary (change C2)
tvd_A <- summ_tvd(dist_argmax)     # thesis definition

# Expected vs observed bracket shares, overall (calibration-in-the-large check)
cat("\n--- Overall bracket shares: observed vs expected (mean prob.) vs argmax ---\n")
print(
  dist_expected$overall$dist %>%
    select(bin, share_true, share_pred_expected = share_pred) %>%
    left_join(dist_argmax$overall$dist %>% select(bin, share_pred_argmax = share_pred),
              by = "bin")
)

# €-space metrics
true_code     <- as.character(model_df$q8a[test_idx])
true_mid      <- as.numeric(MIDPOINTS_MAP[true_code])
mp_vec        <- as.numeric(MIDPOINTS_MAP[y_levels])

if (!is.null(colnames(probs_test)) && !all(colnames(probs_test) == y_levels))
  probs_test <- probs_test[, y_levels, drop = FALSE]

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

tail <- tail_risk_metrics(y_test_int, pred_int, true_mid, pred_mid_snap)

# Euro financing gap metrics
gap_base <- model_df[test_idx, c("permid", "wave_agg")] %>%
  left_join(dataset_A_raw, by = c("permid", "wave_agg")) %>%
  mutate(
    loan_demanded_true = true_mid,
    loan_demanded_pred = pred_mid_snap,
    q7a_a_num      = to_numeric_safe(q7a_a),
    q7b_a_num      = to_numeric_safe(q7b_a),
    q32_num        = to_numeric_safe(q32),
    vulnerable_num = if ("vulnerable" %in% names(.)) to_numeric_safe(vulnerable) else NA_real_,
    q32_clean      = if_else(q7a_a_num == 1 | !is.na(q7b_a_num), NA_real_, q32_num),
    loan_need_factor_a = case_when(
      q7a_a_num == 3 | q7b_a_num == 1 ~ 0,
      q7b_a_num == 5                  ~ 0.2,
      q7b_a_num == 6                  ~ 0.8,
      q7b_a_num %in% c(3, 4)          ~ 1,
      TRUE                            ~ NA_real_
    ),
    loan_need_factor_b = case_when(
      q7a_a_num == 2                                       ~ 1,
      !is.na(q32_clean) & q32_clean >= 1 & q32_clean <= 6 ~ 1,
      TRUE                                                 ~ loan_need_factor_a
    ),
    loan_gap_a_true = case_when(
      loan_need_factor_a == 0 ~ 0,
      loan_need_factor_a >  0 ~ loan_need_factor_a * loan_demanded_true,
      TRUE                    ~ NA_real_
    ),
    loan_gap_a_pred = case_when(
      loan_need_factor_a == 0 ~ 0,
      loan_need_factor_a >  0 ~ loan_need_factor_a * loan_demanded_pred,
      TRUE                    ~ NA_real_
    ),
    loan_gap_b_true = case_when(
      loan_need_factor_b == 0 ~ 0,
      loan_need_factor_b >  0 ~ loan_need_factor_b * loan_demanded_true,
      TRUE                    ~ NA_real_
    ),
    loan_gap_b_pred = case_when(
      loan_need_factor_b == 0 ~ 0,
      loan_need_factor_b >  0 ~ loan_need_factor_b * loan_demanded_pred,
      TRUE                    ~ NA_real_
    )
  )

gap_eval_df <- if ("vulnerable" %in% names(gap_base)) {
  gap_base %>% filter(is.na(vulnerable_num) | vulnerable_num == 0)
} else {
  gap_base
}

relative_error <- function(p, y) {
  rel <- (p - y) / y
  rel[!is.finite(rel)] <- NA_real_
  rel
}

gap_metrics <- function(p, y) {
  rel <- relative_error(p, y)
  list(
    rmse          = rmse(p, y),
    mae           = mae(p, y),
    medae         = medae(p, y),
    spearman      = spearman_fun(p, y),
    r2            = r2(p, y),
    bias_mean     = bias_mean(p, y),
    bias_median   = bias_median(p, y),
    relerr_mean   = if (all(is.na(rel))) NA_real_ else mean(rel,   na.rm = TRUE),
    relerr_median = if (all(is.na(rel))) NA_real_ else median(rel, na.rm = TRUE),
    n             = sum(is.finite(p) & is.finite(y))
  )
}

metrics_gap_a <- gap_metrics(gap_eval_df$loan_gap_a_pred, gap_eval_df$loan_gap_a_true)
metrics_gap_b <- gap_metrics(gap_eval_df$loan_gap_b_pred, gap_eval_df$loan_gap_b_true)

cat("\n================ TABLE: Out-of-time performance on ordered brackets (bin-space) ================\n")

cat("\n--- Accuracy metrics ---\n")
cat(sprintf("Exact accuracy              : %.4f\n", acc_exact))
cat(sprintf("Within-1 accuracy           : %.4f\n", acc_within1))
cat(sprintf("Quadratic weighted kappa    : %.4f\n", qwk))
cat(sprintf("Macro F1                    : %.4f\n", macro_f1))

cat(sprintf("\nClass-weight mode           : %s\n", CLASS_WEIGHT_MODE))
cat("\n--- Distributional credibility: expected shares (mean prob.) [primary] ---\n")
cat(sprintf("Overall TVD                 : %.4f\n", tvd_E$overall_tvd))
cat(sprintf("Avg TVD by wave             : %.4f\n", tvd_E$avg_tvd_wave))
cat(sprintf("Avg TVD by country          : %.4f\n", tvd_E$avg_tvd_country))
cat(sprintf("Avg TVD by sector (%s)  : %.4f\n", SECTOR_VAR, tvd_E$avg_tvd_sector))
cat(sprintf("Avg TVD by size (%s)    : %.4f\n", SIZE_VAR,   tvd_E$avg_tvd_size))
cat(sprintf("Overall mean shift (bins)   : %.4f\n", tvd_E$overall_shift))
cat(sprintf("Avg mean shift by country   : %.4f\n", tvd_E$avg_shift_country))
cat(sprintf("Avg mean shift by sector    : %.4f\n", tvd_E$avg_shift_sector))
cat(sprintf("Avg mean shift by size      : %.4f\n", tvd_E$avg_shift_size))

cat("\n--- Distributional credibility: argmax labels [thesis definition] ---\n")
cat(sprintf("Overall TVD                 : %.4f\n", tvd_A$overall_tvd))
cat(sprintf("Avg TVD by wave             : %.4f\n", tvd_A$avg_tvd_wave))
cat(sprintf("Avg TVD by country          : %.4f\n", tvd_A$avg_tvd_country))
cat(sprintf("Avg TVD by sector           : %.4f\n", tvd_A$avg_tvd_sector))
cat(sprintf("Avg TVD by size             : %.4f\n", tvd_A$avg_tvd_size))
cat(sprintf("Overall mean shift (bins)   : %.4f\n", tvd_A$overall_shift))

cat("\n--- Tail-risk diagnostics ---\n")
cat(sprintf("Big mistake rate (|Δbin|≥2) : %.6f\n", tail$big_mistake_rate))
cat(sprintf("Top-bin underprediction rate: %.6f\n", tail$topbin_under_rate))
cat(sprintf("Catastrophic flip rate      : %.6f\n", tail$catastrophic_rate))

cat("\n================ TABLE: Out-of-time performance on euro-valued midpoints ================\n")

cat("\n--- €-space (snapped midpoints) ---\n")
cat(sprintf("RMSE (€)                    : %.1f\n",  metrics_eur$rmse_eur))
cat(sprintf("MAE (€)                     : %.1f\n",  metrics_eur$mae_eur))
cat(sprintf("Median AE (€)               : %.1f\n",  metrics_eur$medae_eur))
cat(sprintf("Spearman (€)                : %.4f\n", metrics_eur$spearman_eur))
cat(sprintf("R2 (€)                      : %.4f\n", metrics_eur$r2_eur))
cat(sprintf("Bias mean (€)               : %.1f\n",  metrics_eur$bias_mean_eur))
cat(sprintf("Bias median (€)             : %.1f\n",  metrics_eur$bias_median_eur))

cat("\n--- log1p(€)-space ---\n")
cat(sprintf("RMSE (log1p)                : %.4f\n", metrics_log$rmse_log1p))
cat(sprintf("MAE (log1p)                 : %.4f\n", metrics_log$mae_log1p))
cat(sprintf("Median AE (log1p)           : %.4f\n", metrics_log$medae_log1p))
cat(sprintf("Spearman (log1p)            : %.4f\n", metrics_log$spearman_log1p))
cat(sprintf("R2 (log1p)                  : %.4f\n", metrics_log$r2_log1p))
cat(sprintf("Bias mean (log1p)           : %.4f\n", metrics_log$bias_mean_log1p))
cat(sprintf("Bias median (log1p)         : %.4f\n", metrics_log$bias_median_log1p))

cat("\n================ TABLE: Out-of-time performance on financing gap (gap-space) ================\n")
cat(sprintf("N used (gap eval)           : %d\n", metrics_gap_a$n))

cat("\n--- Gap A (loan_financing_need_a) ---\n")
cat(sprintf("RMSE (€)        : %.1f\n", metrics_gap_a$rmse))
cat(sprintf("MAE (€)         : %.1f\n", metrics_gap_a$mae))
cat(sprintf("Median AE (€)   : %.1f\n", metrics_gap_a$medae))
cat(sprintf("Spearman        : %.4f\n", metrics_gap_a$spearman))
cat(sprintf("R2              : %.4f\n", metrics_gap_a$r2))
cat(sprintf("Bias mean (€)   : %.1f\n", metrics_gap_a$bias_mean))
cat(sprintf("Bias median (€) : %.1f\n", metrics_gap_a$bias_median))

cat("\n--- Gap B (loan_financing_need_b) ---\n")
cat(sprintf("RMSE (€)        : %.1f\n", metrics_gap_b$rmse))
cat(sprintf("MAE (€)         : %.1f\n", metrics_gap_b$mae))
cat(sprintf("Median AE (€)   : %.1f\n", metrics_gap_b$medae))
cat(sprintf("Spearman        : %.4f\n", metrics_gap_b$spearman))
cat(sprintf("R2              : %.4f\n", metrics_gap_b$r2))
cat(sprintf("Bias mean (€)   : %.1f\n", metrics_gap_b$bias_mean))
cat(sprintf("Bias median (€) : %.1f\n", metrics_gap_b$bias_median))

# ==============================================================================
# BLOCK 9 — Probabilistic imputation on FULL dataset_A
# ==============================================================================
q8a_was_missing <- is.na(dataset_A$q8a)
miss_idx        <- which(q8a_was_missing)
q8a_imputed_chr <- as.character(dataset_A$q8a)

if (length(miss_idx) > 0) {
  X_miss    <- dataset_A[miss_idx, feature_names, drop = FALSE]
  miss_pool <- catboost.load_pool(X_miss, cat_features = cat_features)
  
  probs_miss <- catboost.predict(final_model, miss_pool, prediction_type = "Probability")
  if (is.null(dim(probs_miss))) probs_miss <- matrix(probs_miss, ncol = n_classes, byrow = TRUE)
  # Change C1: the draw must use probabilities on the original class frequencies
  probs_miss <- correct_probs(probs_miss, w_map_all, n_classes)
  
  set.seed(123)
  q8a_imputed_chr[miss_idx] <- apply(probs_miss, 1, function(p) sample(y_levels, 1L, prob = p))
}

dataset_A <- dataset_A %>%
  mutate(
    q8a          = factor(as.character(q8a), levels = y_levels, ordered = TRUE),
    q8a_imputed  = factor(q8a_imputed_chr,   levels = y_levels, ordered = TRUE),
    q8a_imp_flag = ifelse(q8a_was_missing, 1L, 0L)
  )

# ==============================================================================
# BLOCK 10 — Merge imputed midpoints back to the general SAFE dataset 
# ==============================================================================
q8a_mid_tbl <- dataset_A %>%
  transmute(
    permid,
    wave_agg,
    q8a_midpoints = as.numeric(MIDPOINTS_MAP[as.character(q8a_imputed)]),
    q8a_imp_flag  = q8a_imp_flag
  )

dataset_general <- dat_agg %>%
  filter(
    d0 %in% c("BE","DK","DE","EE","IE","GR","ES","FR","HR","IT",
              "CY","LV","LT","LU","MT","NL","AT","PL","PT","RO",
              "SI","SK","FI","SE","BG","CZ","HU"),
    wave_agg > 10
  ) %>%
  arrange(permid, wave_agg) %>%
  left_join(q8a_mid_tbl, by = c("permid", "wave_agg"))

write_dta(dataset_general, path_out_data)
cat("\nSaved merged dataset:", path_out_data, "\n")

# ==============================================================================
# BLOCK 11 — Export Table 10 metrics (CSV + LaTeX)
# ==============================================================================
test_waves_label <- paste(test_waves, collapse = ", ")

metrics_tbl <- tribble(
  ~section, ~metric, ~value,
  "Setup", "Holdout waves",          test_waves_label,
  "Setup", "N holdout observations", as.character(length(test_idx)),
  
  "Bin-space", "Exact accuracy",           sprintf("%.4f", acc_exact),
  "Bin-space", "Within-1 accuracy",        sprintf("%.4f", acc_within1),
  "Bin-space", "Quadratic weighted kappa", sprintf("%.4f", qwk),
  "Bin-space", "Macro F1",                 sprintf("%.4f", macro_f1),
  
  "Setup", "Class-weight mode",      CLASS_WEIGHT_MODE,
  
  "Distribution (TVD, expected shares)", "Overall TVD",         sprintf("%.4f", tvd_E$overall_tvd),
  "Distribution (TVD, expected shares)", "Avg. TVD by wave",    sprintf("%.4f", tvd_E$avg_tvd_wave),
  "Distribution (TVD, expected shares)", "Avg. TVD by country", sprintf("%.4f", tvd_E$avg_tvd_country),
  "Distribution (TVD, expected shares)", "Avg. TVD by sector",  sprintf("%.4f", tvd_E$avg_tvd_sector),
  "Distribution (TVD, expected shares)", "Avg. TVD by size",    sprintf("%.4f", tvd_E$avg_tvd_size),
  
  "Distribution (Mean shift, expected shares)", "Overall mean shift",         sprintf("%.4f", tvd_E$overall_shift),
  "Distribution (Mean shift, expected shares)", "Avg. mean shift by country", sprintf("%.4f", tvd_E$avg_shift_country),
  "Distribution (Mean shift, expected shares)", "Avg. mean shift by sector",  sprintf("%.4f", tvd_E$avg_shift_sector),
  "Distribution (Mean shift, expected shares)", "Avg. mean shift by size",    sprintf("%.4f", tvd_E$avg_shift_size),
  
  "Distribution (TVD, argmax - thesis def.)", "Overall TVD",         sprintf("%.4f", tvd_A$overall_tvd),
  "Distribution (TVD, argmax - thesis def.)", "Avg. TVD by wave",    sprintf("%.4f", tvd_A$avg_tvd_wave),
  "Distribution (TVD, argmax - thesis def.)", "Avg. TVD by country", sprintf("%.4f", tvd_A$avg_tvd_country),
  "Distribution (TVD, argmax - thesis def.)", "Avg. TVD by sector",  sprintf("%.4f", tvd_A$avg_tvd_sector),
  "Distribution (TVD, argmax - thesis def.)", "Avg. TVD by size",    sprintf("%.4f", tvd_A$avg_tvd_size),
  "Distribution (TVD, argmax - thesis def.)", "Overall mean shift",  sprintf("%.4f", tvd_A$overall_shift),
  
  "€-space", "RMSE (€)",      sprintf("%.1f",  metrics_eur$rmse_eur),
  "€-space", "MAE (€)",       sprintf("%.1f",  metrics_eur$mae_eur),
  "€-space", "Median AE (€)", sprintf("%.1f",  metrics_eur$medae_eur),
  "€-space", "Spearman",      sprintf("%.4f", metrics_eur$spearman_eur),
  "€-space", "R2",            sprintf("%.4f", metrics_eur$r2_eur),
  "€-space", "Bias mean",     sprintf("%.1f",  metrics_eur$bias_mean_eur),
  "€-space", "Bias median",   sprintf("%.1f",  metrics_eur$bias_median_eur),
  
  "log1p(€)", "RMSE",        sprintf("%.4f", metrics_log$rmse_log1p),
  "log1p(€)", "MAE",         sprintf("%.4f", metrics_log$mae_log1p),
  "log1p(€)", "Median AE",   sprintf("%.4f", metrics_log$medae_log1p),
  "log1p(€)", "Spearman",    sprintf("%.4f", metrics_log$spearman_log1p),
  "log1p(€)", "R2",          sprintf("%.4f", metrics_log$r2_log1p),
  "log1p(€)", "Bias mean",   sprintf("%.4f", metrics_log$bias_mean_log1p),
  "log1p(€)", "Bias median", sprintf("%.4f", metrics_log$bias_median_log1p),
  
  "Tail-risk", "Big mistake rate (|Δbin|≥2)",        sprintf("%.8f", tail$big_mistake_rate),
  "Tail-risk", "AE p95 (€)",                         sprintf("%.1f",  tail$ae_p95_eur),
  "Tail-risk", "Top-bin underprediction rate",        sprintf("%.8f", tail$topbin_under_rate),
  "Tail-risk", "Catastrophic flip rate (bottom↔top)", sprintf("%.8f", tail$catastrophic_rate),
  
  "€-gap A", "Median AE (€)",   sprintf("%.1f",  metrics_gap_a$medae),
  "€-gap A", "Bias mean (€)",   sprintf("%.1f",  metrics_gap_a$bias_mean),
  "€-gap A", "Bias median (€)", sprintf("%.1f",  metrics_gap_a$bias_median),
  
  "€-gap B", "Median AE (€)",   sprintf("%.1f",  metrics_gap_b$medae),
  "€-gap B", "Bias mean (€)",   sprintf("%.1f",  metrics_gap_b$bias_mean),
  "€-gap B", "Bias median (€)", sprintf("%.1f",  metrics_gap_b$bias_median),
  
  "€-gap (vulnerable==0)", "N used",               as.character(metrics_gap_a$n),
  "€-gap A", "RMSE (€)",      sprintf("%.1f",  metrics_gap_a$rmse),
  "€-gap A", "MAE (€)",       sprintf("%.1f",  metrics_gap_a$mae),
  "€-gap A", "Spearman",      sprintf("%.4f", metrics_gap_a$spearman),
  "€-gap A", "R2",            sprintf("%.4f", metrics_gap_a$r2),
  "€-gap A", "Rel. err mean", sprintf("%.6f", metrics_gap_a$relerr_mean),
  "€-gap B", "RMSE (€)",      sprintf("%.1f",  metrics_gap_b$rmse),
  "€-gap B", "MAE (€)",       sprintf("%.1f",  metrics_gap_b$mae),
  "€-gap B", "Spearman",      sprintf("%.4f", metrics_gap_b$spearman),
  "€-gap B", "R2",            sprintf("%.4f", metrics_gap_b$r2),
  "€-gap B", "Rel. err mean", sprintf("%.6f", metrics_gap_b$relerr_mean)
)

readr::write_csv(metrics_tbl, path_out_metrics_csv)

tex_lines <- c(
  "\\begin{table}[!htbp]\\centering",
  "\\caption{Holdout comparison metrics — Method B, Model I (CatBoost MultiClass)}",
  "\\label{tab:methodB_table10_modelI}",
  "\\begin{tabular}{lll}",
  "\\hline",
  "Section & Metric & Value \\\\",
  "\\hline",
  paste0(metrics_tbl$section, " & ", metrics_tbl$metric, " & ", metrics_tbl$value, " \\\\"),
  "\\hline",
  "\\end{tabular}",
  "\\end{table}"
)
writeLines(tex_lines, path_out_metrics_tex)

cat("\nSaved Table 10 metrics:\n - CSV:", path_out_metrics_csv,
    "\n - TEX:", path_out_metrics_tex, "\n")

t1 <- Sys.time()
cat("\nTOTAL runtime:\n"); print(t1 - t0)

# ==============================================================================
# BLOCK 12 — Heatmap: log(1+gap) error by country × sector/size
# ==============================================================================
heatmap_df <- gap_eval_df %>%
  mutate(
    d0 = as.character(d0),
    across(any_of(c("d1_rec", "d3_rec")), ~ as.character(to_numeric_safe(.x)))
  ) %>%
  mutate(log_err_b = log1p(loan_gap_b_pred) - log1p(loan_gap_b_true))

agg_fun <- function(df, y_var) {
  df %>%
    filter(!is.na(.data[[y_var]]), !is.na(d0), is.finite(log_err_b)) %>%
    group_by(d0, y_cat = .data[[y_var]]) %>%
    summarise(mean_log_err = mean(log_err_b, na.rm = TRUE), n = n(), .groups = "drop") %>%
    mutate(small_n = n < 10)
}

has_d3      <- "d3_rec" %in% names(heatmap_df)
has_d1_heat <- "d1_rec" %in% names(heatmap_df)

if (has_d3)      hs <- agg_fun(heatmap_df, "d3_rec")
if (has_d1_heat) hz <- agg_fun(heatmap_df, "d1_rec")

all_v <- c(if (has_d3) hs$mean_log_err else NULL,
           if (has_d1_heat) hz$mean_log_err else NULL)
sl    <- ceiling(max(abs(all_v), na.rm = TRUE) * 10) / 10
bk    <- quantile(all_v[is.finite(all_v)], c(0.10, 0.90), na.rm = TRUE)

sector_labels <- c("1" = "Industry", "2" = "Construction", "3" = "Trade",  "4" = "Services")
size_labels   <- c("1" = "Micro",    "2" = "Small",        "3" = "Medium", "4" = "Large")

make_panel2 <- function(df, y_label, panel_title, label_map = NULL) {
  if (!is.null(label_map)) df <- df %>% mutate(y_cat = recode(y_cat, !!!label_map))
  
  abs_vals <- abs(df$mean_log_err[is.finite(df$mean_log_err)])
  
  df <- df %>%
    mutate(
      abs_err   = abs(mean_log_err),
      sign_char = if_else(mean_log_err >= 0, "+", "\u2212"),
      txt_col   = if_else(abs(mean_log_err) > 0.45 * sl, "white", "grey10")
    )
  
  bk_bw <- scales::rescale(
    quantile(abs_vals, c(0, 0.33, 0.66, 1), na.rm = TRUE),
    from = range(abs_vals)
  )
  
  ggplot(df, aes(x = d0, y = y_cat, fill = abs_err)) +
    geom_tile(colour = "grey50", linewidth = 0.3) +
    geom_text(
      aes(label = sign_char, colour = txt_col),
      size = 3.2, fontface = "plain", family = "sans"
    ) +
    scale_colour_identity() +
    geom_text(
      data = filter(df, small_n),
      aes(x = d0, y = y_cat, label = "*"),
      inherit.aes = FALSE, colour = "grey50", size = 3.0,
      nudge_x = 0.33, nudge_y = -0.33
    ) +
    scale_fill_gradientn(
      colours = c("grey97", "grey60", "grey25", "grey5"),
      values  = bk_bw,
      limits  = c(0, sl),
      oob     = scales::squish,
      name    = "log(1+gap)\n|error|"
    ) +
    labs(title = panel_title, x = "Country (d0)", y = y_label) +
    theme_minimal(base_size = 9) +
    theme(
      axis.text.x       = element_text(angle = 45, hjust = 1, size = 7),
      axis.text.y       = element_text(size = 8),
      panel.grid        = element_blank(),
      plot.title        = element_text(face = "bold", size = 10),
      legend.key.height = unit(1.8, "cm")
    )
}

panels <- list()
if (has_d3)      panels <- c(panels, list(make_panel2(hs, "Sector (d3)",    "Panel A \u2014 Country \u00d7 Sector",     sector_labels)))
if (has_d1_heat) panels <- c(panels, list(make_panel2(hz, "Firm size (d1)", "Panel B \u2014 Country \u00d7 Size Class", size_labels)))

ggsave(path_out_heatmap, Reduce(`/`, panels), width = 14, height = 10, units = "in", device = "pdf")
cat("\nSaved:", path_out_heatmap, "\n")
