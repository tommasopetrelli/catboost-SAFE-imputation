# ===============================================================
# Load required packages
# ===============================================================

if (!requireNamespace("missRanger")) install.packages("missRanger")
if (!requireNamespace("ranger")) install.packages("ranger")
if (!requireNamespace("dplyr")) install.packages("dplyr")
if (!requireNamespace("haven")) install.packages("haven")
if (!requireNamespace("forcats")) install.packages("forcats")

library(missRanger)
library(ranger)
library(dplyr)
library(haven)
library(forcats)


# ===============================================================
# BLOCK A: Data cleaning and creation of lag/gap variables
# ===============================================================

cat("Starting BLOCK A:", as.character(Sys.time()), "\n")

# ---------------------------------------------------------------
# A.1 LOAD & FILTER
# ---------------------------------------------------------------
cat("A.1 Loading dataset and filtering:", as.character(Sys.time()), "\n")

dataset_A <- read_dta("C:/Users/RDS-petreto/Desktop/Financing Gap/safepanel_allrounds.dta") %>%
  filter(
    d0 %in% c("BE","DK","DE","EE","IE","GR","ES","FR","HR","IT",
              "CY","LV","LT","LU","MT","NL","AT","PL","PT","RO",
              "SI","SK","FI","SE","BG","CZ","HU")
  ) %>%
  select(
    -c(wgtcommon, wgtentr, wgtoldentr, wgtoldcommon, ida,
       experiment_1, experiment_2, intdate)
  ) %>%
  arrange(permid, wave)


# ---------------------------------------------------------------
# A.2 IDENTIFY VARIABLES FOR LAGS
# ---------------------------------------------------------------
cat("A.2 Identifying lag-eligible variables:", as.character(Sys.time()), "\n")

static_vars <- c(
  "permid", "wave", "d0",
  "d1_rec","d2","d3_rec","d4",
  "d5_rec","d6","d6_rec","d6b","d7","d7_rec"
)

all_vars       <- names(dataset_A)
candidate_vars <- setdiff(all_vars, static_vars)

change_share <- sapply(candidate_vars, function(v) {
  x     <- dataset_A[[v]]
  by_id <- split(x, dataset_A$permid)
  mean(sapply(by_id, function(z) dplyr::n_distinct(na.omit(z)) > 1))
})

lag_vars <- names(change_share[change_share > 0.05])


# ---------------------------------------------------------------
# A.3 CREATE LAG & GAP VARIABLES
# ---------------------------------------------------------------
cat("A.3 Creating lag & gap variables:", as.character(Sys.time()), "\n")

last_non_na_index <- function(x) {
  idx       = seq_along(x)
  last_incl = cummax(replace(idx, is.na(x), 0L))
  prev_idx  = dplyr::lag(last_incl)
  prev_idx[prev_idx == 0L] <- NA_integer_
  prev_idx
}

dataset_A <- dataset_A %>%
  group_by(permid) %>%
  arrange(wave, .by_group = TRUE) %>%
  mutate(
    across(
      all_of(lag_vars),
      ~ { prev_idx <- last_non_na_index(.); .[prev_idx] },
      .names = "{.col}_lag"
    ),
    across(
      all_of(lag_vars),
      ~ {
        prev_idx <- last_non_na_index(.)
        idx      <- seq_along(.)
        as.numeric(ifelse(is.na(prev_idx), NA, idx - prev_idx))
      },
      .names = "{.col}_gap"
    )
  ) %>%
  ungroup()


# ---------------------------------------------------------------
# A.4 RESTRICT WAVES
# ---------------------------------------------------------------
cat("A.4 Restricting waves:", as.character(Sys.time()), "\n")

dataset_A <- dataset_A %>%
  filter(wave > 10, !wave %in% c(31,33,35))


# ---------------------------------------------------------------
# A.5 TYPE ASSIGNMENTS & BUILD PREDICTORS
# ---------------------------------------------------------------
cat("A.5 Assigning types and constructing predictor matrix:", as.character(Sys.time()), "\n")

dataset_A <- dataset_A %>%
  mutate(
    permid     = as.numeric(permid),
    wave       = as.numeric(wave),
    d0         = as.factor(d0),
    Q8A_target = ordered(haven::as_factor(q8a))
  )

predictors_A <- dataset_A %>% 
  select(-q8a, -Q8A_target)

gap_cols     <- grep("_gap$", names(predictors_A), value = TRUE)
id_time_vars <- c("permid", "wave")

# Convert *non-gap*, non-ID variables to factors
predictors_A <- predictors_A %>%
  mutate(
    across(
      .cols = -c(all_of(gap_cols), all_of(id_time_vars), d0),
      .fns  = ~ as.factor(as.character(.))
    )
  )


# ---------------------------------------------------------------
# A.6 FILTERING: missingness, low variance, lag-gap triplets
# ---------------------------------------------------------------
cat("A.6 Filtering predictors:", as.character(Sys.time()), "\n")

non_missing_threshold <- 0.60
non_na_pct <- colMeans(!is.na(predictors_A))

predictors_A <- predictors_A[, non_na_pct >= non_missing_threshold, drop = FALSE]
cat("Variables kept after missingness filtering:", ncol(predictors_A), "\n")

low_var_threshold <- 0.97

is_low_var <- sapply(predictors_A, function(x) {
  if (!is.factor(x)) return(FALSE)
  x_non_na <- x[!is.na(x)]
  if (!length(x_non_na)) return(TRUE)
  tab <- table(x_non_na)
  (max(tab) / sum(tab)) > low_var_threshold
})

predictors_A <- predictors_A[, !is_low_var, drop = FALSE]
cat("Variables kept after low-variance filtering:", ncol(predictors_A), "\n")

# Identify *_lag / *_gap / original triplets
lag_vars_in_data <- grep("_lag$", names(predictors_A), value = TRUE)
gap_vars_in_data <- grep("_gap$", names(predictors_A), value = TRUE)

lag_bases <- sub("_lag$", "", lag_vars_in_data)
gap_bases <- sub("_gap$", "", gap_vars_in_data)

bases_with_lag_and_gap <- intersect(lag_bases, gap_bases)
original_vars <- bases_with_lag_and_gap[bases_with_lag_and_gap %in% names(predictors_A)]

complete_bases <- original_vars

vars_to_keep <- c(
  unlist(lapply(complete_bases, function(b) c(b, paste0(b,"_lag"), paste0(b,"_gap")))),
  c("permid", "wave", "d0")
)

vars_to_keep <- intersect(vars_to_keep, names(predictors_A))

predictors_A <- predictors_A[, vars_to_keep, drop = FALSE]

cat("Variables kept after triplet filtering:", ncol(predictors_A), "\n")


# ---------------------------------------------------------------
# SAVE
# ---------------------------------------------------------------
write_dta(dataset_A,  "Dataset_BlockA_V1.dta")
write_dta(predictors_A, "Dataset_BlockA_Predictors_V1.dta")

cat("Block A saved successfully:", as.character(Sys.time()), "\n")


# ===============================================================
# BLOCK B: Imputation of predictors with missRanger
# ===============================================================

cat("\nStarting BLOCK B:", as.character(Sys.time()), "\n")

# B.1 Row weights -----------------------------------------------------------
cat("B.1 Computing row weights:", as.character(Sys.time()), "\n")

cw <- rowSums(!is.na(predictors_A))
cw <- cw / mean(cw)

# B.2 Run missRanger --------------------------------------------------------
cat("B.2 Running missRanger (this may take time):", as.character(Sys.time()), "\n")

set.seed(123)
X_A <- missRanger(
  predictors_A,
  num.trees    = 10,
  num.threads  = parallel::detectCores(),
  pmm.k        = 3,
  case.weights = cw,
  returnOOB    = TRUE,
  verbose      = 0
)

oob_err <- attr(X_A, "oob")

# B.3 Build full dataset ----------------------------------------------------
cat("B.3 Reconstructing dataset with imputed predictors:", as.character(Sys.time()), "\n")

dataset_B <- dataset_A %>%
  select(permid, wave, d0, Q8A_target) %>%
  bind_cols(as.data.frame(X_A))

write_dta(dataset_B, "Dataset_BlockB_V1.dta")
cat("Block B saved: Dataset_BlockB_V1.dta –", as.character(Sys.time()), "\n")

# ===============================================================
# BLOCK C: Random Forest Imputation of missing Q8A_target
#           (with Grid Search)
# ===============================================================

cat("\nStarting BLOCK C:", as.character(Sys.time()), "\n")

target_var <- "Q8A_target"
has_target <- !is.na(dataset_B[[target_var]])
no_target  <- !has_target

train_df <- dataset_B[has_target, , drop = FALSE]

predictor_vars <- setdiff(
  names(dataset_B),
  c("Q8A_target", "q8a")         # ensure no raw q8a leakage
)

# C.1 Build RF formula -----------------------------------------

cat("C.1 Building RF formula:", as.character(Sys.time()), "\n")

rf_formula <- as.formula(
  paste(target_var, "~", paste(predictor_vars, collapse = " + "))
)


# C.2 GRID SEARCH for best hyperparameters ----------------------

cat("C.2 Starting GRID SEARCH:", as.character(Sys.time()), "\n")

p <- length(predictor_vars)

grid_C <- expand.grid(
  mtry                      = c(floor(sqrt(p)), floor(p/3)),
  min.node.size             = c(5, 20, 50),
  sample.fraction           = c(0.6, 0.8, 1.0),
  respect.unordered.factors = c("ignore", "order"),
  stringsAsFactors          = FALSE
)

oob_errors <- numeric(nrow(grid_C))

for (i in seq_len(nrow(grid_C))) {
  
  cat("  Testing config", i, "of", nrow(grid_C), ":", as.character(Sys.time()), "\n")
  
  cfg <- grid_C[i, ]
  
  model_tmp <- ranger(
    formula                  = rf_formula,
    data                     = train_df,
    num.trees                = 200,
    mtry                     = cfg$mtry,
    min.node.size            = cfg$min.node.size,
    sample.fraction          = cfg$sample.fraction,
    classification           = TRUE,
    probability              = FALSE,
    importance               = "impurity_corrected",
    respect.unordered.factors = cfg$respect.unordered.factors,
    num.threads              = parallel::detectCores()
  )
  
  oob_errors[i] <- model_tmp$prediction.error
}

grid_C$oob <- oob_errors

cat("Grid search completed:", as.character(Sys.time()), "\n")
cat("Best configuration:\n")
print(grid_C[which.min(grid_C$oob), ])

# C.3 Fit FINAL RF MODEL with best hyperparameters  -------------------

cat("C.3 Fitting FINAL Random Forest:", as.character(Sys.time()), "\n")

best_cfg <- grid_C[which.min(grid_C$oob), ]

rf_model_final <- ranger(
  formula                  = rf_formula,
  data                     = train_df,
  num.trees                = 500,
  mtry                     = best_cfg$mtry,
  min.node.size            = best_cfg$min.node.size,
  sample.fraction          = best_cfg$sample.fraction,
  classification           = TRUE,
  probability              = FALSE,
  importance               = "impurity_corrected",
  respect.unordered.factors = best_cfg$respect.unordered.factors,
  num.threads              = parallel::detectCores()
)

cat("Final RF OOB error:", rf_model_final$prediction.error, "\n")

# C.4 Predict missing target values ------------------------------------

cat("C.4 Predicting missing Q8A_target:", as.character(Sys.time()), "\n")

if (any(no_target)) {
  newdata <- dataset_B[no_target, predictor_vars, drop = FALSE]
  pred_RF <- predict(rf_model_final, data = newdata)
  
  Q8A_RF <- dataset_B[[target_var]]
  Q8A_RF[no_target] <- pred_RF$predictions
} else {
  Q8A_RF <- dataset_B[[target_var]]
}

dataset_C <- dataset_B
dataset_C$Q8A_target_RF <- Q8A_RF


# C.5 Save final dataset ------------------------------------------------

write_dta(dataset_C, "Imputed_dataset_Method_A.dta")
cat("Final dataset saved as Imputed_dataset_Method_A.dta –", as.character(Sys.time()), "\n")

cat("\nBLOCK C completed.", as.character(Sys.time()), "\n")

