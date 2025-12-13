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
# BLOCK A: Data cleaning and creation of lag and gap variables to add temporal structure
# ===============================================================

# A.1 Load & initial filtering ----------------------------------------------
dataset_A <- read_dta("C:/Users/RDS-petreto/Desktop/Financing Gap/safepanel_allrounds.dta") %>%
  filter(
    d0 %in% c("BE", "DK", "DE", "EE", "IE","GR", "ES", "FR", "HR", "IT",
              "CY", "LV", "LT", "LU", "MT", "NL", "AT", "PL", "PT", "RO",
              "SI", "SK", "FI", "SE", "BG", "CZ", "HU" )
  ) %>%
  select(
    -c(wgtcommon, wgtentr, wgtoldentr, wgtoldcommon, ida,
       experiment_1, experiment_2, intdate)
  )


# A.2

dataset_A <- dataset_A %>% 
  arrange(permid, wave)

## Identify variables to lag (exclude static fields) ---------
static_vars <- c(
  "permid", "wave", "d0")   # ID, time, country

all_vars <- names(dataset_A)
candidate_vars <- setdiff(all_vars, static_vars)

# Detect which vars actually change within firms
change_share <- sapply(candidate_vars, function(v) {
  x <- dataset_A[[v]]
  by_id <- split(x, dataset_A$permid)
  mean(sapply(by_id, function(z) dplyr::n_distinct(na.omit(z)) > 1))
})

# Keep only variables where >5% of firms show change over time
lag_vars <- names(change_share[change_share > 0.05])

## Create lag and gap variables ------------------------------
dataset_A <- dataset_A %>%
  group_by(permid) %>%
  arrange(wave, .by_group = TRUE) %>%
  group_modify(~{
    g <- .
    n <- nrow(g)
    idx <- seq_len(n)
    
    for (v in lag_vars) {
      x <- g[[v]]
      
      # track previous non-NA index
      last_non_na_idx <- integer(n)
      last_seen <- NA_integer_
      
      for (i in seq_len(n)) {
        last_non_na_idx[i] <- last_seen
        if (!is.na(x[i])) last_seen <- i
      }
      
      # lag value
      lag_val <- ifelse(is.na(last_non_na_idx), NA, x[last_non_na_idx])
      
      # gap value = how many waves since last observation
      gap_val <- ifelse(is.na(last_non_na_idx), NA, idx - last_non_na_idx)
      
      g[[paste0(v, "_lag")]] <- lag_val
      g[[paste0(v, "_gap")]] <- gap_val
    }
    
    g
  }) %>%
  ungroup()

dataset_A <- dataset_A %>%
  filter(wave > 10, !wave %in% c(31, 33, 35))

write_dta(dataset_A, "Dataset_test_V12.dta")

# A.3 Define target variable Q8A --------------------------------------------

dataset_A <- dataset_A %>%
  mutate(Q8A_target = as.factor(haven::as_factor(q8a)))

# A.4 Build predictors: gaps numeric, everything else factor -----------------

predictors_A <- dataset_A %>% 
  select(-q8a, -Q8A_target)

gap_cols <- grep("_gap$", names(predictors_A), value = TRUE)

predictors_A <- predictors_A %>%
  mutate(
    # all non-gap vars as factors
    across(!all_of(gap_cols), ~ as.factor(as.character(.)))
    # gap vars stay as they are (numeric)
  )

# A.5 Filters: high missingness, low variance -------

non_missing_threshold <- 0.60               
non_na_pct <- colMeans(!is.na(predictors_A))

predictors_A <- predictors_A[, non_na_pct >= non_missing_threshold, drop = FALSE]
cat("Variables kept after filtering for missingness:", ncol(predictors_A), "\n")

low_var_threshold <- 0.97   # max proportion allowed for the most common category

is_low_var <- sapply(predictors_A, function(x) {
  if (!is.factor(x)) return(FALSE)   # <- numeric vars (gaps) are never dropped here
  x_non_na <- x[!is.na(x)]
  if (!length(x_non_na)) return(TRUE)
  tab <- table(x_non_na)
  (max(tab) / sum(tab)) > low_var_threshold
})

predictors_A <- predictors_A[, !is_low_var, drop = FALSE]

cat("Variables kept after low variance filtering:", ncol(predictors_A), "\n")

# Enforce lag/gap pairing (drop unmatched members) 

lag_gap_cols <- grep("_(lag|gap)$", names(predictors_A), value = TRUE)
bases        <- sub("_(lag|gap)$", "", lag_gap_cols)
tab_bases    <- table(bases)

# bases that appear only once among lag/gap cols → unmatched
unmatched_bases <- names(tab_bases)[tab_bases == 1]

to_drop <- unlist(
  lapply(unmatched_bases, function(b) {
    grep(paste0("^", b, "_(lag|gap)$"), names(predictors_A), value = TRUE)
  })
)

if (length(to_drop) > 0) {
  message("Dropping unmatched lag/gap variables: ", paste(to_drop, collapse = ", "))
  predictors_A <- predictors_A[, !names(predictors_A) %in% to_drop, drop = FALSE]
}

cat("Variables kept after pairs filtering:", ncol(predictors_A), "\n")

# ===============================================================
# BLOCK B: Predictors imputation
# ===============================================================
       
# B.1 Setting the weight based on missingness -------------------

## Case weights = how informative each respondent is
cw <- rowSums(!is.na(predictors_A))  # number of answered questions per obs
cw <- cw / mean(cw)                  # rescaling

# B.2 Impute predictors with missRanger ------------------------
       
set.seed(123)

X_A <- missRanger(
  predictors_A,
  num.trees    = 10,
  num.threads  = parallel::detectCores(),
  pmm.k        = 3,
  case.weights = cw,       # <- weight each *row* based on the missingness
  returnOOB    = TRUE,
  verbose      = 1
  )
       
## OOB misclassification error per variable
oob_err <- attr(X_A, "oob")
head(oob_err)
sort(oob_err, decreasing = FALSE)[1:20]  

# ===============================================================
# BLOCK C: Target Variable Prediction (q8a)
# ===============================================================


       
