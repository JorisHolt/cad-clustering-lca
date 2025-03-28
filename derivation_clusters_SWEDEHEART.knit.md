---
title: "Derivation of Clusters in SWEDEHEART"
output: 
  github_document:
    toc: true
    toc_depth: 2
---

# Overview

This markdown file contains the full code used to derive patient clusters in the SWEDEHEART registry using Latent Class Analysis (LCA). The aim was to identify clinically meaningful phenotypes in patients with established coronary artery disease (CAD), based on routinely collected clinical variables.

> **Note**: This script uses a single imputation dataset and follows the complete LCA modelling procedure, from variable selection to model fitting and posterior probability extraction. The goal is transparency; this script is not meant to be run as-is without prior preparation.

---

# 📦 1. Package Setup

LCA and associated analysis steps rely on multiple packages. The following block ensures all required packages are installed and loaded.

```r
install_and_load_packages <- function(packages) {
  missing_packages <- packages[!sapply(packages, function(x) requireNamespace(x, quietly = TRUE))]
  if (length(missing_packages) > 0) {
    install.packages(missing_packages, dependencies = TRUE)
  }
  lapply(packages, library, character.only = TRUE)
}

my_packages <- c("broom", "data.table", "dplyr", "ggplot2", "cmprsk", "riskRegression",
                 "ggpubr", "magrittr", "purrr", "readxl", "stringr", "ggalluvial", 
                 "stringi", "tidyr", "tidyselect", "tidyverse", "lubridate", "cluster",
                 "haven", "pamr", "moments", "visdat", "mice", "foreign", "Hmisc", "poLCA",
                 "tableone", "survminer", "rms")

install_and_load_packages(my_packages)
```

```r
# Set reproducibility seed
set.seed(1337)
```

---

# 📂 2. Data Preparation

LCA requires complete categorical data. All variables were recoded into **factors**, using clinically meaningful cut-offs (e.g., KDIGO eGFR stages). Variables with values `== 2` generally indicate *presence of condition* (e.g., diabetes == 2 means diabetes is present).

> **Missing data were handled using single imputation**, which is assumed to be done prior to running this script.

```r
data <- data %>% 
  mutate(
    # Sociodemographics
    sex_nr = factor(ifelse(sex == 1, 2, 1)),
    LVEF_class_nr = factor(ifelse(LVEF_class == "HFrEF <40%", 1,
                           ifelse(LVEF_class == "HFmrEF 40-50%", 2, 3)), ordered = TRUE),
    multivessel_nr = factor(ifelse(multivessel == "Ja", 2, 1)),
    age_nr = factor(case_when(age <= 55 ~ 1, age <= 70 ~ 2, TRUE ~ 3), ordered = TRUE),
    isSmoking_nr = factor(ifelse(isSmoking == 1, 2, 1)),

    # Measurements
    bmi_nr = factor(case_when(bmi <= 25 ~ 1, bmi <= 30 ~ 2, TRUE ~ 3), ordered = TRUE),
    dbp_nr = factor(case_when(dbp <= 70 ~ 1, dbp <= 80 ~ 2, TRUE ~ 3), ordered = TRUE),
    nonHDL_nr = factor(ifelse(nhdl > 2.59, 2, 1)),
    KDIGO_egfr_nr = factor(case_when(gfr > 90 ~ 1, gfr > 60 ~ 2, gfr > 45 ~ 3, TRUE ~ 4)),
    crp_nr = factor(ifelse(crp > 3, 2, 1)),

    # History
    history_afib_nr = factor(ifelse(history_afib == 1, 2, 1)),
    diabetesDiagnosis_nr = factor(ifelse(diabetesDiagnosis == 1, 2, 1)),
    poly_vasc_nr = factor(ifelse(cerebrovascularDisease == 1 | aorticAneurysm == 1 | peripheralArteryDisease == 1, 2, 1))
  )
```

---

# 📊 3. Variable Correlation Inspection

Before clustering, we visualise the **Spearman correlation** between candidate variables. Since most variables are ordinal/binary, this rank-based method is more appropriate than Pearson.

## 3.1 Select Variables

```r
correlation_data.full.SWinSW <- data %>%
  dplyr::select(contains("_nr")) %>%
  mutate_all(as.numeric) %>%
  dplyr::select(sex_nr, age_nr, isSmoking_nr, bmi_nr, dbp_nr, nonHDL_nr, KDIGO_egfr_nr,
                multivessel_nr, LVEF_class_nr, crp_nr, poly_vasc_nr, diabetesDiagnosis_nr, history_afib_nr)
```

## 3.2 Compute Spearman Correlation Matrix

```r
cor_matrices.SWinSW <- cor(as.matrix(correlation_data.full.SWinSW), method = "spearman")
```

## 3.3 Correlation Plot

Pretty label mapping for readability:

```r
rename_vector <- c(
  "poly_vasc_nr" = "Polyvascular disease",
  "age_nr" = "Age",
  "isSmoking_nr" = "Current smoking",
  "nonHDL_nr" = "non-HDL cholesterol",
  "bmi_nr" = "BMI",
  "KDIGO_egfr_nr" = "eGFR",
  "history_afib_nr" = "Atrial fibrillation",
  "crp_nr" = "C-reactive protein",
  "sex_nr" = "Sex",
  "diabetesDiagnosis_nr" = "History of diabetes",
  "dbp_nr" = "Diastolic blood pressure",
  "multivessel_nr" = "Multivessel disease",
  "LVEF_class_nr" = "LVEF class")
```

Create the full correlation plot:

```r
correlation_plot.SWinSW <- ggcorrplot::ggcorrplot(cor_matrices.SWinSW, hc.order = FALSE, type = "full",
  outline.color = "black", lab = TRUE, lab_size = 3.5, method = "square",
  colors = c("red", "white", "blue")) +
  scale_x_discrete(labels = rename_vector) +
  scale_y_discrete(labels = rename_vector) +
  labs(title = "Correlation of variables in dataset") +
  theme_minimal()
```

---

# 🔍 4. Latent Class Model Setup

Define the model formulas, one with and one without LVEF/multivessel:

```r
y_lvef <- cbind(sex_nr, age_nr, education_levels_nr, isSmoking_nr, bmi_nr, dbp_nr, nonHDL_nr, KDIGO_egfr_nr,
                crp_nr, poly_vasc_nr, multivessel_nr, diabetesDiagnosis_nr, history_afib_nr, LVEF_class_nr) ~ 1

y_nolvef <- cbind(sex_nr, age_nr, education_levels_nr, isSmoking_nr, bmi_nr, dbp_nr, nonHDL_nr, KDIGO_egfr_nr,
                  crp_nr, poly_vasc_nr, diabetesDiagnosis_nr, history_afib_nr) ~ 1
```

---

# ⚙️ 5. Model Fitting: LCA across 2–10 Classes

Prepare lists to store model objects and diagnostics:

```r
list_listfull.lvef.SW <- list()
list_listfull.nolvef.SW <- list()
BIC_start.nolvef.SW <- list()

ncla_range <- 2:10

for (ncla in ncla_range) {
  model <- poLCA(y_nolvef, xx2.SWinSW, nclass = ncla, maxiter = 10000, nrep = 10, verbose = FALSE)
  list_listfull.nolvef.SW[[paste0("ncla_", ncla)]] <- model
  BIC_start.nolvef.SW[[paste0("ncla_", ncla)]] <- model$bic
}
```

---

# 📈 6. Posterior Probability Extraction

Once models are fitted, we extract posterior class membership probabilities.

```r
list_posterior.nolvef.SWinSW <- list()

for (l in ncla_range) {
  post <- poLCA.posterior(lc = list_listfull.nolvef.SW[[paste0("ncla_", l)]], y = xx2.SWinSW)
  class_max <- as.factor(apply(post, 1, which.max))
  post <- as.data.frame(post)
  colnames(post) <- paste0("C", 1:l)

  summarised <- cbind(class_max, post) %>%
    group_by(class_max) %>%
    summarise_all(~paste0(round(quantile(.x, probs = 0.5), 4),
                          " [", round(quantile(.x, probs = 0.25), 4), "-",
                          round(quantile(.x, probs = 0.75), 4), "]")) %>%
    arrange(as.numeric(as.character(class_max)))

  list_posterior.nolvef.SWinSW[[paste0("ncla_", l)]] <- summarised
}
```
