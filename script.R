# This script accompanies the mini-project "A machine learning approach to map and predict 
# inequalities in COVID-19 vaccine uptake in Scotland". Please find the final report on 
# the project's GitHub repository: 
# https://github.com/sohailferdous/covid-vaccine-uptake_ml_project
# 
# Please ensure the necessary data files are available in your working folder. Data is
# available at Public Health Scotland's Scottish Health and Social Care Open Data 
# repository: https://www.opendata.nhs.scot/dataset/flu-covid-vaccinations
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1. DATA AND PACKAGE IMPORT ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
#* 1.1 Loading libraries =================================================================
library(tidyverse)
library(tidymodels)
library(ranger)
library(janitor)
library(lubridate)

#* 1.2 Importing data ====================================================================
data_simd <- read.csv("covid_simd_hb_20231806.csv")
data_ethnicity <- read.csv("covid_ethnicity_hb_20231806.csv")

#* 1.3 Checking imported data ============================================================
view(data_simd)
glimpse(data_simd)

view(data_ethnicity)
glimpse(data_ethnicity)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2. DATA PREPARATION ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#* 2.1 Cleaning column names =============================================================
data2_simd <- clean_names(data_simd)
data2_ethnicity <- clean_names(data_ethnicity)

#* 2.2 Date is reported in both yyyymmdd and yyyyddmm formats - standardising to yyyymmdd format and changing variable type to date ====
data2_simd <- data2_simd %>% 
  mutate(date_new = parse_date_time(date, orders = c("Ymd", "Ydm")),
         date_new = as.Date(date_new))

data2_ethnicity <- data2_ethnicity %>% 
  mutate(date_new = parse_date_time(date, orders = c("Ymd", "Ydm")),
         date_new = as.Date(date_new))

#* 2.3 Removing rows for incompatible and unknown health board, ethnicity, and simd variables ====
clean_simd <- data2_simd %>% 
  filter(hb_name != "Scotland",
         hb_name != "Unknown",
         simd_decile != "Not Known")

clean_ethnicity <- data2_ethnicity %>% 
  filter(hb_name != "Scotland",
         hb_name != "Unknown",
         ethnicity != "Not Known")

#* 2.4 Factorising variables and calculating vaccination counts ==========================
clean_simd <- clean_simd %>%
  mutate(hbname = factor(hb_name),
         simd_decile = factor(simd_decile, levels = c("10","9","8","7","6","5","4","3",
                                                      "2","1"),ordered = TRUE),
         n_vaccinated   = count,
         n_unvaccinated = population - count,
         prop_vaccinated = n_vaccinated / population)

clean_ethnicity <- clean_ethnicity %>%
  mutate(hbname = factor(hb_name),
         ethnicity = 
           factor(ethnicity, 
                  levels = c("African, African Scottish or African British",
                             "Other African",
                             "Bangladeshi, Bangladeshi Scottish or Bangladeshi British",
                             "Chinese, Chinese Scottish or Chinese British",
                             "Indian, Indian Scottish or Indian British",
                             "Other Asian, Asian Scottish or Asian British",
                             "Pakistani, Pakistani Scottish or Pakistani British",
                             "Black, Black Scottish or Black British",
                             "Caribbean, Caribbean Scottish or Caribbean British",
                             "Other Caribbean or Black",
                             "Any mixed or multiple ethnic groups",
                             "Arab, Arab Scottish or Arab British",
                             "Other ethnic group",
                             "Gypsy/Traveller",
                             "Irish",
                             "Other British",
                             "Other white ethnic group",
                             "Polish",
                             "Scottish")),
         n_vaccinated   = count,
         n_unvaccinated = population - count,
         prop_vaccinated = n_vaccinated / population)

#* 2.5 Creating subsets only with required variables and renaming the date_new variable ==
subset_simd <- subset(clean_simd, select = c(date_new, hb_name, simd_decile, population, 
                                             n_vaccinated, n_unvaccinated, 
                                             prop_vaccinated))
subset_simd <- subset_simd %>% 
  rename(date = date_new)

subset_ethnicity <- subset(clean_ethnicity, select = c(date_new, hb_name, ethnicity, 
                                                       population, n_vaccinated, 
                                                       n_unvaccinated, prop_vaccinated))

subset_ethnicity <- subset_ethnicity %>% 
  rename(date = date_new)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3. PATTERNS OF VACCINE UPTAKE (AS OF JUNE 18, 2023) ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
#* 3.1 Plot - Vaccine uptake by SIMD deciles/ethnicity ===================================
#
#** 3.1.1 Preparing data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subset_simd_plot <- subset_simd %>%
  group_by(simd_decile) %>% 
  summarise(total_vaccinated   = sum(n_vaccinated),
            total_population   = sum(population),
            prop_vaccinated    = total_vaccinated / total_population,
            .groups = "drop")

subset_ethnicity_plot <- subset_ethnicity %>%
  group_by(ethnicity) %>% 
  summarise(total_vaccinated   = sum(n_vaccinated),
            total_population   = sum(population),
            prop_vaccinated    = total_vaccinated / total_population,
            .groups = "drop")

#** 3.1.2 Plot code ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
simd_plot <- ggplot(subset_simd_plot, aes(x = simd_decile, y = prop_vaccinated)) +
  geom_col() +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(title = "COVID-19 vaccine uptake by SIMD decile (as of June 18, 2023)",
       x = "SIMD decile (10 = least deprived)",
       y = "Uptake") +
  theme_minimal()

print(simd_plot)

ethnicity_plot <- ggplot(subset_ethnicity_plot, 
                         aes(x = reorder(ethnicity, -prop_vaccinated), 
                             y = prop_vaccinated)) +
  geom_col() +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(title = "COVID-19 vaccine uptake by ethnicity (as of June 18, 2023)",
       x = "Ethnicity",
       y = "Uptake") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(ethnicity_plot)

#* 3.2 ANOVA testing for mean vaccine uptake and simd deciles/ethnicity ==================
simd_aov <- aov(prop_vaccinated ~ simd_decile, data = subset_simd)
print(simd_aov)

ethnicity_aov <- aov(prop_vaccinated ~ ethnicity, data = subset_ethnicity)
print(ethnicity_aov)

# Checking for assumption of normality by plotting histogram of ANOVA test residuals
simd_resid <- residuals(simd_aov)
hist(simd_resid)

ethnicity_resid <- residuals(ethnicity_aov)
hist(ethnicity_resid)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 4. MACHINE LEARNING ANALYSIS ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
#* 4.1 Calculating low coverage threshold ================================================
threshold_simd <- subset_simd %>%
  summarise(q25 = quantile(prop_vaccinated, 0.25, na.rm = TRUE)) %>%
  pull(q25)

print(threshold_simd)  # 71.07%

threshold_ethnicity <- subset_ethnicity %>%
  summarise(q25 = quantile(prop_vaccinated, 0.25, na.rm = TRUE)) %>%
  pull(q25)

print(threshold_ethnicity)  # 44.18%

#* 4.2 Labelling low uptake rows =========================================================
data_ml_simd <- subset_simd %>%
  mutate(low_uptake = prop_vaccinated < threshold_simd,
         low_uptake = factor(low_uptake, levels = c(TRUE, FALSE), labels = c("yes","no")))

data_ml_ethnicity <- subset_ethnicity %>%
  mutate(low_uptake = prop_vaccinated < threshold_ethnicity,
         low_uptake = factor(low_uptake, levels = c(TRUE, FALSE), labels = c("yes","no")))

#* 4.3 Splitting datasets based on date ==================================================
set.seed(123)              # Specifying the random seed to make the analysis reproducible

simd_train <- data_ml_simd %>%
  filter(date == "2023-01-29")

simd_test <- data_ml_simd %>%
  filter(date == "2023-06-18")

ethnicity_train <- data_ml_ethnicity %>%
  filter(date == "2023-01-29")

ethnicity_test <- data_ml_ethnicity %>%
  filter(date == "2023-06-18")

#* 4.4 Creating a recipe which clearly defines the explanatory and outcome variables =====
simd_recipe <- recipe(low_uptake ~ simd_decile + hb_name + population, data = simd_train)

ethnicity_recipe <- recipe(low_uptake ~ ethnicity + hb_name + population,
                           data = ethnicity_train)

#* 4.5 Specifying algorithm specifications ===============================================
algorithm_spec <- rand_forest(trees = 500) %>%
  set_engine("ranger") %>%
  set_mode("classification")

#* 4.6 Creating workflows (could consider these untrained algorithms) by combining recipes and algorithm specifications ====
simd_ml_workflow <- workflow() %>%
  add_model(algorithm_spec) %>%
  add_recipe(simd_recipe)

ethnicity_ml_workflow <- workflow() %>%
  add_model(algorithm_spec) %>%
  add_recipe(ethnicity_recipe)

#* 4.7 Training the algorithms ===========================================================
simd_fit <- simd_ml_workflow %>%
  fit(data = simd_train)

ethnicity_fit <- ethnicity_ml_workflow %>%
  fit(data = ethnicity_train)

# 4.8 Testing trained algorithms =========================================================
# 
#* 4.8.1 Calculating probabilities for low_uptake ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
simd_prob <- predict(simd_fit, simd_test, type = "prob")
ethnicity_prob <- predict(ethnicity_fit, ethnicity_test, type = "prob")

#* 4.8.2 Calculating binary yes/no predictions for low_uptake ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
simd_cate <- predict(simd_fit, simd_test, type = "class")
ethnicity_cate <- predict(ethnicity_fit, ethnicity_test, type = "class")

# 4.9 Combining the predictions with test data to view the output (this might appear a bit complicated as I am trying to rearrange and bind columns in one step) ====
simd_prediction <- simd_test %>% select(hb_name, simd_decile, low_uptake) %>% 
  bind_cols(simd_cate) %>%
  bind_cols(simd_prob) %>%
  bind_cols(simd_test %>% select(date))

view(simd_prediction)

ethnicity_prediction <- ethnicity_test %>% select(hb_name, ethnicity, low_uptake) %>% 
  bind_cols(ethnicity_cate) %>%
  bind_cols(ethnicity_prob) %>%
  bind_cols(ethnicity_test %>% select(date))

view(ethnicity_prediction)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 5. ALGORITHM PERFORMANCE TESTS ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
#* 5.1 Confusion matrices ================================================================
simd_confusion <-simd_prediction %>%
  conf_mat(truth = low_uptake, estimate = .pred_class)

print(simd_confusion)

ethnicity_confusion <- ethnicity_prediction %>%
  conf_mat(truth = low_uptake, estimate = .pred_class)

print(ethnicity_confusion)

#* 5.2 Accuracy and kappa values =========================================================
simd_perf <- simd_prediction %>%
  metrics(truth = low_uptake, estimate = .pred_class)

print(simd_perf)

ethnicity_perf <- ethnicity_prediction %>%
  metrics(truth = low_uptake, estimate = .pred_class)

print(ethnicity_perf)

#* 5.3 Area under curve calculation and ROC curves =======================================
simd_auc <- roc_auc(simd_prediction,truth = low_uptake,.pred_yes)
simd_roc <- simd_prediction %>%
  roc_curve(truth = low_uptake, .pred_yes) %>%
  autoplot() +
  labs(title = "ROC curve: Random forest prediction of low uptake for SIMD decile",
       subtitle = paste0("AUC = ", round(simd_auc$.estimate, 3)),
       x = "1 - specificity",
       y = "Sensitivity") +
  theme_minimal() +
  theme(plot.title = element_text(size = 10),
        plot.subtitle = element_text(size = 10))

print(simd_roc)

ethnicity_auc <- roc_auc(ethnicity_prediction,truth = low_uptake,.pred_yes)
ethnicity_roc <- ethnicity_prediction %>%
  roc_curve(truth = low_uptake, .pred_yes) %>%
  autoplot() +
  labs(title = "ROC curve: Random forest prediction of low uptake for ethnicity",
       subtitle = paste0("AUC = ", round(ethnicity_auc$.estimate, 3)),
       x = "1 - specificity",
       y = "Sensitivity") +
  theme_minimal() +
  theme(plot.title = element_text(size = 10),
        plot.subtitle = element_text(size = 10))

print(ethnicity_roc) 

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ END OF SCRIPT ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
