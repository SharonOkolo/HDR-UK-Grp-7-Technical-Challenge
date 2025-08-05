#############################################
# Project: Breast Cancer (C50) Data Exploration
# Dataset: Simulacrum v2.1.0
# Purpose:
# - Filter for breast cancer (ICD10 C50) tumours only
# - Join tumour data to patient demographics
# - Prepare clean cohort for further linking to SACT, RTDS, and Genomic data
# Author: Elsie Osifeso
# Date: 08/07/2025
#############################################

## Set up -----------------

# Clear environment
rm(list = ls())

# Set working directory
setwd("/Users/elsieosifeso/Desktop/Internship 2025/HDRUK/Technical Challenge/simulacrum_v2.1.0/Data")

# --------------------------------------------
# Install required packages
# --------------------------------------------

install.packages("tidyverse")
install.packages("readr")
install.packages("data.table")
install.packages("officer")
install.packages("flextable")
install.packages("dplyr")
install.packages("gtsummary")
install.packages("stringr")
install.packages("ggpubr")

# --------------------------------------------
# Load libraries
library(tidyverse)
library(readr)
library(officer)
library(flextable)
library(dplyr)
library(gtsummary)
library(stringr)
library(ggpubr)
library(ggplot2)

# --------------------------------------------
# Read tumour and patient CSVs
tumour <- read_csv("sim_av_tumour.csv")
patient <- read_csv("sim_av_patient.csv")

# Filter for breast cancer tumours (C50)
c50_tumour <- tumour %>%
  filter(SITE_ICD10_O2_3CHAR == "C50")

# Join with patient data
c50_joined <- c50_tumour %>%
  left_join(patient, by = "PATIENTID")

# Save cleaned cohort
write_csv(c50_joined, "C50_patients_joined.csv")

# Count gender codes and assign labels
gender_counts <- c50_joined %>%
  count(GENDER.x) %>%
  arrange(desc(n))

print(gender_counts)

most_frequent_code <- gender_counts$GENDER.x[1]
least_frequent_code <- gender_counts$GENDER.x[2]

c50_joined <- c50_joined %>%
  mutate(GENDER.x = factor(GENDER.x,
                           levels = c(most_frequent_code, least_frequent_code),
                           labels = c("Female", "Male")))

write_csv(c50_joined, "C50_patients_joined_w_gender.csv")

# -----------------------------------------------
# Breast Cancer Variable Summaries & Visuals
file_path <- "C50_patients_joined_w_gender.csv"
data_dir <- dirname(file_path)
C50_patients_joined <- read_csv(file_path)

# Clean and standardise QUINTILE_2019
C50_patients_joined <- C50_patients_joined %>%
  mutate(
    QUINTILE_2019_clean = as.numeric(str_extract(QUINTILE_2019, "^[1-5]"))
  ) %>%
  mutate(
    QUINTILE_2019_clean = factor(
      QUINTILE_2019_clean,
      levels = 1:5,
      labels = c("Most deprived", "Q2", "Q3", "Q4", "Least deprived")
    ),
    STAGE_BEST = as.factor(STAGE_BEST)
  )

# Create summary table by deprivation quintile - Social Differences Analysis: C50 Patients Only
table1_social <- C50_patients_joined %>%
  select(QUINTILE_2019_clean, GENDER.x, AGE, STAGE_BEST) %>%
  tbl_summary(
    by = QUINTILE_2019_clean,
    missing = "no",
    type = list(STAGE_BEST ~ "categorical")
  ) %>%
  add_p(test = everything() ~ "kruskal.test") %>%
  modify_header(label = "**Variable**") %>%
  bold_labels()

# View summary table
print(table1_social)

# Boxplot: Age at diagnosis by deprivation quintile
ggplot(C50_patients_joined, aes(x = QUINTILE_2019_clean, y = AGE)) +
  geom_boxplot() +
  labs(title = "Age at Diagnosis by Deprivation Quintile",
       x = "Deprivation Quintile",
       y = "Age at Diagnosis") +
  theme_minimal()

ggsave("age_by_deprivation.png", width = 8, height = 6)

# Bar Plot: Stage at diagnosis by deprivation quintile
ggplot(C50_patients_joined, aes(x = QUINTILE_2019_clean, fill = STAGE_BEST)) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  labs(title = "Stage at Diagnosis by Deprivation Quintile",
       x = "Deprivation Quintile",
       y = "Proportion",
       fill = "Stage") +
  theme_minimal()

ggsave("stage_by_deprivation.png", width = 8, height = 6)


# Statistical tests with gender labels in dataset
kruskal_age <- kruskal.test(AGE ~ QUINTILE_2019_clean, data = C50_patients_joined)
print(kruskal_age)

stage_table <- table(C50_patients_joined$STAGE_BEST, C50_patients_joined$QUINTILE_2019_clean)
chisq_stage <- chisq.test(stage_table)
print(chisq_stage)

# -----------------------------------------------
# Further data prep for diagnosis & treatment time analyses

breast_cancer_df <- C50_patients_joined %>%
  filter(SITE_ICD10_O2_3CHAR == "C50") %>%
  mutate(
    DIAGNOSISDATEBEST = as.Date(DIAGNOSISDATEBEST),
    DATE_FIRST_SURGERY = as.Date(DATE_FIRST_SURGERY),
    TIME_TO_SURGERY_DAYS = as.numeric(DATE_FIRST_SURGERY - DIAGNOSISDATEBEST),
    
    # New age stratification
    AGE_GROUP = case_when(
      AGE >= 18 & AGE <= 24 ~ "18-24",
      AGE >= 25 & AGE <= 34 ~ "25-34",
      AGE >= 35 & AGE <= 54 ~ "35-54",
      AGE >= 55 ~ "55+",
      TRUE ~ NA_character_
    )
  )

# ---- Summary statistics (overall) ----
summary_stats <- breast_cancer_df %>%
  summarise(
    Total_Cases = n(),
    Mean_Age = mean(AGE, na.rm = TRUE),
    Median_Age = median(AGE, na.rm = TRUE),
    Mean_Time_To_Surgery = mean(TIME_TO_SURGERY_DAYS, na.rm = TRUE),
    Median_Time_To_Surgery = median(TIME_TO_SURGERY_DAYS, na.rm = TRUE)
  )
print(summary_stats)

# ---- Frequency tables (overall) ----
table_stage <- as.data.frame(table(Stage = breast_cancer_df$STAGE_BEST))
table_er <- as.data.frame(table(ER_Status = breast_cancer_df$ER_STATUS))
table_pr <- as.data.frame(table(PR_Status = breast_cancer_df$PR_STATUS))
table_her2 <- as.data.frame(table(HER2_Status = breast_cancer_df$HER2_STATUS))

print(table_stage)
print(table_er)
print(table_pr)
print(table_her2)

# ---- Stratified summaries by new age group ----

treatment_by_age <- breast_cancer_df %>%
  filter(!is.na(AGE_GROUP)) %>%
  group_by(AGE_GROUP) %>%
  summarise(
    n = n(),
    mean_time_to_surgery = mean(TIME_TO_SURGERY_DAYS, na.rm = TRUE),
    median_time_to_surgery = median(TIME_TO_SURGERY_DAYS, na.rm = TRUE),
    sd_time_to_surgery = sd(TIME_TO_SURGERY_DAYS, na.rm = TRUE),
    min_time_to_surgery = min(TIME_TO_SURGERY_DAYS, na.rm = TRUE),
    max_time_to_surgery = max(TIME_TO_SURGERY_DAYS, na.rm = TRUE)
  ) %>%
  arrange(AGE_GROUP)
print(treatment_by_age)

stage_by_age <- breast_cancer_df %>%
  filter(!is.na(AGE_GROUP)) %>%
  count(AGE_GROUP, STAGE_BEST) %>%
  group_by(AGE_GROUP) %>%
  mutate(prop = round(100 * n / sum(n), 1)) %>%
  arrange(AGE_GROUP)
print(stage_by_age)

diag_by_age <- breast_cancer_df %>%
  filter(!is.na(AGE_GROUP)) %>%
  count(AGE_GROUP, SCREENINGSTATUSFULL_CODE) %>%
  group_by(AGE_GROUP) %>%
  mutate(prop = round(100 * n / sum(n), 1)) %>%
  arrange(AGE_GROUP)
print(diag_by_age)

outcome_by_age <- breast_cancer_df %>%
  filter(!is.na(AGE_GROUP)) %>%
  count(AGE_GROUP, VITALSTATUS) %>%
  group_by(AGE_GROUP) %>%
  mutate(prop = round(100 * n / sum(n), 1)) %>%
  arrange(AGE_GROUP)
print(outcome_by_age)

outcome_perf_by_age <- breast_cancer_df %>%
  filter(!is.na(AGE_GROUP)) %>%
  count(AGE_GROUP, PERFORMANCESTATUS) %>%
  group_by(AGE_GROUP) %>%
  mutate(prop = round(100 * n / sum(n), 1)) %>%
  arrange(AGE_GROUP)
print(outcome_perf_by_age)

treatment_by_diag <- breast_cancer_df %>%
  group_by(SCREENINGSTATUSFULL_CODE) %>%
  summarise(
    n = n(),
    mean_time = mean(TIME_TO_SURGERY_DAYS, na.rm = TRUE),
    median_time = median(TIME_TO_SURGERY_DAYS, na.rm = TRUE),
    sd_time = sd(TIME_TO_SURGERY_DAYS, na.rm = TRUE)
  )
print(treatment_by_diag)

# --------------------------------------------
# Create Word document ----
doc <- read_docx()

doc <- doc %>%
  body_add_par("Breast Cancer Data Summary", style = "heading 1") %>%
  
  body_add_par("Summary Statistics (Overall)", style = "heading 2") %>%
  body_add_flextable(flextable(summary_stats)) %>%
  
  body_add_par("Tumour Stage Frequencies (Overall)", style = "heading 2") %>%
  body_add_flextable(flextable(table_stage)) %>%
  
  body_add_par("ER Status Frequencies (Overall)", style = "heading 2") %>%
  body_add_flextable(flextable(table_er)) %>%
  
  body_add_par("PR Status Frequencies (Overall)", style = "heading 2") %>%
  body_add_flextable(flextable(table_pr)) %>%
  
  body_add_par("HER2 Status Frequencies (Overall)", style = "heading 2") %>%
  body_add_flextable(flextable(table_her2)) %>%
  
  body_add_par("Treatment Time by Age Group", style = "heading 2") %>%
  body_add_flextable(flextable(treatment_by_age)) %>%
  
  body_add_par("Stage at Diagnosis by Age Group", style = "heading 2") %>%
  body_add_flextable(flextable(stage_by_age)) %>%
  
  body_add_par("Diagnostic Methods by Age Group", style = "heading 2") %>%
  body_add_flextable(flextable(diag_by_age)) %>%
  
  body_add_par("Treatment Outcomes by Age Group (Vital Status)", style = "heading 2") %>%
  body_add_flextable(flextable(outcome_by_age)) %>%
  
  body_add_par("Treatment Outcomes by Age Group (Performance Status)", style = "heading 2") %>%
  body_add_flextable(flextable(outcome_perf_by_age))

print(doc, target = file.path(data_dir, "breast_cancer_summary.docx"))

# --------------------------------------------
# Save stratified summaries as CSV ----
write_csv(treatment_by_age, file.path(data_dir, "treatment_time_by_agegroup.csv"))
write_csv(stage_by_age, file.path(data_dir, "stage_distribution_by_agegroup.csv"))
write_csv(diag_by_age, file.path(data_dir, "diagnostic_methods_by_agegroup.csv"))
write_csv(outcome_by_age, file.path(data_dir, "treatment_outcomes_by_agegroup_vitalstatus.csv"))
write_csv(outcome_perf_by_age, file.path(data_dir, "treatment_outcomes_by_agegroup_performance.csv"))

# --------------------------------------------
# Plots ----

# Age distribution
p1 <- ggplot(breast_cancer_df, aes(x = AGE)) +
  geom_histogram(binwidth = 5, fill = "lightblue", colour = "black") +
  theme_minimal() +
  labs(title = "Age Distribution of Breast Cancer Patients", x = "Age", y = "Count")
ggsave(filename = file.path(data_dir, "age_distribution.png"), plot = p1)

# Stage distribution
p2 <- ggplot(breast_cancer_df, aes(x = STAGE_BEST)) +
  geom_bar(fill = "lightpink", colour = "black") +
  theme_minimal() +
  labs(title = "Stage at Diagnosis", x = "Stage", y = "Count")
ggsave(filename = file.path(data_dir, "stage_distribution.png"), plot = p2)

# Time to surgery distribution
p3 <- ggplot(breast_cancer_df, aes(x = TIME_TO_SURGERY_DAYS)) +
  geom_histogram(binwidth = 10, fill = "lightgreen", colour = "black") +
  theme_minimal() +
  labs(title = "Time from Diagnosis to Surgery", x = "Days", y = "Count")
ggsave(filename = file.path(data_dir, "time_to_surgery_distribution.png"), plot = p3)

# Receptor status - ER status
p4 <- ggplot(breast_cancer_df, aes(x = ER_STATUS)) +
  geom_bar(fill = "orchid", colour = "black") +
  theme_minimal() +
  labs(title = "ER Status", x = "ER Status", y = "Count")
ggsave(filename = file.path(data_dir, "er_status_distribution.png"), plot = p4)

# Boxplot: Treatment time by age group
p5 <- ggplot(breast_cancer_df, aes(x = AGE_GROUP, y = TIME_TO_SURGERY_DAYS)) +
  geom_boxplot(fill = "#74c476") +
  theme_minimal() +
  labs(
    title = "Time to Surgery by Age Group",
    x = "Age Group",
    y = "Days between diagnosis and surgery"
  )
ggsave(filename = file.path(data_dir, "time_to_surgery_by_agegroup_boxplot.png"), plot = p5)

# Stacked bar chart: Stage by age group
p6 <- ggplot(stage_by_age, aes(x = AGE_GROUP, y = prop, fill = STAGE_BEST)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_minimal() +
  labs(
    title = "Stage at Diagnosis by Age Group",
    x = "Age Group",
    y = "Percentage",
    fill = "Stage"
  )
ggsave(filename = file.path(data_dir, "stage_by_agegroup_stacked_bar.png"), plot = p6)

