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
# Install and load packages

install.packages(c("tidyverse", "readr", "data.table", "officer", "flextable", "dplyr", "gtsummary", "stringr", "ggpubr", "scales"))

library(tidyverse)
library(readr)
library(officer)
library(flextable)
library(dplyr)
library(gtsummary)
library(stringr)
library(ggpubr)
library(ggplot2)
library(scales)


# Set colour palettes
okabe_ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")

colour_scheme <- c("#285238",  # Dark green
                   "#94C9A9",   # Light green
                   "#B7C3F3",   # Light lavender
                   "#453F78",   # Deep navy
                   "#007E8B"    # Teal
)

# Create outputs directory if it doesn't exist
if (!dir.exists("outputs")) {
  dir.create("outputs")
}

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

# Create summary table by deprivation quintile
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

# Boxplot: Age at diagnosis by deprivation quintile
ggplot(C50_patients_joined, aes(x = QUINTILE_2019_clean, y = AGE, fill = QUINTILE_2019_clean)) +
  geom_boxplot() +
  scale_fill_manual(values = colour_scheme) +  # Using the defined scheme
  labs(title = "Age at Diagnosis by Deprivation Quintile",
       x = "Deprivation Quintile",
       y = "Age at Diagnosis") +
  theme_minimal() +
  theme(legend.position = "none")  # Remove legend for cleaner look

ggsave("age_by_deprivation.png", width = 8, height = 6)


# Statistical tests
kruskal_age <- kruskal.test(AGE ~ QUINTILE_2019_clean, data = C50_patients_joined)
stage_table <- table(C50_patients_joined$STAGE_BEST, C50_patients_joined$QUINTILE_2019_clean)
chisq_stage <- chisq.test(stage_table)

# Create breast cancer dataframe with additional variables
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

# Clean STAGE_BEST 
stage_by_age_clean <- breast_cancer_df %>%
  filter(!is.na(AGE_GROUP)) %>%
  mutate(
    STAGE_BEST_CLEAN = case_when(
      STAGE_BEST %in% c("0", "1", "2", "3", "4") ~ STAGE_BEST,
      TRUE ~ NA_character_  # Group all others as NA
    )
  ) %>%
  count(AGE_GROUP, STAGE_BEST_CLEAN) %>%
  group_by(AGE_GROUP) %>%
  mutate(prop = round(100 * n / sum(n), 1)) %>%
  arrange(AGE_GROUP, STAGE_BEST_CLEAN)

# ---- Summary statistics (overall) ----
summary_stats <- breast_cancer_df %>%
  summarise(
    Total_Cases = n(),
    Mean_Age = mean(AGE, na.rm = TRUE),
    Median_Age = median(AGE, na.rm = TRUE),
    Mean_Time_To_Surgery = mean(TIME_TO_SURGERY_DAYS, na.rm = TRUE),
    Median_Time_To_Surgery = median(TIME_TO_SURGERY_DAYS, na.rm = TRUE)
  )

# ---- Frequency tables (overall) ----
table_stage <- as.data.frame(table(Stage = breast_cancer_df$STAGE_BEST))
table_er <- as.data.frame(table(ER_Status = breast_cancer_df$ER_STATUS))
table_pr <- as.data.frame(table(PR_Status = breast_cancer_df$PR_STATUS))
table_her2 <- as.data.frame(table(HER2_Status = breast_cancer_df$HER2_STATUS))

# Age group analyses
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

stage_by_age <- breast_cancer_df %>%
  filter(!is.na(AGE_GROUP)) %>%
  count(AGE_GROUP, STAGE_BEST) %>%
  group_by(AGE_GROUP) %>%
  mutate(prop = round(100 * n / sum(n), 1)) %>%
  arrange(AGE_GROUP)

diag_by_age <- breast_cancer_df %>%
  filter(!is.na(AGE_GROUP)) %>%
  count(AGE_GROUP, SCREENINGSTATUSFULL_CODE) %>%
  group_by(AGE_GROUP) %>%
  mutate(prop = round(100 * n / sum(n), 1)) %>%
  arrange(AGE_GROUP)

outcome_by_age <- breast_cancer_df %>%
  filter(!is.na(AGE_GROUP)) %>%
  count(AGE_GROUP, VITALSTATUS) %>%
  group_by(AGE_GROUP) %>%
  mutate(prop = round(100 * n / sum(n), 1)) %>%
  arrange(AGE_GROUP)

outcome_perf_by_age <- breast_cancer_df %>%
  filter(!is.na(AGE_GROUP)) %>%
  count(AGE_GROUP, PERFORMANCESTATUS) %>%
  group_by(AGE_GROUP) %>%
  mutate(prop = round(100 * n / sum(n), 1)) %>%
  arrange(AGE_GROUP)

treatment_by_diag <- breast_cancer_df %>%
  group_by(SCREENINGSTATUSFULL_CODE) %>%
  summarise(
    n = n(),
    mean_time = mean(TIME_TO_SURGERY_DAYS, na.rm = TRUE),
    median_time = median(TIME_TO_SURGERY_DAYS, na.rm = TRUE),
    sd_time = sd(TIME_TO_SURGERY_DAYS, na.rm = TRUE)
  )

# --------------------------------------------
# Create Word document
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

print(doc, target = file.path("outputs", "breast_cancer_summary.docx"))

# Create visualisations

# Common theme for all plots
my_theme <- function() {
  theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white"),
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_blank(),
      legend.background = element_rect(fill = "white", color = NA),
      legend.box.background = element_rect(fill = "white", color = NA)
    )
}

# ----------------------------------------------------------

# PLOT 1: Diagnostic tools by age group - Color by diagnostic method
p_diag_methods <- ggplot(
  breast_cancer_df %>% 
    filter(!is.na(AGE_GROUP)) %>% 
    count(AGE_GROUP, SCREENING_CLEAN) %>% 
    group_by(AGE_GROUP) %>% 
    mutate(prop = n/sum(n)), 
  aes(x = AGE_GROUP, y = prop, fill = SCREENING_CLEAN)) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  scale_fill_manual(values = colour_scheme[1:3]) +  # Colors for diagnostic methods
  scale_y_continuous(labels = scales::percent_format()) +
  labs(title = "Diagnostic Methods by Age Group",
       x = "Age Group",
       y = "Proportion of Cases",
       fill = "Diagnostic Method") +
  my_theme()

# Time to diagnosis - Color by diagnostic method
p_diag_time <- ggplot(
  breast_cancer_df %>% 
    filter(!is.na(SCREENING_CLEAN) & !is.na(AGE_GROUP)),
  aes(x = AGE_GROUP, y = TIME_TO_SURGERY_DAYS, fill = SCREENING_CLEAN)) +
  geom_boxplot(outlier.shape = NA) +
  coord_cartesian(ylim = quantile(breast_cancer_df$TIME_TO_SURGERY_DAYS, c(0.1, 0.9), na.rm = TRUE)) +
  scale_fill_manual(values = colour_scheme[1:3]) +  # Colors for diagnostic methods
  labs(title = "Time to Surgery by Diagnostic Method and Age Group",
       x = "Age Group",
       y = "Days from Diagnosis to Surgery",
       fill = "Diagnostic Method") +
  my_theme()

# PLOT 2: Social group variation plots

# Time to surgery by deprivation - Color by quintile
p_deprivation_time <- ggplot(
  breast_cancer_df %>% 
    filter(!is.na(QUINTILE_2019_clean)),
  aes(x = QUINTILE_2019_clean, y = TIME_TO_SURGERY_DAYS, fill = QUINTILE_2019_clean)) +
  geom_boxplot(outlier.shape = NA) +
  coord_cartesian(ylim = quantile(breast_cancer_df$TIME_TO_SURGERY_DAYS, c(0.1, 0.9), na.rm = TRUE)) +
  scale_fill_manual(values = colour_scheme[1:5]) +  # Different color per quintile
  labs(title = "Time to Surgery by Deprivation Quintile",
       x = "Deprivation Quintile (1 = Most deprived)",
       y = "Days from Diagnosis to Surgery",
       fill = "Deprivation Quintile") +
  my_theme()

# Stage distribution by deprivation - Color by stage
p_deprivation_stage <- ggplot(
  breast_cancer_df %>% 
    filter(!is.na(QUINTILE_2019_clean)) %>% 
    mutate(
      STAGE_CLEAN = case_when(
        STAGE_BEST %in% c("0","0A") ~ "0", 
        STAGE_BEST %in% c("1","1A","1A1","1A2","1A3","1B","1B1","1C","1C1","1C2","1C3") ~ "1",
        STAGE_BEST %in% c("2","2A","2A1","2A2","2B","2C") ~ "2", 
        STAGE_BEST %in% c("3","3A","3A1","3A1ii","3A2","3B","3C") ~ "3",
        STAGE_BEST %in% c("4","4A","4B","4C") ~ "4", 
        STAGE_BEST %in% c("?","U","A", "Unknown") ~ "Other/Unknown", 
        TRUE ~ "Other/Unknown"
      )) %>% 
    count(QUINTILE_2019_clean, STAGE_CLEAN) %>% 
    group_by(QUINTILE_2019_clean) %>% 
    mutate(prop = n/sum(n)),
  aes(x = QUINTILE_2019_clean, y = prop, fill = STAGE_CLEAN)) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  scale_fill_manual(values = c(colour_scheme[1:5],"grey70")) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(title = "Stage Distribution by Deprivation Quintile", 
       x = "Deprivation Quintile (1 = Most deprived)",
       y = "Proportion of Cases",
       fill = "Stage") +
  my_theme() +
  theme(
    axis.text.x = element_text(size = 12, color = "black", hjust = 0.5),
    axis.text.y = element_text(size = 12, color = "black"),
    plot.title = element_text(size = 14, face = "bold")
  )

# Age distribution by deprivation - Color by quintile
p_deprivation_age <- ggplot(
  breast_cancer_df %>% 
    filter(!is.na(QUINTILE_2019_clean)),
  aes(x = QUINTILE_2019_clean, y = AGE, fill = QUINTILE_2019_clean)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.2, fill = "white", outlier.shape = NA) +
  scale_fill_manual(values = colour_scheme[1:5]) +
  labs(title = "Age Distribution by Deprivation Quintile",
       x = "Deprivation Quintile (1 = Most deprived)",
       y = "Age at Diagnosis",
       fill = "Deprivation Quintile") +
  my_theme() +
  theme(legend.position = "none")

# Stage distribution by age group - Color by stage
p_stage_age <- ggplot(
  breast_cancer_df %>% 
    filter(!is.na(AGE_GROUP)) %>% 
    mutate(
      STAGE_CLEAN = case_when(
        STAGE_BEST %in% c("0","0A") ~ "0", 
        STAGE_BEST %in% c("1","1A","1A1","1A2","1A3","1B","1B1","1C","1C1","1C2","1C3") ~ "1",
        STAGE_BEST %in% c("2","2A","2A1","2A2","2B","2C") ~ "2", 
        STAGE_BEST %in% c("3","3A","3A1","3A1ii","3A2","3B","3C") ~ "3",
        STAGE_BEST %in% c("4","4A","4B","4C") ~ "4", 
        STAGE_BEST %in% c("?","U","A", "Other") ~ "Other/Unknown", 
        TRUE ~ "Other/Unknown"
      )) %>% 
    count(AGE_GROUP, STAGE_CLEAN) %>% 
    group_by(AGE_GROUP) %>% 
    mutate(prop = n/sum(n)),
  aes(x = AGE_GROUP, y = prop, fill = STAGE_CLEAN)) +
  geom_bar(stat = "identity", position = "fill", colour = "black", width = 0.5) +  # Control bar width (default is 0.9)) 
  scale_fill_manual(values = c(colour_scheme[1:5], "grey70")) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(title = "Stage Distribution by Age Group", 
       x = "Age Group",
       y = "Percentage of Cases",
       fill = "Stage") +
  my_theme() +
  theme(
    axis.text.x = element_text(size = 12, color = "black", hjust = 0.5),
    axis.text.y = element_text(size = 12, color = "black"),
    plot.title = element_text(size = 14, face = "bold")
  )

# ----------------------------------------------------------
# Save all plots
ggsave("outputs/diagnostic_methods_by_age.png", p_diag_methods, width = 8, height = 6, bg = "white")
ggsave("outputs/diagnostic_time_by_method_age.png", p_diag_time, width = 8, height = 6, bg = "white")
ggsave("outputs/time_by_deprivation.png", p_deprivation_time, width = 8, height = 6, bg = "white")
ggsave("outputs/stage_by_deprivation.png", p_deprivation_stage, width = 8, height = 6, bg = "white")
ggsave("outputs/age_by_deprivation.png", p_deprivation_age, width = 8, height = 6, bg = "white")
ggsave("outputs/stage_distribution_by_age.png", p_stage_age, width = 9, height = 6, bg = "white")
