library(gtsummary)
library(tidyverse)
library(tidyr)
library(ggplot2)
library(lubridate)
library(dplyr)
getwd()
setwd(dir = "C:/Users/BlossomOkparocha/Downloads")

patient_data <- read.csv(file = "simulacrum_v2.1.0/simulacrum_v2.1.0/Data/sim_av_patient.csv")
tumour_data <- read.csv(file = "simulacrum_v2.1.0/simulacrum_v2.1.0/Data/sim_av_tumour.csv")
tumour_data <- filter(tumour_data, SITE_ICD10_O2_3CHAR == "C50", .by = NULL, .preserve = TRUE)
tumour_data <- tumour_data[,-c(6:12)]
tumour_data <- tumour_data[,-c(6:10)]
first_tumour_data <- tumour_data[,-c(9, 10, 12:15, 18, 20, 22:25)]
first_tumour_data$DIAGNOSISDATEBEST <- ymd(first_tumour_data$DIAGNOSISDATEBEST)
class(first_tumour_data$DIAGNOSISDATEBEST)
first_tumour_data$DATE_FIRST_SURGERY <- ymd(first_tumour_data$DATE_FIRST_SURGERY)
class(first_tumour_data$DATE_FIRST_SURGERY)
date_diff <- first_tumour_data$DATE_FIRST_SURGERY - first_tumour_data$DIAGNOSISDATEBEST
head(date_diff)  
first_tumour_data$date_diff <- date_diff


first_tumour_data$PATIENTID <- trimws(as.character(first_tumour_data$PATIENTID))
patient_data$PATIENTID <- trimws(as.character(patient_data$PATIENTID))
completely_fixed_data <- first_tumour_data %>%
  left_join(patient_data %>% select(PATIENTID, ETHNICITY), by = "PATIENTID")

class(completely_fixed_data$ETHNICITY)
completely_fixed_data <- completely_fixed_data[completely_fixed_data$ETHNICITY != "X", ]
completely_fixed_data <- completely_fixed_data[completely_fixed_data$ETHNICITY != "Z", ]
completely_fixed_data$date_diff <- as.numeric(completely_fixed_data$date_diff, units = "days")
class(completely_fixed_data$date_diff)
