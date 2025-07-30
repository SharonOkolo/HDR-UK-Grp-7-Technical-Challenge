# R packages required
library(gtsummary)
library(tidyverse)
library(tidyr)
library(ggplot2)
library(lubridate)
library(dplyr)

# Set directory
getwd()
setwd(dir = "C:/Users/BlossomOkparocha/Downloads")

# uploaded both the tumour and patient data sets
patient_data <- read.csv(file = "simulacrum_v2.1.0/simulacrum_v2.1.0/Data/sim_av_patient.csv")
tumour_data <- read.csv(file = "simulacrum_v2.1.0/simulacrum_v2.1.0/Data/sim_av_tumour.csv")

# only tasked to look at C50 breast cancer patients
# code removes patients that are not C50
tumour_data <- filter(tumour_data, SITE_ICD10_O2_3CHAR == "C50", .by = NULL, .preserve = TRUE)

# removes unnecessary columns 
tumour_data <- tumour_data[,-c(6:12)]
tumour_data <- tumour_data[,-c(6:10)]
first_tumour_data <- tumour_data[,-c(9, 10, 12:15, 18, 20, 22:25)]

# changed the class of diagnosis date and first treatment date from character to date
first_tumour_data$DIAGNOSISDATEBEST <- ymd(first_tumour_data$DIAGNOSISDATEBEST)
first_tumour_data$DATE_FIRST_SURGERY <- ymd(first_tumour_data$DATE_FIRST_SURGERY)

# created a new variable for time until treatment, and inputed the variable column back into the data frame
date_diff <- first_tumour_data$DATE_FIRST_SURGERY - first_tumour_data$DIAGNOSISDATEBEST
head(date_diff)  
first_tumour_data$date_diff <- date_diff 

# changed patientID from numeric to character
first_tumour_data$PATIENTID <- trimws(as.character(first_tumour_data$PATIENTID))
patient_data$PATIENTID <- trimws(as.character(patient_data$PATIENTID))

# new data frame only has patientIDs in both fixed tumour frame and patient ID frame
completely_fixed_data <- first_tumour_data %>%
  left_join(patient_data %>% select(PATIENTID, ETHNICITY), by = "PATIENTID")

# optional: removed unknown ethnicities
class(completely_fixed_data$ETHNICITY)
completely_fixed_data <- completely_fixed_data[completely_fixed_data$ETHNICITY != "X", ]
completely_fixed_data <- completely_fixed_data[completely_fixed_data$ETHNICITY != "Z", ]

# changed time until treatment to a numeric value
completely_fixed_data$date_diff <- as.numeric(completely_fixed_data$date_diff, units = "days")
class(completely_fixed_data$date_diff)
