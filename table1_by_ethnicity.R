completely_fixed_data$date_diff <- as.numeric(completely_fixed_data$date_diff, units = "days")
colnames(completely_fixed_data)
completely_fixed_data %>% select(c(date_diff, GENDER)) %>%
  tbl_summary(by = GENDER,
              type = list(date_diff ~ "continuous"),
    statistic = list(date_diff ~ "{mean} (±{sd})")) %>% 
  add_overall() %>%
  add_p(test = list(date_diff ~ "t.test"))

completely_fixed_data %>% select(c(date_diff, ETHNICITY)) %>%
  tbl_summary(by = ETHNICITY,
              type = list(date_diff ~ "continuous"),
              statistic = list(date_diff ~ "{mean} (±{sd})")) %>% 
  add_overall() %>%
  add_p(test = list(date_diff ~ "aov"))

completely_fixed_data %>% select(c(date_diff, QUINTILE_2019)) %>%
  tbl_summary(by = QUINTILE_2019,
              type = list(date_diff ~ "continuous"),
              statistic = list(date_diff ~ "{mean} (±{sd})")) %>% 
  add_overall() %>%
  add_p(test = list(date_diff ~ "aov"))

completely_fixed_data %>% select(c(date_diff, PERFORMANCESTATUS)) %>%
  tbl_summary(by = PERFORMANCESTATUS,
              type = list(date_diff ~ "continuous"),
              statistic = list(date_diff ~ "{mean} (±{sd})")) %>% 
  add_overall() %>%
  add_p(test = list(date_diff ~ "aov"))
