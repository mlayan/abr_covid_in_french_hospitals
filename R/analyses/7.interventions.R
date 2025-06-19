##################################################
## CREATE DATAFRAME WITH INTERVENTIONS
## PERIODS AGAINST SARS-COV-2 IN THE 
## COMMUNITY
##################################################
rm(list = ls())
library(tidyverse)
source("R/helper/helper_functions.R")
source("R/helper/dictionaries.R")

##################################################
# Periods in Juliette's paper
# https://doi.org/10.1186/s12879-023-08106-1
##################################################
# Data frame of the measures before and after the study period of Juliette's paper
unstudied_periods = expand.grid(
  date = c(seq(min(all_dates), as.Date("2019-12-31"), 1), seq(as.Date("2021-05-19"), max(all_dates), 1)),
  dep = corres$dep,
  Confinement = "none",
  `Couvre.feu` = "none"
) 

# Load raw data from Juliette's paper 
int = read.csv("data-raw/interventions/Data_dep_explanatory_model_20240226_Maylis.csv", header = T) %>%
  filter(!dep %in% c("2A", "2B")) %>%
  mutate(date = as.Date(date, "%Y-%m-%d")) %>%
  bind_rows(unstudied_periods) %>%
  arrange(date, dep)

# Verifications
  # Only one intervention per day x department
  int %>% group_by(dep, date) %>% mutate(n = n()) %>% filter(n > 1)
  # No missing day
  int %>% group_by(dep) %>% mutate(d = as.numeric(difftime(date, lag(date)))) %>% filter(d != 1 & !is.na(d))
  # All types of lockdowns
  sort(unique(int$Confinement))
  # All types of curfew
  sort(unique(int$Couvre.feu))
  # Start and end of the intervention dataframe  
  c(max(int$date), max(all_dates))
  c(min(int$date), min(all_dates))

# Classification into the 4 restriction periods in our study
# See excel sheet with correspondences from Paireau's paper (https://doi.org/10.1186/s12879-023-08106-1)
first_wave = c("first", "first_light", "first_light2")
high_restrictions = c("second", "second_light", "third", "third_light")

int_dep = int %>%
  left_join(., corres, by = "dep") %>%
  mutate(
    p_first_wave = ifelse(Confinement %in% first_wave, 1 , 0),
    p_strong_res = ifelse(Confinement %in% high_restrictions, 1 , 0), 
    p_mild_res = ifelse((!Confinement %in% c(first_wave, high_restrictions, "none") & date <= as.Date("2021-05-18")) |
                          (date >= as.Date("2021-05-19") & date <= as.Date("2021-06-29")) |
                          (date >= as.Date("2021-12-10") & date <= as.Date("2022-01-24"))
                          , 1, 0),
    p_no_res = ifelse((date >= as.Date("2021-06-30") & date <= as.Date("2021-12-09")) |
                        date >= as.Date("2022-01-25")
                        , 1, 0)
  )

# https://www.associations.gouv.fr/les-restrictions-sanitaires-en-janvier-2022-adaptation-des-activites-associatives.html

##################################################
# Interventions by region
##################################################
# We consider that the strength of the intervention is determined by the 
# maximal strength across the departments of the region
int_region = int_dep %>%
  group_by(date, region) %>%
  summarise(
    p_first_wave = max(p_first_wave),
    p_strong_res = max(p_strong_res),
    p_mild_res = max(p_mild_res),
    p_no_res = max(p_no_res), 
    .groups = 'drop'
    ) %>%
  mutate(
    Date_week = as.Date(cut(date, "week")),
    p_strong_res = ifelse(p_first_wave > 0 , 0, p_strong_res),
    p_mild_res = ifelse(p_strong_res == 1 | p_first_wave > 0 , 0, p_mild_res),
    p_no_res = ifelse(p_first_wave == 1 | p_strong_res == 1 | p_mild_res == 1, 0, p_no_res)
    ) %>%
  group_by(Date_week, region) %>%
  summarise(
    p_first_wave = mean(p_first_wave),
    p_strong_res = mean(p_strong_res),
    p_mild_res = mean(p_mild_res),
    p_no_res = mean(p_no_res),
    .groups = 'drop'
    )

# Save the data 
save(int_region, file = "data/int_region.rda")

# Verifications
  # Dates of the prepandemic period
  int_region %>%
    filter(p_first_wave == 0, p_strong_res == 0, p_mild_res == 0, p_no_res == 0) %>%
    group_by(region) %>%
    summarise(m = min(Date_week), M =max(Date_week))
  # Regions with pandemic periods all at one
  int_region %>%
    filter(p_first_wave > 0, p_strong_res > 0, p_mild_res > 0, p_no_res > 0) 
  # Sum per week of pandemic periods should be equal to 1
  # Should be null except for week 2020-03-16 as the first lockdown started on 
  # Tuesday the 17th of March, 2020
  int_region %>%
    mutate(s = p_first_wave + p_strong_res + p_mild_res + p_no_res) %>%
    filter(s != 1, p_first_wave > 0 | p_strong_res > 0 | p_mild_res > 0 | p_no_res > 0, Date_week != as.Date("2020-03-16"))

##################################################
# Interventions at the national level
##################################################
# Create database 
int_national = int_dep %>%
  group_by(date) %>%
  summarise(
    p_first_wave = max(p_first_wave),
    p_strong_res = max(p_strong_res),
    p_mild_res = max(p_mild_res),
    p_no_res = max(p_no_res), 
    .groups = 'drop'
  ) %>%
  mutate(
    Date_week = as.Date(cut(date, "week")),
    p_strong_res = ifelse(p_first_wave == 1, 0, p_strong_res),
    p_mild_res = ifelse(p_first_wave == 1 | p_strong_res == 1, 0, p_mild_res),
    p_no_res = ifelse(p_first_wave == 1 | p_strong_res == 1 | p_mild_res == 1, 0, p_no_res)
  ) %>%
  group_by(Date_week) %>%
  summarise(
    p_first_wave = mean(p_first_wave),
    p_strong_res = mean(p_strong_res),
    p_mild_res = mean(p_mild_res),
    p_no_res = mean(p_no_res),
    .groups = 'drop'
  )

# Verifications
  # Dates of the prepandemic period
  int_national %>%
    filter(p_first_wave == 0, p_strong_res == 0, p_mild_res == 0, p_no_res == 0) %>%
    group_by() %>%
    summarise(m = min(Date_week), M =max(Date_week))
  # Dates with all four pandemic periods
  int_national %>%
    filter(p_first_wave > 0, p_strong_res > 0, p_mild_res > 0, p_no_res > 0) 
  # Dates with pandemic periods not equal to 1
  int_national %>%
    mutate(s = p_first_wave + p_strong_res + p_mild_res + p_no_res) %>%
    filter(s != 1, p_first_wave > 0 | p_strong_res > 0 | p_mild_res > 0 | p_no_res > 0, Date_week != as.Date("2020-03-16"))

# Modify int_national to make it more explicit for the regression analysis
int_national = int_national %>%
  mutate(
    periods = case_when(
      p_strong_res + p_mild_res + p_no_res + p_first_wave == 0 ~ "pre-pandemic",
      p_first_wave > 0.5 ~ "first wave", #  & Date_week <= as.Date("2020-04-13") 
      p_strong_res > 0.5 ~ "strong res", # | (p_first > 0.5 & Date_week > as.Date("2020-04-13")) 
      p_mild_res > 0.5 ~ "mild res",
      p_no_res > 0.5 ~ "low to no res",
      .default = NA
    ) 
  ) %>%
  dplyr::select(Date_week, periods) %>%
  mutate(periods = factor(periods, c("pre-pandemic", "first wave", "strong res", "mild res", "low to no res")))

# Save the data
save(int_national, file = "data/int_national.rda")
  
# Get dates of period changes
# We consider the most restrictive intervention to be dominant even though it 
# does not apply to all departments
int_national_start_end = int_dep %>%
  group_by(date) %>%
  summarise(
    p_first_wave = max(p_first_wave), 
    p_strong_res = max(p_strong_res),
    p_mild_res = max(p_mild_res),
    p_no_res = max(p_no_res), 
    .groups = 'drop'
  ) %>%
  mutate(
    p_strong_res = ifelse(p_first_wave == 1, 0, p_strong_res),
    p_mild_res = ifelse(p_first_wave == 1 | p_strong_res == 1, 0, p_mild_res),
    p_no_res = ifelse(p_first_wave == 1 | p_strong_res == 1 | p_mild_res == 1, 0, p_no_res),
  ) %>%
  arrange(date) %>%
  mutate(
    p_first_wave = p_first_wave - lag(p_first_wave),
    p_strong_res = p_strong_res - lag(p_strong_res),
    p_mild_res = p_mild_res - lag(p_mild_res),
    p_no_res = p_no_res - lag(p_no_res)
  ) %>%
  filter(p_first_wave != 0 | p_strong_res != 0 | p_mild_res != 0 | p_no_res != 0) %>%
  mutate(
    p_first_wave = case_when(p_first_wave == 1 ~ "start", p_first_wave == -1 ~ "end", .default = NA), 
    p_strong_res = case_when(p_strong_res == 1 ~ "start", p_strong_res == -1 ~ "end", .default = NA), 
    p_mild_res = case_when(p_mild_res == 1 ~ "start", p_mild_res == -1 ~ "end", .default = NA), 
    p_no_res = case_when(p_no_res == 1 ~ "start", p_no_res == -1 ~ "end", .default = NA)
    ) %>%
  pivot_longer(matches("^p_"), names_to = "restrictions", values_to = "date_type") %>%
  filter(!is.na(date_type)) %>%
  arrange(restrictions, date) %>%
  group_by(restrictions, date_type) %>%
  mutate(n = 1:n()) %>%
  ungroup() %>%
  pivot_wider(names_from = "date_type", values_from = "date") %>%
  arrange(start)

int_national_start_end$end[is.na(int_national_start_end$end)] = as.Date("2022-12-19")
save(int_national_start_end, file = "data/int_national_start_end.rda")

##################################################
# COVID-19 periods
##################################################
# # Dataframe of Covid-19 periods - by region
# intervention_dates_corrected = expand.grid(Date_week = all_dates, region = all_regions) %>%
#   mutate(
#     #p1 = case_when(Date_week <= as.Date("2020-03-15") ~ 1, Date_week == as.Date("2020-03-16") ~ 1.5/7, .default = 0),
#     # 1st lockdown: 17th March 2020 (tuesday) - 10th May 2020 (sunday)
#     p1 = case_when(Date_week >= as.Date("2020-03-23") & Date_week < as.Date("2020-10-26") ~ 1, Date_week == as.Date("2020-10-26") ~ 4/7, Date_week == as.Date("2020-03-16") ~ 5.5/7, .default = 0),
#     # 2nd lockdown: 30th October 2020 (friday) - 14th December 2020 (monday)
#     p2 = case_when(
#       Date_week == as.Date("2020-10-26") ~ 4/7, 
#       
#       Date_week > as.Date("2020-10-26") & Date_week < as.Date("2021-03-15") & region %in% c("Île-de-France", "Hauts-de-France") ~ 1, 
#       Date_week == as.Date("2021-03-15") & region %in% c("Île-de-France", "Hauts-de-France") ~ 4/7, 
#       
#       Date_week > as.Date("2020-10-26") & Date_week < as.Date("2021-03-29") & !region %in% c("Île-de-France", "Hauts-de-France") ~ 1,       
#       Date_week == as.Date("2021-03-29") & !region %in% c("Île-de-France", "Hauts-de-France") ~ 5/7, 
#       .default = 0
#     ),
#     # 3rd lockdown: 
#     p3 = case_when(
#       #   - 19th March 2021 (friday) - 2nd May 2021 (sunday): Île-de-France, Hauts-France, Eure (Normandie 1/5), Seine-Maritime (Normandie 1/5), Alpes-Maritimes (PACA 1/6)
#       Date_week > as.Date("2021-03-15") & Date_week < as.Date("2021-06-28") & region %in% c("Île-de-France", "Hauts-de-France") ~ 1, 
#       Date_week == as.Date("2021-03-15") & region %in% c("Île-de-France", "Hauts-de-France") ~ 3/7, 
#       
#       #   - 25th March 2021 (thursday) - 2nd May 2021 (sunday): Aube (Grand Est 1/10), Nièvre (Bourgogne-Franche-Comté 1/8), Rhône (Auvergne-Rhône-Alpes 1/12)
#       #   - 3rd April 2021 (saturday) - 2nd May 2021 (sunday): tous les départements hexagonaux
#       Date_week > as.Date("2021-03-29") & Date_week < as.Date("2021-06-28") & !region %in% c("Île-de-France", "Hauts-de-France") ~ 1, 
#       Date_week == as.Date("2021-03-29") & !region %in% c("Île-de-France", "Hauts-de-France") ~ 2/7, 
#       
#       Date_week == as.Date("2021-06-28") ~ 2/7, 
#       .default = 0
#     ),
#     # post-interventions: 30th June 2021 (wednesday)
#     p4 = case_when(Date_week == as.Date("2021-06-28") ~ 5/7, Date_week >= as.Date("2021-07-05") ~ 1, .default = 0)
#   ) 
# 
# 
# # Dataframe of Covid-19 periods - national 
# intervention_dates = data.frame(Date_week = all_dates) %>%
#   mutate(
#     #p1 = case_when(Date_week <= as.Date("2020-03-15") ~ 1, Date_week == as.Date("2020-03-16") ~ 1.5/7, .default = 0),
#     # 1st lockdown: 17th March 2020 (tuesday) - 10th May 2020 (sunday)
#     p1 = case_when(Date_week >= as.Date("2020-03-23") & Date_week < as.Date("2020-10-26") ~ 1, Date_week == as.Date("2020-10-26") ~ 4/7, Date_week == as.Date("2020-03-16") ~ 5.5/7, .default = 0),
#     # 2nd lockdown: 30th October 2020 (friday) - 14th December 2020 (monday)
#     p2 = case_when(
#       Date_week == as.Date("2020-10-26") ~ 4/7, 
#       Date_week > as.Date("2020-10-26") & Date_week < as.Date("2021-03-15") ~ 1, 
#       Date_week == as.Date("2021-03-15") ~ 4/7, 
#       .default = 0
#     ),
#     # 3rd lockdown: 
#     p3 = case_when(
#       #   - 19th March 2021 (friday) - 2nd May 2021 (sunday): Île-de-France, Hauts-France, Eure (Normandie 1/5), Seine-Maritime (Normandie 1/5), Alpes-Maritimes (PACA 1/6)
#       Date_week > as.Date("2021-03-15") & Date_week < as.Date("2021-06-28") ~ 1, 
#       Date_week == as.Date("2021-03-15") ~ 3/7, 
#       Date_week == as.Date("2021-06-28") ~ 2/7, 
#       .default = 0
#     ),
#     # post-interventions: 30th June 2021 (wednesday)
#     p4 = case_when(Date_week == as.Date("2021-06-28") ~ 5/7, Date_week >= as.Date("2021-07-05") ~ 1, .default = 0)
#   ) 
# 
# # Data frame of Covid-19 restriction type periods - national
# intervention_dates = data.frame(Date_week = all_dates) %>%
#   mutate(
#     # 1st lockdown: 17th March 2020 (tuesday) - 10th May 2020 (sunday)
#     # 2nd lockdown: 30th October 2020 (friday) - 14th December 2020 (monday)
#     # 3rd lockdown: 19th March 2021 (friday) - 2nd May 2021 (sunday): Île-de-France, Hauts-France, Eure (Normandie 1/5), Seine-Maritime (Normandie 1/5), Alpes-Maritimes (PACA 1/6)
#     p_very_restrictive = case_when(
#       Date_week >= as.Date("2020-03-23") & Date_week < as.Date("2020-10-26") ~ 1, 
#       Date_week == as.Date("2020-10-26") ~ 4/7, 
#       Date_week == as.Date("2020-03-16") ~ 5.5/7, 
#       Date_week == as.Date("2020-10-26") ~ 4/7, 
#       Date_week > as.Date("2020-10-26") & Date_week < as.Date("2021-03-15") ~ 1, 
#       Date_week == as.Date("2021-03-15") ~ 4/7, 
#       Date_week > as.Date("2021-03-15") & Date_week < as.Date("2021-06-28") ~ 1, 
#       Date_week == as.Date("2021-03-15") ~ 3/7, 
#       Date_week == as.Date("2021-06-28") ~ 2/7, 
#       .default = 0
#     ),
#     
#     # After the 1st lockdowns, cinemas and theaters reopened on the 22nd of June 2020 
#     # After the 2nd lockdown, moderate restrictions
#     # After the 3rd lockdown, cinemas, theaters, shops and restaurants reopened on the 19 of May, 2021
#     p_int_restrictive = case_when(
#       Date_week <= 
#         
#         .default = 0
#     ),
#     
#     # post-interventions: 30th June 2021 (wednesday)
#     p_not_restrictive = case_when(Date_week == as.Date("2021-06-28") ~ 5/7, Date_week >= as.Date("2021-07-05") ~ 1, .default = 0)
#   ) 
