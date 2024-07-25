rm(list = ls())
library(tidyverse)
library(haven)
library(sf)
library(geofacet)
library(RColorBrewer)

load("data/metadata_admin_espic.rda")
load("data/cohort_final.rda")
load("data/icu_cohort_final.rda")
load("data/france.rda")
load("data/my_regional_grid.rda")

my_grid = subset(my_regional_grid, name != "Corse")

################################################################################
# Load data
################################################################################
# Load intubation data
intub = read_sas("data-raw/atih/mco_intubation_weekly.sas7bdat") %>%
  rename(finess = FinessGeo) %>%
  left_join(., metadata_admin_espic %>% dplyr::select(finess, region) %>% distinct(), by = "finess") %>%
  filter(!is.na(region)) %>%
  group_by(Date_week, region, status) %>%
  summarise(nintub = sum(nbjh), .groups = "drop")

# Load number of bed-days
bd = read_sas("data-raw/atih/mco_sejours_weekly.sas7bdat") %>%
  rename(finess = FinessGeo) %>%
  left_join(., metadata_admin_espic %>% dplyr::select(finess, region) %>% distinct(), by = "finess") %>%
  filter(!is.na(region)) %>%
  group_by(Date_week, region) %>%
  summarise(nbjh = sum(nbjh), .groups = "drop")

# Load number of COVID-19 DP
covid_dp = read_sas("data-raw/atih/mco_covid_weekly.sas7bdat") %>%
  rename(finess = FinessGeo) %>%
  left_join(., metadata_admin_espic %>% dplyr::select(finess, region) %>% distinct(), by = "finess") %>%
  filter(diag_type %in% "DP", !is.na(region)) %>%
  group_by(Date_week, region) %>%
  summarise(ncovid = sum(covid_jh), .groups = "drop")

# Load number of COVID-19 
covid = read_sas("data-raw/atih/mco_covid_weekly.sas7bdat") %>%
  rename(finess = FinessGeo) %>%
  left_join(., metadata_admin_espic %>% dplyr::select(finess, region) %>% distinct(), by = "finess") %>%
  filter(!is.na(region)) %>%
  group_by(Date_week, region) %>%
  summarise(ncovid = sum(covid_jh), .groups = "drop")

################################################################################
# Plots 
################################################################################
# Plot - Counts
intub %>%
  ggplot(., aes(x = as.Date(Date_week), y = nintub, col = status)) +
  annotate("rect", xmin = as.Date("2020-03-16"), xmax = as.Date("2020-05-10"), ymin = 0, ymax = Inf, fill = "gray90") +
  annotate("rect", xmin = as.Date("2020-10-30"), xmax = as.Date("2020-12-15"), ymin = 0, ymax = Inf, fill = "gray90") +
  geom_line() +
  facet_geo(facets = vars(region), grid = my_grid) +
  scale_color_manual(values = c("coral", "darkorchid")) +
  scale_y_continuous(trans = "log10") +
  theme_bw() +
  labs(x = "", y = "Counts of intubated patients by week", col = "")
ggsave("plots/pmsi_data_cleaning/intubation_counts.png", height = 8,
       width = 10)  

# Plot - Prevalence 
intub %>%
  left_join(., bd, by = c("region", "Date_week")) %>%
  mutate(prevalence = nintub / nbjh * 1000) %>%
  ggplot(., aes(x = as.Date(Date_week), y = prevalence, col = status)) +
  annotate("rect", xmin = as.Date("2020-03-16"), xmax = as.Date("2020-05-10"), ymin = 0, ymax = Inf, fill = "gray90") +
  annotate("rect", xmin = as.Date("2020-10-30"), xmax = as.Date("2020-12-15"), ymin = 0, ymax = Inf, fill = "gray90") +
  geom_line() +
  facet_geo(facets = vars(region), grid = my_grid) +
  scale_color_manual(values = c("coral", "darkorchid")) +
  theme_bw() +
  labs(x = "", y = "Prevalence of intubated patients\nby week and 1,000 bed-days", col = "")
ggsave("plots/pmsi_data_cleaning/intubation_prevalence.png", height = 7,
       width = 12)  


# Plot total intubation versus COVID-19 DP
intub %>%
  pivot_wider(names_from = status, values_from = nintub) %>%
  left_join(., covid_dp, by = c("region", "Date_week")) %>%
  left_join(., bd, by = c("region", "Date_week")) %>%
  mutate(
    ncovid = ifelse(is.na(ncovid), 0, ncovid),
    covid = ifelse(is.na(covid), 0, covid)
    ) %>%
  mutate(
    pintub_other = other / nbjh * 1000,
    pintub_covid = covid / nbjh * 1000,
    pintub = (other+covid) / nbjh * 1000,
    pcovid = ncovid / nbjh * 1000
    ) %>%
  dplyr::select(Date_week, region, matches("^p")) %>%
  pivot_longer(matches("^p"), names_to = "prevalence", values_to = "values") %>%
  ggplot(., aes(x = as.Date(Date_week), y = values, col = prevalence)) +
  annotate("rect", xmin = as.Date("2020-03-16"), xmax = as.Date("2020-05-10"), ymin = 0, ymax = Inf, fill = "gray90") +
  annotate("rect", xmin = as.Date("2020-10-30"), xmax = as.Date("2020-12-15"), ymin = 0, ymax = Inf, fill = "gray90") +
  geom_line() +
  facet_geo(facets = vars(region), grid = my_grid) +
  theme_bw() +
  labs(x = "", y = "Prevalence of intubated patients\nby week and 1,000 bed-days", col = "")
ggsave("plots/pmsi_data_cleaning/intubation_covid_prevalence.png", height = 7,
       width = 12)  

# Compare COVID-19 cases and intubated COVID-19 cases at the national level
bd_nat = bd %>%
  group_by(Date_week) %>%
  summarise(nbjh = sum(nbjh), .groups = "drop")

covid_nat = covid_dp %>%
  group_by(Date_week) %>%
  summarise(ncovid = sum(ncovid), .groups = "drop")

intub %>%
  pivot_wider(values_from = nintub, names_from = status) %>%
  mutate(covid = ifelse(is.na(covid), 0, covid), 
         nintub = case_when(is.na(covid) ~ other, 
                            is.na(other) ~ covid, 
                            .default = covid+other)) %>%
  group_by(Date_week) %>%
  summarise(covid_intub = sum(covid), intub = sum(nintub), .groups = "drop") %>%
  left_join(., covid_nat, by = "Date_week") %>%
  left_join(., bd_nat, by = "Date_week") %>%
  mutate(
    ncovid = ifelse(is.na(ncovid), 0, ncovid),
    covid_intub = ifelse(is.na(covid_intub), 0, covid_intub),
    intub = ifelse(is.na(intub), 0, intub)
  ) %>%
  mutate(
    pintub_covid = covid_intub / nbjh * 1000,
    pintub = intub / nbjh * 1000,
    pcovid = ncovid / nbjh * 1000
  ) %>%
  dplyr::select(Date_week, pintub_covid, pintub, pcovid) %>%
  pivot_longer(matches("^p"), names_to = "prevalence", values_to = "values") %>%
  ggplot(., aes(x = as.Date(Date_week), y = values, col = prevalence)) +
  annotate("rect", xmin = as.Date("2020-03-16"), xmax = as.Date("2020-05-10"), ymin = 0, ymax = Inf, fill = "gray90") +
  annotate("rect", xmin = as.Date("2020-10-30"), xmax = as.Date("2020-12-15"), ymin = 0, ymax = Inf, fill = "gray90") +
  geom_line() +
  theme_bw() +
  labs(x = "", y = "Prevalence of intubated patients\nby week and 1,000 bed-days", col = "")
ggsave("plots/pmsi_data_cleaning/intubation_covid_prevalence.png", height = 8,
       width = 10) 

# Plot proportion of intubated patients among COVID-19 patients
intub %>%
  filter(status == "covid") %>%
  full_join(., covid, by = c("Date_week", "region")) %>%
  mutate(
    nintub = ifelse(is.na(nintub), 0, nintub), 
    ncovid = ifelse(is.na(ncovid), 0, ncovid)
  ) %>%
  mutate(prop = nintub/ncovid) %>%
  filter(region != "Corse") %>%
  ggplot(., aes(x = as.Date(Date_week), y = prop)) +
  geom_line() +
  facet_geo(facets = vars(region), grid = my_grid) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 7)) +
  labs(x = "", y = "Proportion of COVID-19 bed-days with intubation")
ggsave("plots/pmsi_data_cleaning/intubation_proportion.png", height = 7, width = 10)

