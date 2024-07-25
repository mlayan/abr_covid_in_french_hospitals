rm(list = ls())

# Data management and visualisation
library(tidyverse)
library(haven)
library(ggpubr)

# Helper functions
source("R/helper_functions.R")
source("R/helper_plots.R")

# Basic data and dictionaries 
load("data/metadata_admin_espic.rda")
load("data/int_national.rda")
load("data/cohort_final.rda")

##################################################
# Load data 
##################################################
# Resistances
res = read.table("data-raw/spares/combined/resistance_cohortfinal.txt", header = T, sep = "\t") %>%
  mutate(Date_week = as.Date(Date_week))

# Antibiotic consumption
atb = read.table("data-raw/spares/combined/antibiotics_cohortfinal.txt", header = T, sep = "\t") %>%
  filter(secteur == "Total établissement") %>%
  group_by(Date_year, region) %>%
  summarise(
    Penicillins  = sum(Penicillins),
    Third_generation_Cephalosporins  = sum(Third_generation_Cephalosporins),
    Carbapenems  = sum(Carbapenems),
    Nbhosp  = sum(Nbhosp),
    .groups = "drop"
  ) 

# PMSI JH
pmsi_jh = read.table("data-raw/atih/pmsi_sejours.txt", header = T, sep = "\t") %>%
  group_by(Date_week, region) %>%
  summarise(nbjh = sum(nbjh), .groups = "drop") %>%
  mutate(Date_week = as.Date(Date_week))

pmsi_jh_icu = read.table("data-raw/atih/pmsi_sejours.txt", header = T, sep = "\t") %>%
  filter(secteur == "Réanimation") %>%
  group_by(Date_week, region) %>%
  summarise(nbjh = sum(nbjh), .groups = "drop") %>%
  mutate(Date_week = as.Date(Date_week))

# PMSI - COVID
pmsi_covid = read.table("data-raw/atih/pmsi_sejours_covid.txt", header = T, sep = "\t") %>%
  group_by(Date_week, region) %>%
  summarise(covid_jh  = sum(covid_jh), .groups = "drop")  

pmsi_covid_icu = read.table("data-raw/atih/pmsi_sejours_covid.txt", header = T, sep = "\t") %>%
  filter(secteur == "Réanimation") %>%
  group_by(Date_week, region) %>%
  summarise(covid_jh  = sum(covid_jh), .groups = "drop") 

# PMSI - Intubation
pmsi_intubation = read_sas("data-raw/atih/mco_intubation_weekly.sas7bdat") %>%
  rename(finess = FinessGeo) %>%
  inner_join(., metadata_admin_espic %>% filter(code %in% cohort_final) %>% dplyr::select(finess, region) %>% distinct(), by = "finess") %>%
  filter(!is.na(region), status == "covid") %>%
  group_by(Date_week, region) %>%
  summarise(nintub = sum(nbjh), .groups = "drop")

pmsi_intubation_icu = read_sas("data-raw/atih/mco_intubation_weekly.sas7bdat") %>%
  rename(finess = FinessGeo) %>%
  inner_join(., metadata_admin_espic %>% filter(code %in% cohort_final) %>% dplyr::select(finess, region) %>% distinct(), by = "finess") %>%
  filter(!is.na(region), status == "covid", secteur == "Réanimation") %>%
  group_by(Date_week, region) %>%
  summarise(nintub = sum(nbjh), .groups = "drop")

##################################################
# Correlation between Covid-19 data and incidence 
# of resistant infections
##################################################
# Proportion of Covid-19 patients that are intubated in hospitals
df_corr_etab = res %>%
  group_by(bacterie, Date_week) %>%
  summarise(n_res = sum(n_res), .groups = "drop") %>%
  left_join(., pmsi_intubation %>% group_by(Date_week) %>% summarise(nintub = sum(nintub), .groups = "drop"), by = c("Date_week")) %>%
  left_join(., pmsi_jh %>% group_by(Date_week) %>% summarise(nbjh = sum(nbjh), .groups = "drop"), by = c("Date_week")) %>%
  filter(!is.na(nintub)) %>%
  mutate(covid_intub_prev = nintub / nbjh * 1000, res_i = n_res / nbjh * 1000, secteur = "Hospital") %>%
  dplyr::select(bacterie, secteur, Date_week, res_i, covid_intub_prev)

# Data with incidence of resistant infections 
# and proportion of Covid-19 patients that are intubated in ICUs
df_corr_icu = res %>%
  filter(secteur == "ICU") %>%
  group_by(bacterie, Date_week) %>%
  summarise(n_res = sum(n_res), .groups = "drop") %>%
  left_join(., pmsi_intubation_icu %>% group_by(Date_week) %>% summarise(nintub = sum(nintub), .groups = "drop"), by = c("Date_week")) %>%
  left_join(., pmsi_jh_icu %>% group_by(Date_week) %>% summarise(nbjh = sum(nbjh), .groups = "drop"), by = c("Date_week")) %>%
  filter(!is.na(nintub)) %>%
  mutate(covid_intub_prev = nintub / nbjh * 1000, res_i = n_res / nbjh * 1000, secteur = "ICU") %>%
  dplyr::select(bacterie, secteur, Date_week, res_i, covid_intub_prev)

# Dataframes of Pearson correlation
df_corr = bind_rows(df_corr_etab, df_corr_icu) %>%
  group_by(secteur, bacterie) %>%
  summarise(corr = paste0(
    'tau = ',
    round(cor.test(res_i, covid_intub_prev, method = "kendall")$estimate, 3), 
    "\np = ",
    round(cor.test(res_i, covid_intub_prev, method = "kendall")$p.value, 3)
  ), 
  .groups = "drop"
  ) %>%
  mutate(ymax = ifelse(secteur == "Hospital", 0.3, 2.50),
         xmax = ifelse(secteur == "Hospital", 15, 250))

# Covid-19 prevalence vs incidence of resistant infections
p1 = bind_rows(df_corr_etab, df_corr_icu) %>%
  ggplot(., aes(y = res_i, x = covid_intub_prev)) +
  geom_point(col = "grey70") +
  geom_text(data = df_corr, aes(x = xmax, y = ymax, label = corr)) +
  theme_bw() +
  ggh4x::facet_grid2(cols = vars(bacterie), rows = vars(secteur), scales = "free", independent  ="x") +
  labs(x = "Weekly number of intubated Covid-19 bed-days (for 1,000 bed-days)", 
       y = "Weekly incidence of resistant infections\n(for 1,000 bed-days)")
p1

# Covid-19 periods
int_national = int_national %>%
  mutate(p_first = ifelse(Date_week >= as.Date("2020-03-16") & Date_week <= as.Date("2020-06-15"), p_strong_res, 0)) %>%
  mutate(p_strong_res = ifelse(Date_week >= as.Date("2020-03-16") & Date_week <= as.Date("2020-06-15"), 0, p_strong_res)) %>%
  mutate(
    periods = case_when(
      p_strong_res + p_mild_res + p_no_res + p_first == 0 ~ "pre-pandemic",
      p_first > 0.5 ~ "first wave", #  & Date_week <= as.Date("2020-04-13") 
      p_strong_res > 0.5 ~ "strong res", # | (p_first > 0.5 & Date_week > as.Date("2020-04-13")) 
      p_mild_res > 0.5 ~ "mild res",
      p_no_res > 0.5 ~ "low to no res",
      .default = NA
    ) 
  ) %>%
  dplyr::select(Date_week, periods) %>%
  mutate(periods = factor(periods, c("pre-pandemic", "first wave", "strong res", "mild res", "low to no res")))

# Correlation with dummy variables
df_int_icu = res %>%
  filter(secteur == "ICU") %>%
  group_by(bacterie, Date_week, secteur) %>%
  summarise(n_res = sum(n_res), .groups = "drop") %>%
  left_join(., pmsi_jh_icu %>% group_by(Date_week) %>% summarise(nbjh = sum(nbjh), .groups = "drop"), by = c("Date_week"))

df_int = res %>%
  group_by(bacterie, Date_week) %>%
  summarise(n_res = sum(n_res), .groups = "drop") %>%
  left_join(., pmsi_jh %>% group_by(Date_week) %>% summarise(nbjh = sum(nbjh), .groups = "drop"), by = c("Date_week")) %>%
  mutate(secteur = "Hospital") %>%
  bind_rows(., df_int_icu) %>%
  mutate(Date_week = as.Date(Date_week), res_i = n_res/nbjh*1000) %>%
  left_join(., int_national, by = "Date_week") 

df_anova = df_int %>%
  group_by(bacterie, secteur) %>%
  nest() %>%
  mutate(aov_res = map(data, function(x) kruskal.test(x$res_i ~ x$periods)$p.value)) %>%
  dplyr::select(-data) %>% 
  unnest(cols = aov_res) %>%
  mutate(y = ifelse(secteur == "ICU", 3, 0.35))

p2 = ggplot(data = df_int, aes(x = periods, y = res_i)) +
  geom_boxplot() +
  geom_text(data = df_anova, aes(x = 2.5, y = y, label = paste0("p=",round(aov_res,3)))) +
  facet_grid(cols = vars(bacterie), rows = vars(secteur), scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), axis.title.x = element_blank()) +
  labs(x = "", y = "Weekly incidence of resistant infections\n(for 1,000 bed-days)") +
  expand_limits(y=0)
p2

# Final correlation figure
p = ggarrange(p1, p2, nrow = 2, labels = c("A", "B"))
p
ggsave("../Paper/Supplementary/amr_covid_correlation.png", p, height = 9, width = 11)


# ##################################################
# # Other figures
# ##################################################
# # Add figures on the composite variables  
# df_full = bind_rows(df_corr_etab, df_corr_icu) %>%
#   mutate(Date_week = as.Date(Date_week)) %>%
#   left_join(., int_national, by = "Date_week") %>%
#   filter(!(p_strong_res==0 & p_mild_res==0 & p_no_res==0)) %>%
#   mutate(
#     p_strong_res = ifelse(p_strong_res == 0, -1, p_strong_res),
#     p_mild_res = ifelse(p_mild_res == 0, -1, p_mild_res),
#     p_no_res = ifelse(p_no_res == 0, -1, p_no_res)
#   ) %>%
#   mutate(
#     p_strong_res = p_strong_res*covid_prev,
#     p_mild_res = p_mild_res*covid_prev,
#     p_no_res = p_no_res*covid_prev
#   ) %>%
#   pivot_longer(c(p_strong_res, p_mild_res, p_no_res), names_to = "period", values_to = "prev") %>%
#   filter(prev >= 0) %>%
#   mutate(period = factor(case_when(period=="p_strong_res" ~ "Strong restrictions",
#                                    period=="p_mild_res" ~ "Mild restrictions",
#                                    period=="p_no_res" ~ "Low to no restrictions"
#   ),
#   c("Strong restrictions", "Mild restrictions", "Low to no restrictions")
#   ))
# 
# df_cor_test = df_full %>%
#   group_by(bacterie, secteur, period) %>%
#   summarise(
#     corr = paste0(
#       'tau = ',
#       round(cor.test(res_i, covid_prev, method = "kendall", exact=F)$estimate, 3), 
#       "\np = ",
#       round(cor.test(res_i, covid_prev, method = "kendall", exact=F)$p.value, 3)
#     ), 
#     .groups = "drop"
#   )
# 
# figureS5 = df_full %>%
#   filter(secteur == "Hospital") %>%
#   ggplot(., aes(x = prev, y = res_i)) +
#   geom_point(col = "grey70") +
#   geom_text(data = subset(df_cor_test, secteur == "Hospital"), aes(label = corr, x = 0.2, y = 0.3)) +
#   facet_grid(cols = vars(bacterie), rows = vars(period)) +
#   theme_bw() +
#   labs(x = "Weekly Covid-19 prevalence", y = "Weekly incidence of resistant infections")
# ggsave("../Paper/Supplementary/FigureS5.png", figureS5, height = 6, width = 10)
# 
# 
# figureS6 = df_full %>%
#   filter(secteur == "ICU") %>%
#   ggplot(., aes(x = prev, y = res_i)) +
#   geom_point(col = "grey70") +
#   geom_text(data = subset(df_cor_test, secteur == "ICU"), aes(label = corr, x = 0.55, y = 2.5)) +
#   facet_grid(cols = vars(bacterie), rows = vars(period)) +
#   theme_bw() +
#   labs(x = "Weekly Covid-19 prevalence", y = "Weekly incidence of resistant infections")
# ggsave("../Paper/Supplementary/FigureS6.png", figureS6, height = 6, width = 10)
# 
