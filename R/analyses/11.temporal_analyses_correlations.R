##################################################
## PLOTS OF CORRELATION BETWEEN INCIDENCE OF 
## RESISTANT BACTERIA AND COVID-19 RELATED 
## VARIABLES 
##################################################
rm(list = ls())
library(tidyverse)
library(ggpubr)

source("R/helper/helper_functions.R")
source("R/helper/helper_plots.R")

load("data/int_national.rda")
load("data/covid_intub_icu.rda")
load("data/covid_intub_hospital.rda")
load("data/res_hospital.rda")
load("data/res_icu.rda")

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

