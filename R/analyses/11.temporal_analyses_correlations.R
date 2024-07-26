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

load("data/bd_pmsi_icu.rda")
load("data/bd_pmsi_hospital.rda")

load("data/res_hospital.rda")
load("data/res_icu.rda")

##################################################
# Correlation between Covid-19 data and incidence 
# of resistant infections
##################################################
# Incidence at the hospital level
df_corr_etab = res_hospital %>%
  left_join(., covid_intub_hospital %>% mutate(Date_week = as.character(Date_week)), by = c("Date_year", "Date_week")) %>%
  left_join(., bd_pmsi_hospital, by = c("Date_year", "Date_week")) %>%
  filter(!is.na(covid_intub)) %>%
  mutate(covid_intub_prev = covid_intub / nbjh * 1000, 
         res_i = n_res / nbjh * 1000, 
         setting = "Hospital") %>%
  dplyr::select(bacterie, setting, Date_week, res_i, covid_intub_prev)

# Incidence at the ICU level
df_corr_icu = res_icu %>%
  left_join(., covid_intub_icu %>% mutate(Date_week = as.character(Date_week)), by = c("Date_year", "Date_week")) %>%
  left_join(., bd_pmsi_icu, by = c("Date_year", "Date_week")) %>%
  filter(!is.na(covid_intub)) %>%
  mutate(covid_intub_prev = covid_intub / nbjh * 1000, 
         res_i = n_res / nbjh * 1000, 
         setting = "ICU") %>%
  dplyr::select(bacterie, setting, Date_week, res_i, covid_intub_prev)

# Dataframes of Pearson correlation
df_corr = bind_rows(df_corr_etab, df_corr_icu) %>%
  group_by(setting, bacterie) %>%
  summarise(corr = paste0(
    'tau = ',
    round(cor.test(res_i, covid_intub_prev, method = "kendall")$estimate, 3), 
    "\np = ",
    round(cor.test(res_i, covid_intub_prev, method = "kendall")$p.value, 3)
  ), 
  .groups = "drop"
  ) %>%
  mutate(ymax = ifelse(setting == "Hospital", 0.35, 3.50),
         xmax = ifelse(setting == "Hospital", 20, 400))

# Covid-19 prevalence vs incidence of resistant infections
p1 = bind_rows(df_corr_etab, df_corr_icu) %>%
  ggplot(., aes(y = res_i, x = covid_intub_prev)) +
  geom_point(col = "grey70") +
  geom_text(data = df_corr, aes(x = xmax, y = ymax, label = corr)) +
  theme_bw() +
  ggh4x::facet_grid2(cols = vars(bacterie), rows = vars(setting), scales = "free", independent  ="x") +
  labs(x = "Weekly number of intubated Covid-19 bed-days (for 1,000 bed-days)", 
       y = "Weekly incidence of resistant infections\n(for 1,000 bed-days)")

# Correlation with dummy variables
df_int = bind_rows(
    res_icu %>% mutate(setting = "ICU") %>% left_join(., bd_pmsi_icu, by = c("Date_year", "Date_week")),
    res_hospital %>% mutate(setting = "Hospital") %>% left_join(., bd_pmsi_hospital, by = c("Date_year", "Date_week"))
  ) %>%
  mutate(Date_week = as.Date(Date_week), res_incidence = n_res/nbjh*1000) %>%
  left_join(., int_national, by = "Date_week") 

df_anova = df_int %>%
  group_by(bacterie, setting) %>%
  nest() %>%
  mutate(aov_res = map(data, function(x) kruskal.test(x$res_incidence ~ x$periods)$p.value)) %>%
  dplyr::select(-data) %>% 
  unnest(cols = aov_res) %>%
  mutate(y = ifelse(setting == "ICU", 4, 0.4))

p2 = ggplot(data = df_int, aes(x = periods, y = res_incidence)) +
  geom_boxplot() +
  geom_text(data = df_anova, aes(x = 4.5, y = y, label = paste0("p=",round(aov_res,3)))) +
  facet_grid(cols = vars(bacterie), rows = vars(setting), scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), axis.title.x = element_blank()) +
  labs(x = "", y = "Weekly incidence of resistant infections\n(for 1,000 bed-days)") +
  expand_limits(y=0)
p2

# Final correlation figure
p = ggarrange(p1, p2, nrow = 2, labels = c("A", "B"))
ggsave("../Paper/Supplementary/amr_covid_correlation.png", p, height = 9, width = 11)
ggsave("plots/antibiotic_resistance/amr_covid_correlation.png", p, height = 9, width = 11)

