rm(list = ls())

# Data management and visualisation
library(readxl)
library(tidyverse)
library(qqplotr)
library(geofacet)
library(ggpubr)
library(cowplot)
library(gt)
library(DescTools)
library(haven)

# Count regression models
library(performance)
library(MASS)

# Autoregressive models
library(forecast)

# Mixed effects models
library(lme4)
library(glmmTMB)
library(ggeffects)

# Helper functions
source("R/helper_functions.R")
source("R/helper_plots.R")

# Basic data and dictionaries 
load("data/metadata_admin_espic.rda")

# Variables of intervention levels - national level
load("data/int_national.rda")
load("data/int_national_start_end.rda")
load("data/cohort_final.rda")

# Model names
model_names = c(
  "model0" = "No Covid-19 variable",
  # "model0bis" = "Year effect", 
  "model1" = "Pandemic periods w",
  "model2" = "Pandemic periods w-1",
  "model3" = "Pandemic periods w-2",
  "model4" = "Covid-19 intubation prevalence w", 
  "model5" = "Covid-19 intubation prevalence w-1", 
  "model6" = "Covid-19 intubation prevalence w-2"#,
  # "model7" = "Prop. intubation w", 
  # "model8" = "Prop. intubation w-1", 
  # "model9" = "Prop. intubation w-2"
)

##################################################
# Load data 
##################################################
# Resistances
res = read.table("data-raw/spares/combined/resistance_cohortfinal.txt", header = T, sep = "\t")

# Antibiotic consumption
atb = read.table("data-raw/spares/combined/antibiotics_cohortfinal.txt", header = T, sep = "\t") %>%
  filter(secteur %in% c("Hospital", "ICU")) %>%
  group_by(Date_year, secteur) %>%
  summarise(
    Penicillins  = sum(Penicillins),
    Third_generation_Cephalosporins = sum(Third_generation_Cephalosporins),
    Carbapenems  = sum(Carbapenems),
    nbjh_spares  = sum(Nbhosp),
    .groups = "drop"
  )

# PMSI JH
pmsi_jh = read.table("data-raw/atih/pmsi_sejours.txt", header = T, sep = "\t") %>%
  group_by(Date_week) %>%
  summarise(nbjh = sum(nbjh), .groups = "drop") 

pmsi_jh_icu = read.table("data-raw/atih/pmsi_sejours.txt", header = T, sep = "\t") %>%
  filter(secteur == "Réanimation") %>%
  group_by(Date_week) %>%
  summarise(nbjh = sum(nbjh), .groups = "drop") 

# PMSI - COVID
pmsi_covid = read.table("data-raw/atih/pmsi_sejours_covid.txt", header = T, sep = "\t") %>%
  group_by(Date_week) %>%
  summarise(covid_jh  = sum(covid_jh), .groups = "drop") 

pmsi_covid_icu = read.table("data-raw/atih/pmsi_sejours_covid.txt", header = T, sep = "\t") %>%
  filter(secteur == "Réanimation") %>%
  group_by(Date_week) %>%
  summarise(covid_jh  = sum(covid_jh), .groups = "drop") 

# PMSI - Intubation
pmsi_intubation = read_sas("data-raw/atih/mco_intubation_weekly.sas7bdat") %>%
  rename(finess = FinessGeo) %>%
  left_join(., metadata_admin_espic %>% filter(code %in% cohort_final) %>% dplyr::select(finess, region) %>% distinct(), by = "finess") %>%
  filter(!is.na(region), status == "covid") %>%
  group_by(Date_week) %>%
  summarise(nintub = sum(nbjh), .groups = "drop")

pmsi_intubation_icu = read_sas("data-raw/atih/mco_intubation_weekly.sas7bdat") %>%
  rename(finess = FinessGeo) %>%
  left_join(., metadata_admin_espic %>% filter(code %in% cohort_final) %>% dplyr::select(finess, region) %>% distinct(), by = "finess") %>%
  filter(!is.na(region), status == "covid", secteur == "Réanimation") %>%
  group_by(Date_week) %>%
  summarise(nintub = sum(nbjh), .groups = "drop")

##################################################
# Merge databases
##################################################
# Create an additional variable for the 1st wave
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

# Merge national hospital data bases 
res_nat = res %>%
  group_by(Date_year, Date_week, bacterie) %>%
  summarise(n_res = sum(n_res), n_tot = sum(n_tot), .groups = "drop") %>%
  left_join(., atb %>% filter(secteur == "Hospital"), by = "Date_year") %>%
  left_join(., pmsi_jh, by = "Date_week") %>%
  left_join(., pmsi_covid, by = "Date_week") %>%
  mutate(Date_week = as.Date(Date_week)) %>%
  left_join(., pmsi_intubation, by = "Date_week") %>%
  mutate(
    covid_jh = ifelse(is.na(covid_jh), 0, covid_jh),
    nintub = ifelse(is.na(nintub), 0, nintub)
  ) %>%
  left_join(., int_national, by = "Date_week") %>%
  mutate(
    covid_intub_prop = ifelse(covid_jh == 0, 0, nintub/covid_jh),
    covid_jh = covid_jh / nbjh * 1000,
    Penicillins = Penicillins / nbjh_spares * 1000,
    Third_generation_Cephalosporins = Third_generation_Cephalosporins / nbjh_spares * 1000,
    Carbapenems = Carbapenems / nbjh_spares * 1000,
    nintub = nintub / nbjh * 1000 
  )

# National ICU databases
res_nat_icu = res %>%
  filter(secteur == "ICU") %>%
  group_by(Date_year, Date_week, bacterie) %>%
  summarise(n_res = sum(n_res), n_tot = sum(n_tot), .groups = "drop") %>%
  left_join(., atb %>% filter(secteur == "ICU"), by = "Date_year") %>%
  left_join(., pmsi_jh_icu, by = "Date_week") %>%
  left_join(., pmsi_covid_icu, by = "Date_week") %>%
  mutate(Date_week = as.Date(Date_week)) %>%
  left_join(., pmsi_intubation_icu, by = "Date_week") %>%
  mutate(
    covid_jh = ifelse(is.na(covid_jh), 0, covid_jh),
    nintub = ifelse(is.na(nintub), 0, nintub)
  ) %>%
  left_join(., int_national, by = "Date_week") %>%
  mutate(
    covid_intub_prop = ifelse(covid_jh == 0, 0, nintub/covid_jh),
    covid_jh = covid_jh / nbjh * 1000,
    Penicillins = Penicillins / nbjh_spares * 1000,
    Third_generation_Cephalosporins = Third_generation_Cephalosporins / nbjh_spares * 1000,
    Carbapenems = Carbapenems / nbjh_spares * 1000,
    nintub = nintub / nbjh * 1000 
  )

##################################################
# Correlation plots
##################################################
# Comparison of prevalence of Covid-19 patients and 
# prevalence of intubated Covid-19 in hospitals
p1 = res_nat %>%
  dplyr::select(Date_week, covid_jh, nintub) %>%
  distinct() %>%
  ggplot(., aes(x = Date_week)) +
  geom_line(aes(y = covid_jh)) +
  geom_line(aes(y = nintub), col = "red") +
  theme_bw()

p2 = res_nat %>%
  dplyr::select(Date_week, covid_intub_prop) %>%
  distinct() %>%
  ggplot(., aes(x = Date_week)) +
  geom_line(aes(y = covid_intub_prop)) +
  theme_bw()

ggarrange(p1, p2, ncol = 2)

# Covid-19 intubated patients in ICUs and hospitals 
bind_rows(
  res_nat %>%
    dplyr::select(Date_week, nintub, secteur) %>%
    distinct(),
  res_nat_icu %>%
    dplyr::select(Date_week, nintub, secteur) %>%
    distinct()
  ) %>%
  ggplot(., aes(x = Date_week, y = nintub, col = secteur)) +
  geom_line() +
  theme_bw()

##################################################
# Regression models
##################################################
# Bacterias and models
bacterias=c("CR P. aeruginosa", "ESBL K. pneumoniae", "ESBL E. cloacae", "ESBL E. coli", "MRSA")
n_bacterias = length(bacterias)

models = names(model_names)
n_models = length(models)

# Outputs
results = expand.grid(setting = c("icu", "hospital"), model = models, bacteria = bacterias) %>%
  mutate(
    poisson_OD = NA, 
    negbin_OD = NA, 
    theta = NA, 
    theta_min = NA,
    theta_max = NA,
    pseudoR2 = NA, 
    aic = NA
    ) %>%
  arrange(setting, model, bacteria)

all_estimates = data.frame()
all_residuals = data.frame()
all_fits = data.frame()
all_vif = vector("list", n_bacterias*n_models)
all_models = vector("list", n_bacterias*n_models)
k=1

# Test regression models 
for (db in c("icu", "hospital")) {
  for (mod in models) {
    for (i in seq_along(bacterias)) {
      
      # Create final database 
      if (db == "hospital") final_db = res_nat
      if (db == "icu") final_db = res_nat_icu
      
      final_db = final_db %>% 
        filter(bacterie == bacterias[i]) %>%
        arrange(Date_week) %>%
        mutate(
          lag1_i_res = lag(n_res / nbjh * 1000),
          lag1_prop_intub = lag(covid_intub_prop),
          lag2_prop_intub = lag(covid_intub_prop, 2),
          lag1_covid_jh = lag(covid_jh),
          lag2_covid_jh = lag(covid_jh,2),
          lag1_nintub = lag(nintub),
          lag2_nintub = lag(nintub,2),
          lag1_periods = lag(periods, 1),
          lag2_periods = lag(periods, 2),
          # lag1_p_first = lag(p_first),
          # lag2_p_first = lag(p_first, 2),
          # lag1_p_strong_res = lag(p_strong_res),
          # lag2_p_strong_res = lag(p_strong_res,2),
          # lag1_p_mild_res = lag(p_mild_res),
          # lag2_p_mild_res = lag(p_mild_res,2),
          # lag1_p_no_res = lag(p_no_res),
          # lag2_p_no_res = lag(p_no_res,2),
          nbjh = nbjh/1000
        ) %>% 
        filter(!is.na(lag2_covid_jh)) %>%
        mutate(
          Carbapenems = (Carbapenems - mean(Carbapenems)) / sd(Carbapenems),
          Penicillins = (Penicillins - mean(Penicillins)) / sd(Penicillins),
          Third_generation_Cephalosporins = (Third_generation_Cephalosporins - mean(Third_generation_Cephalosporins)) / sd(Third_generation_Cephalosporins),
          lag1_i_res = (lag1_i_res - mean(lag1_i_res)) / sd(lag1_i_res),

          covid_jh = (covid_jh - mean(covid_jh)) / sd(covid_jh),
          lag1_covid_jh = (lag1_covid_jh - mean(lag1_covid_jh)) / sd(lag1_covid_jh),
          lag2_covid_jh = (lag2_covid_jh - mean(lag2_covid_jh)) / sd(lag2_covid_jh),

          nintub = (nintub - mean(nintub)) / sd(nintub),
          lag1_nintub = (lag1_nintub - mean(lag1_nintub)) / sd(lag1_nintub),
          lag2_nintub = (lag2_nintub - mean(lag2_nintub)) / sd(lag2_nintub),

          covid_intub_prop  = (covid_intub_prop - mean(covid_intub_prop )) / sd(covid_intub_prop),
          lag1_prop_intub = (lag1_prop_intub - mean(lag1_prop_intub)) / sd(lag1_prop_intub),
          lag2_prop_intub = (lag2_prop_intub - mean(lag2_prop_intub)) / sd(lag2_prop_intub),
          
          Date_year = as.character(Date_year)
        )
      
      # Multivariate model
      if (mod == "model0") {
        m_eq = "n_res ~ lag1_i_res + Third_generation_Cephalosporins + offset(log(nbjh))"
        if (bacterias[i] == "CR P. aeruginosa") m_eq = "n_res ~ lag1_i_res + Carbapenems + offset(log(nbjh))"
        if (bacterias[i] == "MRSA") m_eq = "n_res ~ lag1_i_res + Penicillins + offset(log(nbjh))"
      } 
      
      if (mod == "model0bis") {
        m_eq = "n_res ~ lag1_i_res + Date_year + Third_generation_Cephalosporins + offset(log(nbjh))"
        if (bacterias[i] == "CR P. aeruginosa") m_eq = "n_res ~ lag1_i_res + Date_year + Carbapenems + offset(log(nbjh))"
        if (bacterias[i] == "MRSA") m_eq = "n_res ~ lag1_i_res + Date_year + Penicillins + offset(log(nbjh))"
      } 
      
      if (mod == "model1") {
        m_eq = "n_res ~ lag1_i_res + periods + Third_generation_Cephalosporins + offset(log(nbjh))"
        if (bacterias[i] == "CR P. aeruginosa") m_eq = "n_res ~ lag1_i_res + periods + Carbapenems + offset(log(nbjh))"
        if (bacterias[i] == "MRSA") m_eq = "n_res ~ lag1_i_res + periods + Penicillins + offset(log(nbjh))"
      }
      
      if (mod == "model2") {
        m_eq = "n_res ~ lag1_i_res + lag1_periods + Third_generation_Cephalosporins + offset(log(nbjh))"
        if (bacterias[i] == "CR P. aeruginosa") m_eq = "n_res ~ lag1_i_res + lag1_periods + Carbapenems + offset(log(nbjh))"
        if (bacterias[i] == "MRSA") m_eq = "n_res ~ lag1_i_res + lag1_periods + Penicillins + offset(log(nbjh))"
      }
      
      if (mod == "model3") {
        m_eq = "n_res ~ lag1_i_res + lag2_periods + Third_generation_Cephalosporins + offset(log(nbjh))"
        if (bacterias[i] == "CR P. aeruginosa") m_eq = "n_res ~ lag1_i_res + lag2_periods + Carbapenems + offset(log(nbjh))"
        if (bacterias[i] == "MRSA") m_eq = "n_res ~ lag1_i_res + lag2_periods + Penicillins + offset(log(nbjh))"
      }
      
      if (mod == "model4") {
        m_eq = "n_res ~ lag1_i_res + nintub + Third_generation_Cephalosporins + offset(log(nbjh))"
        if (bacterias[i] == "CR P. aeruginosa") m_eq = "n_res ~ lag1_i_res + nintub + Carbapenems + offset(log(nbjh))"
        if (bacterias[i] == "MRSA") m_eq = "n_res ~ lag1_i_res + nintub + Penicillins + offset(log(nbjh))"
      }
      
      if (mod == "model5") {
        m_eq = "n_res ~ lag1_i_res + lag1_nintub + Third_generation_Cephalosporins + offset(log(nbjh))"
        if (bacterias[i] == "CR P. aeruginosa") m_eq = "n_res ~ lag1_i_res + lag1_nintub + Carbapenems + offset(log(nbjh))"
        if (bacterias[i] == "MRSA") m_eq = "n_res ~ lag1_i_res + lag1_nintub + Penicillins + offset(log(nbjh))"
      }
      
      if (mod == "model6") {
        m_eq = "n_res ~ lag1_i_res + lag2_nintub + Third_generation_Cephalosporins + offset(log(nbjh))"
        if (bacterias[i] == "CR P. aeruginosa") m_eq = "n_res ~ lag1_i_res + lag2_nintub + Carbapenems + offset(log(nbjh))"
        if (bacterias[i] == "MRSA") m_eq = "n_res ~ lag1_i_res + lag2_nintub + Penicillins + offset(log(nbjh))"
      }
      
      results$poisson_OD[k] = check_overdispersion(glm(m_eq, data = final_db, family = poisson))$p_value
      
      # Negative binomial regression model
      m = glm.nb(m_eq, data = final_db, link = log)
      
      # m_eq = n_res ~ lag1_i_res + Penicillins + lag2_nintub + (1|Date_year) + offset(log(nbjh))
      # m = glmer.nb(m_eq, data=final_db, control=glmerControl(optimizer = "bobyqa"))
      # summary(m)
      # check_singularity(m)
      # check_overdispersion(m)
      # cbind(c(NA, exp(fixef(m))), exp(confint(m)))
      
      results$negbin_OD[k] = check_overdispersion(m)$p_value
      results$theta[k] = m$theta 
      results$theta_min[k] = m$theta + qnorm(0.025) * m$SE.theta
      results$theta_max[k] = m$theta + qnorm(0.975) * m$SE.theta
      results$pseudoR2[k] = DescTools::PseudoR2(m)
      
      # Save the models
      all_models[[k]] = m
      
      # Colinearity of predictors
      all_vif[[k]] = VIF(m)
      
      # AIC
      results$aic[k] = AIC(m)
      
      # Residuals
      mod_residuals = data.frame(
        setting = db,
        model = mod,
        bacteria = bacterias[i],
        residuals = residuals(m)
      )
      all_residuals = bind_rows(all_residuals, mod_residuals)
      
      # Estimates and p-values
      mod_estimates = data.frame(cbind(summary(m)$coefficients, confint(m))) %>%
        rownames_to_column(var = 'variable') %>%
        rename(q2_5 = `X2.5..`, q97_5 = `X97.5..`, p = `Pr...z..`, z = `z.value`, sd = `Std..Error`) %>%
        mutate(bacteria = bacterias[i], model = mod, setting = db)
      all_estimates = bind_rows(all_estimates, mod_estimates)
      
      # Model fit
      mod_fit = data.frame(
        model = mod,
        bacteria = bacterias[i],
        setting = db,
        Date_week = final_db$Date_week, 
        n_res = final_db$n_res, 
        nbjh = final_db$nbjh,
        fit = predict(m, type = "link", se.fit = T)$fit,
        se.fit = predict(m, type = "link", se.fit = T)$se.fit
      )
      all_fits = bind_rows(all_fits, mod_fit)
      
      k = k+1
    } 
  } 
}

# Best models
results %>%
  dplyr::select(setting, model, bacteria, aic) %>%
  arrange(desc(setting), bacteria) %>%
  mutate(aic = round(aic)) %>%
  pivot_wider(names_from = model, values_from = aic)

best_models = data.frame(
  model = c("model6", "model0", "model0", "model1", "model1", "model6", "model1", "model0", "model1", "model1"),
  bacteria = rep(bacterias, 2), 
  setting = rep(c("hospital", "icu"), each = n_bacterias)
)

# Pseudo R2
results %>% 
  filter(bacteria == "CR P. aeruginosa") %>%
  View()

# # VIF of model with all Covid-19 variables
# vif_tab = data.frame()
# for (b in bacterias) {
#   for (db in c("hospital", "icu")) {
#     i = which(results$setting == db & results$bacteria == b & results$model == "model6")
#     vif_mod_tab = all_estimates %>%
#       filter(model == "model6", bacteria == b, setting == db, variable != "(Intercept)") %>%
#       dplyr::select(model, bacteria, setting, variable) %>%
#       left_join(., data.frame(vif = all_vif[[i]]) %>% rownames_to_column(var = "variable"), by = "variable")
#     vif_tab = bind_rows(vif_tab, vif_mod_tab)
#   }
# }
# 
# vif_tab %>%
#   mutate(
#     setting = ifelse(setting == "icu", "ICU", "Hospital"),
#     variable = case_when(
#       variable == "lag_i_res" ~ "Incidence lag1",
#       variable == "covid_jh" ~ "Covid-19 prevalence",
#       variable == "p_no_res" ~ "Low to no restriction",
#       variable == "p_mild_res" ~ "Mild restrictions",
#       variable == "p_strong_res" ~ "Strong restrictions",
#       variable %in% c("Carbapenems", "Broad_Penicillins", "Narrow_Penicillins") ~ "Target antibiotic class"
#     ),
#     vif = round(vif, 1)
#   ) %>%
#   pivot_wider(names_from = bacteria, values_from = vif) %>%
#   group_by(setting) %>%
#   gt() %>%
#   cols_label(
#     variable = ""
#   ) %>%
#   tab_style(
#     style = list(cell_text(weight = "bold")),
#     locations = list(cells_row_groups(), 
#                      cells_column_labels())
#   ) %>%
#   gtsave("../Paper/Supplementary/national_vif.png")

# GT table of AIC
results %>%
  mutate(
    setting = ifelse(setting == "icu", "ICU", "Hospital"),
    aic = round(aic),
    model = recode(model, !!!model_names)
  ) %>%
  dplyr::select(setting, model, bacteria, aic) %>%
  pivot_wider(names_from = bacteria, values_from = aic) %>%
  arrange(setting) %>%
  group_by(setting) %>%
  gt() %>%
  cols_label(
    model = ""
  ) %>%
  tab_style(
    style = list(cell_text(weight = "bold")),
    locations = list(cells_row_groups(), 
                     cells_column_labels(), 
                     cells_body(rows = model == "No Covid-19 variable" & setting == "Hospital", columns = c(`ESBL K. pneumoniae`, `ESBL E. cloacae`)),
                     cells_body(rows = model == "No Covid-19 variable" & setting == "ICU", columns = `ESBL E. cloacae`),
                     cells_body(rows = model == "Pandemic periods w", columns = c(`ESBL E. coli`, `MRSA`)),
                     cells_body(rows = model == "Pandemic periods w" & setting == "ICU", columns = `ESBL K. pneumoniae`),
                     cells_body(rows = model == "Covid-19 intubation prevalence w-2", columns = `CR P. aeruginosa`)
    )
  ) %>%
  cols_align(
    align = "left",
    columns = model
  ) %>%
  gtsave(filename = "../Paper/Supplementary/national_aic.png")

# Verify that negative binomial models are better than poisson models
results %>%
  inner_join(., best_models, by = c("model", "setting", "bacteria")) %>%
  filter(poisson_OD > 0.05)

results %>%
  inner_join(., best_models, by = c("model", "setting", "bacteria")) %>%
  filter(negbin_OD <= 0.05)

results %>%
  inner_join(., best_models, by = c("model", "setting", "bacteria")) %>%
  ggplot(., aes(x = bacteria, y = theta, ymin = theta_min, ymax = theta_max, col = model)) +
  geom_hline(yintercept = 1) + 
  geom_pointrange() +
  facet_grid(cols = vars(setting)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  labs(x = "", y = "Dispersion parameter (95% Wald CI)", col = "")

# Verify normality of residuals
sw_p_df = all_residuals %>%
  inner_join(., best_models, by = c("model", "setting", "bacteria")) %>%
  group_by(setting, bacteria) %>%
  summarise(p = round(shapiro.test(residuals)$p.value, 4), .groups = "drop")
  
all_residuals %>%
  inner_join(., best_models, by = c("model", "setting", "bacteria")) %>%
  ggplot(., aes(sample = residuals)) +
  qqplotr::stat_qq_band(bandType = "boot", mapping = aes(fill = "Bootstrap"), alpha = 0.5, fill = "grey80") +
  stat_qq_line() +
  stat_qq() +
  geom_text(data = sw_p_df, aes(sample = NULL, label = p, x = 0, y = 2)) +
  facet_grid(cols = vars(setting), rows = vars(bacteria)) +
  theme_bw()

# Verify absence of autocorrelation
par(mfrow = c(n_bacterias,2))
for (i in 1:nrow(best_models)) {
  t = paste0(best_models$bacteria[i], " in ", best_models$setting[i])
  r = all_residuals %>%
    filter(
      bacteria == best_models$bacteria[i], 
      setting == best_models$setting[i],
      model == best_models$model[i]
      ) %>%
    .$residuals
  
  acf(r, main = t)
  pacf(r, main = t)
}

# Verify absence of multi-colinearity
best_models_vifs = all_vif[which(apply(results[, c("model", "bacteria", "setting")], 1, paste0, collapse = "_") %in% 
          apply(best_models, 1, paste0, collapse = "_"))]

names(best_models_vifs) = apply(results[, c("model", "bacteria", "setting")], 1, paste0, collapse = "_")[
  which(apply(results[, c("model", "bacteria", "setting")], 1, paste0, collapse = "_") %in% 
          apply(best_models, 1, paste0, collapse = "_"))]
best_models_vifs

##################################################
# Final plots
##################################################
# Autocorrelation of residuals
png("../Paper/Supplementary/national_acf_residuals.png", width = 15, height = 30, 
    res = 300, units = "cm")
par(mfrow = c(n_bacterias,2))
for (i in 1:nrow(best_models)) {
  t = paste0(best_models$bacteria[i], " in ", best_models$setting[i])
  r = all_residuals %>%
    filter(
      bacteria == best_models$bacteria[i], 
      setting == best_models$setting[i],
      model == best_models$model[i]
    ) %>%
    .$residuals
  
  acf(r, main = t)
}
dev.off()
  
# Normality of residuals 
sw_p_df = all_residuals %>%
  inner_join(., best_models, by = c("model", "setting", "bacteria")) %>%
  group_by(setting, bacteria) %>%
  summarise(p = round(shapiro.test(residuals)$p.value, 4), .groups = "drop") %>%
  mutate(setting = ifelse(setting == "hospital", "Hospital", "ICU"))

all_residuals %>%
  inner_join(., best_models, by = c("model", "setting", "bacteria")) %>%
  mutate(setting = ifelse(setting == "hospital", "Hospital", "ICU")) %>%
  ggplot(., aes(sample = residuals)) +
  stat_qq_band(bandType = "boot", mapping = aes(fill = "Bootstrap"), alpha = 0.5, fill = "grey60") +
  stat_qq_line() +
  stat_qq() +
  geom_text(data = sw_p_df, aes(sample = NULL, label = p, x = 0, y = 2)) +
  facet_grid(cols = vars(setting), rows = vars(bacteria)) +
  theme_bw()
ggsave("../Paper/Supplementary/national_normality_residuals.png", height = 8, width = 6)

# Plot fits
int_national_start_end$restrictions[int_national_start_end$start == as.Date("2020-03-17")] = "p_first"

all_fits %>%
  inner_join(., best_models, by = c("bacteria", "model", "setting")) %>%
  mutate(pred = exp(fit), #/nbjh*1000, 
         lw = exp(fit -1.96*se.fit), #/nbjh*1000,
         ur = exp(fit+1.96*se.fit), #/nbjh*1000,
         setting = ifelse(setting == "hospital", "Hospital", "ICU")
         ) %>%
  ggplot(., aes(x = Date_week)) +
  ggh4x::facet_grid2(cols = vars(setting), rows = vars(bacteria), scales = "free_y", independent  ="y") +
  expand_limits(y = 0) +
  geom_rect(data = int_national_start_end, 
            aes(NULL, NULL, xmin=start, xmax=end, fill=restrictions, ymin=-Inf, ymax=Inf)) +
  geom_line(aes(y = n_res, col = "Data")) +
  geom_ribbon(aes(ymin = lw, ymax = ur), fill = "red", alpha = 0.4) +
  geom_line(aes(y = pred, col = "Fit")) +
  scale_color_manual(name = "", values = c("Data" = "black", "Fit" = "red")) +
  scale_fill_manual(
    name = "Anti-Covid-19 interventions",
    labels = c("First wave", "Strong", "Intermediary", "Light to none"),
    breaks = c("p_first", "p_strong_res", "p_mild_res", "p_no_res"),
    values = c("p_first" = col_interventions(1), "p_strong_res" = col_interventions(2), "p_mild_res" = col_interventions(3), "p_no_res" = col_interventions(4))
  ) +
  guides(fill=guide_legend(order=1, override.aes = list(col = 'black'))) +
  theme_bw() +
  theme(legend.title = element_text(hjust=0.5), legend.position = "bottom") +
  labs(x = "", y = "Weekly no. of resistant acquisitions")
ggsave("../Paper/Supplementary/national_fits.png", height = 8, width = 8)


# Estimates of the association with Covid-19 variables
covid_var_names = c(
  # "p_first" = "First wave w",
  # "p_strong_res" = "Strong w",
  # "p_mild_res" = "Mild w",
  # "p_no_res" = "Low to none w",
  
  'periodsfirst wave' = "First wave w",
  'periodsstrong res' = "Strong w",
  'periodsmild res' = "Mild w",
  'periodslow to no res' = "Low to none w",
  
  "lag2_covid_jh" = "Covid-19 prev. w-2",
  "lag2_nintub" = "Covid-19 intubation\nprev. w-2"
)

covid_estimates = all_estimates %>%
  mutate(Estimate = exp(Estimate), q2_5 = exp(q2_5), q97_5 = exp(q97_5)) %>%
  inner_join(., best_models, by = c("bacteria", "model", "setting")) %>%
  filter(!variable %in% c("(Intercept)", "Penicillins", "Third_generation_Cephalosporins", "Carbapenems", "lag1_i_res")) %>%
  mutate(variable = factor(recode(variable, !!!covid_var_names), c("Covid-19 intubation\nprev. w-2", "Low to none w", "Mild w", "Strong w", "First wave w")),
         setting = ifelse(setting == "hospital", "Hospital", "ICU")) %>%
  ggplot(., aes(x = variable, y = Estimate, ymin = q2_5, ymax = q97_5, col = fct_rev(bacteria))) +
  geom_hline(yintercept = 1, linetype = "dashed", col = "grey50") +
  geom_pointrange(position = position_dodge(width = 0.5), fatten = 1) +
  coord_flip() +
  facet_grid(cols = vars(setting)) +
  expand_limits(y = c(0.9,1.2)) +
  scale_color_manual(values = c("MRSA" = "#CC79A7", "ESBL K. pneumoniae" = "#0072B2", 
                                "ESBL E. coli" = "#56B4E9", "ESBL E. cloacae" = "#E69F00", 
                                "CR P. aeruginosa" = "#F0E442")) +
  theme_classic() +
  theme(legend.position = "bottom",
        axis.text.y = element_text(hjust = 1),
        axis.line.y = element_blank(),
        axis.ticks.y= element_blank(),
        axis.title.y= element_blank()) +
  guides(col = guide_legend(reverse = T)) +
  labs(y = "Incidence Rate Ratio (95% CI)", col = "")
  
# Estimates of the association with antibiotic consumption and
# incidence of the preceding week
m = all_estimates %>%
  mutate(Estimate = exp(Estimate), q2_5 = exp(q2_5), q97_5 = exp(q97_5)) %>%
  inner_join(., best_models, by = c("bacteria", "model", "setting")) %>%
  filter(variable != "(Intercept)") %>%
  summarise(m = min(q2_5)) %>% 
  .$m
M = all_estimates %>%
  mutate(Estimate = exp(Estimate), q2_5 = exp(q2_5), q97_5 = exp(q97_5)) %>%
  inner_join(., best_models, by = c("bacteria", "model", "setting")) %>%
  filter(variable != "(Intercept)") %>%
  summarise(M = max(q97_5)) %>% 
  .$M

other_estimates = all_estimates %>%
  mutate(Estimate = exp(Estimate), q2_5 = exp(q2_5), q97_5 = exp(q97_5)) %>%
  inner_join(., best_models, by = c("bacteria", "model", "setting")) %>%
  filter(variable %in% c("Third_generation_Cephalosporins", "Penicillins", "Carbapenems", "lag1_i_res")) %>%
  mutate(variable = factor(ifelse(variable == "lag1_i_res", "Incidence w-1", "Target antibiotic"), c("Target antibiotic", "Incidence w-1")),
         setting = ifelse(setting == "hospital", "Hospital", "ICU")
  ) %>%
  ggplot(., aes(x = variable, y = Estimate, ymin = q2_5, ymax = q97_5, col = fct_rev(bacteria))) +
  geom_hline(yintercept = 1, linetype = "dashed", col = "grey50") +
  geom_pointrange(position = position_dodge(width = 0.6), fatten = 1) +
  coord_flip() +
  facet_grid(cols = vars(setting)) +
  expand_limits(y = c(0.9,1.2)) +
  scale_color_manual(values = c("MRSA" = "#CC79A7", "ESBL K. pneumoniae" = "#0072B2", 
                                "ESBL E. coli" = "#56B4E9", "ESBL E. cloacae" = "#E69F00", 
                                "CR P. aeruginosa" = "#F0E442")) +
  theme_classic() +
  theme(legend.position = "bottom",
        axis.text.y = element_text(hjust = 1),
        axis.line.y = element_blank(),
        axis.ticks.y= element_blank(),
        axis.title.y= element_blank()) +
  guides(col = guide_legend(reverse = T)) +
  labs(y = "Incidence Rate Ratio (95% CI)", col = "") +
  ylim(c(m, M))

# Final figure
figure5 = plot_grid(
  other_estimates + theme(legend.position = "none"), 
  covid_estimates + theme(legend.position = "none"),
  ggpubr::get_legend(other_estimates),
  nrow = 3, 
  rel_heights = c(0.8, 1, 0.1), 
  labels = c("A", "B", "")
  )
figure5
ggsave("../Paper/Figures/Figure5.png", figure5, height = 8, width = 7)

# Gt table of results
var_names = c(
  "(Intercept)" = "Intercept",
  "lag1_i_res" = "Incidence w-1",

  "periodsfirst wave" = "Pandemic periods",
  "periodsstrong res" = "Pandemic periods",
  "periodsmild res" = "Pandemic periods",
  "periodslow to no res" = "Pandemic periods",

  "lag2_nintub" = "Covid-19 intubation prev. w-2",
  
  "Penicillins" = "Penicillins",
  "Third_generation_Cephalosporins" = "3rd generation Cephalosporins",
  "Carbapenems" = "Imipenem + Meropenem"
)

all_estimates %>%
  inner_join(., best_models, by = c("model", "setting", "bacteria")) %>%
  mutate(
    estimate = paste0(round(exp(Estimate),2), " (", round(exp(q2_5),2), "-", round(exp(q97_5), 2), ")"),
    p = round(p,3),
    variable_precision = case_when(
      grepl("first", variable) ~ "First wave",
      grepl("strong", variable) ~ "Strong",
      grepl("mild", variable) ~ "Mild",
      grepl("no", variable) ~ "Low to none",
      .default = ""
    ),
    variable = factor(recode(variable, !!!var_names), unique(var_names)),
    setting = ifelse(setting == "icu", "ICU", "Hospital")
  ) %>%
  dplyr::select(setting, bacteria, variable, variable_precision, estimate, p) %>%
  pivot_wider(names_from = bacteria, values_from = c(estimate, p)) %>%
  arrange(setting, variable) %>%
  group_by(setting) %>%
  gt(.) %>%
  tab_spanner(
    label = "CR P. aeruginosa",
    columns = c(`estimate_CR P. aeruginosa`, `p_CR P. aeruginosa`)
  ) %>%
  tab_spanner(
    label = "ESBL E. coli",
    columns = c(`estimate_ESBL E. coli`, `p_ESBL E. coli`)
  ) %>%
  tab_spanner(
    label = "ESBL K. pneumoniae",
    columns = c(`estimate_ESBL K. pneumoniae`, `p_ESBL K. pneumoniae`)
  ) %>%
  tab_spanner(
    label = "ESBL E. cloacae",
    columns = c(`estimate_ESBL E. cloacae`, `p_ESBL E. cloacae`)
  ) %>%
  tab_spanner(
    label = "MRSA",
    columns = c(`estimate_MRSA`, `p_MRSA`)
  ) %>%
  cols_label(
    `estimate_CR P. aeruginosa` = "IRR",  
    `estimate_ESBL K. pneumoniae` = "IRR",
    `estimate_ESBL E. coli` = "IRR",
    `estimate_ESBL E. cloacae` = "IRR",
    estimate_MRSA = "IRR",
    `p_CR P. aeruginosa` = "p-value",
    `p_ESBL K. pneumoniae` = "p-value",
    `p_ESBL E. coli` = "p-value",
    `p_ESBL E. cloacae` = "p-value",
    p_MRSA = "p-value",
    variable = "",
    variable_precision = ""
  ) %>%
  cols_align(
    align = "left",
    columns = variable
  ) %>%
  tab_style(
    style = list(cell_text(weight = "bold")),
    locations = list(
      cells_body(rows = `p_CR P. aeruginosa` <= 0.05, columns = `estimate_CR P. aeruginosa`),
      cells_body(rows = `p_ESBL K. pneumoniae` <= 0.05, columns = `estimate_ESBL K. pneumoniae`),
      cells_body(rows = `p_ESBL E. coli` <= 0.05, columns = `estimate_ESBL E. coli`),
      cells_body(rows = `p_ESBL E. cloacae` <= 0.05, columns = `estimate_ESBL E. cloacae`),
      cells_body(rows = p_MRSA <= 0.05, columns = estimate_MRSA),
      cells_row_groups(),
      cells_column_labels(),
      cells_column_spanners()
    )
  ) %>%
  sub_missing(missing_text = "-") %>%
  tab_options(table.font.size = 11) %>%
  sub_values(
    columns = variable,
    rows = variable_precision %in% c("Strong", "Mild", "Low to none"),
    replacement = "",
    values = c("Pandemic periods")
  ) %>%
  gtsave("../Paper/Supplementary/national_estimates.png")

# Value of theta 
results %>%
  inner_join(., best_models, by = c("model", "setting", "bacteria")) %>%
  mutate(
    setting = ifelse(setting == "hospital", "Hospital", "ICU"),
    model = recode(model, !!!model_names), 
    theta = paste0(round(theta, 1), " (", round(theta_min), "-", round(theta_max), ")")
    ) %>%
  dplyr::select(setting, bacteria, model, theta) %>%
  arrange(setting, bacteria) %>%
  group_by(setting) %>%
  gt() %>%
  cols_label(
    bacteria = "",  
    model = "Best model",
    theta = "Overdispersion (95% CI)"
  ) %>%
  tab_style(
    style = list(cell_text(weight = "bold")),
    locations = list(cells_row_groups(),cells_column_labels())
  ) %>%
  sub_missing(missing_text = "-") %>%
  cols_align(align = "left", columns = c(bacteria, model)) %>%
  gtsave(., "../Paper/Supplementary/national_overdispersion.png")

# Fit comparison for model1 and model0 applied to CR P. aeruginosa in 
# ICUs
all_fits %>% 
  filter(bacteria == "CR P. aeruginosa", setting == "icu", model %in% c("model0", "model1")) %>%
  mutate(model = ifelse(model == "model0", "No Covid-19 variable", "Pandemic periods w")) %>%
  ggplot(., aes(x = Date_week)) +
  facet_grid(rows = vars(model)) +
  geom_rect(data = int_national_start_end, 
            aes(NULL, NULL, xmin=start, xmax=end, fill=restrictions, ymin=-Inf, ymax=Inf)) +
  geom_line(aes(y = n_res, col = "Data")) +
  geom_ribbon(aes(ymin = exp(fit + qnorm(0.025) * se.fit), ymax = exp(fit + qnorm(0.975) * se.fit)), fill = "red", alpha = 0.4) +
  geom_line(aes(y = exp(fit), col = "Fit")) + 
  scale_color_manual(name = "", values = c("Data" = "black", "Fit" = "red")) +
  scale_fill_manual(
    name = "Anti-Covid-19 interventions",
    labels = c("First wave", "Strong", "Intermediary", "Light to none"),
    breaks = c("p_first", "p_strong_res", "p_mild_res", "p_no_res"),
    values = c("p_first" = col_interventions(1), "p_strong_res" = col_interventions(2), "p_mild_res" = col_interventions(3), "p_no_res" = col_interventions(4))
  ) +
  guides(fill=guide_legend(order=1, override.aes = list(col = 'black'))) +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(x = "", y = "Weekly no. of incident ESBL K. pneumoniae in ICUs")
ggsave("../Paper/Supplementary/kpneumoniae_icus_fits.png", height = 6, width = 7)

# Ljung-Box autocorrelation test
all_residuals %>%
  inner_join(., best_models, by = c("model", "setting", "bacteria")) %>%
  group_by(bacteria, setting) %>%
  nest() %>%
  mutate(test = map(data, box_test)) %>%
  dplyr::select(-data) %>%
  unnest(test) %>%
  mutate(setting = ifelse(setting == "icu", "ICU", "Hospital")) %>%
  arrange(setting, bacteria) %>%
  group_by(setting) %>%
  gt() %>%
  cols_label(bacteria = "") %>%
  tab_style(
    style = list(cell_text(weight = "bold")),
    locations = list(cells_row_groups(),cells_column_labels())
  ) %>%
  gtsave(., "../Paper/Supplementary/national_box_ljung.png")

# Save best models
usethis::use_data(best_models, overwrite = TRUE)

##################################################
# Sensitivity analysis
# Multicolinearity between Covid-19 variables in 
# model 7
##################################################
# Bacterias
bacterias=c("CR P. aeruginosa", "ESBL K. pneumoniae", "ESBL E. cloacae", "ESBL E. coli", "MRSA")
n_bacterias = length(bacterias)

# Outputs
results_vif = expand.grid(setting = c("icu", "hospital"), bacteria = bacterias) %>%
  mutate(
    poisson_OD = NA, 
    negbin_OD = NA, 
    theta = NA, 
    theta_min = NA,
    theta_max = NA,
    aic = NA
  ) %>%
  arrange(setting, bacteria)
all_vif = data.frame()
k=1

# Test regression models 
for (db in c("icu", "hospital")) {
  for (i in seq_along(bacterias)) {
    
    # Create final database 
    if (db == "hospital") final_db = res_nat
    if (db == "icu") final_db = res_nat_icu
    
    final_db = final_db %>% 
      filter(bacterie == bacterias[i]) %>%
      arrange(Date_week) %>%
      mutate(
        lag1_i_res = lag(n_res / nbjh * 1000),
        lag1_prop_intub = lag(covid_intub_prop),
        lag2_prop_intub = lag(covid_intub_prop, 2),
        lag1_covid_jh = lag(covid_jh),
        lag2_covid_jh = lag(covid_jh,2),
        lag1_nintub = lag(nintub),
        lag2_nintub = lag(nintub,2),
        lag1_periods = lag(periods, 1),
        lag2_periods = lag(periods, 2),
        nbjh = nbjh/1000
      ) %>% 
      filter(!is.na(lag2_covid_jh)) %>%
      mutate(
        Carbapenems = (Carbapenems - mean(Carbapenems)) / sd(Carbapenems),
        Penicillins = (Penicillins - mean(Penicillins)) / sd(Penicillins),
        Third_generation_Cephalosporins = (Third_generation_Cephalosporins - mean(Third_generation_Cephalosporins)) / sd(Third_generation_Cephalosporins),
        lag1_i_res = (lag1_i_res - mean(lag1_i_res)) / sd(lag1_i_res),
        
        covid_jh = (covid_jh - mean(covid_jh)) / sd(covid_jh),
        lag1_covid_jh = (lag1_covid_jh - mean(lag1_covid_jh)) / sd(lag1_covid_jh),
        lag2_covid_jh = (lag2_covid_jh - mean(lag2_covid_jh)) / sd(lag2_covid_jh),
        
        nintub = (nintub - mean(nintub)) / sd(nintub),
        lag1_nintub = (lag1_nintub - mean(lag1_nintub)) / sd(lag1_nintub),
        lag2_nintub = (lag2_nintub - mean(lag2_nintub)) / sd(lag2_nintub),
        
        covid_intub_prop  = (covid_intub_prop - mean(covid_intub_prop )) / sd(covid_intub_prop),
        lag1_prop_intub = (lag1_prop_intub - mean(lag1_prop_intub)) / sd(lag1_prop_intub),
        lag2_prop_intub = (lag2_prop_intub - mean(lag2_prop_intub)) / sd(lag2_prop_intub)
      )
    
    # Multivariate model
    m_eq = "n_res ~ lag1_i_res + nintub + periods + Third_generation_Cephalosporins + offset(log(nbjh))"
    if (bacterias[i] == "CR P. aeruginosa") m_eq = "n_res ~ lag1_i_res + nintub + periods + Carbapenems + offset(log(nbjh))"
    if (bacterias[i] == "MRSA") m_eq = "n_res ~ lag1_i_res + nintub + periods + Penicillins + offset(log(nbjh))"
    
    results_vif$poisson_OD[k] = check_overdispersion(glm(m_eq, data = final_db, family = poisson))$p_value
    
    m = glm.nb(m_eq, data = final_db, link = log)
    results_vif$theta[k] = m$theta 
    results_vif$theta_min[k] = m$theta + qnorm(0.025) * m$SE.theta
    results_vif$theta_max[k] = m$theta + qnorm(0.975) * m$SE.theta
    results_vif$negbin_OD[k] = check_overdispersion(m)$p_value
    
    # Colinearity of predictors
    all_vif = bind_rows(
      all_vif, 
      data.frame(VIF(m)) %>%
        rownames_to_column(var = "variable") %>%
        rename(GVIF_df = GVIF..1..2.Df..) %>%
        mutate(bacteria = bacterias[i], setting = db)
    )
    
    # AIC
    results_vif$aic[k] = AIC(m)
    
    k = k+1
  } 
}

# AIC
results_vif %>%
  dplyr::select(setting, bacteria, aic) %>%
  arrange(desc(setting), bacteria) %>%
  mutate(aic = round(aic))

# VIFs
var_names2 = c(
  "lag1_i_res" = "Incidence w-1",
  "periods" = "Pandemic periods",
  "nintub" = "Covid-19 intubation prev. w",
  "Penicillins" = "Penicillins",
  "Third_generation_Cephalosporins" = "3rd generation Cephalosporins",
  "Carbapenems" = "Imipenem + Meropenem"
)

all_vif %>%
  mutate(
    setting = ifelse(setting == "icu", "ICU", "Hospital"),
    variable = recode(variable, !!!var_names2),
    GVIF = round(GVIF, 1)
    ) %>%
  dplyr::select(setting, variable, bacteria, GVIF) %>%
  pivot_wider(names_from = bacteria, values_from = GVIF) %>%
  arrange(setting) %>%
  group_by(setting) %>%
  gt() %>%
  cols_label(variable = "") %>%
  tab_style(
    style = list(cell_text(weight = "bold")),
    locations = list(cells_row_groups(), cells_column_labels(),
                     cells_body(rows = `CR P. aeruginosa` > 5, columns = `CR P. aeruginosa`),
                     cells_body(rows = `ESBL K. pneumoniae` > 5, columns = `ESBL K. pneumoniae`),
                     cells_body(rows = `ESBL E. cloacae` > 5, columns = `ESBL E. cloacae`),
                     cells_body(rows = `ESBL E. coli` > 5, columns = `ESBL E. coli`),
                     cells_body(rows = MRSA > 5, columns = MRSA)
                     )
  ) %>%
  sub_missing(missing_text = "-") %>%
  gtsave("../Paper/Supplementary/national_vif.png")

##################################################
# Comparison of regression models with Covid-19 
# variables of w+1 or w+2
##################################################
# Residuals of best models of w-2
all_aic = data.frame()
all_estimates_leads = data.frame()
all_fits_leads = data.frame()
models = paste0("model", 7:8)
b = "CR P. aeruginosa"

for (db in c("icu", "hospital")) {
  for (mod in models) {
    # Create final database 
    if (db == "hospital") final_db = res_nat
    if (db == "icu") final_db = res_nat_icu
    
    final_db = final_db %>% 
      filter(bacterie == b) %>%
      arrange(Date_week) %>%
      mutate(
        lag1_i_res = lag(n_res / nbjh * 1000, 2),
        lag2_nintub = lag(nintub, 2),
        lead1_nintub = lead(nintub, 1),
        lead2_nintub = lead(nintub,2),
        lead3_nintub = lead(nintub,3),
        nbjh = nbjh/1000
      ) %>% 
      filter(!is.na(lead3_nintub), !is.na(lag1_i_res)) %>%
      mutate(
        Carbapenems = (Carbapenems - mean(Carbapenems)) / sd(Carbapenems),
        lag1_i_res = (lag1_i_res - mean(lag1_i_res)) / sd(lag1_i_res),
        
        lag2_nintub = (lag2_nintub - mean(lag2_nintub)) / sd(lag2_nintub),
        lead1_nintub = (lead1_nintub - mean(lead1_nintub)) / sd(lead1_nintub),
        lead2_nintub = (lead2_nintub - mean(lead2_nintub)) / sd(lead2_nintub),
        lead3_nintub = (lead3_nintub - mean(lead3_nintub)) / sd(lead3_nintub)
      )
    
    # Multivariate model
    if (mod == "model7") {
      m_eq = "n_res ~ lag1_i_res + lead1_nintub + Carbapenems + offset(log(nbjh))"
    }
    
    if (mod == "model8") {
      m_eq = "n_res ~ lag1_i_res + lead2_nintub + Carbapenems + offset(log(nbjh))"
    }
    
    m1 = glm.nb(n_res ~ lag1_i_res + Carbapenems + offset(log(nbjh)), data = final_db, link = log)
    m2 = glm.nb(n_res ~ lag1_i_res + lead1_nintub + Carbapenems + offset(log(nbjh)), data = final_db, link = log)
    m3 = glm.nb(n_res ~ lag1_i_res + lead2_nintub + Carbapenems + offset(log(nbjh)), data = final_db, link = log)
    m4 = glm.nb(n_res ~ lag1_i_res + lead3_nintub + Carbapenems + offset(log(nbjh)), data = final_db, link = log)
    m5 = glm.nb(n_res ~ lag1_i_res + lag2_nintub + Carbapenems + offset(log(nbjh)), data = final_db, link = log)
    anova(m1, m2, m3, m4, m5)
    
    # AIC
    all_aic = bind_rows(all_aic, data.frame(
      bacteria = b,
      setting = db,
      model = mod,
      aic = AIC(m)
    )) 
    
    # Estimates
    mod_estimates = data.frame(cbind(summary(m)$coefficients, confint(m))) %>%
      rownames_to_column(var = 'variable') %>%
      rename(q2_5 = `X2.5..`, q97_5 = `X97.5..`, p = `Pr...z..`, z = `z.value`, sd = `Std..Error`) %>%
      mutate(bacteria = b, model = mod, setting = db)
    all_estimates_leads = bind_rows(all_estimates_leads, mod_estimates)
    
    # Fits
    mod_fit = data.frame(
      model = mod,
      setting = db,
      Date_week = final_db$Date_week, 
      n_res = final_db$n_res, 
      nbjh = final_db$nbjh,
      fit = predict(m, type = "link", se.fit = T)$fit,
      se.fit = predict(m, type = "link", se.fit = T)$se.fit
    )
    all_fits_leads = bind_rows(all_fits_leads, mod_fit)
  }
} 

# Comparison of AICs
p1=bind_rows(
  results %>%
    inner_join(., best_models, by = c("bacteria", "model", "setting")) %>%
    dplyr::select(model, bacteria, setting, aic) %>%
    mutate(model = "best_model"),
  all_aic
  ) %>%
  filter(bacteria == "CR P. aeruginosa") %>%
  dplyr::select(setting, model, aic) %>%
  mutate(
    aic = round(aic), 
    setting = ifelse(setting == "icu", "ICU", "Hospital"),
    model = recode(model, !!!c("best_model" = "Covid intubation w-2 (best model)", 
                               "model7" = "Covid intubation w+1",
                               "model8" = "Covid intubation w+2"))
    ) %>%
  arrange(setting) %>%
  rename(Setting = setting, Model = model, AIC = aic) %>%
  ggtexttable(., rows = NULL, theme = ttheme("light",base_size = 9))

# Estimates
p2=all_estimates_leads %>%
  filter(variable != "(Intercept)") %>%
  mutate(
    setting = ifelse(setting == "icu", "ICU", "Hospital"),
    model = ifelse(model == "model7", "Covid-19 intubation w+1", "Covid-19 intubation w+2"),
    variable = recode(variable, !!!c("lag1_i_res" = "Incidence w-1", "lead1_nintub" = "Covid-19 intubation", "lead2_nintub" = "Covid-19 intubation"))
  ) %>%
  ggplot(., aes(x = variable, y = exp(Estimate), ymin = exp(q2_5), ymax = exp(q97_5), col = setting)) +
  geom_hline(yintercept = 1) +
  geom_pointrange(position = position_dodge(width = 0.5)) +
  facet_grid(cols = vars(model)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  labs(x = "", y = "IRRs (95% CI)", col = "")
# ggsave("../Paper/Supplementary/crpa_leads_estimates.png", height = 4, width = 6)

# Fits
p3=all_fits %>%
  inner_join(., best_models, by = c("setting", "bacteria", "model")) %>%
  filter(bacteria == "CR P. aeruginosa") %>%
  bind_rows(., all_fits_leads) %>%
  mutate(
    setting = ifelse(setting == "icu", "ICU", "Hospital"),
    model = recode(model, !!!c("model6" = "Covid-19 intubation w-2\n(best model)", "model7" = "Covid-19 intubation w+1", "model8" = "Covid-19 intubation w+2"))
  ) %>%
  ggplot(., aes(x = Date_week)) +
  facet_grid(cols = vars(setting), rows = vars(model)) +
  geom_line(aes(y = n_res, col = "Data")) +
  geom_line(aes(y = exp(fit), col = "Fit")) +
  geom_ribbon(aes(ymin = exp(fit + qnorm(0.025) * se.fit), ymax = exp(fit + qnorm(0.975) * se.fit)), fill = "red", alpha = 0.4) +
  scale_color_manual(values = c("Data" = "black", "Fit" = "red")) +
  theme_bw() +
  labs(x = "", y = "Weekly no. of incident CR P. aeruginosa", col = "")
#ggsave("../Paper/Supplementary/crpa_leads_fits.png", height = 6, width = 10)

# Final Supplementary figure 
fig = plot_grid(
  plot_grid(
    p1, p2, ncol = 2, labels = c("A", "B"), rel_widths = c(0.6, 1)
    ),
  p3, nrow = 2, rel_heights = c(0.6,1), labels = c("", "C")
)
ggsave("../Paper/Supplementary/crpa_leads.png", fig, height = 8, width = 8)

