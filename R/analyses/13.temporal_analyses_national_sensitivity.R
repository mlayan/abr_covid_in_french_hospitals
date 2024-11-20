##################################################
## COUNT REGRESSION ANALYSIS AT THE NATIONAL
## LEVEL IN HOSPITALS AND ICUs
##################################################
rm(list = ls())
library(tidyverse)
library(MASS)
library(performance)
library(DescTools)
library(qqplotr)
library(ggpubr)
library(cowplot)
library(gt)

# Helper functions and dictionaries
source("R/helper/helper_functions.R")
source("R/helper/helper_plots.R")
source("R/helper/dictionaries.R")

##################################################
# Load data 
##################################################
# Variables of intervention levels - national level
load("data/int_national.rda")
load("data/int_national_start_end.rda")

# Resistance data
load("data/res_hospital.rda")
load("data/res_icu.rda")

# Covid-19 intubation data
load("data/covid_intub_hospital.rda")
load("data/covid_intub_icu.rda")

# Bed-days data
load("data/bd_pmsi_hospital.rda")
load("data/bd_pmsi_icu.rda")

##################################################
# Merge databases
##################################################
# Merge national hospital data bases 
res_national = res_hospital %>%
  left_join(., bd_pmsi_hospital, by = "Date_week") %>%
  left_join(., covid_intub_hospital, by = "Date_week") %>%
  mutate(
    covid_intub = ifelse(is.na(covid_intub), 0, covid_intub),
    Date_year = lubridate::year(Date_week)
    ) %>%
  mutate(covid_intub_prev = covid_intub / nbjh * 1000) %>%
  left_join(., int_national, by = "Date_week")

sum(all_dates %in% res_national$Date_week)
length(all_dates)

# Merge national ICU databases
res_national_icu = res_icu %>%
  left_join(., bd_pmsi_icu, by = "Date_week") %>%
  left_join(., covid_intub_icu, by = "Date_week") %>%
  mutate(
    covid_intub = ifelse(is.na(covid_intub), 0, covid_intub),
    Date_year = lubridate::year(Date_week)
  ) %>%
  mutate(covid_intub_prev = covid_intub / nbjh * 1000) %>%
  left_join(., int_national, by = "Date_week")
  
sum(all_dates %in% res_national_icu$Date_week)
length(all_dates)

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

# Regression models 
for (db in c("icu", "hospital")) {
  for (mod in models) {
    for (i in seq_along(bacterias)) {
      
      # Create final database 
      if (db == "hospital") final_db = res_national
      if (db == "icu") final_db = res_national_icu
      
      final_db = final_db %>% 
        filter(bacterie == bacterias[i]) %>%
        arrange(Date_week) %>%
        mutate(
          lag1_i_res = lag(n_res / nbjh * 1000),
          
          lag1_covid_intub_prev = lag(covid_intub_prev),
          lag2_covid_intub_prev = lag(covid_intub_prev,2),
          
          lag1_periods = lag(periods, 1),
          lag2_periods = lag(periods, 2),

          nbjh = nbjh/1000
        ) %>% 
        filter(!is.na(lag2_periods)) %>%
        mutate(
          lag1_i_res = (lag1_i_res - mean(lag1_i_res)) / sd(lag1_i_res),

          covid_intub_prev = (covid_intub_prev - mean(covid_intub_prev)) / sd(covid_intub_prev),
          lag1_covid_intub_prev = (lag1_covid_intub_prev - mean(lag1_covid_intub_prev)) / sd(lag1_covid_intub_prev),
          lag2_covid_intub_prev = (lag2_covid_intub_prev - mean(lag2_covid_intub_prev)) / sd(lag2_covid_intub_prev),
          
          Date_year = as.character(Date_year)
        )
      
      # Multivariate model
      if (mod == "model0") {
        m_eq = "n_res ~ lag1_i_res + offset(log(nbjh))"
        if (bacterias[i] == "CR P. aeruginosa") m_eq = "n_res ~ lag1_i_res + offset(log(nbjh))"
        if (bacterias[i] == "MRSA") m_eq = "n_res ~ lag1_i_res + offset(log(nbjh))"
      } 
      
      if (mod == "model1") {
        m_eq = "n_res ~ lag1_i_res + periods  + offset(log(nbjh))"
        if (bacterias[i] == "CR P. aeruginosa") m_eq = "n_res ~ lag1_i_res + periods + offset(log(nbjh))"
        if (bacterias[i] == "MRSA") m_eq = "n_res ~ lag1_i_res + periods + offset(log(nbjh))"
      }
      
      if (mod == "model2") {
        m_eq = "n_res ~ lag1_i_res + lag1_periods + offset(log(nbjh))"
        if (bacterias[i] == "CR P. aeruginosa") m_eq = "n_res ~ lag1_i_res + lag1_periods + offset(log(nbjh))"
        if (bacterias[i] == "MRSA") m_eq = "n_res ~ lag1_i_res + lag1_periods + offset(log(nbjh))"
      }
      
      if (mod == "model3") {
        m_eq = "n_res ~ lag1_i_res + lag2_periods + offset(log(nbjh))"
        if (bacterias[i] == "CR P. aeruginosa") m_eq = "n_res ~ lag1_i_res + lag2_periods + offset(log(nbjh))"
        if (bacterias[i] == "MRSA") m_eq = "n_res ~ lag1_i_res + lag2_periods + offset(log(nbjh))"
      }
      
      if (mod == "model4") {
        m_eq = "n_res ~ lag1_i_res + covid_intub_prev + offset(log(nbjh))"
        if (bacterias[i] == "CR P. aeruginosa") m_eq = "n_res ~ lag1_i_res + covid_intub_prev + offset(log(nbjh))"
        if (bacterias[i] == "MRSA") m_eq = "n_res ~ lag1_i_res + covid_intub_prev + offset(log(nbjh))"
      }
      
      if (mod == "model5") {
        m_eq = "n_res ~ lag1_i_res + lag1_covid_intub_prev + offset(log(nbjh))"
        if (bacterias[i] == "CR P. aeruginosa") m_eq = "n_res ~ lag1_i_res + lag1_covid_intub_prev + offset(log(nbjh))"
        if (bacterias[i] == "MRSA") m_eq = "n_res ~ lag1_i_res + lag1_covid_intub_prev + offset(log(nbjh))"
      }
      
      if (mod == "model6") {
        m_eq = "n_res ~ lag1_i_res + lag2_covid_intub_prev + offset(log(nbjh))"
        if (bacterias[i] == "CR P. aeruginosa") m_eq = "n_res ~ lag1_i_res + lag2_covid_intub_prev + offset(log(nbjh))"
        if (bacterias[i] == "MRSA") m_eq = "n_res ~ lag1_i_res + lag2_covid_intub_prev + offset(log(nbjh))"
      }
      
      # Overdispersion in Poisson regression model
      results$poisson_OD[k] = check_overdispersion(glm(m_eq, data = final_db, family = poisson))$p_value
      
      # # Negative binomial regression model
      m = glm.nb(m_eq, data = final_db, link = log)
      
      results$negbin_OD[k] = check_overdispersion(m)$p_value
      results$theta[k] = m$theta 
      results$theta_min[k] = m$theta + qnorm(0.025) * m$SE.theta
      results$theta_max[k] = m$theta + qnorm(0.975) * m$SE.theta
      results$pseudoR2[k] = DescTools::PseudoR2(m)
      
      # Save the models
      all_models[[k]] = m
      
      # Colinearity of predictors
      if (mod != "model0") all_vif[[k]] = VIF(m)
      
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

best_models_sensitivity = data.frame(
  model = c("model2", "model0", "model0", "model1", "model1", "model6", "model1", "model0", "model1", "model1"),
  bacteria = rep(bacterias, 2), 
  setting = rep(c("hospital", "icu"), each = n_bacterias)
)
save(best_models_sensitivity, file = "data/best_models_sensitivity.rda")

# For CR-PA at the hospital level comparison with LRT
# of models 1 and 6
models = which(results$bacteria =="CR P. aeruginosa" & 
                 results$model %in% c("model1", "model2", "model6") &
                 results$setting == "hospital")
anova(all_models[[models[1]]], all_models[[models[2]]], all_models[[models[3]]])

# Table of AIC
aic_tab = results %>%
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
                     cells_body(rows = model == "Covid-19 intubation prevalence w-2" & setting == "ICU", columns = `CR P. aeruginosa`),
                     cells_body(rows = model == "Pandemic periods w-1" & setting == "Hospital", columns = `CR P. aeruginosa`)
    )
  ) %>%
  cols_align(
    align = "left",
    columns = model
  )
gtsave(aic_tab, filename = "../Paper/Supplementary/national_aic_sensitivity.docx")
gtsave(aic_tab, filename = "tables/national_aic_sensitivity.docx")

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
png("../Paper/Supplementary/national_acf_residuals_sensitivity.png", width = 15, height = 30, 
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

# Plot fits
all_fits %>%
  inner_join(., best_models, by = c("bacteria", "model", "setting")) %>%
  mutate(pred = exp(fit), 
         lw = exp(fit -1.96*se.fit), 
         ur = exp(fit+1.96*se.fit), 
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
    breaks = c("p_first_wave", "p_strong_res", "p_mild_res", "p_no_res"),
    values = c("p_first_wave" = col_interventions(1), "p_strong_res" = col_interventions(2), "p_mild_res" = col_interventions(3), "p_no_res" = col_interventions(4))
  ) +
  guides(fill=guide_legend(order=1, override.aes = list(col = 'black'))) +
  theme_bw() +
  theme(legend.title = element_text(hjust=0.5), legend.position = "bottom") +
  labs(x = "", y = "Weekly no. of resistant isolates")
ggsave("../Paper/Supplementary/national_fits_sensitivity.png", height = 8, width = 8)
ggsave("plots/regressions/national_fits_sensitivity.png", height = 8, width = 8)

# Estimates of the association with Covid-19 variables
covid_var_names = c(
  'periodsfirst wave' = "First wave w",
  'periodsstrong res' = "Strong w",
  'periodsmild res' = "Mild w",
  'periodslow to no res' = "Low to none w",
  
  'lag1_periodsfirst wave' = "First wave w-1",
  'lag1_periodsstrong res' = "Strong w-1",
  'lag1_periodsmild res' = "Mild w-1",
  'lag1_periodslow to no res' = "Low to none w-1",
  
  "lag2_covid_intub_prev" = "Covid-19 intubation\nprev. w-2"
)

covid_estimates = all_estimates %>%
  mutate(Estimate = exp(Estimate), q2_5 = exp(q2_5), q97_5 = exp(q97_5)) %>%
  inner_join(., best_models, by = c("bacteria", "model", "setting")) %>%
  filter(!variable %in% c("(Intercept)", "Penicillins", "TGC", "Carbapenems", "lag1_i_res")) %>%
  mutate(variable = factor(recode(variable, !!!covid_var_names), c("Covid-19 intubation\nprev. w-2", 
                                                                   "Low to none w", "Mild w", "Strong w", "First wave w", 
                                                                   "Low to none w-1", "Mild w-1", "Strong w-1", "First wave w-1")),
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
  filter(variable %in% c("TGC", "Penicillins", "Carbapenems", "lag1_i_res")) %>%
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
figure6 = plot_grid(
  covid_estimates + theme(legend.position = "none"),
  other_estimates + theme(legend.position = "none"), 
  ggpubr::get_legend(other_estimates),
  nrow = 3, 
  rel_heights = c(1, 0.5, 0.1), 
  labels = c("A", "B", "")
  )
figure6
ggsave("../Paper/Supplementary/national_estimates_sensitivity.png", figure6, height = 7, width = 7)
ggsave("plots/regressions/national_estimates_sensitivity.png", figure6, height = 7, width = 7)

# Value of theta 
overdispersion_tab = results %>%
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
  cols_align(align = "left", columns = c(bacteria, model)) 
gtsave(overdispersion_tab, "../Paper/Supplementary/national_overdispersion_sensitivity.docx")
gtsave(overdispersion_tab, "tables/national_overdispersion_sensitivity.docx")

# Ljung-Box autocorrelation test
ljung_box_tab = all_residuals %>%
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
  ) 
gtsave(ljung_box_tab, "../Paper/Supplementary/national_box_ljung_sensitivity.docx")
gtsave(ljung_box_tab, "tables/national_box_ljung_sensitivity.docx")



