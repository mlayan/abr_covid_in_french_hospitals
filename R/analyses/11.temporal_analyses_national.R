##################################################
## COUNT REGRESSION ANALYSIS AT THE NATIONAL
## LEVEL IN HOSPITALS AND ICUs
##################################################
rm(list = ls())
library(tidyverse)
library(MASS)
library(performance)
library(DescTools)
library(ggpubr)
library(cowplot)
library(gt)

# Helper functions and dictionaries
source("R/helper/helper_functions.R")
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
load("data/res_not_icu.rda")

# Covid-19 data
load("data/covid_hospital.rda")
load("data/covid_icu.rda")
load("data/covid_not_icu.rda")

# Bed-days data
load("data/bd_pmsi_hospital.rda")
load("data/bd_pmsi_icu.rda")
load("data/bd_pmsi_not_icu.rda")

# Antibiotic consumption data
load("data/atb_use.rda")
atb_use = atb_use %>%
  filter(atb_class %in% c("Imipenem + Meropenem", "Penicillins", "Third generation Cephalosporins")) %>%
  mutate(consumption = molDDD/Nbhosp*1000, 
         atb_class = case_when(
           atb_class == "Imipenem + Meropenem" ~ "Carbapenems", 
           atb_class == "Third generation Cephalosporins" ~ "TGC",
           .default = atb_class
    )) %>%
  dplyr::select(Date_year, atb_class, consumption) %>%
  pivot_wider(names_from = atb_class, values_from = consumption)

load("data/atb_use_icu.rda")
atb_use_icu = atb_use_icu %>%
  filter(atb_class %in% c("Imipenem + Meropenem", "Penicillins", "Third generation Cephalosporins")) %>%
  mutate(consumption = molDDD/Nbhosp*1000, 
         atb_class = case_when(
           atb_class == "Imipenem + Meropenem" ~ "Carbapenems", 
           atb_class == "Third generation Cephalosporins" ~ "TGC",
           .default = atb_class
         )) %>%
  dplyr::select(Date_year, atb_class, consumption) %>%
  pivot_wider(names_from = atb_class, values_from = consumption)

load("data/atb_use_not_icu.rda")
atb_use_not_icu = atb_use_not_icu %>%
  filter(atb_class %in% c("Imipenem + Meropenem", "Penicillins", "Third generation Cephalosporins")) %>%
  mutate(consumption = molDDD/Nbhosp*1000, 
         atb_class = case_when(
           atb_class == "Imipenem + Meropenem" ~ "Carbapenems", 
           atb_class == "Third generation Cephalosporins" ~ "TGC",
           .default = atb_class
         )) %>%
  dplyr::select(Date_year, atb_class, consumption) %>%
  pivot_wider(names_from = atb_class, values_from = consumption)

##################################################
# Merge databases
##################################################
# Merge national hospital data bases 
res_national = res_hospital %>%
  left_join(., bd_pmsi_hospital, by = "Date_week") %>%
  left_join(., covid_hospital, by = "Date_week") %>%
  mutate(
    covid = ifelse(is.na(covid), 0, covid),
    Date_year = lubridate::year(Date_week)
    ) %>%
  left_join(., atb_use, by = "Date_year") %>%
  mutate(covid_prev = covid / nbjh * 1000) %>%
  left_join(., int_national, by = "Date_week")

sum(all_dates %in% res_national$Date_week)
length(all_dates)

# Merge national ICU databases
res_national_icu = res_icu %>%
  left_join(., bd_pmsi_icu, by = "Date_week") %>%
  left_join(., covid_icu, by = "Date_week") %>%
  mutate(
    covid = ifelse(is.na(covid), 0, covid),
    Date_year = lubridate::year(Date_week)
  ) %>%
  left_join(., atb_use_icu, by = "Date_year") %>%
  mutate(covid_prev = covid / nbjh * 1000) %>%
  left_join(., int_national, by = "Date_week")
  
sum(all_dates %in% res_national_icu$Date_week)
length(all_dates)

# Merge national not ICU databases
res_national_not_icu = res_not_icu %>%
  left_join(., bd_pmsi_not_icu, by = "Date_week") %>%
  left_join(., covid_not_icu, by = "Date_week") %>%
  mutate(
    covid = ifelse(is.na(covid), 0, covid),
    Date_year = lubridate::year(Date_week)
  ) %>%
  left_join(., atb_use_not_icu, by = "Date_year") %>%
  mutate(covid_prev = covid / nbjh * 1000) %>%
  left_join(., int_national, by = "Date_week")

sum(all_dates %in% res_national_not_icu$Date_week)
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
    pseudoR2_mcfadden = NA, 
    pseudoR2_efron = NA, 
    pseudoR2_veall = NA, 
    aic = NA,
    ll_null = NA,
    ll = NA
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
          
          lag1_covid_prev = lag(covid_prev),
          lag2_covid_prev = lag(covid_prev,2),
          
          lag1_periods = lag(periods, 1),
          lag2_periods = lag(periods, 2),

          nbjh = nbjh/1000
        ) %>% 
        filter(!is.na(lag2_periods)) %>%
        mutate(
          Carbapenems = (Carbapenems - mean(Carbapenems)) / sd(Carbapenems),
          Penicillins = (Penicillins - mean(Penicillins)) / sd(Penicillins),
          TGC = (TGC - mean(TGC)) / sd(TGC),
          
          lag1_i_res = (lag1_i_res - mean(lag1_i_res)) / sd(lag1_i_res),

          covid_prev = (covid_prev - mean(covid_prev)) / sd(covid_prev),
          lag1_covid_prev = (lag1_covid_prev - mean(lag1_covid_prev)) / sd(lag1_covid_prev),
          lag2_covid_prev = (lag2_covid_prev - mean(lag2_covid_prev)) / sd(lag2_covid_prev),
          
          Date_year = as.character(Date_year)
        )
      
      # Multivariate model
      if (mod == "model0") {
        m_eq = "n_res ~ lag1_i_res + TGC + offset(log(nbjh))"
        if (bacterias[i] == "CR P. aeruginosa") m_eq = "n_res ~ lag1_i_res + Carbapenems + offset(log(nbjh))"
        if (bacterias[i] == "MRSA") m_eq = "n_res ~ lag1_i_res + Penicillins + offset(log(nbjh))"
      } 
      
      if (mod == "model1") {
        m_eq = "n_res ~ lag1_i_res + periods + TGC + offset(log(nbjh))"
        if (bacterias[i] == "CR P. aeruginosa") m_eq = "n_res ~ lag1_i_res + periods + Carbapenems + offset(log(nbjh))"
        if (bacterias[i] == "MRSA") m_eq = "n_res ~ lag1_i_res + periods + Penicillins + offset(log(nbjh))"
      }
      
      if (mod == "model2") {
        m_eq = "n_res ~ lag1_i_res + lag1_periods + TGC + offset(log(nbjh))"
        if (bacterias[i] == "CR P. aeruginosa") m_eq = "n_res ~ lag1_i_res + lag1_periods + Carbapenems + offset(log(nbjh))"
        if (bacterias[i] == "MRSA") m_eq = "n_res ~ lag1_i_res + lag1_periods + Penicillins + offset(log(nbjh))"
      }
      
      if (mod == "model3") {
        m_eq = "n_res ~ lag1_i_res + lag2_periods + TGC + offset(log(nbjh))"
        if (bacterias[i] == "CR P. aeruginosa") m_eq = "n_res ~ lag1_i_res + lag2_periods + Carbapenems + offset(log(nbjh))"
        if (bacterias[i] == "MRSA") m_eq = "n_res ~ lag1_i_res + lag2_periods + Penicillins + offset(log(nbjh))"
      }
      
      if (mod == "model4") {
        m_eq = "n_res ~ lag1_i_res + covid_prev + TGC + offset(log(nbjh))"
        if (bacterias[i] == "CR P. aeruginosa") m_eq = "n_res ~ lag1_i_res + covid_prev + Carbapenems + offset(log(nbjh))"
        if (bacterias[i] == "MRSA") m_eq = "n_res ~ lag1_i_res + covid_prev + Penicillins + offset(log(nbjh))"
      }
      
      if (mod == "model5") {
        m_eq = "n_res ~ lag1_i_res + lag1_covid_prev + TGC + offset(log(nbjh))"
        if (bacterias[i] == "CR P. aeruginosa") m_eq = "n_res ~ lag1_i_res + lag1_covid_prev + Carbapenems + offset(log(nbjh))"
        if (bacterias[i] == "MRSA") m_eq = "n_res ~ lag1_i_res + lag1_covid_prev + Penicillins + offset(log(nbjh))"
      }
      
      if (mod == "model6") {
        m_eq = "n_res ~ lag1_i_res + lag2_covid_prev + TGC + offset(log(nbjh))"
        if (bacterias[i] == "CR P. aeruginosa") m_eq = "n_res ~ lag1_i_res + lag2_covid_prev + Carbapenems + offset(log(nbjh))"
        if (bacterias[i] == "MRSA") m_eq = "n_res ~ lag1_i_res + lag2_covid_prev + Penicillins + offset(log(nbjh))"
      }
      
      # Overdispersion in Poisson regression model
      results$poisson_OD[k] = check_overdispersion(glm(m_eq, data = final_db, family = poisson))$p_value
      
      # Loglikelihood of the null model
      m_null = glm.nb(n_res ~ 1 + offset(log(nbjh)), data = final_db, link = log)
      results$ll_null[k] = m_null$twologlik
      
      # # Negative binomial regression model
      m = glm.nb(m_eq, data = final_db, link = log)
      
      results$ll[k] = m$twologlik
      results$negbin_OD[k] = check_overdispersion(m)$p_value
      results$theta[k] = m$theta 
      results$theta_min[k] = m$theta + qnorm(0.025) * m$SE.theta
      results$theta_max[k] = m$theta + qnorm(0.975) * m$SE.theta
      results$pseudoR2_mcfadden[k] = DescTools::PseudoR2(m, which = "McFadden")
      results$pseudoR2_efron[k] = DescTools::PseudoR2(m, which = "Efron")
      results$pseudoR2_veall[k] = DescTools::PseudoR2(m, which = "VeallZimmermann")
      
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

##################################################
# Model comparison and selection
##################################################
# Select and save best models
results %>%
  dplyr::select(setting, model, bacteria, aic) %>%
  arrange(desc(setting), bacteria) %>%
  mutate(aic = round(aic)) %>%
  pivot_wider(names_from = model, values_from = aic)

best_models = data.frame(
  model = c("model6", "model0", "model0", "model1", "model1", "model6", "model2", "model0", "model1", "model1"),
  bacteria = rep(bacterias, 2), 
  setting = rep(c("hospital", "icu"), each = n_bacterias)
)
save(best_models, file = "data/best_models.rda")

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
                     cells_body(rows = model == "Pandemic periods w-1" & setting == "ICU", columns = `ESBL K. pneumoniae`),
                     cells_body(rows = model == "Covid-19 prevalence w-2", columns = `CR P. aeruginosa`)
    )
  ) %>%
  cols_align(
    align = "left",
    columns = model
  )
aic_tab
gtsave(aic_tab, filename = "../Paper/Supplementary/national_aic.docx")
gtsave(aic_tab, filename = "tables/national_aic.docx")

##################################################
# Final plots and tables with best model
##################################################
# Table with regression results (regression coefficients, standard errors, z-values)
results_tab = left_join(best_models, all_estimates, by = c("model", "bacteria", "setting")) %>%
  mutate(
    setting = ifelse(setting == "icu", "ICU", "Hospital"), 
    model = recode(model, !!!model_names),
    variable = recode(variable, !!!c(covid_var_names, "(Intercept)" = "Intercept", "lag1_i_res" = "Incidence w-1", "TGC" = "3GC")),
    Estimate_ex = paste0(sprintf("%.3f", round(Estimate, 2)), " (", sprintf("%.3f", round(q2_5, 3)), ", ", sprintf("%.3f", round(q97_5, 3)), ")"),
    sd = sprintf("%.3f", round(sd, 3)),
    z = sprintf("%.3f", round(z, 1))
  ) %>%
  dplyr::select(setting, bacteria, model, variable, Estimate_ex, sd, z) %>%
  group_by(setting) %>%
  gt() %>%
  cols_label(
    bacteria = "", 
    model = "Best model",
    variable = "Variable",
    Estimate_ex = 'Estimate (95% CI)',
    sd = "Standard deviation",
    z = "Z value"
  ) %>%
  tab_style(
    style = list(cell_text(weight = "bold")),
    locations = list(
      cells_row_groups(),
      cells_column_labels(),
      cells_column_spanners()
    )
  ) %>% 
  tab_options(table.font.size = 11)
results_tab
gtsave(results_tab, "../Paper/Supplementary/national_regression_results.docx")
gtsave(results_tab, "tables/national_regression_results.docx")

# Table of estimates
all_estimates_tab = all_estimates %>%
  inner_join(., best_models, by = c("model", "setting", "bacteria")) %>%
  mutate(
    estimate = paste0(round(exp(Estimate),2), " (", round(exp(q2_5),2), ", ", round(exp(q97_5), 2), ")"),
    p = ifelse(p < 0.0001, "<0.0001", format(round(p,4), nsmall = 4)),
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
    `estimate_CR P. aeruginosa` = "IRR (95% CI)",  
    `estimate_ESBL K. pneumoniae` = "IRR (95% CI)",
    `estimate_ESBL E. coli` = "IRR (95% CI)",
    `estimate_ESBL E. cloacae` = "IRR (95% CI)",
    estimate_MRSA = "IRR (95% CI)",
    `p_CR P. aeruginosa` = "p",
    `p_ESBL K. pneumoniae` = "p",
    `p_ESBL E. coli` = "p",
    `p_ESBL E. cloacae` = "p",
    p_MRSA = "p",
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
  )
gtsave(all_estimates_tab, "../Paper/Supplementary/national_estimates.docx")
gtsave(all_estimates_tab, "tables/national_estimates.docx")

# Plot of the estimates of the association with Covid-19 variables
covid_estimates = all_estimates %>%
  mutate(Estimate = exp(Estimate), q2_5 = exp(q2_5), q97_5 = exp(q97_5)) %>%
  inner_join(., best_models, by = c("bacteria", "model", "setting")) %>%
  filter(!variable %in% c("(Intercept)", "Penicillins", "TGC", "Carbapenems", "lag1_i_res")) %>%
  mutate(variable = factor(recode(variable, !!!covid_var_names), c("Covid-19 prev. w-2", 
                                                                   "Low to none w-1", "Mild w-1", "Strong w-1", "First wave w-1",
                                                                   "Low to none w", "Mild w", "Strong w", "First wave w")),
         setting = ifelse(setting == "hospital", "Hospital", "ICU")) %>%
  ggplot(., aes(x = variable, y = Estimate, ymin = q2_5, ymax = q97_5, col = fct_rev(bacteria))) +
  geom_hline(yintercept = 1, linetype = "dashed", col = "grey50") +
  geom_vline(aes(xintercept = 5.5), linewidth = 0.1) + 
  geom_vline(aes(xintercept = 1.5), linewidth = 0.1) +
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

# Save panels to make figure in script on regional data
save(covid_estimates, other_estimates, file = "plots/figure_composite.rda")

# Metrics for the abstract
all_estimates %>% 
  filter(bacteria == "CR P. aeruginosa", model == "model6", variable == "lag2_covid_prev") %>%
  dplyr::select(bacteria, model, variable, setting, p)

all_estimates %>%
  filter(bacteria %in% c("MRSA", "ESBL E. coli"), model == "model1", grepl("periods", variable),
         exp(q97_5) < 1) %>%
  dplyr::select(bacteria, model, variable, setting, Estimate, q2_5, q97_5) %>%
  arrange(bacteria, setting) %>%
  mutate(Estimate = round(1-exp(Estimate),2), q2_5 = round(1-exp(q2_5),2), 
         q97_5 = round(1-exp(q97_5), 2))

##################################################
# Model validity 
##################################################
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
ggsave("../Paper/Supplementary/national_fits.png", height = 8, width = 8)
ggsave("plots/regressions/national_fits.png", height = 8, width = 8)

# Get pseudo-R2 from McFadden, Efron, and Veall & Zimmermann 
r_square_tab = left_join(best_models, results, by =c("model", "bacteria", "setting")) %>%
  mutate(
    setting = ifelse(setting == "icu", "ICU", "Hospital"), 
    model = recode(model, !!!model_names),
    pseudoR2_mcfadden = sprintf("%.3f", round(pseudoR2_mcfadden, 3)),
    pseudoR2_efron = sprintf("%.3f", round(pseudoR2_efron, 3)),
    pseudoR2_veall = sprintf("%.3f", round(pseudoR2_veall, 3))
    ) %>%
  dplyr::select(setting, bacteria, model, contains("pseudoR2")) %>%
  group_by(setting) %>%
  gt(.) %>%
  fmt_markdown() %>%
  cols_label(
    bacteria = "", 
    model = "Best model",
    pseudoR2_mcfadden = md("McFadden's $R^{2}$"),
    pseudoR2_efron = md("Efron's $R^{2}$"),
    pseudoR2_veall = md("Veall & Zimmermann's $R^{2}$")
  ) %>%
  tab_style(
    style = list(cell_text(weight = "bold")),
    locations = list(
      cells_row_groups(),
      cells_column_labels(),
      cells_column_spanners()
    )
  ) %>%
  cols_align(
    align = "center",
    columns = contains("pseudoR2")
  ) %>%
  tab_options(table.font.size = 11)
gtsave(r_square_tab, "../Paper/Supplementary/national_pseudoR2.docx")
gtsave(r_square_tab, "tables/national_pseudoR2.docx")

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

# Value of theta and overdispersion test
overdispersion_tab = results %>%
  inner_join(., best_models, by = c("model", "setting", "bacteria")) %>%
  mutate(
    setting = ifelse(setting == "hospital", "Hospital", "ICU"),
    model = recode(model, !!!model_names), 
    theta = paste0(round(theta), " (", round(theta_min), ", ", round(theta_max), ")"),
    poisson_OD = ifelse(poisson_OD<0.0001, "<0.0001", format(round(poisson_OD, 4), nsmall = 4))
    ) %>%
  dplyr::select(setting, bacteria, model, poisson_OD, theta) %>%
  arrange(setting, bacteria) %>%
  group_by(setting) %>%
  gt() %>%
  cols_label(
    bacteria = "",  
    model = "Best model",
    poisson_OD = "Overdispersion test",
    theta = "Overdispersion (95% CI)"
  ) %>%
  tab_style(
    style = list(cell_text(weight = "bold")),
    locations = list(cells_row_groups(),cells_column_labels())
  ) %>%
  sub_missing(missing_text = "-") %>%
  cols_align(align = "left", columns = c(bacteria, model))
overdispersion_tab
gtsave(overdispersion_tab, "../Paper/Supplementary/national_overdispersion.docx")
gtsave(overdispersion_tab, "tables/national_overdispersion.docx")

# Ljung-Box autocorrelation test
ljung_box_tab = all_residuals %>%
  inner_join(., best_models, by = c("model", "setting", "bacteria")) %>%
  group_by(bacteria, setting) %>%
  nest() %>%
  mutate(test = map(data, box_test)) %>%
  dplyr::select(-data) %>%
  unnest(test) %>%
  mutate(
    setting = ifelse(setting == "icu", "ICU", "Hospital"),
    p.value = ifelse(p.value < 0.0001, "<0.0001", format(round(p.value, 4), nsmall = 4))
    ) %>%
  arrange(setting, bacteria) %>%
  group_by(setting) %>%
  gt() %>%
  cols_label(bacteria = "") %>%
  tab_style(
    style = list(cell_text(weight = "bold")),
    locations = list(cells_row_groups(),cells_column_labels())
  ) 
ljung_box_tab
gtsave(ljung_box_tab, "../Paper/Supplementary/national_box_ljung.docx")
gtsave(ljung_box_tab, "tables/national_box_ljung.docx")

##################################################
# Table with unadjusted estimates
##################################################
# Outputs
univariate_results = data.frame()
for (i in 1:nrow(best_models)) {

  out = data.frame()
  b = best_models$bacteria[i]
  s = best_models$setting[i]

  if (s == "hospital") final_db = res_national
  if (s == "icu") final_db = res_national_icu

  final_db = final_db %>%
    filter(bacterie == b) %>%
    arrange(Date_week) %>%
    mutate(
      lag1_i_res = lag(n_res / nbjh * 1000),

      lag1_covid_prev = lag(covid_prev),
      lag2_covid_prev = lag(covid_prev,2),

      lag1_periods = lag(periods, 1),
      lag2_periods = lag(periods, 2),

      nbjh = nbjh/1000
    ) %>%
    filter(!is.na(lag2_periods)) %>%
    mutate(
      Carbapenems = (Carbapenems - mean(Carbapenems)) / sd(Carbapenems),
      Penicillins = (Penicillins - mean(Penicillins)) / sd(Penicillins),
      TGC = (TGC - mean(TGC)) / sd(TGC),

      lag1_i_res = (lag1_i_res - mean(lag1_i_res)) / sd(lag1_i_res),

      covid_prev = (covid_prev - mean(covid_prev)) / sd(covid_prev),
      lag1_covid_prev = (lag1_covid_prev - mean(lag1_covid_prev)) / sd(lag1_covid_prev),
      lag2_covid_prev = (lag2_covid_prev - mean(lag2_covid_prev)) / sd(lag2_covid_prev),

      Date_year = as.character(Date_year)
    )

  # Univariate negative binomial regression model
  for (v in c("periods", "lag1_periods", "lag2_periods")) {
    m_eq = paste0("n_res ~ ", v, "+ offset(log(nbjh))")
    m_out = glm.nb(m_eq, data = final_db, link = log)
    m0 = glm.nb(n_res ~ offset(log(nbjh)), data = final_db, link = log)

    out_temp = data.frame(
      var_name = rownames(confint(m_out))[-1],
      p = summary(m_out)$coefficients[-1, "Pr(>|z|)"],
      estimate = exp(summary(m_out)$coefficients[-1, "Estimate"]),
      estimate_lw = exp(confint(m_out)[-1, "2.5 %"]),
      estimate_up = exp(confint(m_out)[-1, "97.5 %"])
    )
    out_temp2 = data.frame(
      var_name = v,
      p = anova(m_out, m0)$`Pr(Chi)`[2],
      estimate = NA,
      estimate_lw = NA,
      estimate_up = NA
    )

    out = bind_rows(out, out_temp2, out_temp)
  }

  for (v in c("covid_prev", "lag1_covid_prev", "lag2_covid_prev")) {
    m_eq = paste0("n_res ~ ", v, "+ offset(log(nbjh))")
    m_out = glm.nb(m_eq, data = final_db, link = log)

    out_temp = data.frame(
      var_name = v,
      p = summary(m_out)$coefficients[v,"Pr(>|z|)"],
      estimate = exp(summary(m_out)$coefficients[v,"Estimate"]),
      estimate_lw = exp(confint(m_out)[v,"2.5 %"]),
      estimate_up = exp(confint(m_out)[v,"97.5 %"])
    )
    out = bind_rows(out, out_temp)
  }

  out$bacteria = b
  out$setting = s
  univariate_results = bind_rows(univariate_results, out)
}

# Tables for each drug-bacterium pair
covid_var_names = c("covid_prev" = "Week w",
                    "lag1_covid_prev" = "Week w-1",
                    "lag2_covid_prev" = "Week w-2",

                    "periods" = "Pandemic periods w",
                    "periodsfirst wave" = "First wave",
                    "periodslow to no res" = "Low to no res",
                    "periodsmild res" = "Mild res",
                    "periodsstrong res" = "Strong res",

                    "lag1_periods" = "Pandemic periods w-1",
                    "lag1_periodsfirst wave" = "First wave",
                    "lag1_periodslow to no res" = "Low to no res",
                    "lag1_periodsmild res" = "Mild res",
                    "lag1_periodsstrong res" = "Strong res",

                    "lag2_periods" = "Pandemic periods w-2",
                    "lag2_periodsfirst wave" = "First wave",
                    "lag2_periodslow to no res" = "Low to no res",
                    "lag2_periodsmild res" = "Mild res",
                    "lag2_periodsstrong res" = "Strong res"
                    )

for (b in bacterias) {

  univariate_tab = univariate_results %>%
    filter(bacteria == b) %>%
    mutate(
      IRR = ifelse(is.na(estimate), "", paste0(round(estimate, 2), " (", round(estimate_lw, 2), ", ", round(estimate_up,2), ")")),
      setting = ifelse(setting == "icu", "ICU", "Hospital"),
      SubVariable = recode(var_name, !!!covid_var_names),
      Variable = case_when(
        grepl("covid_prev", var_name) ~ "COVID-19 prev.",
        grepl("^periods", var_name) ~ "Pandemic periods w",
        grepl("^lag1_periods", var_name) ~ "Pandemic periods w-1",
        grepl("^lag2_periods", var_name) ~ "Pandemic periods w-2"
      ),
      p = ifelse(p<0.0001, "<0.0001", format(round(p, 4), nsmall = 4))
    ) %>%
    mutate(SubVariable = ifelse(Variable == SubVariable, "", SubVariable)) %>%
    dplyr::select(setting, Variable, SubVariable, IRR, p) %>%
    arrange(setting) %>%
    group_by(setting, Variable) %>%
    mutate(id = 1:n()) %>%
    ungroup() %>%
    group_by(setting) %>%
    gt(.) %>%
    cols_hide(columns = id) %>%
    cols_label(Variable = b, SubVariable = "", p = "p-value") %>%
    tab_style(
      style = list(cell_text(weight = "bold")),
      locations = list(cells_row_groups(), cells_column_labels(), cells_body(columns = p, rows = p <= 0.05))) %>%
    # tab_options(table.font.size = 11) %>%
    sub_values(
      columns = Variable,
      rows = id > 1,
      replacement = "",
      pattern = ".*"
    ) %>%
    tab_footnote(
      locations = cells_body(rows = grepl("periods", Variable) & SubVariable == "",
                             columns = p),
      footnote = "Global p-value",
      placement = "right"
      )
  gtsave(univariate_tab, paste0("../Paper/Supplementary/national_estimates_univariate_", tolower(paste0(b, collapse = "")),".docx"))
  gtsave(univariate_tab, paste0("tables/national_estimates_univariate_", tolower(paste0(b, collapse = "")),".docx"))
}

##################################################
# Sensitivity analysis
# Multicolinearity between Covid-19 variables in
# model 7
##################################################
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
    if (db == "hospital") final_db = res_national
    if (db == "icu") final_db = res_national_icu

    final_db = final_db %>%
      filter(bacterie == bacterias[i]) %>%
      arrange(Date_week) %>%
      mutate(
        lag1_i_res = lag(n_res / nbjh * 1000),

        lag1_covid_prev = lag(covid_prev),
        lag2_covid_prev = lag(covid_prev,2),

        lag1_periods = lag(periods, 1),
        lag2_periods = lag(periods, 2),

        nbjh = nbjh/1000
      ) %>%
      filter(!is.na(lag2_periods)) %>%
      mutate(
        Carbapenems = (Carbapenems - mean(Carbapenems)) / sd(Carbapenems),
        Penicillins = (Penicillins - mean(Penicillins)) / sd(Penicillins),
        TGC = (TGC - mean(TGC)) / sd(TGC),
        lag1_i_res = (lag1_i_res - mean(lag1_i_res)) / sd(lag1_i_res),

        covid_prev  = (covid_prev - mean(covid_prev )) / sd(covid_prev),
        lag1_covid_prev = (lag1_covid_prev - mean(lag1_covid_prev)) / sd(lag1_covid_prev),
        lag2_covid_prev = (lag2_covid_prev - mean(lag2_covid_prev)) / sd(lag2_covid_prev)
      )

    # Multivariate model
    m_eq = "n_res ~ lag1_i_res + covid_prev + periods + TGC + offset(log(nbjh))"
    if (bacterias[i] == "CR P. aeruginosa") m_eq = "n_res ~ lag1_i_res + covid_prev + periods + Carbapenems + offset(log(nbjh))"
    if (bacterias[i] == "MRSA") m_eq = "n_res ~ lag1_i_res + covid_prev + periods + Penicillins + offset(log(nbjh))"

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
  "covid_prev" = "Covid-19 prev. w",
  "Penicillins" = "Penicillins",
  "TGC" = "3rd generation Cephalosporins",
  "Carbapenems" = "Imipenem + Meropenem"
)

vif_sensitivity_tab = all_vif %>%
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
  sub_missing(missing_text = "-")
vif_sensitivity_tab
gtsave(vif_sensitivity_tab, "../Paper/Supplementary/national_vif.docx")
gtsave(vif_sensitivity_tab, "tables/national_vif.docx")

##################################################
# Comparison of regression models with Covid-19 
# variables of w+1, w+2 or w+3
##################################################
all_aic = data.frame()
all_estimates_leads = data.frame()
all_fits_leads = data.frame()
models = paste0("model", 7:9)
names(models) = paste0("Covid-19 prevalence w+", 1:3)
b = "CR P. aeruginosa"

for (db in c("icu", "hospital")) {
  for (mod in models) {
    # Create final database 
    if (db == "hospital") final_db = res_national
    if (db == "icu") final_db = res_national_icu
    
    final_db = final_db %>% 
      filter(bacterie == b) %>%
      arrange(Date_week) %>%
      mutate(
        lag1_i_res = lag(n_res / nbjh * 1000, 2),
        
        lead1_covid_prev = lead(covid_prev, 1),
        lead2_covid_prev = lead(covid_prev,2),
        lead3_covid_prev = lead(covid_prev,3),
        
        lag1_covid_prev = lag(covid_prev), 
        lag2_covid_prev = lag(covid_prev, 2),
        
        nbjh = nbjh/1000
      ) %>% 
      filter(!is.na(lag2_covid_prev), !is.na(lead3_covid_prev)) %>%
      mutate(
        Carbapenems = (Carbapenems - mean(Carbapenems)) / sd(Carbapenems),
        lag1_i_res = (lag1_i_res - mean(lag1_i_res)) / sd(lag1_i_res),
        
        covid_prev = (covid_prev -mean(covid_prev))/sd(covid_prev),
        lag1_covid_prev = (lag1_covid_prev -mean(lag1_covid_prev))/sd(lag1_covid_prev),
        lag2_covid_prev = (lag2_covid_prev -mean(lag2_covid_prev))/sd(lag2_covid_prev),
        
        lead1_covid_prev = (lead1_covid_prev - mean(lead1_covid_prev)) / sd(lead1_covid_prev),
        lead2_covid_prev = (lead2_covid_prev - mean(lead2_covid_prev)) / sd(lead2_covid_prev),
        lead3_covid_prev = (lead3_covid_prev - mean(lead3_covid_prev)) / sd(lead3_covid_prev)
      )
    
    # Multivariate model
    if (mod == "model7") {
      m_eq = "n_res ~ lag1_i_res + lead1_covid_prev + Carbapenems + offset(log(nbjh))"
    }
    
    if (mod == "model8") {
      m_eq = "n_res ~ lag1_i_res + lead2_covid_prev + Carbapenems + offset(log(nbjh))"
    }
    
    if (mod == "model9") {
      m_eq = "n_res ~ lag1_i_res + lead3_covid_prev + Carbapenems + offset(log(nbjh))"
    }
    
    m = glm.nb(m_eq, data=final_db, link=log)
    
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
  results %>%
    filter(model == "model0") %>%
    dplyr::select(model, bacteria, setting, aic) %>%
    mutate(model = "Baseline model"),
  all_aic
  ) %>%
  filter(bacteria == "CR P. aeruginosa") %>%
  dplyr::select(setting, model, aic) %>%
  mutate(
    aic = round(aic), 
    setting = ifelse(setting == "icu", "ICU", "Hospital"),
    model = recode(model, !!!c("best_model" = "Covid prev. w-2 (best model)", 
                               "model7" = "Covid prev. w+1",
                               "model8" = "Covid prev. w+2",
                               "model9" = "Covid prev. w+3"
                               ))
    ) %>%
  arrange(setting) %>%
  rename(Setting = setting, Model = model, AIC = aic) %>%
  ggtexttable(., rows = NULL, theme = ttheme("light",base_size = 9))

# Estimates
p2=all_estimates_leads %>%
  filter(!variable %in% c("(Intercept)", "Carbapenems", "lag1_i_res")) %>%
  mutate(
    setting = ifelse(setting == "icu", "ICU", "Hospital"),
    variable = case_when(model == "model7" ~ "Covid-19 prev. w+1", 
                      model == "model8" ~ "Covid-19 prev. w+2",
                      model == "model9" ~ "Covid-19 prev. w+3")
  ) %>%
  ggplot(., aes(x = variable, y = exp(Estimate), ymin = exp(q2_5), ymax = exp(q97_5), col = setting)) +
  geom_hline(yintercept = 1) +
  geom_pointrange(position = position_dodge(width = 0.5)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  labs(x = "", y = "IRRs (95% CI)", col = "")

# Fits
p3=all_fits %>%
  inner_join(., best_models, by = c("setting", "bacteria", "model")) %>%
  filter(bacteria == "CR P. aeruginosa") %>%
  bind_rows(., all_fits_leads) %>%
  mutate(
    setting = ifelse(setting == "icu", "ICU", "Hospital"),
    model = recode(model, !!!c("model6" = "Covid-19 prev. w-2\n(best model)", 
                               "model7" = "Covid-19 prev. w+1", 
                               "model8" = "Covid-19 prev. w+2",
                               "model9" = "Covid-19 prev. w+3"))
  ) %>%
  ggplot(., aes(x = Date_week)) +
  geom_rect(data = int_national_start_end, 
            aes(NULL, NULL, xmin=start, xmax=end, fill=restrictions, ymin=-Inf, ymax=Inf)) +
  facet_grid(cols = vars(setting), rows = vars(model)) +
  geom_line(aes(y = n_res, col = "Data")) +
  geom_line(aes(y = exp(fit), col = "Fit")) +
  geom_ribbon(aes(ymin = exp(fit + qnorm(0.025) * se.fit), ymax = exp(fit + qnorm(0.975) * se.fit)), fill = "red", alpha = 0.4) +
  scale_color_manual(values = c("Data" = "black", "Fit" = "red")) +
  scale_fill_manual(
    name = "Anti-Covid-19 interventions",
    labels = c("First wave", "Strong", "Intermediary", "Light to none"),
    breaks = c("p_first_wave", "p_strong_res", "p_mild_res", "p_no_res"),
    values = c("p_first_wave" = col_interventions(1), "p_strong_res" = col_interventions(2), "p_mild_res" = col_interventions(3), "p_no_res" = col_interventions(4))
  ) +
  guides(fill=guide_legend(order=1, override.aes = list(col = 'black'))) +
  theme_bw() +
  labs(x = "", y = "Weekly no. of incident CR P. aeruginosa isolates", col = "")

# Final Supplementary figure 
fig = plot_grid(
  plot_grid(
    p1, p2, ncol = 2, labels = c("A", "B")#, rel_widths = c(0.6, 1)
    ),
  p3, nrow = 2, rel_heights = c(0.5,1), labels = c("", "C")
)
fig
ggsave("../Paper/Supplementary/crpa_leads.png", fig, height = 9, width = 8)
ggsave("plots/regressions/crpa_leads.png", fig, height = 9, width = 8)


##################################################
# Hospital without ICUs at the national level
# for CR-PA
##################################################
# Plot resistance incidence data
res_national_not_icu %>%
  filter(bacterie=="CR P. aeruginosa") %>%
  ggplot(., aes(x = Date_week, y = n_res)) +
  geom_line() +
  theme_bw()

# Best regression model
results_not_icu = data.frame(
    setting = "not icu",
    model = names(model_names),
    bacteria = "CR P. aeruginosa",
    poisson_OD = NA,
    negbin_OD = NA,
    theta = NA,
    theta_min = NA,
    theta_max = NA,
    pseudoR2 = NA,
    aic = NA
  ) %>%
  arrange(model)

all_estimates_not_icu = data.frame()
all_residuals_not_icu = data.frame()
all_vif_not_icu = vector("list", n_models)
all_models_not_icu = vector("list", n_models)

# Regression models
b = "CR P. aeruginosa"
db = "not icu"
for (i in 1:nrow(results_not_icu)) {
  mod =results_not_icu$model[i]

  # Create final database
  final_db = res_national_not_icu

  final_db = final_db %>%
    filter(bacterie == b) %>%
    arrange(Date_week) %>%
    mutate(
      lag1_i_res = lag(n_res / nbjh * 1000),

      lag1_covid_prev = lag(covid_prev),
      lag2_covid_prev = lag(covid_prev,2),

      lag1_periods = lag(periods, 1),
      lag2_periods = lag(periods, 2),

      nbjh = nbjh/1000
    ) %>%
    filter(!is.na(lag2_periods)) %>%
    mutate(
      Carbapenems = (Carbapenems - mean(Carbapenems)) / sd(Carbapenems),

      lag1_i_res = (lag1_i_res - mean(lag1_i_res)) / sd(lag1_i_res),

      covid_prev = (covid_prev - mean(covid_prev)) / sd(covid_prev),
      lag1_covid_prev = (lag1_covid_prev - mean(lag1_covid_prev)) / sd(lag1_covid_prev),
      lag2_covid_prev = (lag2_covid_prev - mean(lag2_covid_prev)) / sd(lag2_covid_prev),

      Date_year = as.character(Date_year)
    )

  # Multivariate model
  if (mod == "model0") m_eq = "n_res ~ lag1_i_res + Carbapenems + offset(log(nbjh))"
  if (mod == "model1") m_eq = "n_res ~ lag1_i_res + periods + Carbapenems + offset(log(nbjh))"
  if (mod == "model2") m_eq = "n_res ~ lag1_i_res + lag1_periods + Carbapenems + offset(log(nbjh))"
  if (mod == "model3") m_eq = "n_res ~ lag1_i_res + lag2_periods + Carbapenems + offset(log(nbjh))"
  if (mod == "model4") m_eq = "n_res ~ lag1_i_res + covid_prev + Carbapenems + offset(log(nbjh))"
  if (mod == "model5") m_eq = "n_res ~ lag1_i_res + lag1_covid_prev + Carbapenems + offset(log(nbjh))"
  if (mod == "model6") m_eq = "n_res ~ lag1_i_res + lag2_covid_prev + Carbapenems + offset(log(nbjh))"

  # Overdispersion in Poisson regression model
  results_not_icu$poisson_OD[i] = check_overdispersion(glm(m_eq, data = final_db, family = poisson))$p_value

  # # Negative binomial regression model
  # m1 = glm.nb(m_eq1, data = final_db, link = log)
  m = glm.nb(m_eq, data = final_db, link = log)
  # anova(m, m1)

  results_not_icu$negbin_OD[i] = check_overdispersion(m)$p_value
  results_not_icu$theta[i] = m$theta
  results_not_icu$theta_min[i] = m$theta + qnorm(0.025) * m$SE.theta
  results_not_icu$theta_max[i] = m$theta + qnorm(0.975) * m$SE.theta
  results_not_icu$pseudoR2[i] = DescTools::PseudoR2(m)

  # Save the models
  all_models_not_icu[[i]] = m

  # Colinearity of predictors
  if (mod != "model0") all_vif_not_icu[[i]] = VIF(m)

  # AIC
  results_not_icu$aic[i] = AIC(m)

  # Residuals
  mod_residuals = data.frame(
    setting = db,
    model = mod,
    bacteria = b,
    residuals = residuals(m)
  )
  all_residuals_not_icu = bind_rows(all_residuals_not_icu, mod_residuals)

  # Estimates and p-values
  mod_estimates = data.frame(cbind(summary(m)$coefficients, confint(m))) %>%
    rownames_to_column(var = 'variable') %>%
    rename(q2_5 = `X2.5..`, q97_5 = `X97.5..`, p = `Pr...z..`, z = `z.value`, sd = `Std..Error`) %>%
    mutate(bacteria = b, model = mod, setting = db)
  all_estimates_not_icu = bind_rows(all_estimates_not_icu, mod_estimates)
}

# Best model
best_not_icu = results_not_icu %>%
  mutate(aic = round(aic), setting = "Hospital (except ICU)") %>%
  dplyr::select(setting, model, aic) %>%
  group_by(setting) %>%
  gt(.) %>%
  cols_label(model = "", aic = "CR P. aeruginosa") %>%
  tab_style(
    style = list(cell_text(weight = "bold")),
    locations = list(cells_row_groups(), cells_column_labels(), cells_body(columns = aic, rows = aic == min(aic))))
best_not_icu
gtsave(best_not_icu, "../Paper/Supplementary/best_model_not_icu_crpa.docx")
gtsave(best_not_icu, "tables/best_model_not_icu_crpa.docx")
