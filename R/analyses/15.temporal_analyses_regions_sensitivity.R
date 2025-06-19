##################################################
## COUNT REGRESSION ANALYSIS AT THE REGIONAL
## LEVEL IN HOSPITALS 
##################################################
rm(list = ls())
library(tidyverse)
library(performance)
library(MASS)
library(DescTools)
library(ggpubr)
library(gt)

# Helper functions
source("R/helper/helper_functions.R")
source("R/helper/dictionaries.R")

##################################################
# Load data 
##################################################
# Best regression models at the national level
load("data/best_models_sensitivity.rda")

# Interventions
load("data/int_region.rda")
int_region = int_region %>%
  mutate(
    periods = case_when(
      p_strong_res + p_mild_res + p_no_res + p_first_wave == 0 ~ "pre-pandemic",
      p_first_wave > 0.5 ~ "first wave",
      p_strong_res > 0.5 ~ "strong res", 
      p_mild_res > 0.5 ~ "mild res",
      p_no_res > 0.5 ~ "low to no res",
      .default = NA
    ) 
  ) %>%
  dplyr::select(Date_week, region, periods) %>%
  mutate(periods = factor(periods, c("pre-pandemic", "first wave", "strong res", "mild res", "low to no res")))

# Resistances
load("data/res_regions.rda")

# Bed-days
load("data/bd_pmsi_regions.rda")

# Covid-19 intubation data
load("data/covid_region.rda")

##################################################
# Combine data 
##################################################
regional_df = res_regions %>%
  filter(bacterie %in% c("CR P. aeruginosa", "ESBL E. coli", "MRSA")) %>%
  left_join(., int_region, by = c("Date_week", "region")) %>%
  left_join(., bd_pmsi_regions, by = c("Date_week", "region")) %>%
  left_join(., covid_region, by = c("Date_week", "region")) %>%
  mutate(
    Date_year = lubridate::year(Date_week),
    covid = ifelse(is.na(covid), 0, covid)
  ) %>%
  mutate(covid_prev = covid / nbjh * 1000)

##################################################
# Regression analysis
##################################################
# French regions
regions=sort(unique(regional_df$region))
n_regions = length(regions)

# Initialize objects to store results of regression analyses
results = expand.grid(region = regions, bacterie = c("MRSA", "ESBL E. coli", "CR P. aeruginosa")) %>%
  mutate(
    poisson_OD = NA, 
    negbin_OD = NA, 
    theta_OD = NA, 
    aic = NA
  )

all_vif = vector("list", nrow(results))
all_estimates = data.frame()
all_residuals = data.frame()
all_fits = data.frame()
k=1

# Regression 
for (i in 1:nrow(results)) {
  
  r = as.character(results$region[i])
  b = as.character(results$bacterie[i])
    
  # Create final database 
  final_db = regional_df %>% 
    filter(bacterie == b, region == r) %>%
    arrange(Date_week) %>%
    mutate(
      lag1_i_res = lag(n_res / nbjh * 1000),
      lag1_periods = lag(periods, 1),
      lag2_covid_prev = lag(covid_prev, 2)
    ) %>% 
    filter(!is.na(lag2_covid_prev)) %>%
    mutate(
      lag1_i_res = (lag1_i_res - mean(lag1_i_res)) / sd(lag1_i_res),
      lag2_covid_prev = (lag2_covid_prev - mean(lag2_covid_prev)) / sd(lag2_covid_prev)
      )
  
  # Multivariate model
  mod = best_models_sensitivity$model[best_models_sensitivity$setting == "hospital" & best_models_sensitivity$bacteria == b]
  
  if (mod == "model1") {
    m_eq = "n_res ~ lag1_i_res + periods  + offset(log(nbjh))"
  }
  
  if (mod == "model2") {
    m_eq = "n_res ~ lag1_i_res + lag1_periods + offset(log(nbjh))"
  }
  
  if (check_overdispersion(glm(m_eq, family = poisson, data = final_db))$p_value > 0.05) {
    m = glm(m_eq, family = poisson, data = final_db)
    results$poisson_OD[k] = check_overdispersion(m)$p_value

  } else {
    m = glm.nb(m_eq, data = final_db, link = log)
    results$poisson_OD[k] = check_overdispersion(glm(m_eq, family = poisson, data = final_db))$p_value
    results$negbin_OD[k] = check_overdispersion(glm.nb(m_eq, data = final_db, link = log))$p_value
    results$theta_OD[k] = m$theta
  }

  # Model checks
  all_vif[[k]] = VIF(m)
  results$aic[k] = AIC(m)
  
  # Residuals
  mod_residuals = data.frame(
    region = r,
    model = mod,
    bacteria = b,
    residuals = residuals(m)
  )
  all_residuals = bind_rows(all_residuals, mod_residuals)
  
  # Estimates and p-values
  mod_estimates = data.frame(cbind(summary(m)$coefficients, confint(m))) %>%
    rownames_to_column(var = 'variable') %>%
    rename(q2_5 = `X2.5..`, q97_5 = `X97.5..`, p = `Pr...z..`, z = `z.value`, sd = `Std..Error`) %>%
    mutate(bacteria = b, region = r)
  all_estimates = bind_rows(all_estimates, mod_estimates)

  # Model fit
  mod_fit = data.frame(
    region = r,
    bacteria = b,
    Date_week = final_db$Date_week, 
    n_res = final_db$n_res, 
    nbjh = final_db$nbjh,
    fit = predict(m, type = "link", se.fit = T)$fit,
    se.fit = predict(m, type = "link", se.fit = T)$se.fit
  )
  all_fits = bind_rows(all_fits, mod_fit)
  
  k = k+1
} 

##################################################
# Final plot
##################################################
p1 = all_estimates %>%
  filter(grepl("lag1_periods", variable)) %>%
  mutate(
    variable = factor(recode(variable, !!!covid_var_names), as.character(covid_var_names)),
    region = recode(region, !!!dict_regions),
    alpha_level = ifelse(p <= 0.05, 1, 0.1)
  ) %>%
  ggplot(., aes(x = fct_rev(region), y = exp(Estimate), ymin = exp(q2_5), ymax = exp(q97_5), col = fct_rev(variable))) +
  facet_grid(cols = vars(bacteria)) +
  coord_flip() +
  geom_hline(yintercept = 1, linetype = "dashed", col = "grey70") +
  geom_pointrange(aes(alpha = alpha_level), position = position_dodge(width = 0.5), fatten = 1) +
  scale_color_manual(values = rev(c("dodgerblue3", "deeppink", "black", "goldenrod1"))) +
  scale_alpha_identity() +
  theme_bw() +
  theme(legend.position = "bottom", legend.title = element_blank(), legend.direction = "vertical") +
  guides(alpha = "none", color = guide_legend(reverse = T)) +
  labs(x = "", y = "Incidence rate ratio (95% CI)") 


p2 = all_estimates %>%
  filter(grepl("^periods", variable)) %>%
  mutate(
    variable = factor(recode(variable, !!!covid_var_names), as.character(covid_var_names)),
    region = recode(region, !!!dict_regions),
    alpha_level = ifelse(p <= 0.05, 1, 0.1)
  ) %>%
  ggplot(., aes(x = fct_rev(region), y = exp(Estimate), ymin = exp(q2_5), ymax = exp(q97_5), col = fct_rev(variable))) +
  facet_grid(cols = vars(bacteria)) +
  coord_flip() +
  geom_hline(yintercept = 1, linetype = "dashed", col = "grey70") +
  geom_pointrange(aes(alpha = alpha_level), position = position_dodge(width = 0.5), fatten = 1) +
  scale_color_manual(values = rev(c("dodgerblue3", "deeppink", "black", "goldenrod1"))) +
  scale_alpha_identity() +
  theme_bw() +
  theme(legend.position = "bottom", legend.title = element_blank(), legend.direction = "vertical") +
  expand_limits(y = c(0,5)) +
  guides(alpha = "none", color = guide_legend(reverse = T)) +
  labs(x = "", y = "Incidence rate ratio (95% CI)") 
p = ggarrange(p1, p2, labels = c("A", "B"), widths = c(1.8,1))
p
ggsave("../Paper/Supplementary/estimates_regions_sensitivity.png", p, height = 7, width = 8)
ggsave("plots/regressions/estimates_regions_sensitivity.png", p, height = 7, width = 8)

# P-values COVID-19-related variable
all_estimates %>%
  filter(grepl("lag1_periods", variable),
         region %in% c("Auvergne-Rhône-Alpes", "Grand-Est", "Provence-Alpes-Côte d'Azur", "Île-de-France")) %>%
  dplyr::select(region, variable, Estimate, q2_5, q97_5, p) %>%
  mutate(Estimate = exp(Estimate), q2_5 = exp(q2_5), q97_5 = exp(q97_5))
