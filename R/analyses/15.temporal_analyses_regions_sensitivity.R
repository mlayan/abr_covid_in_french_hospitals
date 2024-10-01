##################################################
## COUNT REGRESSION ANALYSIS AT THE REGIONAL
## LEVEL IN HOSPITALS 
##################################################
rm(list = ls())
library(tidyverse)
library(performance)
library(MASS)
library(DescTools)
library(qqplotr)
library(geofacet)
library(ggpubr)
library(cowplot)
library(gt)

# Helper functions
source("R/helper/helper_functions.R")
source("R/helper/helper_plots.R")
source("R/helper/dictionaries.R")

##################################################
# Load data 
##################################################
# Best regression models at the national level
load("data/best_models.rda")

# Interventions
load("data/int_region.rda")
int_region = int_region %>%
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
  dplyr::select(Date_week, region, periods) %>%
  mutate(periods = factor(periods, c("pre-pandemic", "first wave", "strong res", "mild res", "low to no res")))

# Resistances
load("data/res_regions.rda")

# Bed-days
load("data/bd_pmsi_regions.rda")

# Covid-19 intubation data
load("data/covid_intub_region.rda")

##################################################
# Combine data 
##################################################
regional_df = res_regions %>%
  left_join(., int_region, by = c("Date_week", "region")) %>%
  left_join(., bd_pmsi_regions, by = c("Date_week", "region")) %>%
  left_join(., covid_intub_region, by = c("Date_week", "region")) %>%
  mutate(
    Date_year = lubridate::year(Date_week),
    covid_intub = ifelse(is.na(covid_intub), 0, covid_intub)
  ) %>%
  mutate(covid_intub_prev = covid_intub / nbjh * 1000)

##################################################
# Regression analysis
##################################################
# French regions
regions=sort(unique(regional_df$region))
n_regions = length(regions)

# Initialize objects to store results of regression analyses
results = data.frame(
  region = regions,
  poisson_OD = NA, 
  negbin_OD = NA, 
  theta_OD = NA, 
  aic = NA
  )

all_vif = vector("list", n_regions)
all_estimates = data.frame()
all_residuals = data.frame()
all_fits = data.frame()
k=1

# Regression 
b = "CR P. aeruginosa"
for (r in regions) {
    
    # Create final database 
    final_db = regional_df %>% 
      filter(bacterie == b, region == r) %>%
      arrange(Date_week) %>%
      mutate(
        lag1_i_res = lag(n_res / nbjh * 1000),
        lag1_periods = lag(periods, 1),
        lag2_covid_intub_prev = lag(covid_intub_prev, 2)
      ) %>% 
      filter(!is.na(lag2_covid_intub_prev)) %>%
      mutate(
        lag1_i_res = (lag1_i_res - mean(lag1_i_res)) / sd(lag1_i_res),
        lag2_covid_intub_prev = (lag2_covid_intub_prev - mean(lag2_covid_intub_prev)) / sd(lag2_covid_intub_prev)
        )
    
    # Multivariate model
    mod = "model2"
    
    if (mod == "model2") {
      m_eq = "n_res ~ lag1_i_res + lag1_periods + offset(log(nbjh))"
    }
    
    if (mod == "model6") {
      m_eq = "n_res ~ lag1_i_res + lag2_covid_intub_prev + offset(log(nbjh))"
    }
    
    if (r == "Hauts-de-France") {
      m = glm(m_eq, family = poisson, data = final_db)
      results$poisson_OD[k] = check_overdispersion(glm(m_eq, family = poisson, data = final_db))$p_value

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

# Is poisson necessary
results %>%
  filter(poisson_OD > 0.05)

# P-values COVID-19-related variable
all_estimates %>%
  filter(grepl("lag1_periods", variable),
         region %in% c("Auvergne-Rhône-Alpes", "Grand-Est", "Provence-Alpes-Côte d'Azur", "Île-de-France")) %>%
  dplyr::select(region, variable, Estimate, q2_5, q97_5, p) %>%
  mutate(Estimate = exp(Estimate), q2_5 = exp(q2_5), q97_5 = exp(q97_5))

##################################################
# Final plot
##################################################
covid_var_names = c(
  "lag1_periodsfirst wave" = "First wave w-1",
  "lag1_periodsstrong res" = "Strong w-1",
  "lag1_periodsmild res" = "Mild w-1",
  "lag1_periodslow to no res" = "Low to none w-1"
)

model2_estimates = all_estimates %>%
  filter(grepl("lag1_periods", variable)) %>%
  mutate(
    variable = factor(recode(variable, !!!covid_var_names), as.character(covid_var_names)),
    region = recode(region, !!!dict_regions),
    alpha_level = ifelse(p <= 0.05, 1, 0.2)
  ) %>%
  ggplot(., aes(x = fct_rev(region), y = exp(Estimate), ymin = exp(q2_5), ymax = exp(q97_5), col = fct_rev(variable))) +
  facet_grid(cols = vars(bacteria)) +
  coord_flip() +
  geom_hline(yintercept = 1, linetype = "dashed", col = "grey70") +
  geom_pointrange(aes(alpha = alpha_level), position = position_dodge(width = 0.5), fatten = 1) +
  scale_color_manual(values = rev(colorRampPalette(c("dodgerblue3", "lightblue1"))(4))) +
  scale_alpha_identity() +
  theme_bw() +
  theme(legend.position = "right", legend.title = element_blank()) +
  guides(alpha = "none") +
  labs(x = "", y = "Incidence rate ratio (95% CI)") 
ggsave("../Paper/Supplementary/estimates_regions_sensitivity.png", model2_estimates, height = 7, width = 6)
ggsave("plots/regressions/estimates_regions_sensitivity.png", model2_estimates, height = 7, width = 6)

