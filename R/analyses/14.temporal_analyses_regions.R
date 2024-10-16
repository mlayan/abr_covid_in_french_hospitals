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

# Antibiotic consumption
load("data/atb_use_regions.rda")
atb_use_regions = atb_use_regions %>%
  filter(atb_class %in% c("Imipenem + Meropenem", "Penicillins", "Third generation Cephalosporins")) %>%
  mutate(consumption = molDDD/Nbhosp*1000, 
         atb_class = case_when(
           atb_class == "Imipenem + Meropenem" ~ "Carbapenems", 
           atb_class == "Third generation Cephalosporins" ~ "TGC",
           .default = atb_class
         )) %>%
  dplyr::select(Date_year, region, atb_class, consumption) %>%
  pivot_wider(names_from = atb_class, values_from = consumption)

# Bed-days
load("data/bd_pmsi_regions.rda")

# Covid-19 intubation data
load("data/covid_intub_region.rda")

##################################################
# Combine data 
##################################################
regional_df = res_regions %>%
  filter(bacterie %in% c("CR P. aeruginosa", "ESBL E. coli", "MRSA")) %>%
  left_join(., int_region, by = c("Date_week", "region")) %>%
  left_join(., bd_pmsi_regions, by = c("Date_week", "region")) %>%
  left_join(., covid_intub_region, by = c("Date_week", "region")) %>%
  mutate(
    Date_year = lubridate::year(Date_week),
    covid_intub = ifelse(is.na(covid_intub), 0, covid_intub)
  ) %>%
  left_join(., atb_use_regions, by = c("Date_year", "region")) %>%
  mutate(covid_intub_prev = covid_intub / nbjh * 1000)

##################################################
# Regions most affected by the pandemic
##################################################
df = bd_pmsi_regions %>%
  left_join(., covid_intub_region, by = c("Date_week", "region")) %>%
  mutate(Date_year = lubridate::year(Date_week),
         region = recode(region, !!!dict_regions),
         covid_intub = ifelse(is.na(covid_intub), 0, covid_intub)) %>%
  group_by(Date_year, region) %>%
  summarise(p = sum(covid_intub)/sum(nbjh)*1000, .groups = "drop") %>%
  arrange(Date_year, desc(p)) %>%
  group_by(Date_year) %>%
  mutate(o = 1:n()) %>%
  ungroup() 

p2020 = df %>% 
  filter(Date_year == 2020) %>%
  ggplot(., aes(x = reorder(region, desc(o)), y = p)) +
  geom_bar(stat = "identity", width = 0.5) +
  facet_wrap(facets = vars(Date_year), ncol = 3, scales = "free_y") +
  coord_flip() +
  theme_bw() + 
  theme(axis.title = element_blank()) 

p2021 = df %>% 
  filter(Date_year == 2021) %>%
  ggplot(., aes(x = reorder(region, desc(o)), y = p)) +
  geom_bar(stat = "identity", width = 0.5) +
  facet_wrap(facets = vars(Date_year), ncol = 3, scales = "free_y") +
  coord_flip() +
  theme_bw() + 
  theme(axis.title = element_blank()) +
  labs(y = "", x = "") 

p2022 = df %>% 
  filter(Date_year == 2022) %>%
  ggplot(., aes(x = reorder(region, desc(o)), y = p)) +
  geom_bar(stat = "identity", width = 0.5) +
  facet_wrap(facets = vars(Date_year), ncol = 3, scales = "free_y") +
  coord_flip() +
  theme_bw() + 
  theme(axis.title = element_blank()) +
  labs(y = "", x = "") 

p_title = title <- ggdraw() + 
  draw_label(
    "Intubated COVID-19 patients prevalence (per 1,000 bed-days)",
    hjust = 0.5,
    size = 11
  )

p_covid = plot_grid(
  plot_grid(p2020, p2021, p2022, ncol = 3),
  p_title, labels = c("C", "", ""),
  nrow = 2, rel_heights = c(1, 0.05)
)

##################################################
# Plot intubated COVID-19 patient prevalence 
##################################################
regional_df %>% 
  filter(bacterie == "CR P. aeruginosa") %>%
  mutate(region = recode(region, !!!dict_regions)) %>%
  ggplot(., aes(x = Date_week, y = covid_intub_prev)) +
  geom_line() +
  facet_wrap(facets = vars(region), ncol = 4) +
  theme_bw() +
  labs(x = "", y = "Intubated COVID-19 patients prevalence (per 1,000 bed-days)")

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
      lag2_covid_intub_prev = lag(covid_intub_prev, 2)
    ) %>% 
    filter(!is.na(lag2_covid_intub_prev)) %>%
    mutate(
      lag1_i_res = (lag1_i_res - mean(lag1_i_res)) / sd(lag1_i_res),
      lag2_covid_intub_prev = (lag2_covid_intub_prev - mean(lag2_covid_intub_prev)) / sd(lag2_covid_intub_prev),
      
      Carbapenems = (Carbapenems - mean(Carbapenems)) / sd(Carbapenems),      
      Penicillins = (Penicillins - mean(Penicillins)) / sd(Penicillins),
      TGC = (TGC - mean(TGC)) / sd(TGC),
      
      lag1_i_res = (lag1_i_res - mean(lag1_i_res)) / sd(lag1_i_res)
    )
  
  # Multivariate model
  mod = best_models$model[best_models$bacteria == b & best_models$setting == "hospital"]
  
  if (mod == "model1") {
    m_eq = "n_res ~ lag1_i_res + periods + TGC + offset(log(nbjh))"
    if (b == "CR P. aeruginosa") m_eq = "n_res ~ lag1_i_res + periods + Carbapenems + offset(log(nbjh))"
    if (b == "MRSA") m_eq = "n_res ~ lag1_i_res + periods + Penicillins + offset(log(nbjh))"
  }
  
  if (mod == "model6") {
    m_eq = "n_res ~ lag1_i_res + lag2_covid_intub_prev + Carbapenems + offset(log(nbjh))"
  }
  
  if (check_overdispersion(glm(m_eq, family = poisson, data = final_db))$p_value > 0.05) {
    m = glm(m_eq, family = poisson, data = final_db)
    results$poisson_OD[i] = check_overdispersion(m)$p_value
  } else {
    m = glm.nb(m_eq, data = final_db, link = log)
    results$poisson_OD[i] = check_overdispersion(glm(m_eq, family = poisson, data = final_db))$p_value
    results$negbin_OD[i] = check_overdispersion(glm.nb(m_eq, data = final_db, link = log))$p_value
    results$theta_OD[i] = m$theta    
  }
    
  # Model checks
  all_vif[[i]] = VIF(m)
  results$aic[i] = AIC(m)
  
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
}


# Is poisson necessary
results %>%
  filter(poisson_OD > 0.05)

# P-values COVID-19-related variable
all_estimates %>%
  filter(variable == "lag2_covid_intub_prev",
         region %in% c("Auvergne-Rhône-Alpes", "Grand-Est", "Provence-Alpes-Côte d'Azur", "Île-de-France")) %>%
  dplyr::select(region, Estimate, q2_5, q97_5, p)


##################################################
# Final plot
##################################################
var_names = c(
  "lag2_covid_intub_prev" = "COVID-19 intub. prev. w-2", 
  "periodsfirst wave" = "First wave w", 
  "periodsstrong res" = "Strong w", 
  "periodsmild res" = "Mild w", 
  "periodslow to no res" = "Low to none w"
)

# Specific case of CR-PA (not on the same scale compared to 
# ESBL E. coli and MRSA + not the same COVID-19 related variables)
model6_estimates_crpa = all_estimates %>%
  filter(variable == "lag2_covid_intub_prev", bacteria == "CR P. aeruginosa") %>%
  mutate(
    variable = factor(recode(variable, !!!var_names), unname(var_names)),
    region = recode(region, !!!dict_regions),
    alpha_level = ifelse(p <= 0.05, 1, 0.2)
  ) %>%
  ggplot(., aes(x = fct_rev(region), y = exp(Estimate), ymin = exp(q2_5), ymax = exp(q97_5), col = variable)) +
  facet_grid(cols = vars(bacteria)) +
  coord_flip() +
  geom_hline(yintercept = 1, linetype = "dashed", col = "grey70") +
  geom_pointrange(aes(alpha = alpha_level), position = position_dodge(width = 0.5), fatten = 1) +
  scale_color_manual(values = "black") +
  scale_alpha_identity() +
  theme_bw() +
  theme(axis.title.y = element_blank(), legend.position = "bottom", legend.title = element_blank()) +
  guides(alpha = "none") +
  labs(y = "Incidence rate ratio (95% CI)") 

  
# ESBL-E. coli and MRSA estimates 
model6_estimates_other = all_estimates %>%
    filter(variable %in% names(var_names), bacteria %in% c("MRSA", "ESBL E. coli")) %>%
    mutate(
      variable = factor(recode(variable, !!!var_names), unname(var_names)),
      region = recode(region, !!!dict_regions),
      alpha_level = ifelse(p <= 0.05, 1, 0.2)
    ) %>%
    ggplot(., aes(x = fct_rev(region), y = exp(Estimate), ymin = exp(q2_5), ymax = exp(q97_5), 
                  col = fct_rev(variable))) +
    facet_grid(cols = vars(bacteria)) +
    coord_flip() +
    geom_hline(yintercept = 1, linetype = "dashed", col = "grey70") +
    geom_pointrange(aes(alpha = alpha_level), position = position_dodge(width = 0.7), fatten = 1) +
    scale_color_manual(values = c(
      "First wave w" = "#313695", 
      "Strong w" = colorRampPalette(c("#313695", "lightskyblue"))(4)[2], 
      "Mild w" = colorRampPalette(c("#313695", "lightskyblue"))(4)[3], 
      "Low to none w" = "lightskyblue"
    )) +
    scale_alpha_identity() +
    theme_bw() +
    theme(axis.title.y = element_blank(), legend.position = "bottom", legend.title = element_blank()) +
    guides(alpha = "none", color = guide_legend(reverse=T)) +
    labs(y = "Incidence rate ratio (95% CI)") 
  

# Save figure 7
figure7 = plot_grid(
  plot_grid(model6_estimates_crpa, model6_estimates_other, 
            ncol = 2, rel_widths = c(1, 1.8), labels = c("A", "B")), 
  plot_grid(p2020, p2021, p2022, ncol = 3, labels = c("C", "", "")),
  p_title,
  nrow = 3,
  rel_heights = c(1.5, 1, 0.05)
)
figure7
ggsave("../Paper/Figures/Figure7.png", figure7, height = 9, width = 7)
ggsave("plots/Figure7.png", figure7, height = 9, width = 7)

