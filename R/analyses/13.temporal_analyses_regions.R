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
load("data/my_regional_grid.rda")
load("data/dict_regions.rda")
load("data/int_region.rda")
load("data/best_models.rda")

##################################################
# Load data 
##################################################
# Resistances
res = read.table("data-raw/spares/combined/resistance_cohortfinal.txt", header = T, sep = "\t")

# Antibiotic consumption
atb = read.table("data-raw/spares/combined/antibiotics_cohortfinal.txt", header = T, sep = "\t") %>%
  filter(secteur == "Hospital") %>%
  group_by(Date_year, region) %>%
  summarise(
    Penicillins  = sum(Penicillins),
    Third_generation_Cephalosporins  = sum(Third_generation_Cephalosporins),
    Carbapenems  = sum(Carbapenems),
    nbjh_spares  = sum(Nbhosp),
    .groups = "drop"
  ) 

# PMSI JH
pmsi_jh = read.table("data-raw/atih/pmsi_sejours.txt", header = T, sep = "\t") %>%
  group_by(Date_week, region) %>%
  summarise(nbjh = sum(nbjh), .groups = "drop") 

# PMSI - COVID
pmsi_covid = read.table("data-raw/atih/pmsi_sejours_covid.txt", header = T, sep = "\t") %>%
  group_by(Date_week, region) %>%
  summarise(covid_jh  = sum(covid_jh), .groups = "drop") 

# PMSI - Intubation
pmsi_intubation = read_sas("data-raw/atih/mco_intubation_weekly.sas7bdat") %>%
  rename(finess = FinessGeo) %>%
  left_join(., metadata_admin_espic %>% dplyr::select(finess, region) %>% distinct(), by = "finess") %>%
  filter(!is.na(region), status == "covid") %>%
  group_by(Date_week, region) %>%
  summarise(nintub = sum(nbjh), .groups = "drop")

##################################################
# Combine data only at the hospital level 
# (too few ICUs by region in our final cohort)
##################################################
# Regional database
regional_df = res %>%
  group_by(Date_year, Date_week, bacterie, region) %>%
  summarise(n_res = sum(n_res), .groups = "drop") %>%
  left_join(., pmsi_covid, by = c("region", "Date_week")) %>%
  left_join(., pmsi_jh, by = c("region", "Date_week")) %>%
  left_join(., atb, by = c("region", "Date_year")) %>%
  mutate(Date_week = as.Date(Date_week)) %>%
  left_join(., pmsi_intubation, by = c("region", "Date_week")) %>%
  mutate(
    covid_jh = ifelse(is.na(covid_jh), 0, covid_jh), 
    nintub = ifelse(is.na(nintub), 0, nintub)
  ) %>%
  mutate(
    covid_jh = covid_jh / nbjh * 1000,
    covid_intub_prop = ifelse(covid_jh == 0, 0, nintub/covid_jh),
    Penicillins = Penicillins / nbjh_spares * 1000,
    Third_generation_Cephalosporins = Third_generation_Cephalosporins / nbjh_spares * 1000,
    Carbapenems = Carbapenems / nbjh_spares * 1000
  ) %>%
  left_join(., int_region, by = c("Date_week", "region")) %>%
  mutate(p_first = ifelse(Date_week >= as.Date("2020-03-16") & Date_week <= as.Date("2020-06-15"), p_strong_res, 0)) %>%
  mutate(p_strong_res = ifelse(Date_week >= as.Date("2020-03-16") & Date_week <= as.Date("2020-06-15"), 0, p_strong_res))

##################################################
# Regions most affected by the pandemic
##################################################
df = pmsi_intubation %>%
  left_join(., pmsi_jh %>% mutate(Date_week = as.Date(Date_week)), by = c("Date_week", "region")) %>%
  mutate(Date_year = as.numeric(format(Date_week,'%Y')),
         region = recode(region, !!!dict_regions)) %>%
  group_by(Date_year, region) %>%
  summarise(p = sum(nintub)/sum(nbjh)*1000, .groups = "drop") %>%
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
  theme(axis.title = element_blank()) +
  ylim(c(0,13))

p2021 = df %>% 
  filter(Date_year == 2021) %>%
  ggplot(., aes(x = reorder(region, desc(o)), y = p)) +
  geom_bar(stat = "identity", width = 0.5) +
  facet_wrap(facets = vars(Date_year), ncol = 3, scales = "free_y") +
  coord_flip() +
  theme_bw() + 
  theme(axis.title = element_blank()) +
  labs(y = "", x = "") +
  ylim(c(0,13))

p2022 = df %>% 
  filter(Date_year == 2022) %>%
  ggplot(., aes(x = reorder(region, desc(o)), y = p)) +
  geom_bar(stat = "identity", width = 0.5) +
  facet_wrap(facets = vars(Date_year), ncol = 3, scales = "free_y") +
  coord_flip() +
  theme_bw() + 
  theme(axis.title = element_blank()) +
  labs(y = "", x = "") +
  ylim(c(0,13))

p_title = title <- ggdraw() + 
  draw_label(
    "Intubated COVID-19 patients prevalence (per 1,000 bed-days)",
    hjust = 0.5,
    size = 11
  )

p_covid = plot_grid(
  plot_grid(p2020, p2021, p2022, ncol = 3),
  p_title, labels = c("B", "", ""),
  nrow = 2, rel_heights = c(1, 0.05)
)

##################################################
# Regression analysis
##################################################
# Basic data
bacterias=c("CR P. aeruginosa", "ESBL E. coli", "MRSA")
regions=sort(unique(regional_df$region))
n_bacterias = length(bacterias)
n_regions = length(regions)

# Initialize objects to store results of regression analyses
results = expand.grid(region = regions, bacteria = bacterias) %>%
  mutate(
    poisson_OD = NA, 
    negbin_OD = NA, 
    theta_OD = NA, 
    aic = NA
  ) %>%
  arrange(bacteria, region)

all_vif = vector("list", n_bacterias*n_regions)
all_estimates = data.frame()
all_residuals = data.frame()
all_fits = data.frame()
k=1

# Regression 
for (i in seq_along(bacterias)) {
  for (r in regions) {
    
    print(paste0(bacterias[i], " - ", r))
    
    # Create final database 
    final_db = regional_df %>% 
      filter(bacterie == bacterias[i], region == r) %>%
      arrange(Date_week) %>%
      mutate(
        lag1_i_res = lag(n_res / nbjh * 1000),
        lag2_nintub = lag(nintub, 2)
      ) %>% 
      filter(!is.na(lag2_nintub)) %>%
      mutate(
        lag1_i_res = (lag1_i_res - mean(lag1_i_res)) / sd(lag1_i_res),
        lag2_nintub = (lag2_nintub - mean(lag2_nintub)) / sd(lag2_nintub),
        Penicillins = (Penicillins - mean(Penicillins)) / sd(Penicillins),
        Carbapenems = (Carbapenems - mean(Carbapenems)) / sd(Carbapenems),
        Third_generation_Cephalosporins = (Third_generation_Cephalosporins - mean(Third_generation_Cephalosporins)) / sd(Third_generation_Cephalosporins)
      )
    
    # Multivariate model
    mod = best_models$model[best_models$bacteria == bacterias[i] & best_models$setting == "hospital"]
    
    if (mod == "model1") {
      m_eq = "n_res ~ lag1_i_res + p_first + p_strong_res + p_mild_res + p_no_res + Third_generation_Cephalosporins + offset(log(nbjh))"
      if (bacterias[i] == "MRSA") m_eq = "n_res ~ lag1_i_res + p_first + p_strong_res + p_mild_res + p_no_res + Penicillins + offset(log(nbjh))"
    }
    
    if (mod == "model6") {
      m_eq = "n_res ~ lag1_i_res + lag2_nintub + Carbapenems + offset(log(nbjh))"
    }
    
    m = glm.nb(m_eq, data = final_db, link = log)
    results$poisson_OD[k] = check_overdispersion(glm(m_eq, family = poisson, data = final_db))$p_value
    results$negbin_OD[k] = check_overdispersion(glm.nb(m_eq, data = final_db, link = log))$p_value
    results$theta_OD[k] = m$theta
  
    if (results$poisson_OD[k] > 0.05) {
      m = glm(m_eq, family = poisson, data = final_db)
    }
    
    # Model checks
    all_vif[[k]] = VIF(m)
    results$aic[k] = AIC(m)
    
    # Residuals
    mod_residuals = data.frame(
      region = r,
      model = mod,
      bacteria = bacterias[i],
      residuals = residuals(m)
    )
    all_residuals = bind_rows(all_residuals, mod_residuals)
    
    # Estimates and p-values
    mod_estimates = data.frame(cbind(summary(m)$coefficients, confint(m))) %>%
      rownames_to_column(var = 'variable') %>%
      rename(q2_5 = `X2.5..`, q97_5 = `X97.5..`, p = `Pr...z..`, z = `z.value`, sd = `Std..Error`) %>%
      mutate(bacteria = bacterias[i], region = r)
    all_estimates = bind_rows(all_estimates, mod_estimates)

    # Model fit
    mod_fit = data.frame(
      region = r,
      bacteria = bacterias[i],
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

# Is poisson necessary
results %>%
  filter(poisson_OD > 0.05)

##################################################
# Final plot
##################################################
# Plot estimates for Covid-19 variables by resistant
# pathogens, except for ESBL E. coli
model1_estimates = all_estimates %>%
  filter(!variable %in% c("(Intercept)", "Penicillins", "Carbapenems", "Third_generation_Cephalosporins", "lag1_i_res"),
         bacteria %in% c("ESBL E. coli", "MRSA")) %>%
  mutate(variable = factor(case_when(
    variable == "p_first" ~ "First wave w",
    variable == "p_mild_res" ~ "Mild w", 
    variable == "p_no_res" ~ "Low to none w", 
    variable == "p_strong_res" ~ "Strong w",
  ),
  c("First wave w", "Strong w", "Mild w", "Low to none w")
  ),
  region = recode(region, !!!dict_regions),
  alpha_level = ifelse(p <= 0.05, 1, 0.2)
  ) %>%
  ggplot(., aes(x = fct_rev(region), y = exp(Estimate), ymin = exp(q2_5), ymax = exp(q97_5), col = fct_rev(variable))) +
  facet_grid(cols = vars(bacteria)) +
  coord_flip() +
  geom_hline(yintercept = 1, linetype = "dashed", col = "grey70") +
  geom_pointrange(aes(alpha = alpha_level), position = position_dodge(width = 0.7), fatten = 1) + 
  scale_color_manual(values = c("First wave w" = "black", "Strong w" = "#56B4E9", "Mild w" = "#CC79A7", "Low to none w"="#E69F00")) +
  # scale_color_manual(values = c("First wave w" = "#181b4a", "Strong w" = "#313695", "Mild w" = "#4575b4", "Low to none w"="#a6bddb")) +
  scale_alpha_identity() +
  theme_bw() +
  guides(alpha = "none", col = guide_legend(reverse = T)) +
  theme(legend.position = "bottom", legend.title = element_blank()) +
  labs(x = "", y = "Incidence rate ratio\n(95% CI)")

model6_estimates = all_estimates %>%
  filter(variable == "lag2_nintub", bacteria == "CR P. aeruginosa") %>%
  mutate(
    variable = "Covid-19 intub.\nprev. w-2",
    region = recode(region, !!!dict_regions),
    alpha_level = ifelse(p <= 0.05, 1, 0.2)
  ) %>%
  ggplot(., aes(x = fct_rev(region), y = exp(Estimate), ymin = exp(q2_5), ymax = exp(q97_5), col = variable)) +
  facet_grid(cols = vars(bacteria)) +
  coord_flip() +
  geom_hline(yintercept = 1, linetype = "dashed", col = "grey70") +
  geom_pointrange(aes(alpha = alpha_level), position = position_dodge(width = 0.5), fatten = 1) +
  scale_color_manual(values = "#313695") +
  scale_alpha_identity() +
  theme_bw() +
  theme(legend.position = "right", legend.title = element_blank()) +
  guides(alpha = "none") +
  labs(x = "", y = "Incidence rate ratio\n(95% CI)") +
  ylim(c(0.55, 1.55))

# Save figure 6
figure6 = ggarrange(model1_estimates, model6_estimates, 
                    ncol = 2, 
                    labels = c("A", "B"), 
                    widths = c(1,0.6))
figure6
ggsave("../Paper/Figures/Figure6.png", figure6, height = 6, width = 8)

# Save figure 6bis
figure6bis = plot_grid(
  plot_grid(model6_estimates, ggplot()+theme_void(), ncol = 2, rel_widths = c(1, 0.2), labels = c("A", "")), 
  p_covid,
  nrow = 2
)
ggsave("../Paper/Figures/Figure6bis.png", figure6bis, height = 8, width = 6)

# Plot estimates for non Covid-19 variables 
others_estimates = all_estimates %>%
  filter(variable %in% c("Penicillins", "Carbapenems", "Third_generation_Cephalosporins", "lag1_i_res")) %>%
  mutate(variable = ifelse(variable %in% c("Penicillins", "Carbapenems", "Third_generation_Cephalosporins"),
                           "Target antibiotic", 
                           "Incidence w-1"),
  region = factor(recode(region, !!!dict_regions), rev(dict_regions))
  ) %>%
  ggplot(., aes(x = region, y = exp(Estimate), ymin = exp(q2_5), ymax = exp(q97_5), col = fct_rev(variable))) +
  facet_grid(cols = vars(bacteria)) +
  coord_flip() +
  geom_hline(yintercept = 1, linetype = "dashed", col = "grey70") +
  geom_pointrange(aes(linetype = p > 0.05), position = position_dodge(width = 0.5), fatten = 1) +
  scale_color_manual(values = c("brown1", "brown4")) +
  theme_bw() +
  theme(legend.position = "bottom") +
  guides(col = guide_legend(reverse = T), linetype = "none") +
  labs(x = "", y = "Incidence rate ratio\n(95% CI)", col = "Non-Covid-19 variables") +  
  ylim(c(0.55, 1.55))






all_estimates %>%
  filter(variable != "(Intercept)", bacteria == "MRSA") %>%
  mutate(variable = factor(case_when(
    variable == "lag_i_res" ~ "Incidence w-1",
    variable == "p_strong_res" ~ "Strong restrictions",
    variable == "p_mild_res" ~ "Mild restrictions",
    variable == "p_no_res" ~ "No restriction",
    variable == "Narrow_Penicillins" ~ "Narrow spec. P.",
    .default = variable
  ),
  c("Narrow spec. P.", "Incidence w-1", "Strong restrictions", "Mild restrictions", "No restriction")
  )) %>%
  mutate(
    region = recode(region, !!!dict_regions),
    significant = ifelse(p<=0.05, "p≤0.05", "p>0.05")
  ) %>%
  mutate(region = ifelse(region == "NOR", "NOR*", region)) %>%
  ggplot(., aes(x = variable, y = exp(Estimate), ymin = exp(q2_5), ymax = exp(q97_5), col = significant)) +
  geom_hline(yintercept = 1, linetype = "dashed", col = "grey70") +
  geom_pointrange(position = position_dodge(width = 1), fatten = 1) +
  scale_color_manual(values = c("black", "orange")) +
  facet_wrap(facets = vars(region)) +
  theme_bw() +
  theme(
    legend.position = "bottom", 
    axis.text.x = element_text(angle = 90, hjust = 1),
    legend.title = element_blank()
  ) +
  labs(x = "", y = "Incidence rate ratio (95% CI)", title = "MRSA") +
  expand_limits(y = 0)
ggsave("plots/regressions/regional/regional_estimates_mrsa.png", 
       height = 7, width = 8)

# Regional plots 

(i == 1 & r == "Hauts-de-France") | 
  (i == 2 & r == "Bourgogne-Franche-Comté") |
  (i == 4 & r == "Occitanie") |
  (i == 5 & r == "Normandie")

for (i in 1:5) {
  b= bacterias[i]
  out = regional_df %>%
    filter(bacterie == b) %>%
    ggplot(., aes(x = Date_week, y = n_res)) +
    geom_line() +
    facet_geo(facets = vars(region), grid = my_regional_grid) +
    theme_bw() +
    labs(x = "", y = "Weekly no. of resistant infections",
         title = b) 
  print(out)
}



