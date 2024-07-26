##################################################
## PLOTS OF ANTIBIOTIC RESISTANCE DATA
##################################################
rm(list = ls())
library(tidyverse)
library(ggpubr)
library(DescTools)
library(readxl)

load("data/int_national_start_end.rda")

load("data/bd_pmsi_hospital.rda")
load("data/bd_pmsi_icu.rda")
load("data/bd_pmsi_regions.rda")

load("data/res_hospital.rda")
load("data/res_icu.rda")
load("data/res_regions.rda")

source("R/helper/dictionaries.R")
source("R/helper/helper_functions.R")
source("R/helper/helper_plots.R")

##################################################
# Plot weekly resistance rates 
##################################################
# Resistance rate at the hospital and ICU levels
resistance_prop = bind_rows(
  res_icu %>% mutate(setting = "ICU"),
  res_hospital %>% mutate(setting = "Hospital")
  ) %>%
  group_by(Date_year, Date_week, bacterie, setting) %>%
  nest() %>%
  mutate(cfint = map(data, function(.data) getBinomCI(.data, sides = "two.sided", method = "wilson"))) %>%
  dplyr::select(-data) %>%
  unnest(cfint)

# Plot with anti-Covid-19 interventions periods
ggplot(resistance_prop, aes(x = as.Date(Date_week), y = res_rate, ymin = res_rate_lwr, ymax = res_rate_upr)) +
  geom_rect(data = int_national_start_end, 
            aes(NULL, NULL, xmin=start, xmax=end, fill=restrictions, ymin=-Inf, ymax=Inf)) +
  geom_ribbon(fill = "grey10", alpha = 0.2) +
  geom_line() +
  facet_grid(rows = vars(bacterie), cols = vars(setting), scales = "free_y") +
  scale_fill_manual(
    name = "Anti-Covid-19 interventions",
    labels = c("First wave", "Strong", "Intermediary", "Light to none"),
    breaks = c("p_first_wave", "p_strong_res", "p_mild_res", "p_no_res"),
    values = c("p_first_wave" = col_interventions(1), "p_strong_res" = col_interventions(2), "p_mild_res" = col_interventions(3), "p_no_res" = col_interventions(4))
  ) +
  theme_bw() +
  theme(legend.key = element_rect(colour = "black"), legend.position = "bottom") +
  guides(fill = guide_legend(title.position = "top", title.hjust = 0.5)) +
  expand_limits(y=0) +
  labs(x = "", y = "Weekly resistance proportion (95% CI)")
ggsave("../Paper/Supplementary/weekly_resistance_proportion.png", height = 9, width = 8)
ggsave("plots/antibiotic_resistance/weekly_resistance_proportion.png", height = 9, width = 8)

##################################################
# Description of bacterial samples (Figure 2)
##################################################
# Table with the total no. of episodes (R and not R)
tab_tot = res_hospital %>%
  group_by(Date_year, bacterie) %>%
  summarise(n_tot = sum(n_tot), .groups = "drop") %>%
  mutate(bacterie = case_when(
    bacterie == "CR P. aeruginosa" ~ "P. aeruginosa",
    bacterie == "ESBL E. cloacae" ~ "E. cloacae",
    bacterie == "ESBL E. coli" ~ "E. coli",
    bacterie == "ESBL K. pneumoniae" ~ "K. pneumoniae",
    bacterie == "MRSA" ~ "S. aureus"
  )) %>%
  pivot_wider(names_from = bacterie, values_from = n_tot) %>%
  rename(` ` = Date_year) %>%
  ggtexttable(., rows = NULL, theme = ttheme("light",base_size = 9))

# Plot of annual proportion of resistant infections 
plot_res_p = res_hospital %>%
  group_by(Date_year, bacterie) %>%
  summarise(n_res = sum(n_res), n_tot = sum(n_tot), .groups = "drop") %>%
  group_by(Date_year, bacterie) %>%
  nest() %>%
  mutate(cfint = map(data, function(.data) getBinomCI(.data, sides = "two.sided", method = "wilson"))) %>%
  unnest(cols = c(data, cfint)) %>%
  ggplot(., aes(x = bacterie, y = res_rate, ymin = res_rate_lwr, ymax = res_rate_upr, fill = factor(Date_year))) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.6), width = 0.5) +
  geom_linerange(position = position_dodge(width = 0.6)) +
  scale_fill_manual(values = colorRampPalette(c("dodgerblue3", "lightblue1"))(4)) +
  theme_bw() +
  theme(legend.position = "right", axis.text.x = element_text(angle = 30, hjust = 1)) +
  labs(x = "", y = "Annual proportion of\nresistant infections", fill = "")

# Temporal dynamics of weekly incidence rates 
bd_all = bind_rows(
  bd_pmsi_hospital %>% mutate(setting = "Hospital"),
  bd_pmsi_icu %>% mutate(setting = "ICU")
)

plot_res_i = bind_rows(
    res_hospital %>% mutate(setting = "Hospital"),
    res_icu %>% mutate(setting = "ICU"),
  ) %>%
  group_by(Date_week, setting, bacterie) %>%
  summarise(n_res = sum(n_res), .groups = "drop") %>%
  left_join(., bd_all, by = c("Date_week", "setting")) %>%
  mutate(incidence = n_res / nbjh * 1000) %>%
  ggplot(., aes(x = as.Date(Date_week), y = incidence)) +
  geom_rect(data = int_national_start_end, 
            aes(NULL, NULL, xmin=start, xmax=end, fill=restrictions, ymin=-Inf, ymax=Inf)) +
  geom_line() +
  ggh4x::facet_grid2(cols = vars(setting), rows = vars(bacterie), scales = "free_y", independent = "y") +
  theme_bw() +
  theme(legend.position = "bottom", legend.title = element_text(hjust = 0.5)) +
  scale_fill_manual(
    name = "Anti-Covid-19 interventions",
    labels = c("First wave", "Strong", "Intermediary", "Light to none"),
    breaks = c("p_first_wave", "p_strong_res", "p_mild_res", "p_no_res"),
    values = c("p_first_wave" = col_interventions(1), "p_strong_res" = col_interventions(2), "p_mild_res" = col_interventions(3), "p_no_res" = col_interventions(4))
    ) +
  guides(fill = guide_legend(title.position = "top", order=1, override.aes = list(col = 'black'))) +
  expand_limits(y = 0) +
  labs(x = "", y = "Weekly incidence of resistant infections (per 1,000 bed-days)")
ggsave("plots/antibiotic_resistance/national_weekly_incidence.png", plot_res_i, height = 6, width = 12)

# Final figure
figure2 = ggarrange(
  ggarrange(tab_tot, plot_res_p, nrow = 2, heights = c(1, 0.8), labels = c("A", "B")), 
  plot_res_i,
  ncol = 2,
  labels = c("", 'C'),
  widths = c(0.6,1)
)
figure2
ggsave("plots/Figure2.png", figure2, height = 8, width = 12)
ggsave("../Paper/Figures/Figure2.png", figure2, height = 8, width = 12)

##################################################
# Paper figure 3 - Bed days - Total and Covid-19
##################################################
# Total bed-days
hd_icu = hd %>%
  filter(secteur == "Réanimation") %>%
  group_by(Date_week) %>%
  summarise(nbjh = sum(nbjh), .groups = "drop") %>%
  mutate(secteur = "ICU")

plot_bd = hd %>%
  group_by(Date_week) %>%
  summarise(nbjh = sum(nbjh), .groups = "drop") %>%
  mutate(secteur = "Hospital") %>%
  bind_rows(., hd_icu) %>%
  filter(!Date_week %in% c(as.Date("2018-12-31"), as.Date("2022-12-26"))) %>%
  ggplot(., aes(x = as.Date(Date_week), y = nbjh)) +
  geom_rect(data = int_national_start_end, 
            aes(NULL, NULL, xmin=start, xmax=end, fill=restrictions, ymin=-Inf, ymax=Inf)) +
  geom_line() +
  facet_wrap(facets = vars(secteur), scales = "free_y") +
  theme_bw() +
  theme(legend.position = "bottom",  
        legend.key = element_rect(colour = "black"),
        legend.title.align = 0.5) +
  guides(fill = guide_legend(title.position = "top")) +
  scale_fill_manual(
    name = "Anti-Covid-19 interventions",
    labels = c("First wave", "Strong", "Intermediary", "Light to none"),
    breaks = c("p_first", "p_strong_res", "p_mild_res", "p_no_res"),
    values = c("p_first" = col_interventions(1), "p_strong_res" = col_interventions(2), "p_mild_res" = col_interventions(3), "p_no_res" = col_interventions(4))
  ) +
  scale_y_continuous(labels = scales::label_comma()) +
  expand_limits(y = 0) +
  labs(x = "", y = "Weekly bed-days")

# Covid-19 bed-days
to_add = data.frame(
  Date_week = seq(as.Date("2019-01-07"), as.Date("2019-12-02"), by = 7),
  nintub = 0
)

covid_icu = haven::read_sas("data-raw/atih/mco_intubation_weekly.sas7bdat") %>%
  rename(finess = FinessGeo) %>%
  left_join(., metadata_admin_espic %>% dplyr::select(finess, region) %>% distinct(), by = "finess") %>%
  filter(!is.na(region), status == "covid", secteur == "Réanimation") %>%
  group_by(Date_week) %>%
  summarise(nintub = sum(nbjh), .groups = "drop") %>%
  bind_rows(., to_add) %>%
  mutate(secteur = "ICU")


plot_covid = haven::read_sas("data-raw/atih/mco_intubation_weekly.sas7bdat") %>%
  rename(finess = FinessGeo) %>%
  left_join(., metadata_admin_espic %>% dplyr::select(finess, region) %>% distinct(), by = "finess") %>%
  filter(!is.na(region), status == "covid") %>%
  group_by(Date_week) %>%
  summarise(nintub = sum(nbjh), .groups = "drop") %>%
  bind_rows(., to_add) %>%
  mutate(secteur = "Hospital") %>%
  bind_rows(., covid_icu) %>%
  filter(Date_week != "2022-12-26") %>%
  ggplot(., aes(x = as.Date(Date_week), y = nintub)) +
  geom_rect(data = int_national_start_end, 
            aes(NULL, NULL, xmin=start, xmax=end, fill=restrictions, ymin=-Inf, ymax=Inf)) +
  geom_line() +
  facet_wrap(facets = vars(secteur), scales = "free_y") +
  theme_bw() +
  theme(legend.position = "bottom",  
        legend.key = element_rect(colour = "black"),
        legend.title.align = 0.5) +
  guides(fill = guide_legend(title.position = "top")) +
  scale_fill_manual(
    name = "Anti-Covid-19 interventions",
    labels = c("First wave", "Strong", "Intermediary", "Light to none"),
    breaks = c("p_first", "p_strong_res", "p_mild_res", "p_no_res"),
    values = c("p_first" = col_interventions(1), "p_strong_res" = col_interventions(2), "p_mild_res" = col_interventions(3), "p_no_res" = col_interventions(4))
  ) +
  expand_limits(y = 0) +
  scale_y_continuous(labels = scales::label_comma()) +
  labs(x = "", y = "Weekly intubated Covid-19 bed-days")

# Final figure
figure3 = ggarrange(plot_bd, plot_covid, nrow = 2, 
                    labels = c("A", "B"), hjust = -1.5,
                    common.legend = T, legend = "bottom")
figure3
ggsave("../Paper/Figures/Figure3.png", figure3, height = 5, width = 8)


##################################################
# Comparison with annual incidence rate reported
# by SPARES
##################################################
# Load data
report = read_excel("data-raw/spares/reports/incidence_density_from_reports.xlsx") %>%
  mutate(bacterie = case_when(
    bacterie == "Enterobacter cloacae complex" ~ "ESBL E. cloacae",
    bacterie == "Escherichia coli" ~ "ESBL E. coli",
    bacterie == "Klebsiella pneumoniae" ~ "ESBL K. pneumoniae",
    bacterie == "Staphylococcus aureus" ~ "MRSA"
  )
  ) %>%
  rename(SPARES = di_manual, Date_year = qrt)

# Annual number of hospitalisation days
hd_annual_secteur = hd %>%
  mutate(Date_year = ifelse(year(Date_week) == 2018, 2019, year(Date_week))) %>%
  filter(!(secteur == "Réanimation" & !code %in% icu_cohort_final)) %>%
  group_by(Date_year, secteur) %>%
  summarise(nbjh = sum(nbjh), .groups = "drop")

hd_annual_etab = hd %>%
  mutate(Date_year = ifelse(year(Date_week) == 2018, 2019, year(Date_week))) %>%
  group_by(Date_year, code) %>%
  summarise(nbjh = sum(nbjh), .groups = "drop")

hd_annual = hd %>%
  mutate(Date_year = ifelse(year(Date_week) == 2018, 2019, year(Date_week))) %>%
  group_by(Date_year) %>%
  summarise(nbjh = sum(nbjh), .groups = "drop")

# Annual incidence 
res %>%
  group_by(Date_year, secteur, bacterie) %>%
  summarise(n_res = sum(n_res), .groups = "drop") %>%
  left_join(., hd_annual_secteur, by = c("Date_year", "secteur")) %>%
  mutate(incidence = n_res / nbjh * 1000) %>%
  ggplot(., aes(x = factor(Date_year), y = incidence, fill = secteur)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.5) +
  facet_wrap(facets = vars(bacterie), ncol = 3, scales = "free_y") +
  theme_bw() +
  labs(x = "", y = "Incidence rate (per 1,000 hospitalisation days)",
       fill = "")
ggsave("plots/resistances/incidence_secteur.png", height = 8, width = 12)

# Comparison with SPARES
res %>%
  filter(bacterie %in% c("ESBL E. cloacae", "ESBL E. coli",
                         "ESBL K. pneumoniae", "MRSA")) %>%
  group_by(bacterie, Date_year) %>%
  summarise(n_res = sum(n_res), .groups = "drop") %>%
  left_join(., hd_annual, by = "Date_year") %>%
  mutate(`AMR COVID` = n_res / nbjh * 1000) %>%
  left_join(., report[, c("Date_year", "SPARES", "bacterie")], by = c("Date_year", "bacterie")) %>%
  pivot_longer(c(`AMR COVID`, SPARES), names_to = "database", values_to = "incidence") %>%
  ggplot(., aes(x = factor(Date_year), y = incidence, fill = database)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.5) +
  facet_wrap(facets = vars(bacterie), ncol = 2, scales = "free_y") +
  theme_bw() +
  labs(x = "", y = "Incidence rate (per 1,000 hospitalisation days)",
       fill = "")
ggsave("plots/resistances/comparison_incidence_spares_report.png", 
       height = 5, width = 8)

# Weekly number of hospitalisation days
hd_weekly_etab = hd %>%
  group_by(Date_week, code, region) %>%
  summarise(nbjh = sum(nbjh), .groups = "drop")

# Names of hospitals
hosp_names = metadata_admin_espic %>% 
  filter(code %in% cohort_final) %>% 
  dplyr::select(code, name) %>% 
  distinct() 

##################################################
# Plot national level time series 
##################################################
# Weekly bed-days - Hospital
weekly_bd = hd %>%
  group_by(Date_week) %>%
  summarise(nbjh = sum(nbjh), .groups = "drop") 

# Weekly bed-days ICUs
weekly_bd_icu = hd %>%
  filter(secteur == "Réanimation") %>%
  group_by(Date_week) %>%
  summarise(nbjh = sum(nbjh), .groups = "drop") 

# Function to plot counts, incidence and resistance percentage
plot_national_ts = function(df, plot_name, df_bd) {
  
  out = df %>%
    group_by(Date_week, bacterie) %>%
    summarise(n_res = sum(n_res), n_tot = sum(n_tot), .groups = "drop") %>%
    left_join(., df_bd, by = c("Date_week")) %>%
    mutate(
      `Weekly incidence of infections\n(per 1,000 bed-days)` = n_tot / nbjh * 1000,
      `Weekly incidence of resistant infections\n(per 1,000 bed-days)`= n_res / nbjh * 1000,
      `Weekly resistance percentage` = n_res / n_tot * 100
    ) %>% 
    rename(
      `Weekly number of infections` = n_tot,
      `Weekly number of resistant infections` = n_res,
      `Weekly number of bed-days` = nbjh
    ) %>%
    pivot_longer(matches("^Weekly"), names_to = "metrics", values_to = "values") %>%
    mutate(metrics = factor(metrics, c(
      "Weekly number of bed-days",
      "Weekly number of infections", 
      "Weekly number of resistant infections", 
      "Weekly incidence of infections\n(per 1,000 bed-days)",
      "Weekly incidence of resistant infections\n(per 1,000 bed-days)",
      "Weekly resistance percentage"
    ))) %>%
    ggplot(., aes(x = as.Date(Date_week), y = values)) +
    annotate("rect", xmin = as.Date("2020-03-16"), xmax = as.Date("2020-05-10"), ymin = 0, ymax = Inf, fill = "gray90") +
    annotate("rect", xmin = as.Date("2020-10-30"), xmax = as.Date("2020-12-15"), ymin = 0, ymax = Inf, fill = "gray90") +
    geom_line() +
    ggh4x::facet_grid2(cols = vars(metrics), rows = vars(bacterie), scales = "free_y", 
                       independent = "y") +
    geom_hline(yintercept = 0, col = "#00000000") +
    theme_bw() +
    labs(x = "", y = "")
  
  print(out)
  
  ggsave(paste0("plots/resistances/timeseries_", gsub(" ", "_", plot_name), ".png"), 
         out, height = 12, width = 18)  
}

# Number of resistances, number of tests, resistance rates - All infections 
plot_national_ts(res, "all", weekly_bd)

# Number of resistances, number of tests, resistance rates - Blood-sourced infections
plot_national_ts(res[res$site == "blood-sourced infection", ], "blood_sourced", weekly_bd)

# Number of resistances, number of tests, resistance rates - Respiratory infections
plot_national_ts(res[res$site == "lower respiratory tract infection", ], "lrti", weekly_bd)

# Number of resistances, number of tests, resistance rates - CHU
plot_national_ts(res %>% 
                   left_join(., metadata_admin_espic %>% dplyr::select(code, type) %>% distinct(), by = "code") %>% 
                   filter(type %in% c("CHU", "ESPIC")), 
                 "chu-espic", 
                 weekly_bd)

plot_national_ts(res %>% 
                   left_join(., metadata_admin_espic %>% dplyr::select(code, type) %>% distinct(), by = "code") %>% 
                   filter(type == "CHU"), 
                 "chu", 
                 weekly_bd)

plot_national_ts(res %>% 
                   left_join(., metadata_admin_espic %>% dplyr::select(code, type) %>% distinct(), by = "code") %>% 
                   filter(type == "MCO"), 
                 "mco", 
                 weekly_bd)

plot_national_ts(res %>% 
                   left_join(., metadata_admin_espic %>% dplyr::select(code, type) %>% distinct(), by = "code") %>% 
                   filter(type == "CH"), 
                 "ch", 
                 weekly_bd)

plot_national_ts(res %>% 
                   left_join(., metadata_admin_espic %>% dplyr::select(code, type) %>% distinct(), by = "code") %>% 
                   filter(type == "ESSR"), 
                 "essr", 
                 weekly_bd)

# Number of resistances, number of tests, resistance rates - Blood-sourced infections in ICUs
plot_national_ts(res[res$site == "blood-sourced infection" & res$secteur == "Réanimation", ], 
                 "blood_sourced_icu", weekly_bd_icu)

# Number of resistances, number of tests, resistance rates - Respiratory infections in ICUs
plot_national_ts(res[res$site == "lower respiratory tract infection" & res$secteur == "Réanimation", ], 
                 "lrti_icu", weekly_bd_icu)

# Number of resistances, number of tests, resistance rates - All infections in ICUs
plot_national_ts(res[res$secteur == "Réanimation", ], "all_icu", weekly_bd_icu)

##################################################
# Plot regional level time series 
##################################################
# Weekly bed-days - Hospital
weekly_regional_bd = hd %>%
  group_by(Date_week, region) %>%
  summarise(nbjh = sum(nbjh), .groups = "drop") 

# Weekly bed-days ICUs
weekly_regional_bd_icu = hd %>%
  filter(secteur == "Réanimation") %>%
  group_by(Date_week, region) %>%
  summarise(nbjh = sum(nbjh), .groups = "drop") 

# Function to plot counts, incidence and resistance percentage
plot_regional_ts = function(df, b, plot_name, df_bd, df_region = metadata_admin_espic,
                            hexagonal_france = my_regional_grid) {
  
  out = df %>%
    filter(bacterie == b) %>%
    left_join(., df_region %>% select(code, region) %>% distinct(), by = "code") %>%
    group_by(Date_week, region, bacterie) %>%
    summarise(n_res = sum(n_res), n_tot = sum(n_tot), .groups = "drop") %>%
    left_join(., df_bd, by = c("Date_week", "region")) %>%
    mutate(
      incidence_total = n_tot / nbjh * 1000,
      incidence_resistance = n_res / nbjh * 1000,
      resistance_percentage = n_res / n_tot * 100
    )
    
  p1 = ggplot(out, aes(x = as.Date(Date_week), y = incidence_total)) +
    annotate("rect", xmin = as.Date("2020-03-16"), xmax = as.Date("2020-05-10"), ymin = 0, ymax = Inf, fill = "gray90") +
    annotate("rect", xmin = as.Date("2020-10-30"), xmax = as.Date("2020-12-15"), ymin = 0, ymax = Inf, fill = "gray90") +
    geom_line() +
    facet_geo(facets = vars(region), grid = hexagonal_france) +
    geom_hline(yintercept = 0, col = "#00000000") +
    theme_bw() +
    labs(x = "", y = "", title = paste0("Weekly incidence rate of infections - ", b))
  ggsave(paste0("plots/resistances/timeseries_regional_I_", gsub(" ", "", tolower(b)), "_", gsub(" ", "_", plot_name), ".png"), 
         p1, height = 8, width = 10)  
  
  p2 = ggplot(out, aes(x = as.Date(Date_week), y = incidence_resistance)) +
    annotate("rect", xmin = as.Date("2020-03-16"), xmax = as.Date("2020-05-10"), ymin = 0, ymax = Inf, fill = "gray90") +
    annotate("rect", xmin = as.Date("2020-10-30"), xmax = as.Date("2020-12-15"), ymin = 0, ymax = Inf, fill = "gray90") +
    geom_line() +
    facet_geo(facets = vars(region), grid = hexagonal_france) +
    geom_hline(yintercept = 0, col = "#00000000") +
    theme_bw() +
    labs(x = "", y = "", title = paste0("Weekly incidence rate of resistant infections - ", b))
  ggsave(paste0("plots/resistances/timeseries_regional_RI_", gsub(" ", "", tolower(b)), "_", gsub(" ", "_", plot_name), ".png"), 
         p2, height = 8, width = 10)  
  
  p3 = ggplot(out, aes(x = as.Date(Date_week), y = resistance_percentage)) +
    annotate("rect", xmin = as.Date("2020-03-16"), xmax = as.Date("2020-05-10"), ymin = 0, ymax = Inf, fill = "gray90") +
    annotate("rect", xmin = as.Date("2020-10-30"), xmax = as.Date("2020-12-15"), ymin = 0, ymax = Inf, fill = "gray90") +
    geom_line() +
    facet_geo(facets = vars(region), grid = hexagonal_france) +
    geom_hline(yintercept = 0, col = "#00000000") +
    theme_bw() +
    labs(x = "", y = "", title = paste0("Weekly resistance percentage - ", b))
  ggsave(paste0("plots/resistances/timeseries_regional_RP_", gsub(" ", "", tolower(b)), "_", gsub(" ", "_", plot_name), ".png"), 
         p3, height = 8, width = 10)  
  
}

for (b in c("CR P. aeruginosa", "ESBL E. cloacae", "ESBL K. pneumoniae")) {
  # Number of resistances, number of tests, resistance rates - Blood-sourced infections
  plot_regional_ts(res[res$site == "blood-sourced infection", ],
                   b, "blood_sourced", weekly_regional_bd)
  
  # Number of resistances, number of tests, resistance rates - Respiratory infections
  plot_regional_ts(res[res$site == "lower respiratory tract infection", ], 
                   b, "lrti", weekly_regional_bd)
  
  # Number of resistances, number of tests, resistance rates - All infections 
  plot_regional_ts(res, b, "all", weekly_regional_bd)
}

##################################################
# Plot resistance rate time series at the national
# level
##################################################
# Weekly resistance rate
res %>%
  group_by(Date_week, bacterie) %>%
  summarise(
    n_res = sum(n_res),
    n_tot = sum(n_tot),
    .groups = "drop"
  ) %>%
  complete(Date_week, bacterie, fill = list(n_tot = 0 , n_res = 0)) %>%
  group_by(Date_week, bacterie) %>%
  nest() %>%
  mutate(cfint = map(data, function(.data) getBinomCI(.data, sides = "two.sided", method = "wilson"))) %>%
  unnest(cols = c(data, cfint)) %>%
  ungroup() %>%
  ggplot(., aes(x = as.Date(Date_week), y = res_rate, ymin = res_rate_lwr, ymax = res_rate_upr, 
                col = bacterie, fill = bacterie)) +
  facet_wrap(facets = vars(bacterie), ncol = 3, scales = "free_y") +
  annotate("rect", xmin = as.Date("2020-03-16"), xmax = as.Date("2020-05-10"), ymin = 0, ymax = Inf, fill = "gray90") +
  annotate("rect", xmin = as.Date("2020-10-30"), xmax = as.Date("2020-12-15"), ymin = 0, ymax = Inf, fill = "gray90") +
  geom_line() +
  geom_ribbon(alpha = 0.2, col = "#00000000") +
  theme_bw() +
  theme(legend.position = "none") +
  labs(x = "", y = "Weekly resistance rate (95% CI)", col = "")
ggsave("plots/resistances/weekly_rrate_national.png", height = 7, width = 8)


# Weekly resistance rate in ICUs
res %>%
  filter(secteur == "Réanimation", code %in% icu_cohort_final) %>%
  group_by(Date_week, bacterie) %>%
  summarise(
    n_res = sum(n_res),
    n_tot = sum(n_tot),
    .groups = "drop"
  ) %>%
  complete(Date_week, bacterie, fill = list(n_tot = 0 , n_res = 0)) %>%
  group_by(Date_week, bacterie) %>%
  nest() %>%
  mutate(cfint = map(data, function(.data) getBinomCI(.data, sides = "two.sided", method = "wilson"))) %>%
  unnest(cols = c(data, cfint)) %>%
  ungroup() %>%
  ggplot(., aes(x = as.Date(Date_week), y = res_rate, ymin = res_rate_lwr, ymax = res_rate_upr, 
                col = bacterie, fill = bacterie)) +
  facet_wrap(facets = vars(bacterie), ncol = 3, scales = "free_y") +
  annotate("rect", xmin = as.Date("2020-03-16"), xmax = as.Date("2020-05-10"), ymin = 0, ymax = Inf, fill = "gray90") +
  annotate("rect", xmin = as.Date("2020-10-30"), xmax = as.Date("2020-12-15"), ymin = 0, ymax = Inf, fill = "gray90") +
  geom_line() +
  geom_ribbon(alpha = 0.2, col = "#00000000") +
  theme_bw() +
  theme(legend.position = "none") +
  labs(x = "", y = "Weekly resistance rate (95% CI)", col = "")
ggsave("plots/resistances/weekly_rrate_national_icu.png", height = 7, width = 8)


