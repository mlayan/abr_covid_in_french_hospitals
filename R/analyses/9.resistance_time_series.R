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
  group_by(Date_week, bacterie, setting) %>%
  nest() %>%
  mutate(cfint = map(data, function(.data) getBinomCI(.data, sides = "two.sided", method = "wilson"))) %>%
  dplyr::select(-data) %>%
  unnest(cfint)

# Plot with anti-Covid-19 interventions periods
ggplot(resistance_prop, aes(x = Date_week, y = res_rate, ymin = res_rate_lwr, ymax = res_rate_upr)) +
  geom_rect(data = int_national_start_end, 
            aes(NULL, NULL, xmin=start, xmax=end, fill=restrictions, ymin=-Inf, ymax=Inf)) +
  geom_ribbon(fill = "grey10", alpha = 0.2) +
  geom_line() +
  facet_grid(rows = vars(bacterie), cols = vars(setting), scales = "free_y") +
  scale_fill_manual(
    name = "Anti-COVID-19 interventions",
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
# Description of bacterial samples (Figure 3)
##################################################
# Table with the total no. of episodes (R and not R)
tab_tot = res_hospital %>%
  mutate(Date_year = lubridate::year(Date_week)) %>%
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
  mutate(Date_year = lubridate::year(Date_week)) %>%
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
  labs(x = "", y = "Annual proportion of\nresistant isolates", fill = "")

# Cochrane-Armitage tests for trends
prop_test_trend_df = function(df) {
  df$var = (df$n_res / df$n_tot) * (1-(df$n_res / df$n_tot)) / df$n_tot
  out = prop.trend.test(df$n_res, df$n_tot)
  out = as.data.frame(out[c("statistic", "p.value")])
  rownames(out) = NULL
  s = lm(n_res/n_tot ~ Date_year, data = df) #, weights = df$var)
  out$slope = coef(s)[["Date_year"]]*100
  out$slope_lw = confint(s)["Date_year", "2.5 %"]*100
  out$slope_up = confint(s)["Date_year", "97.5 %"]*100
  return(out)
}

plot_trends = res_hospital %>%
  mutate(Date_year = lubridate::year(Date_week)) %>%
  group_by(Date_year, bacterie) %>%
  summarise(n_res = sum(n_res), n_tot = sum(n_tot), .groups = "drop") %>%
  nest(.by=bacterie) %>%
  mutate(out = map(data, prop_test_trend_df)) %>%
  dplyr::select(-data) %>%
  unnest(out) %>%
  mutate(
    yaxis = "Trend",
    p_stars = case_when(p.value > 0.05 ~ "NS",
                        p.value <= 0.05 & p.value > 0.01 ~ "*",
                        p.value <= 0.01 & p.value > 0.001 ~ "**",
                        p.value <= 0.001 ~ "***")
      ) %>%
  ggplot(., aes(x = bacterie, y = yaxis, fill = slope, label = p_stars)) +
  geom_tile() +
  geom_text(color = "gray100") +
  scale_fill_gradient2(low = "cornflowerblue", mid = "white", high = "red") +
  theme_bw() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1),
        axis.title.x = element_blank(),
        legend.position = "bottom",
        plot.margin = margin(5.5, 5.5, 5.5, 10, "pt")) +
  labs(y = "", fill = "Temporal trend (%)")

# Temporal dynamics of weekly incidence rates 
bd_all = bind_rows(
  bd_pmsi_hospital %>% mutate(setting = "Hospital"),
  bd_pmsi_icu %>% mutate(setting = "ICU")
)

plot_res_i = bind_rows(
    res_hospital %>% mutate(setting = "Hospital"),
    res_icu %>% mutate(setting = "ICU"),
  ) %>%
  left_join(., bd_all, by = c("Date_week", "setting")) %>%
  mutate(incidence = n_res / nbjh * 1000) %>%
  ggplot(., aes(x = Date_week, y = incidence)) +
  geom_rect(data = int_national_start_end, 
            aes(NULL, NULL, xmin=start, xmax=end, fill=restrictions, ymin=-Inf, ymax=Inf)) +
  geom_line() +
  ggh4x::facet_grid2(cols = vars(setting), rows = vars(bacterie), scales = "free_y", independent = "y") +
  theme_bw() +
  theme(legend.position = "bottom", legend.title = element_text(hjust = 0.5)) +
  scale_fill_manual(
    name = "Anti-COVID-19 interventions",
    labels = c("First wave", "Strong", "Intermediary", "Light to none"),
    breaks = c("p_first_wave", "p_strong_res", "p_mild_res", "p_no_res"),
    values = c("p_first_wave" = col_interventions(1), "p_strong_res" = col_interventions(2), "p_mild_res" = col_interventions(3), "p_no_res" = col_interventions(4))
    ) +
  guides(fill = guide_legend(title.position = "top", order=1, override.aes = list(col = 'black'))) +
  expand_limits(y = 0) +
  labs(x = "", y = "Weekly incidence of resistant isolates (for 1,000 bed-days)")
ggsave("plots/antibiotic_resistance/national_weekly_incidence.png", plot_res_i, height = 6, width = 12)

# Final figure
figure5 = ggarrange(
  ggarrange(tab_tot, plot_res_p, plot_trends, nrow = 3, heights = c(0.5, 1, 0.7), labels = c("A", "B", "C")), 
  plot_res_i,
  ncol = 2,
  labels = c("", 'D'),
  widths = c(0.6,1)
)
figure5
ggsave("plots/Figure5.png", figure5, height = 8, width = 12)
ggsave("../Paper/Figures/Figure5.png", figure5, height = 8, width = 12)

##################################################
# Comparison with annual incidence rate reported
# by SPARES
##################################################
# Load SPARES annual incidence data
report = read_excel("data-raw/spares/incidence_density_from_reports.xlsx") %>%
  mutate(bacterie = case_when(
    bacterie == "Enterobacter cloacae complex" ~ "ESBL E. cloacae",
    bacterie == "Escherichia coli" ~ "ESBL E. coli",
    bacterie == "Klebsiella pneumoniae" ~ "ESBL K. pneumoniae",
    bacterie == "Staphylococcus aureus" ~ "MRSA"
  )
  ) %>%
  rename(SPARES = di_manual, Date_year = qrt)

# Annual number of hospitalization days
bd_annual = bd_pmsi_hospital %>%
  mutate(Date_year = lubridate::year(Date_week)) %>%
  group_by(Date_year) %>%
  summarise(nbjh = sum(nbjh), .groups = "drop")

# Comparison with SPARES
res_hospital %>%
  mutate(Date_year = lubridate::year(Date_week)) %>%
  group_by(bacterie, Date_year) %>%
  summarise(n_res = sum(n_res), .groups = "drop") %>%
  left_join(., bd_annual, by = "Date_year") %>%
  mutate(Cohort = n_res / nbjh * 1000) %>%
  left_join(., report[, c("Date_year", "SPARES", "bacterie")], by = c("Date_year", "bacterie")) %>%
  filter(!is.na(SPARES)) %>%
  pivot_longer(c(Cohort, SPARES), names_to = "database", values_to = "incidence") %>%
  ggplot(., aes(x = factor(Date_year), y = incidence, fill = database)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.5) +
  facet_wrap(facets = vars(bacterie), ncol = 2, scales = "free_y") +
  theme_bw() +
  labs(x = "", y = "Incidence rate (per 1,000 hospitalisation days)",
       fill = "")
ggsave("plots/antibiotic_resistance/comparison_incidence_spares_report.png", height = 5, width = 8)

##################################################
# Plot regional level time series 
##################################################
# Incidence 
for (b in unique(res_hospital$bacterie)) {
  p = res_regions %>%
    filter(bacterie == b) %>%
    left_join(., bd_pmsi_regions, by = c("Date_week", "region")) %>%
    mutate(res_incidence = n_res / nbjh * 1000, region = recode(region, !!!dict_regions)) %>%
    ggplot(., aes(x = Date_week, y = res_incidence)) +
    geom_line() +
    facet_wrap(facets = vars(region), ncol = 4) +
    theme_bw() + 
    labs(x = "", y = paste0("Weekly incidence of ", b, " (for 1,000 bed-days)"))
  ggsave(paste0("plots/antibiotic_resistance/", gsub(" ", "", b), "_region.png"), 
         p, height = 6, width = 11)
}


