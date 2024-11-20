##################################################
# ANALYSES OF ANTIBIOTIC USE IN FRENCH HOSPITALS
##################################################
rm(list = ls())
library(tidyverse)
library(gt)
library(rstatix)
library(DescTools)
library(readxl)
library(ggpubr)
library(cowplot)
library(sf)

source("R/helper/dictionaries.R")
source("R/helper/helper_functions.R")
source("R/helper/helper_plots.R")

##################################################
# Load data
##################################################
# Hospital metadata 
load("data/metadata_admin_espic.rda")

# Identifiers of hospital and ICU cohorts
load("data/cohort_final.rda")
load("data/icu_cohort_final.rda")

# Antibiotic use data at the hospital level
load("data/atb_use_hospital_level.rda")

##################################################
# Plot raw data - Figure 3
##################################################
atb_class_names = c("Aminoglycosides", "Carbapenems", "Cephalosporins", "Fosfomycin", 
                    "Glycopeptides", "Lipopeptides", "Macrolides", "Monobactams", 
                    "Oxazolidinones", "Penicillins", "Polymyxins", "Quinolones",
                    "Tetracyclines", "Trimethoprim", "Total")

# Median consumption in Hospitals by antibiotic class
median_hospital = atb_use_hospital_level %>%
  filter(atb_class %in% atb_class_names) %>%
  group_by(code, Date_year, atb_class) %>%
  summarise(consumption = sum(molDDD) / sum(Nbhosp) * 1000, .groups = "drop") %>%
  group_by(Date_year, atb_class) %>%
  summarise(med = median(consumption), code = NA, .groups = "drop") %>%
  mutate(atb_class = factor(atb_class, atb_class_names))

# Hospital level consumption 
consumption_hospital = atb_use_hospital_level %>%
  filter(atb_class %in% atb_class_names) %>%
  group_by(code, Date_year, atb_class) %>%
  summarise(consumption = sum(molDDD) / sum(Nbhosp) * 1000, .groups = "drop") %>%
  mutate(atb_class = factor(atb_class, atb_class_names))

# Trajectories at the hospital level
p1 = ggplot(consumption_hospital, aes(x = Date_year, y = consumption, group = code)) +
  geom_line(col = "grey", linewidth = 0.1) +
  geom_point(data = median_hospital, aes(x = Date_year, y = med), col = "red") +
  facet_wrap(facets = vars(atb_class), scales = "free_y", ncol = 5) +
  theme_bw() +
  # scale_y_continuous(trans = scales::pseudo_log_trans()) +
  labs(x = "", y = "Annual antibiotic consumption (in DDD/1,000 bed-days)",
       title = "Hospitals")

# Median consumption in ICUs by antibiotic class
median_icu = atb_use_hospital_level %>%
  filter(secteur == "ICU", code %in% icu_cohort_final, atb_class %in% atb_class_names) %>%
  mutate(consumption = molDDD / Nbhosp * 1000) %>%
  group_by(Date_year, atb_class) %>%
  summarise(med = median(consumption), code = NA, .groups = "drop") %>%
  mutate(atb_class = factor(atb_class, atb_class_names))

# Hospital level consumption 
consumption_icu = atb_use_hospital_level %>%
  filter(secteur == "ICU", code %in% icu_cohort_final, atb_class %in% atb_class_names) %>%
  mutate(consumption = molDDD / Nbhosp * 1000) %>%
  mutate(atb_class = factor(atb_class, atb_class_names))

# Trajectories at the ICU level
p2 = ggplot(consumption_icu, aes(x = Date_year, y = consumption, group = code)) +
  geom_line(col = "grey", linewidth = 0.1) +
  geom_point(data = median_icu, aes(x = Date_year, y = med), col = "red") +
  facet_wrap(facets = vars(atb_class), scales = "free_y", ncol = 5) +
  theme_bw() +
  # scale_y_continuous(trans = scales::pseudo_log_trans()) +
  labs(x = "", y = "Annual antibiotic consumption (in DDD/1,000 bed-days)",
       title = "ICUs")

# Final figure 
figure3 = ggarrange(p1, p2, nrow = 2, labels = c("A", "B"))
figure3
ggsave("../Paper/Figures/Figure3.png", figure3, height = 10, width = 10)
ggsave("plots/Figure3.png", figure3, height = 10, width = 10)

# All hospitals reported antibiotic consumption by class every year
consumption_hospital %>%
  group_by(code, atb_class) %>%
  mutate(n = n()) %>% 
  filter(n<4) 

# ICUs that did not report at least one antibiotic class at least one year
consumption_icu %>%
  group_by(code, atb_class) %>%
  summarise(n = n(), .groups = "drop") %>% 
  filter(n<4) 

##################################################
# Table 1 of national consumption of antibiotics
# by molecular class 
##################################################
atb_class_names2 = c("Aminoglycosides", "Carbapenems", "Imipenem + Meropenem", "Cephalosporins", 
                     "Fosfomycin", "Glycopeptides", "Vancomycin", "Lipopeptides", "Macrolides", 
                     "Azithromycin", "Monobactams", "Oxazolidinones", "Penicillins", "Polymyxins", 
                     "Quinolones", "Tetracyclines", "Trimethoprim", "Total")

# Hospital consumption
atb_use_hospital_level_agg = atb_use_hospital_level %>%
  filter(atb_class %in% atb_class_names2) %>%
  group_by(code, Date_year, atb_class) %>%
  summarise(molDDD = sum(molDDD), Nbhosp = sum(Nbhosp), .groups = "drop") %>%
  mutate(atb_class = factor(atb_class, atb_class_names2))

national_2019_ref_year = atb_use_hospital_level_agg %>%
  group_by(atb_class) %>%
  nest() %>%
  mutate(
    p_f = map(data, atbpval, level = "national"),
    mul_comp = map(data, atbmulcomp, level = "national")
  ) %>%
  select(-data) %>%
  unnest(c(p_f, mul_comp)) %>%
  ungroup() %>%
  mutate(
    comparison = paste0(group1, "vs", group2),
    conf_int = ifelse( -conf.low > -conf.high,
                       paste0(-round(estimate, 1), " (", -round(conf.high, 1), ", ", -round(conf.low,1), ")"),
                       paste0(-round(estimate, 1), " (", -round(conf.low, 1), ", ", -round(conf.high,1), ")")
    ),
    p.adj.bis = round(p.adj, 3),
    p.adj = ifelse(p.adj < 0.001, "<0.001", round(p.adj, 3)),
    p_f = ifelse(p_f < 0.001, "<0.001", round(p_f, 3))
  ) %>%
  select(atb_class, p_f, comparison, conf_int, p.adj, p.adj.bis) %>%
  pivot_wider(names_from = comparison, values_from = c(conf_int, p.adj, p.adj.bis)) %>%
  arrange(atb_class) %>%
  mutate(atb_class = as.character(atb_class)) %>%
  mutate(DCI = ifelse(atb_class %in% c("Azithromycin", "Imipenem + Meropenem", "Vancomycin"), atb_class, ""),
         .after = 1)

# Hospital consumption table
hospital_tab = national_2019_ref_year %>%
  gt(.) %>%
  cols_hide(
    columns = c(p.adj.bis_2019vs2020, p.adj.bis_2019vs2021, p.adj.bis_2019vs2022)
  ) %>%
  tab_spanner(
    label = "2019 vs 2020",
    columns = c(conf_int_2019vs2020, `p.adj_2019vs2020`)
  ) %>%
  tab_spanner(
    label = "2019 vs 2021",
    columns = c(conf_int_2019vs2021, `p.adj_2019vs2021`)
  ) %>%
  tab_spanner(
    label = "2019 vs 2022",
    columns = c(conf_int_2019vs2022, `p.adj_2019vs2022`)
  ) %>%
  cols_label(
    p_f = "Friedman test p-value",
    conf_int_2019vs2020 = "Estimate (98.3% CI)", 
    `p.adj_2019vs2020` = "p-value",
    conf_int_2019vs2021 = "Estimate (98.3% CI)", 
    `p.adj_2019vs2021` = "p-value",
    conf_int_2019vs2022 = "Estimate (98.3% CI)", 
    `p.adj_2019vs2022` = "p-value",
    atb_class = "",
    DCI = ""
  ) %>%
  tab_style(
    style = list(
      cell_text(weight = "bold")
    ),
    locations = list(
      cells_body(
        rows = p.adj.bis_2019vs2020 <= 0.05,
        columns = conf_int_2019vs2020
      ),
      cells_body(
        rows = p.adj.bis_2019vs2021 <= 0.05,
        columns = conf_int_2019vs2021
      ),
      cells_body(
        rows = p.adj.bis_2019vs2022 <= 0.05,
        columns = conf_int_2019vs2022
      ),
      cells_column_spanners(),
      cells_column_labels()
    )
  ) %>%
  tab_options(table.font.size = 11) %>%
  tab_footnote(
    footnote = "Corresponds to trimethoprim and combinations of sulfanomides",
    locations = cells_body(rows = atb_class == "Trimethoprim", columns = atb_class)
  ) %>%
  tab_footnote(
    footnote = "Adjusted p-value with Bonferroni correction",
    locations = cells_column_labels(columns = c(p.adj_2019vs2020, p.adj_2019vs2021, p.adj_2019vs2022))
  ) %>%
  sub_values(
    columns = atb_class,
    rows = DCI != "",
    replacement = "",
    pattern = "Imipenem \\+ Meropenem|Azithromycin|Vancomycin"
  )
gtsave(hospital_tab, filename = "../Paper/Tables/national_hospital.docx")
gtsave(hospital_tab, filename = "tables/Table2.docx")


# ICU consumption
icu_2019_ref_year = atb_use_hospital_level %>%
  filter(secteur == "ICU", code %in% icu_cohort_final, atb_class %in% atb_class_names2) %>%
  mutate(atb_class = factor(atb_class, atb_class_names2)) %>%
  group_by(atb_class) %>%
  nest() %>%
  mutate(
    p_f = map(data, atbpval, level = "national"),
    mul_comp = map(data, atbmulcomp, level = "national")
  ) %>%
  select(-data) %>%
  unnest(c(p_f, mul_comp)) %>%
  ungroup() %>%
  mutate(
    comparison = paste0(group1, "vs", group2),
    conf_int = ifelse( -conf.low > -conf.high,
                       paste0(-round(estimate, 1), " (", -round(conf.high, 1), ", ", -round(conf.low,1), ")"),
                       paste0(-round(estimate, 1), " (", -round(conf.low, 1), ", ", -round(conf.high,1), ")")
    ),
    p.adj.bis = round(p.adj, 3),
    p.adj = ifelse(p.adj < 0.001, "<0.001", round(p.adj, 3)),
    p_f = ifelse(p_f < 0.001, "<0.001", round(p_f, 3))
  ) %>%
  select(atb_class, p_f, comparison, conf_int, p.adj, p.adj.bis) %>%
  pivot_wider(names_from = comparison, values_from = c(conf_int, p.adj, p.adj.bis)) %>%
  arrange(atb_class) %>%
  mutate(atb_class = as.character(atb_class)) %>%
  mutate(DCI = ifelse(atb_class %in% c("Azithromycin", "Imipenem + Meropenem", "Vancomycin"), atb_class, ""),
         .after = 1)

# ICU
icu_tab = icu_2019_ref_year %>%
  gt(.) %>%
  cols_hide(columns = c(p.adj.bis_2019vs2020, p.adj.bis_2019vs2021, p.adj.bis_2019vs2022)) %>%
  tab_spanner(
    label = "2019 vs 2020",
    columns = c(conf_int_2019vs2020, `p.adj_2019vs2020`)
  ) %>%
  tab_spanner(
    label = "2019 vs 2021",
    columns = c(conf_int_2019vs2021, `p.adj_2019vs2021`)
  ) %>%
  tab_spanner(
    label = "2019 vs 2022",
    columns = c(conf_int_2019vs2022, `p.adj_2019vs2022`)
  ) %>%
  cols_label(
    p_f = "Friedman test p-value",
    conf_int_2019vs2020 = "Estimate (98.3% CI)", 
    `p.adj_2019vs2020` = "p-value",
    conf_int_2019vs2021 = "Estimate (98.3% CI)", 
    `p.adj_2019vs2021` = "p-value",
    conf_int_2019vs2022 = "Estimate (98.3% CI)", 
    `p.adj_2019vs2022` = "p-value",
    atb_class = "",
    DCI = ""
  ) %>%
  tab_style(
    style = list(
      cell_text(weight = "bold")
    ),
    locations = list(
      cells_body(
        rows = p.adj.bis_2019vs2020 <= 0.05,
        columns = conf_int_2019vs2020
      ),
      cells_body(
        rows = p.adj.bis_2019vs2021 <= 0.05,
        columns = conf_int_2019vs2021
      ),
      cells_body(
        rows = p.adj.bis_2019vs2022 <= 0.05,
        columns = conf_int_2019vs2022
      ),
      cells_column_spanners(),
      cells_column_labels()
    )
  ) %>%
  tab_options(table.font.size = 11) %>%
  tab_footnote(
    footnote = "Corresponds to trimethoprim and combinations of sulfanomides",
    locations = cells_body(rows = atb_class == "Trimethoprim", columns = atb_class)
  ) %>%
  tab_footnote(
    footnote = "Adjusted p-value with Bonferroni correction",
    locations = cells_column_labels(columns = c(p.adj_2019vs2020, p.adj_2019vs2021, p.adj_2019vs2022))
  ) %>%
  sub_values(
    columns = atb_class,
    rows = DCI != "",
    replacement = "",
    pattern = "Imipenem \\+ Meropenem|Vancomycin|Azithromycin"
  ) 
gtsave(icu_tab, filename = "../Paper/Tables/national_icu.docx", vwidth = 1500)
gtsave(icu_tab, filename = "tables/Table3.docx", vwidth = 1500)

##################################################
# Plot regional changes - Figure 4
##################################################
# Changes from 2019 for specific antibiotic classes
to_add = atb_use_hospital_level %>%
  group_by(atb_class, code, Date_year, region) %>%
  summarise(molDDD = sum(molDDD), Nbhosp = sum(Nbhosp), .groups = "drop") %>%
  mutate(region = recode(region, !!!dict_regions)) %>%
  group_by(atb_class, region) %>%
  nest() %>%
  mutate(out = map(data, atbunicomp)) %>%
  unnest(cols = out) %>%
  select(-data) %>%
  mutate(p = factor(case_when(
    p <= 0.05 & p> 0.01 ~ "*",
    p <= 0.01 & p> 0.001 ~ "**",
    p <= 0.001 ~ "***",
    .default = NA
  ), c("*","**","***")
  ))

france_region = france %>%
  select(region, geometry) %>%
  mutate(region = recode(region, !!!dict_regions)) %>%
  right_join(., to_add, by = "region") %>%
  filter(!atb_class %in% c("Azithromycin", "Imipenem + Meropenem", "Vancomycin", 
                           "Lipopeptides", "Monobactams", "Polymyxins",
                           "Third generation Cephalosporins")) %>%
  mutate(atb_class = factor(
    ifelse(atb_class == "Trimethoprim", "Trimethoprim and\nsulfanomides", atb_class),
    c("Aminoglycosides", "Carbapenems", "Cephalosporins", "Fosfomycin", "Glycopeptides",
      "Macrolides", "Oxazolidinones", "Penicillins", "Quinolones", "Tetracyclines",
      "Trimethoprim and\nsulfanomides", "Total"))
  )

p1 = ggplot() +
  geom_sf(data = france_region, aes(fill = diff)) +
  stat_sf_coordinates(data = subset(france_region, !is.na(p)), aes(size = p), fill = "white", col = "black", shape=21) +
  facet_wrap(facets = vars(atb_class), ncol = 4) +
  scale_fill_gradient2(
    high = "orange", mid = "ivory", low = "deepskyblue4",
    breaks = seq(-0.5, 0.5, 0.25), labels = scales::percent(seq(-0.50, 0.50, 0.25))
    ) +
  scale_size_manual("Paired Wilcoxon test\np-value", values = c("*" = 1.5, "**" = 2.5, "***" = 4.5)) +#, trans = "log10", range = c(0.5,4)) +
  theme_bw() +
  theme(
    legend.title = element_text(size = 10), 
    legend.key.size = unit(0.5, "cm"),
    axis.text = element_text(size = 8)
    ) +
  labs(x = "", y = "", fill = "Median relative %\nchange 2019-2020")
p1

# Regional changes in azithromycin use at the hospital level
df = atb_use_hospital_level %>%
  filter(atb_class == "Azithromycin") %>%
  group_by(code, Date_year, atb_class, region) %>%
  summarise(molDDD = sum(molDDD), Nbhosp = sum(Nbhosp), .groups = "drop") %>%
  mutate(region = recode(region, !!!dict_regions))

p2 = atbplot(df, "", level = "regional", facet_type = "alphabetical") +
  labs(y = "Annual azithromycin consumption\n(DDD/1,000 bed-days)") +
  theme(plot.title = element_blank(), axis.title.x = element_blank())
# ggsave("../Paper/Figures/Figure4B.png", p2, height = 7, width = 8)
  
# Small map of France with abbreviated names of regions
p3 = france %>%
  mutate(region = recode(region, !!!dict_regions)) %>%
  ggplot(.) +
  geom_sf(fill = "white") +
  geom_sf_text(aes(label = region), col = "black", size = 3) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.text = element_blank(),
        axis.ticks = element_blank(), axis.title = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "cm")) 

# Final figure
figure4 = ggarrange(
  p1,
  ggarrange(p2, p3, ncol = 2, widths = c(1,0.3)),
  nrow = 2, heights = c(1,1), labels = c("A", "B")
  )
figure4
ggsave("../Paper/Figures/Figure4.png", figure4, height = 12, width = 8)
ggsave("plots/Figure4.png", figure4, height = 12, width = 8)

##################################################
# Supplementary plot of absolute change between 
# 2019 and 2020
##################################################
# Plot for legend 
p_legend = ggplot() +
  geom_sf(data = france_region) +
  stat_sf_coordinates(data = subset(france_region, !is.na(p)), aes(size = p), fill = "white", col = "black", shape=21) +
  facet_wrap(facets = vars(atb_class), ncol = 4) +
  scale_size_manual("Paired Wilcoxon test p-value", values = c("*" = 1.5, "**" = 2.5, "***" = 4)) +
  theme_bw() +
  theme(legend.title = element_text(size = 9), legend.direction = "horizontal") +
  labs(x = "", y = "", fill = "")

# Changes from 2019 for antibiotic classes
atb_class_new_names = c("Aminoglycosides", "Carbapenems", "Cephalosporins", "Fosfomycin", 
                        "Glycopeptides", "Lipopeptides", "Macrolides", "Monobactams", 
                        "Oxazolidinones", "Penicillins", "Polymyxins", 
                        "Quinolones", "Tetracyclines", "Trimethoprim and\nsulfanomides", "Total")
france_region = france %>%
  select(region, geometry) %>%
  mutate(region = recode(region, !!!dict_regions)) %>%
  right_join(., to_add, by = "region") %>%
  mutate(atb_class = factor(
    ifelse(atb_class == "Trimethoprim", "Trimethoprim and\nsulfanomides", atb_class),
    atb_class_new_names)
  )

# Get all plots separately
all_plots = vector("list", length(atb_class_new_names))
for (a in seq_along(atb_class_new_names)) {
  all_plots[[a]] = ggplot() +
    geom_sf(data = france_region[france_region$atb_class==atb_class_new_names[a], ], aes(fill = -estimate)) +
    stat_sf_coordinates(data = subset(france_region, !is.na(p) & atb_class == atb_class_new_names[a]), aes(size = p), fill = "white", col = "black", shape=21) +
    scale_fill_gradient2(high = "orange", mid = "ivory", low = "deepskyblue4") +
    scale_size_manual(
      "",
      values = c("*" = 1, "**" = 2, "***" = 3),
      guide = "none"
      ) +
    theme_bw() +
    theme(
      legend.key.size = unit(0.3, "cm"),
      legend.text = element_text(size = 6),
      axis.text = element_text(size = 6),
      plot.title = element_text(hjust = 0.5, size = 10, vjust = 0),
      plot.margin = margin(t=0, r=0, b=0, l=0),
    ) +
    labs(x = "", y = "", fill = "", title = atb_class_new_names[a])
  
  ggsave(paste0("plots/antibiotic_use/", atb_class_names[a], ".png"), all_plots[[a]], width = 3.5, height = 3.5)
}

supp_fig = plot_grid(
  plot_grid(plotlist = all_plots, ncol = 4, nrow = 4, align = "hv"),
  cowplot::get_legend(p_legend),
  nrow = 2,
  rel_heights = c(1,0.025)
  )
supp_fig
ggsave("../Paper/Supplementary/antibiotic_consumption_regions.png", 
       supp_fig, height = 8, width = 10)
ggsave("plots/antibiotic_use/antibiotic_consumption_regions.png", 
       supp_fig, height = 8, width = 10)

##################################################
# Table of multiple comparisons at the regional
# level for all antibiotic classes
##################################################
# # Paired Wilcoxon tests for Total consumption
# total_atb_regional = atb %>%
#   filter(secteur == "Hospital", atb_class != "Anti-MDR GNB") %>%
#   group_by(code, Date_year, region) %>%
#   summarise(consumption = sum(molDDD)/unique(Nbhosp)*1000, .groups = "drop") %>%
#   arrange(region, code, Date_year) %>%
#   group_by(region) %>%
#   wilcox_test(
#     consumption ~ Date_year,
#     comparisons = list(c("2019", "2020"), c("2019", "2021"), c("2019", "2022")),
#     paired = T, 
#     p.adjust.method = "bonferroni", 
#     alternative = "two.sided",
#     conf.level = 1-0.05/3,
#     detailed = T
#   ) %>%
#   mutate(atb_class = "Total")
# 
# # Final table of paired Wilcoxon test results for multiple comparisons 
# atb %>%
#   filter(secteur == "Hospital", atb_class != "Anti-MDR GNB") %>%
#   group_by(code, Date_year, atb_class, region) %>%
#   summarise(consumption = sum(molDDD)/unique(Nbhosp)*1000, .groups = "drop") %>%
#   arrange(atb_class, region, code, Date_year) %>%
#   group_by(atb_class, region) %>%
#   wilcox_test(
#     consumption ~ Date_year, 
#     comparisons = list(c("2019", "2020"), c("2019", "2021"), c("2019", "2022")),
#     paired = T, 
#     p.adjust.method = "bonferroni", 
#     alternative = "two.sided",
#     conf.level = 1-0.05/3,
#     detailed = T
#   ) %>%
#   bind_rows(., anti_mdr_regional, total_atb_regional) %>%
#   mutate(
#     comparison = paste0(group1, "vs", group2),
#     conf_int = ifelse( -conf.low > -conf.high,
#                        paste0("(", -round(conf.high, 2), ", ", -round(conf.low,2), ")"),
#                        paste0("(", -round(conf.low, 2), ", ", -round(conf.high,2), ")")
#     ),
#     estimate = -round(estimate, 2),
#     p.adj = round(p.adj, 4)
#   ) %>%
#   select(atb_class, region, comparison, estimate, conf_int, p.adj) %>%
#   filter(comparison %in% c("2019vs2020", "2019vs2021", "2019vs2022")) %>%
#   pivot_wider(names_from = comparison, values_from = c(estimate, conf_int, p.adj)) %>%
#   dplyr::select(atb_class, region, 
#                 estimate_2019vs2020, conf_int_2019vs2020, p.adj_2019vs2020,
#                 estimate_2019vs2021, conf_int_2019vs2021, p.adj_2019vs2021,
#                 estimate_2019vs2022, conf_int_2019vs2022, p.adj_2019vs2022) %>%
#   write.csv(., "../Paper/Supplementary/antibiotic_consumption_hospitals_regional_heterogeneity.csv", row.names = F)
