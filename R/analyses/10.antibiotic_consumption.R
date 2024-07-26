##################################################
# ANALYSES OF ANTIBIOTIC USE IN FRENCH HOSPITALS
##################################################

rm(list = ls())
library(geofacet)
library(gt)
library(tidyverse)
library(rstatix)
library(DescTools)
library(readxl)
library(factoextra)
library(ggpubr)
library(cowplot)
library(gridExtra)
library(NbClust)
library(cluster)
library(ggsignif)
library(sf)

load("data/metadata_admin_espic.rda")
load("data/cohort_final.rda")

source("R/helper/dictionaries.R")
source("R/helper/helper_functions.R")
source("R/helper/helper_plots.R")

##################################################
# Figure 4
##################################################
# Changes from 2019 for antibiotic classes
to_add = atb %>%
  filter(secteur == "Hospital", atb_class != "Anti-MDR GNB") %>%
  group_by(atb_class, code, Date_year, region, Nbhosp) %>%
  summarise(molDDD = sum(molDDD), .groups = "drop") %>%
  bind_rows(., atb_total %>% filter(secteur == "Hospital") %>% dplyr::select(atb_class, code, Date_year, region, Nbhosp, molDDD)) %>%
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
  filter(!atb_class %in% c("Lipopeptides", "Monobactams", "Fosfomycin", "Polymyxins")) %>%
  mutate(atb_class = factor(
    ifelse(atb_class == "Trimethoprim", "Trimethoprim and\nsulfanomides", atb_class),
    c("Aminoglycosides", "Carbapenems", "Cephalosporins", "Glycopeptides", 
      "Macrolides", "Penicillins", "Oxazolidinones", "Quinolones", "Tetracyclines", 
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
# p1

# Azithromycin
df = atb %>%
  filter(secteur == "Hospital", DCI == "Azithromycine") %>%
  mutate(region = recode(region, !!!dict_regions))
p2 = atbplot(df, "", level = "regional", facet_type = "alphabetical") +
  labs(y = "Annual azithromycin consumption\n(DDD/1,000 bed days)") +
  theme(plot.title = element_blank(), axis.title.x = element_blank())
ggsave("../Paper/Figures/Figure4B.png", p2, height = 7, width = 8)
  
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
  ggarrange(p2, p3, ncol = 2, widths = c(1,0.5)),
  nrow = 2, heights = c(1,0.7), labels = c("A", "B")
  )
# figure4
ggsave("../Paper/Figures/Figure4.png", figure4, height = 11, width = 9)

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
france_region = france %>%
  select(region, geometry) %>%
  mutate(region = recode(region, !!!dict_regions)) %>%
  right_join(., to_add, by = "region") %>%
  mutate(atb_class = factor(
    ifelse(atb_class == "Trimethoprim", "Trimethoprim and\nsulfanomides", atb_class),
    c("Aminoglycosides", "Carbapenems", "Cephalosporins", "Fosfomycin", 
      "Glycopeptides", "Lipopeptides", "Macrolides", "Monobactams", 
      "Penicillins", "Oxazolidinones", "Polymyxins", 
      "Quinolones", "Tetracyclines", "Trimethoprim and\nsulfanomides", "Total"))
  )

# Get all plots separately
atb_classes = unique(france_region$atb_class)
all_plots = vector("list", length(atb_classes))
for (a in seq_along(atb_classes)) {
  all_plots[[a]] = ggplot() +
    geom_sf(data = france_region[france_region$atb_class==atb_classes[a], ], aes(fill = -estimate)) +
    stat_sf_coordinates(data = subset(france_region, !is.na(p) & atb_class == atb_classes[a]), aes(size = p), fill = "white", col = "black", shape=21) +
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
    labs(x = "", y = "", fill = "", title = atb_classes[a])
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

##################################################
# Table 1 of national consumption of antibiotics
# by molecular class 
##################################################
# All classes of antibiotics
atb_all_classes = atb %>%
    filter(secteur %in% c("Hospital", "ICU"), atb_class != "Anti-MDR GNB") %>%
    mutate(DCI = "") %>%
    group_by(code, atb_class, secteur, DCI, Date_year, Nbhosp) %>%
    summarise(molDDD = sum(molDDD), .groups = "drop")

# Total antibiotic consumption
atb_total2 = atb_total %>%
  filter(secteur %in% c("Hospital", "ICU")) %>%
  mutate(DCI = "") %>%
  dplyr::select(-c(Nblits, type, region))
  
# Hospital consumption
national_2019_ref_year = atb %>%
  filter(DCI %in% c("Vancomycine", "Imipenem", "Meropenem", "Azithromycine")) %>%
  mutate(DCI = case_when(DCI %in% c("Imipenem", "Meropenem") ~ "Imipenem + Meropenem", 
                         DCI == "Vancomycine" ~ "Vancomycin",
                         DCI == "Azithromycine" ~ "Azithromycin"
                         )) %>%
  group_by(code, atb_class, secteur, DCI, Date_year, Nbhosp) %>%
  summarise(molDDD = sum(molDDD), .groups = "drop") %>%
  bind_rows(., atb_all_classes, atb_total2) %>%
  group_by(secteur, atb_class, DCI) %>%
  nest() %>%
  mutate(
    p_f = map(data, atbpval, level = "national"),
    mul_comp = map(data, atbmulcomp, level = "national")
  ) %>%
  select(-data) %>%
  unnest(c(p_f, mul_comp)) %>%
  mutate(
    comparison = paste0(group1, "vs", group2),
    conf_int = ifelse( -conf.low > -conf.high,
                       paste0("(", -round(conf.high, 1), ", ", -round(conf.low,1), ")"),
                       paste0("(", -round(conf.low, 1), ", ", -round(conf.high,1), ")")
    ),
    estimate = -round(estimate, 2),
    p.adj = round(p.adj, 3),
    p_f = round(p_f, 3),
    atb_class = factor(atb_class, c("Aminoglycosides", "Carbapenems", "Cephalosporins", 
                                    "Fosfomycin", "Glycopeptides", "Lipopeptides", 
                                    "Macrolides", "Monobactams", "Penicillins", 
                                    "Oxazolidinones", "Polymyxins", "Quinolones", 
                                    "Tetracyclines", "Trimethoprim", "Total"))
  ) %>%
  select(secteur, atb_class, DCI, p_f, comparison, estimate, conf_int, p.adj) %>%
  pivot_wider(names_from = comparison, values_from = c(estimate, conf_int, p.adj)) %>%
  arrange(secteur, atb_class, DCI) %>%
  mutate(atb_class = as.character(atb_class)) 


# Hospital
national_2019_ref_year %>%
  filter(secteur == "Hospital") %>%
  group_by(secteur) %>%
  gt(.) %>%
  tab_spanner(
    label = "2019 vs 2020",
    columns = c(estimate_2019vs2020, conf_int_2019vs2020, `p.adj_2019vs2020`)
  ) %>%
  tab_spanner(
    label = "2019 vs 2021",
    columns = c(estimate_2019vs2021, conf_int_2019vs2021, `p.adj_2019vs2021`)
  ) %>%
  tab_spanner(
    label = "2019 vs 2022",
    columns = c(estimate_2019vs2022, conf_int_2019vs2022, `p.adj_2019vs2022`)
  ) %>%
  cols_label(
    p_f = "Friedman test p-value",
    estimate_2019vs2020 = "Estimate",
    conf_int_2019vs2020 = "98.3% CI", 
    `p.adj_2019vs2020` = "adjusted p",
    estimate_2019vs2021 = "Estimate",
    conf_int_2019vs2021 = "98.3% CI", 
    `p.adj_2019vs2021` = "adjusted p",
    estimate_2019vs2022 = "Estimate",
    conf_int_2019vs2022 = "98.3% CI", 
    `p.adj_2019vs2022` = "adjusted p",
    atb_class = "",
    DCI = ""
  ) %>%
  tab_style(
    style = list(
      cell_text(weight = "bold")
    ),
    locations = list(
      cells_body(
        rows = p.adj_2019vs2020 <= 0.05,
        columns = estimate_2019vs2020
      ),
      cells_body(
        rows = p.adj_2019vs2021 <= 0.05,
        columns = estimate_2019vs2021
      ),
      cells_body(
        rows = p.adj_2019vs2022 <= 0.05,
        columns = estimate_2019vs2022
        ),
      cells_row_groups(), 
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
    pattern = "Carbapenems|Glycopeptides|Macrolides"
  ) %>%
  gtsave(filename = "../Paper/Tables/national_hospital.png", vwidth = 1500)

# Hospital
national_2019_ref_year %>%
  filter(secteur == "ICU") %>%
  group_by(secteur) %>%
  gt(.) %>%
  tab_spanner(
    label = "2019 vs 2020",
    columns = c(estimate_2019vs2020, conf_int_2019vs2020, `p.adj_2019vs2020`)
  ) %>%
  tab_spanner(
    label = "2019 vs 2021",
    columns = c(estimate_2019vs2021, conf_int_2019vs2021, `p.adj_2019vs2021`)
  ) %>%
  tab_spanner(
    label = "2019 vs 2022",
    columns = c(estimate_2019vs2022, conf_int_2019vs2022, `p.adj_2019vs2022`)
  ) %>%
  cols_label(
    p_f = "Friedman test p-value",
    estimate_2019vs2020 = "Estimate",
    conf_int_2019vs2020 = "98.3% CI", 
    `p.adj_2019vs2020` = "adjusted p",
    estimate_2019vs2021 = "Estimate",
    conf_int_2019vs2021 = "98.3% CI", 
    `p.adj_2019vs2021` = "adjusted p",
    estimate_2019vs2022 = "Estimate",
    conf_int_2019vs2022 = "98.3% CI", 
    `p.adj_2019vs2022` = "adjusted p",
    atb_class = "",
    DCI = ""
  ) %>%
  tab_style(
    style = list(
      cell_text(weight = "bold")
    ),
    locations = list(
      cells_body(
        rows = p.adj_2019vs2020 <= 0.05,
        columns = estimate_2019vs2020
      ),
      cells_body(
        rows = p.adj_2019vs2021 <= 0.05,
        columns = estimate_2019vs2021
      ),
      cells_body(
        rows = p.adj_2019vs2022 <= 0.05,
        columns = estimate_2019vs2022
      ),
      cells_row_groups(), 
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
    pattern = "Carbapenems|Glycopeptides|Macrolides"
  ) %>%
  gtsave(filename = "../Paper/Tables/national_icu.png", vwidth = 1500)
  
# # Verify that estimates correspond to the median of the differences
# k = atb %>%
#   filter(
#     secteur %in% c("Total établissement", "Réanimation"), 
#     DCI %in% c("Vancomycine", "Imipenem", "Meropenem", "Azithromycine") 
#   ) %>%
#   mutate(DCI = case_when(DCI %in% c("Imipenem", "Meropenem") ~ "Imipenem + Meropenem", 
#                          DCI == "Vancomycine" ~ "Vancomycin",
#                          DCI == "Azithromycine" ~ "Azythromycin"
#   )) %>%
#   group_by(code, atb_class, secteur, DCI, Date_year, Nbhosp) %>%
#   summarise(molDDD = sum(molDDD), .groups = "drop") %>%
#   bind_rows(., atb_total, atb_all_classes) %>%
#   filter(atb_class == "Aminoglycosides", secteur == "Total établissement", Date_year <= 2020) %>%
#   mutate(consumption = molDDD/Nbhosp*1000) %>%
#   dplyr::select(-Nbhosp, -molDDD) %>%
#   pivot_wider(names_from = Date_year, values_from = consumption)
# 
# wilcox.test(k$`2019`, k$`2020`, paired = T, conf.int = T)
# median(sort(outer(k$`2019`, k$ `2020`, "-")))
# median(k$`2019`-k$`2020`)

##################################################
# Table 2 of multiple comparisons for last line 
# antibiotics
##################################################
atb %>%
  filter(DCI %in% last_line, secteur %in% c("ICU", "Hospital")) %>%
  group_by(code, secteur, Date_year, Nbhosp) %>%
  summarise(molDDD = sum(molDDD), .groups = "drop") %>%
  group_by(secteur) %>%
  nest() %>%
  mutate(
    p_f = map(data, atbpval, level = "national"),
    mul_comp = map(data, atbmulcomp, level = "national")
  ) %>%
  select(-data) %>%
  unnest(c(p_f, mul_comp)) %>%
  mutate(
    comparison = paste0(group1, "vs", group2),
    conf_int = ifelse( -conf.low > -conf.high,
                       paste0("(", -round(conf.high, 1), ", ", -round(conf.low,1), ")"),
                       paste0("(", -round(conf.low, 1), ", ", -round(conf.high,1), ")")
    ),
    estimate = -round(estimate, 2),
    p.adj = round(p.adj, 3),
    p_f = round(p_f, 3)
  ) %>%
  dplyr::select(secteur, p_f, comparison, estimate, conf_int, p.adj) %>%
  pivot_wider(names_from = comparison, values_from = c(estimate, conf_int, p.adj)) %>%
  arrange(secteur) %>%
  ungroup() %>%
  gt(.) %>%
  tab_spanner(
    label = "2019 vs 2020",
    columns = c(estimate_2019vs2020, conf_int_2019vs2020, `p.adj_2019vs2020`)
  ) %>%
  tab_spanner(
    label = "2019 vs 2021",
    columns = c(estimate_2019vs2021, conf_int_2019vs2021, `p.adj_2019vs2021`)
  ) %>%
  tab_spanner(
    label = "2019 vs 2022",
    columns = c(estimate_2019vs2022, conf_int_2019vs2022, `p.adj_2019vs2022`)
  ) %>%
  cols_label(
    p_f = "Friedman test p-value",
    estimate_2019vs2020 = "Estimate",
    conf_int_2019vs2020 = "98.3% CI", 
    `p.adj_2019vs2020` = "adjusted p",
    estimate_2019vs2021 = "Estimate",
    conf_int_2019vs2021 = "98.3% CI", 
    `p.adj_2019vs2021` = "adjusted p",
    estimate_2019vs2022 = "Estimate",
    conf_int_2019vs2022 = "98.3% CI", 
    `p.adj_2019vs2022` = "adjusted p",
    secteur = ""
  ) %>%
  tab_style(
    style = list(
      cell_text(weight = "bold")
    ),
    locations = list(
      cells_body(
        rows = p.adj_2019vs2020 <= 0.05,
        columns = estimate_2019vs2020
      ),
      cells_body(
        rows = p.adj_2019vs2021 <= 0.05,
        columns = estimate_2019vs2021
      ),
      cells_body(
        rows = p.adj_2019vs2022 <= 0.05,
        columns = estimate_2019vs2022
      ),
      cells_column_spanners(),
      cells_column_labels()
    )
  ) %>%
  tab_options(table.font.size = 11) %>%
  gtsave(filename = "../Paper/Tables/national_hospital_lastline.png", vwidth = 1500)

##################################################
# Table of multiple comparisons at the regional
# level for all antibiotic classes, including 
# the functional class of last line antibiotics
##################################################
# Paired Wilcoxon tests for anti-MDR GNB antibiotics
anti_mdr_regional = atb %>%
  filter(secteur == "Hospital", DCI %in% last_line) %>%
  group_by(code, Date_year, region) %>%
  summarise(consumption = sum(molDDD)/unique(Nbhosp)*1000, .groups = "drop") %>%
  arrange(region, code, Date_year) %>%
  group_by(region) %>%
  wilcox_test(
    consumption ~ Date_year, 
    comparisons = list(c("2019", "2020"), c("2019", "2021"), c("2019", "2022")),
    paired = T, 
    p.adjust.method = "bonferroni", 
    alternative = "two.sided",
    conf.level = 1-0.05/3,
    detailed = T
  ) %>%
  mutate(atb_class = "Anti-MDR GNB")

# Paired Wilcoxon tests for Total consumption
total_atb_regional = atb %>%
  filter(secteur == "Hospital", atb_class != "Anti-MDR GNB") %>%
  group_by(code, Date_year, region) %>%
  summarise(consumption = sum(molDDD)/unique(Nbhosp)*1000, .groups = "drop") %>%
  arrange(region, code, Date_year) %>%
  group_by(region) %>%
  wilcox_test(
    consumption ~ Date_year,
    comparisons = list(c("2019", "2020"), c("2019", "2021"), c("2019", "2022")),
    paired = T, 
    p.adjust.method = "bonferroni", 
    alternative = "two.sided",
    conf.level = 1-0.05/3,
    detailed = T
  ) %>%
  mutate(atb_class = "Total")

# Final table of paired Wilcoxon test results for multiple comparisons 
atb %>%
  filter(secteur == "Hospital", atb_class != "Anti-MDR GNB") %>%
  group_by(code, Date_year, atb_class, region) %>%
  summarise(consumption = sum(molDDD)/unique(Nbhosp)*1000, .groups = "drop") %>%
  arrange(atb_class, region, code, Date_year) %>%
  group_by(atb_class, region) %>%
  wilcox_test(
    consumption ~ Date_year, 
    comparisons = list(c("2019", "2020"), c("2019", "2021"), c("2019", "2022")),
    paired = T, 
    p.adjust.method = "bonferroni", 
    alternative = "two.sided",
    conf.level = 1-0.05/3,
    detailed = T
  ) %>%
  bind_rows(., anti_mdr_regional, total_atb_regional) %>%
  mutate(
    comparison = paste0(group1, "vs", group2),
    conf_int = ifelse( -conf.low > -conf.high,
                       paste0("(", -round(conf.high, 2), ", ", -round(conf.low,2), ")"),
                       paste0("(", -round(conf.low, 2), ", ", -round(conf.high,2), ")")
    ),
    estimate = -round(estimate, 2),
    p.adj = round(p.adj, 4)
  ) %>%
  select(atb_class, region, comparison, estimate, conf_int, p.adj) %>%
  filter(comparison %in% c("2019vs2020", "2019vs2021", "2019vs2022")) %>%
  pivot_wider(names_from = comparison, values_from = c(estimate, conf_int, p.adj)) %>%
  dplyr::select(atb_class, region, 
                estimate_2019vs2020, conf_int_2019vs2020, p.adj_2019vs2020,
                estimate_2019vs2021, conf_int_2019vs2021, p.adj_2019vs2021,
                estimate_2019vs2022, conf_int_2019vs2022, p.adj_2019vs2022) %>%
  write.csv(., "../Paper/Supplementary/antibiotic_consumption_hospitals_regional_heterogeneity.csv", row.names = F)

##################################################
# Plot raw data
##################################################
atb_class_names = c("Aminoglycosides", "Carbapenems", "Cephalosporins", "Fosfomycin", 
                    "Glycopeptides", "Lipopeptides", "Macrolides", "Monobactams", 
                    "Oxazolidinones", "Penicillins", "Polymyxins", "Quinolones",
                    "Tetracyclines", "Trimethoprim", "Total")

# Total antibiotics
median_total = atb %>%
  filter(secteur %in% c("Hospital", "ICU"), atb_class != "Anti-MDR GNB") %>%
  group_by(code, Date_year, secteur) %>%
  summarise(consumption = sum(molDDD) / unique(Nbhosp) * 1000, .groups = "drop") %>%
  group_by(Date_year, secteur) %>%
  summarise(med = median(consumption), code = NA, .groups = "drop") %>%
  mutate(atb_class = "Total")

atb_total = atb %>%
  filter(secteur %in% c("Hospital", "ICU"), atb_class != "Anti-MDR GNB") %>%
  group_by(code, Date_year, secteur) %>%
  summarise(consumption = sum(molDDD) / unique(Nbhosp) * 1000, .groups = "drop") %>%
  mutate(atb_class = "Total")

# Hospitals
median_consumption = atb %>%
  filter(secteur == "Hospital", atb_class != "Anti-MDR GNB") %>%
  group_by(code, Date_year, atb_class) %>%
  summarise(consumption = sum(molDDD) / unique(Nbhosp) * 1000, .groups = "drop") %>%
  group_by(Date_year, atb_class) %>%
  summarise(med = median(consumption), code = NA, .groups = "drop") %>%
  bind_rows(., median_total %>% filter(secteur == "Hospital")) %>%
  mutate(atb_class = factor(atb_class, atb_class_names))

p1 = atb %>%
  filter(secteur == "Hospital", atb_class != "Anti-MDR GNB") %>%
  group_by(code, Date_year, atb_class) %>%
  summarise(consumption = sum(molDDD) / unique(Nbhosp) * 1000, .groups = "drop") %>%
  bind_rows(., atb_total %>% filter(secteur == "Hospital")) %>%
  mutate(atb_class = factor(atb_class, atb_class_names)) %>%
  ggplot(., aes(x = Date_year, y = consumption, group = code)) +
  geom_line(col = "grey") +
  geom_point(data = median_consumption, aes(x = Date_year, y = med), col = "red") +
  facet_wrap(facets = vars(atb_class), scales = "free_y", ncol = 5) +
  theme_bw() +
  labs(x = "", y = "Annual antibiotic consumption (in DDD/1,000 bed-days)",
       title = "Hospitals")

# ICUs
median_consumption_icu = atb %>%
  filter(secteur == "ICU", atb_class != "Anti-MDR GNB") %>%
  group_by(code, Date_year, atb_class) %>%
  summarise(consumption = sum(molDDD) / unique(Nbhosp) * 1000, .groups = "drop") %>%
  group_by(Date_year, atb_class) %>%
  summarise(med = median(consumption), code = NA, .groups = "drop") %>%
  bind_rows(., median_total %>% filter(secteur == "ICU")) %>%
  mutate(atb_class = factor(atb_class, atb_class_names))

p2 = atb %>%
  filter(secteur == "ICU", atb_class != "Anti-MDR GNB") %>%
  group_by(code, Date_year, atb_class) %>%
  summarise(consumption = sum(molDDD) / unique(Nbhosp) * 1000, .groups = "drop") %>%
  bind_rows(., atb_total %>% filter(secteur == "ICU")) %>%
  mutate(atb_class = factor(atb_class, atb_class_names)) %>%
  ggplot(., aes(x = Date_year, y = consumption, group = code)) +
  geom_line(col = "grey") +
  geom_point(data = median_consumption_icu, aes(x = Date_year, y = med), col = "red") +
  facet_wrap(facets = vars(atb_class), scales = "free_y", ncol = 5) +
  theme_bw() +
  labs(x = "", y = "Annual antibiotic consumption (in DDD/1,000 bed-days)",
       title = "ICUs")

# Final figure 
p = ggarrange(p1, p2, nrow = 2, labels = c("A", "B"))
ggsave("../Paper/Supplementary/antibiotic_consumption.png", 
       p, height = 10, width = 10)

# All hospitals reported antibiotic consumption by class every year
atb %>%
  filter(secteur == "Hospital", atb_class != "Anti-MDR GNB") %>%
  group_by(code, Date_year, atb_class) %>%
  summarise(consumption = sum(molDDD) / unique(Nbhosp) * 1000, .groups = "drop") %>%
  group_by(code, atb_class) %>%
  mutate(n = n()) %>% 
  filter(n<4) 

# ICUs that did not report at least one antibiotic class at least one year
atb %>%
  filter(secteur == "ICU", atb_class != "Anti-MDR GNB") %>%
  group_by(code, Date_year, atb_class) %>%
  summarise(consumption = sum(molDDD) / unique(Nbhosp) * 1000, .groups = "drop") %>%
  group_by(code, atb_class) %>%
  summarise(n = n(), .groups = "drop") %>% 
  filter(n<4) 

##################################################
# Get DDJ per region and antibiotic class
##################################################
# Compute DDJ for total antibiotics
atb_total = atb %>%
  filter(atb_class != "Anti-MDR GNB") %>%
  group_by(Date_year, region, type, secteur) %>%
  summarise(Total = sum(molDDD, na.rm = T), .groups = "drop")

# Get number of beds
total_beds_days = atb %>%
  dplyr::select(Date_year, region, code, secteur, Nbhosp) %>%
  distinct() %>%
  group_by(Date_year, region) %>%
  summarise(Nbhosp = sum(Nbhosp), .groups = "drop") %>%
  mutate(secteur = "Hospital")
  
secteur_bed_days = atb %>%
  dplyr::select(Date_year, region, code, secteur, Nbhosp) %>%
  distinct() %>%
  group_by(Date_year, region, secteur) %>%
  summarise(Nbhosp = sum(Nbhosp), .groups = "drop") %>%
  bind_rows(., total_beds_days)

# Get consumption of 3GC
atb_3gc = atb %>%
  filter(ATC %in% vec_3gc) %>%
  mutate(atb_class = "Third_generation_Cephalosporins") %>%
  group_by(Date_year, region, atb_class, secteur) %>%
  summarise(ddj = sum(molDDD, na.rm = T), .groups = "drop")

# Compute DDJ
atb = atb %>% 
  mutate(atb_class = ifelse(DCI == "Ertapenem", "Carbapenems_others", atb_class)) %>%
  group_by(Date_year, region, atb_class, secteur) %>%
  summarise(ddj = sum(molDDD, na.rm = T), .groups = "drop") %>%
  bind_rows(., atb_3gc) %>%
  pivot_wider(names_from = atb_class, values_from = ddj) %>%
  rename(Anti_MDR_GNB = `Anti-MDR GNB`) %>%
  left_join(., secteur_bed_days, by = c("Date_year", "region", "secteur"))

# Save antibiotic consumption data
write.table(
  atb, 
  "data-raw/spares/combined/antibiotics_cohortfinal.txt", 
  row.names = F, 
  sep = "\t"
  )
