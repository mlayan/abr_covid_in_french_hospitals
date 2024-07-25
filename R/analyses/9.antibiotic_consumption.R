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
load("data/dict_id_site.rda")
load("data/dict_site_plot.rda")
load("data/dict_antibiotic_class.rda")
load("data/dict_secteur_spares.rda")
load("data/cohort_final.rda")
load("data/my_regional_grid.rda")
load("data/dict_regions.rda")
load("data/france.rda")
load("data/last_line.rda")
load("data/vec_3gc.rda")

source("R/helper_functions.R")
source("R/helper_plots.R")

##################################################
# Load antibiotic consumption data 
##################################################
# Load resistance data
atb = bind_rows( 
  read_excel("data-raw/spares/2019/ATB2019_NAT.xlsx", sheet = "ATB2019") %>% mutate(Date_year = 2019),
  read_excel("data-raw/spares/2020/atb2020_2021.xlsx", sheet = "ATB2020") %>% mutate(Date_year = 2020),
  read_excel("data-raw/spares/2020/atb2020_2021.xlsx", sheet = "ATB2021") %>% mutate(Date_year = 2021),
  read_excel("data-raw/spares/2022/ATB2022.xlsx", sheet = "ATB2022") %>% mutate(Date_year = 2022)
) %>%
  dplyr::select(-c(Type, type, region, Region)) %>%
  filter(
    code %in% cohort_final, 
    !secteur %in% c("Hématologie", "Maladie inf", "Pédiatrie", "Psychiatrie", "SLD", "Total établissement")
  ) %>%
  mutate(
    atb_class = dplyr::recode(ATC, !!!dict_antibiotic_class),
    secteur = recode(secteur, !!!dict_secteur_spares)
    )

##################################################
# Selection of antibiotics
##################################################
# Table of missing reporting years
atb %>%
  dplyr::select(secteur, DCI, ATC, Date_year) %>%
  distinct() %>%
  group_by(secteur, DCI, ATC) %>%
  mutate(n = n()) %>%
  filter(n < 4) %>%
  dplyr::select(-n) %>%
  filter(secteur == "Medicine") %>%
  group_by(DCI, ATC) %>%
  summarise(Missing_years = paste0(c(2019:2022)[!2019:2022 %in% Date_year], collapse = ", "), 
            .groups = "drop") %>%
  write.csv(., "tables/antibiotics_missing_years.csv", row.names = F)

# Verify that DCI whose ATC are not for systemic use have only one ATC code
not_for_systemic_use = atb %>%
  filter(!grepl("^J01", ATC)) %>%
  dplyr::select(ATC, DCI) %>%
  distinct()
not_for_systemic_use

# Antibiotics that are in the 'Other' class
other_class = atb %>%
  filter(grepl("^J01", ATC), atb_class == "Other") %>%
  dplyr::select(ATC, DCI) %>%
  distinct()
other_class

# Antibiotics that are not reported every year of the study period
missing_years = atb %>%
  filter(atb_class != "Other") %>%
  select(secteur, DCI, ATC, Date_year) %>%
  distinct() %>%
  group_by(secteur, DCI, ATC) %>%
  summarise(n = n(), .groups = "drop") %>%
  filter(n < 4, !DCI %in% last_line) %>%
  dplyr::select(-c(n, secteur)) %>%
  distinct() %>%
  arrange(DCI)
missing_years

# Clean antibiotic dataset 
atb = atb %>%
  # Remove ATC classes that are not J01 (antiinfectives for systemic use)
  filter(grepl("^J01", ATC)) %>%
  # Remove antibiotics of the "Other category"
  filter(atb_class != "Other") %>%
  # remove antibiotics not reported every year
  anti_join(., missing_years, by = c("DCI", "ATC"))

##################################################
# Final clean database
##################################################
# Add consumption at the hospital level
hospital_atb = atb %>%
  group_by(code, Date_year, ATC, DCI, atb_class, type, region) %>%
  summarise(
    molDDD = sum(molDDD),
    Nbhosp = sum(Nbhosp),
    Nblits = sum(Nblits),
    .groups = "drop"
  ) %>%
  mutate(secteur = "Hospital")

atb = bind_rows(atb, hospital_atb)

# Add total consumption of antibiotics
atb_total = atb %>%
  filter(atb_class != "Anti-MDR GNB") %>%
  group_by(code, secteur, Nbhosp, Nblits, Date_year, type, region) %>%
  summarise(molDDD = sum(molDDD), .groups = "drop") %>%
  mutate(atb_class = "Total")

# Verify that there is only one Nbhosp and one Nblits per code x secteur x year
atb_total %>% group_by(code, secteur, Date_year) %>% mutate(n = n()) %>% filter(n > 1)

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
# AWARE WHO classification
##################################################
# # Load aware classification
# aware_classification = read_excel("data-raw/who/WHO-MHP-HPS-EML-2023.04-eng.xlsx", 
#                                   sheet = "AWaRe classification 2023", range = "A4:E261") %>%
#   rename(ATC = `ATC code`, DCI = Antibiotic) %>%
#   mutate(DCI_mod = case_when(
#     grepl("_IV$", DCI) ~ gsub("_IV", " I", DCI),
#     grepl("_oral$", DCI) ~ gsub("_oral", " O", DCI),
#     DCI == "Fusidic-acid" ~ "Acide fusidique",
#     .default = DCI
#   )) %>%
#   dplyr::select(ATC, Category, DCI_mod)
# 
# # 
# atb %>% 
#   mutate(DCI_mod = ifelse(DCI %in% c("Colistin I", "Fosfomycine I", "Fosfomycine O"), 
#                           DCI, 
#                           gsub(" I^| O^", "", DCI))) %>%
#   left_join(., aware_classification, by = c("ATC", "DCI_mod")) %>%
#   filter(is.na(Category)) %>%
#   dplyr::select(DCI_mod, ATC, Category) %>%
#   distinct()

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

##################################################
# Check antibiotic consumption at the national
# level by hospital type in a gt table
# ##################################################
# # Total Atb
# atb_total = atb %>%
#   filter(secteur=="Total établissement") %>%
#   group_by(code, type, Date_year, Nbhosp) %>%
#   summarise(molDDD = sum(molDDD), .groups = "drop") %>%
#   mutate(atb_class = "Total")
# 
# # Atb by class
# atb_all_classes = atb %>%
#   filter(secteur=="Total établissement") %>%
#   group_by(code, type, atb_class, Date_year, Nbhosp) %>%
#   summarise(molDDD = sum(molDDD), .groups = "drop")
# 
# # Atb classes
# all_plots = atb %>%
#   filter(
#     secteur == "Total établissement", 
#     DCI %in% c("Vancomycine", "Imipenem", "Meropenem")
#   ) %>%
#   mutate(atb_class = ifelse(DCI %in% c("Imipenem", "Meropenem"), "Imipenem+Meropenem","Vancomycine")) %>%
#   group_by(code, type, atb_class, Date_year, Nbhosp) %>%
#   summarise(molDDD = sum(molDDD), .groups = "drop") %>%
#   bind_rows(., atb_total, atb_all_classes) %>%
#   group_by(atb_class) %>%
#   nest() %>%
#   mutate(plots = map2(data, atb_class, atbplot, level = "type"))
# all_plots$plots
#
##################################################
# Check macrolide consumption
##################################################
# # Regions with macrolide consumption that changes 
# # between 2019 and 2021
# all_plots %>%
#   unnest(p_val) %>%
#   filter(p_friedman <= 0.05, atb_class == "Macrolides")
# 
# for (a in c("Azithromycine", "Clarithromycine", "Clindamycine", "Erythromycine", 
#             "Josamycine", "Pristinamycine", "Roxithromycine", "Spiramycine")) {
#   # Midécamycine: no consumption in the hospital cohort
#   # Télithromycine: uniquement 6 établissements qui en ont consommé
#   # Lincomycine: reported only in 2019 and 2020
#   
#   # Regions with azythromycin consumption that changes 
#   # between 2019 and 2021
#   azy_etab = atb %>% 
#     filter(grepl(a, DCI), secteur == "Total établissement", atb_class == "Macrolides") %>%
#     group_by(code, secteur, region, Nbhosp, Nblits, Date_year, type, atb_class) %>%
#     summarise(molDDD = sum(molDDD), .groups = "drop")
#   atbplot(azy_etab, a)  
#   
#   # Multiple comparison using Wilcoxon paired test
#   azy_etab %>%
#     left_join(., atbpval(azy_etab), by = "region") %>%
#     mutate(consumption = molDDD / Nbhosp * 1000) %>%
#     arrange(code, Date_year) %>%
#     group_by(region, p_friedman) %>%
#     wilcox_test(
#       consumption ~ Date_year, 
#       paired = T, 
#       p.adjust.method = "bonferroni", 
#       alternative = "two.sided",
#       detailed = T
#     ) %>%
#     mutate(
#       comparison = paste0(group1, "vs", group2),
#       conf_int = ifelse( -conf.low > -conf.high,
#                          paste0("(", -round(conf.high, 2), ", ", -round(conf.low,2), ")"),
#                          paste0("(", -round(conf.low, 2), ", ", -round(conf.high,2), ")")
#       ),
#       estimate = -round(estimate, 2),
#       p.adj = round(p.adj, 4)
#     ) %>%
#     select(region, p_friedman, comparison, estimate, conf_int, p.adj) %>%
#     pivot_wider(names_from = comparison, values_from = c(estimate, conf_int, p.adj)) %>%
#     gt() %>%
#     tab_spanner(
#       label = "2019 vs 2020",
#       columns = c(estimate_2019vs2020, conf_int_2019vs2020, `p.adj_2019vs2020`)
#     ) %>%
#     tab_spanner(
#       label = "2019 vs 2021",
#       columns = c(estimate_2019vs2021, conf_int_2019vs2021, `p.adj_2019vs2021`)
#     ) %>%
#     tab_spanner(
#       label = "2019 vs 2022",
#       columns = c(estimate_2019vs2022, conf_int_2019vs2022, `p.adj_2019vs2022`)
#     ) %>%
#     tab_spanner(
#       label = "2020 vs 2021",
#       columns = c(estimate_2020vs2021, conf_int_2020vs2021, `p.adj_2020vs2021`)
#     ) %>%
#     tab_spanner(
#       label = "2020 vs 2022",
#       columns = c(estimate_2020vs2022, conf_int_2020vs2022, `p.adj_2020vs2022`)
#     ) %>%
#     tab_spanner(
#       label = "2021 vs 2022",
#       columns = c(estimate_2021vs2022, conf_int_2021vs2022, `p.adj_2021vs2022`)
#     ) %>%
#     cols_label(
#       estimate_2019vs2020 = "Estimate",
#       conf_int_2019vs2020 = "95% CI", 
#       `p.adj_2019vs2020` = "adjusted p",
#       estimate_2019vs2021 = "Estimate",
#       conf_int_2019vs2021 = "95% CI", 
#       `p.adj_2019vs2021` = "adjusted p",
#       estimate_2019vs2022 = "Estimate",
#       conf_int_2019vs2022 = "95% CI", 
#       `p.adj_2019vs2022` = "adjusted p",
#       estimate_2020vs2021 = "Estimate",
#       conf_int_2020vs2021 = "95% CI", 
#       `p.adj_2020vs2021` = "adjusted p",
#       estimate_2020vs2022 = "Estimate",
#       conf_int_2020vs2022 = "95% CI", 
#       `p.adj_2020vs2022` = "adjusted p",
#       estimate_2021vs2022 = "Estimate",
#       conf_int_2021vs2022 = "95% CI", 
#       `p.adj_2021vs2022` = "adjusted p",
#       region = ""
#     ) %>%
#     tab_style(
#       style = list(
#         cell_text(weight = "bold")
#       ),
#       locations = cells_body(
#         rows = p_friedman <= 0.05
#       )
#     ) %>%
#     tab_options(table.font.size = 11) %>%
#     cols_hide(columns = p_friedman) %>%
#     gtsave(filename = paste0("tables/", tolower(a), "_multi_comparisons.png"))
#   
#   # Regions with macrolide antibiotic consumption that changes 
#   # between 2019 and 2021 in ICUs
#   azy_etab = atb %>% 
#     filter(grepl(a, DCI), secteur == "Réanimation", atb_class == "Macrolides") %>%
#     group_by(code, secteur, region, Nbhosp, Nblits, Date_year, type, atb_class) %>%
#     summarise(molDDD = sum(molDDD), .groups = "drop") %>%
#     group_by(code) %>%
#     mutate(n = n()) %>%
#     ungroup() %>%
#     group_by(region) %>%
#     mutate(n_r = n()) %>%
#     ungroup() %>%
#     filter(n == 4, n_r > 4) 
#   atbplot(azy_etab, paste0(a, " in ICUs"))
#   
#   # National analysis for azithromycin consumption in ICUs
#   azy_etab_icus = atb %>% 
#     filter(grepl(a, DCI), secteur == "Réanimation", atb_class == "Macrolides") %>%
#     group_by(code, secteur, region, Nbhosp, Nblits, Date_year, type, atb_class) %>%
#     summarise(molDDD = sum(molDDD), .groups = "drop") %>%
#     group_by(code) %>%
#     mutate(n = n()) %>%
#     ungroup() %>%
#     filter(n == 4) %>%
#     mutate(consumption = molDDD/Nbhosp*1000) %>%
#     arrange(code, Date_year)
#   m = max(azy_etab_icus$consumption)
#   p_friedman = friedman.test(consumption ~ Date_year | code, data = azy_etab_icus)
#   
#   azy_etab_icus_mean = azy_etab_icus %>%
#     group_by(Date_year) %>%
#     summarise(consumption = mean(consumption), .groups = "drop")
#   
#   azy_etab_icus_comp = azy_etab_icus %>%
#     arrange(code, Date_year) %>%
#     wilcox_test(
#       consumption ~ Date_year, 
#       paired = T, 
#       p.adjust.method = "bonferroni", 
#       alternative = "two.sided",
#       detailed = T
#     ) %>%
#     mutate(group1 = factor(group1), group2 = factor(group2)) %>%
#     filter(p.adj <= 0.05, group1 == "2019") %>%
#     mutate(y = case_when(group2 == "2020" ~ m/2, group2 == "2021" ~ m*2/3, group2 == "2022" ~ m*5/6))
#   
#   ggplot() +
#     geom_line(data = azy_etab_icus, aes(x = factor(Date_year), y = consumption, group = code), alpha = 0.1) +
#     geom_point(data = azy_etab_icus_mean, aes(x = factor(Date_year), y = consumption), col = "red") +
#     geom_signif(data = azy_etab_icus_comp, aes(xmin = group1, xmax = group2, annotations = p.adj, y_position = y), 
#                 manual = T) +
#     theme_bw() +
#     labs(x = "", y = "Annual antibiotic consumption\n(DDD/1,000 hospitalisation days)", 
#          title = paste0(a, " in ICUs"))
#   ggsave(paste0("plots/cons_atb/", a, "inicus_national.png"), height = 4, width = 5) 
# }
##################################################
# Investigate colinearity using clustering method
# Euclidean distance k-means
##################################################
# Prepare data for clustering---------------------
# Annual consumption per region in wide format
# regional_atb = atb %>%
#   filter(secteur == "Total établissement") %>%
#   group_by(code, Date_year, atb_class, region) %>%
#   summarise(
#     consumption = sum(molDDD),
#     Nbhosp = unique(Nbhosp),
#     .groups = "drop"
#   ) %>%
#   group_by(Date_year, region, atb_class) %>%
#   summarise(consumption = sum(consumption)/sum(Nbhosp)*1000, .groups = "drop") %>%
#   pivot_wider(names_from = Date_year, values_from = consumption, names_prefix = "year_")
# 
# # Verify that there is no missing value
# regional_atb %>%
#   filter(if_any(matches("year_"), is.na))
#
# # K-means clustering on unnormalized data-------------------------
# # Plot of all data - regional consumption by antibiotic class
# p1 = regional_atb %>%
#   pivot_longer(matches("year_"), names_to = "year", values_to = "consumption") %>%
#   mutate(year = gsub("year_", "", year)) %>%
#   ggplot(., aes(x = year, y = consumption, group = interaction(region, atb_class),
#                 col = region, linetype = atb_class)) +
#   geom_line() +
#   theme_bw() +
#   labs(x = "", y = "Antibiotic consumption in DDD/1,000 bed-days", 
#        col = "Region",
#        linetype = "Antibiotic class", title = "All data")
# 
# # Define number of clusters
# wss = map_dbl(1:10, ~{kmeans(dplyr::select(regional_atb, -c(region, atb_class)), ., nstart=50,iter.max = 15 )$tot.withinss})
# n_clust = 1:10
# p2 = as.data.frame(cbind("n_clust" = n_clust, "wss" = wss)) %>%
#   ggplot(., aes(y = wss, x = n_clust), colour = "#82518c") +
#   geom_line() +
#   theme_bw() +
#   scale_x_continuous(labels = scales::number_format(accuracy = 1)) +
#   labs(x = "Number of clusters", y = "Within cluster variation",
#        title = "Define number of clusters")
# nc = NbClust(data = dplyr::select(regional_atb, -c(region, atb_class)), 
#              distance = "euclidean", min.nc = 2, max.nc = 10, method = "kmeans")
# 
# # Add cluster assignation  
# clusters = kmeans(dplyr::select(regional_atb, -c(region, atb_class)), centers = 3)
# regional_atb_long <- regional_atb %>% 
#   mutate(cluster = clusters$cluster) %>%
#   pivot_longer(matches("year_"), names_to = "year", values_to = "consumption") %>%
#   mutate(year = gsub("year_", "", year))
# 
# centers_long = data.frame(clusters$centers) %>%
#   rownames_to_column(var = "cluster") %>%
#   pivot_longer(cols = -cluster, names_to = "year", values_to = "consumption") %>%  
#   mutate(year = gsub("year_", "", year))
# 
# p3 = ggplot() +
#   geom_line(data = regional_atb_long, 
#             aes(y = consumption, x = year, group = interaction(region, atb_class),
#                 col = region, linetype = atb_class)) +
#   facet_wrap(facets = vars(cluster), nrow = 1) + 
#   geom_line(data = centers_long, aes(y = consumption, x = year, group = cluster), 
#             col = "black", linewidth = 2) +
#   theme_bw() +
#   labs(x = "", y = "Antibiotic consumption in DDD/1,000 bed-days", 
#        col = "Region",
#        linetype = "Antibiotic class", title = "K-means clustering")
# 
# # Plot all data at the same time
# p = grid.arrange(
#   arrangeGrob(
#     arrangeGrob(p1+theme(legend.position = "hidden"),p2,NULL,ncol = 3),
#     p3+theme(legend.position="hidden"),nrow = 2),
#   get_legend(p3), ncol = 2, widths = c(4,1)
# )
# # ggsave("plots/cons_atb/clustering_unnormalized.png", p, 
# #        height = 8, width = 10)
# 
# K-medoids clustering on normalized data aggregated at the regional level-------------------------
# # Normalized data
# regional_atb_norm = regional_atb %>%
#   pivot_longer(matches("year_"), names_to = "year", values_to = "consumption") %>%
#   group_by(region, atb_class) %>%
#   mutate(consumption = (consumption - mean(consumption)) / sd(consumption)) %>%
#   ungroup() %>%
#   pivot_wider(names_from = year, values_from = consumption)
# 
# # Plot of all data - regional consumption by antibiotic class
# p1 = regional_atb_norm %>%
#   pivot_longer(matches("year_"), names_to = "year", values_to = "consumption") %>%
#   mutate(year= gsub("year_", "", year)) %>%
#   ggplot(., aes(x = year, y = consumption, group = interaction(region, atb_class),
#                 col = region, linetype = atb_class)) +
#   geom_line() +
#   theme_bw() +
#   labs(x = "", y = "Antibiotic consumption\n(DDD/1,000 bed-days)",
#        col = "Region",
#        linetype = "Antibiotic class", title = "All data")
# 
# # Define number of clusters
# df = select(regional_atb_norm, -c(region, atb_class))
# summary(df)
# p21 = fviz_nbclust(df, pam, method = "wss", medoids = "random",
#                    nstart = 100, metric = "euclidean", k.max = 13) +
#   labs(title = "")
# p22 = fviz_nbclust(df, pam, method = "silhouette", medoids = "random",
#                    nstart = 100, metric = "euclidean", k.max = 13) +
#   labs(title = "")
# p23 = fviz_gap_stat(cluster::clusGap(df, FUN = pam, nstart = 100, K.max = 13,
#                                      B = 100)) +
#   labs(title = "")
# ggarrange(p21, p22, p23, ncol = 3)
# 
# # Add cluster assignation
# clusters = pam(df, 5, medoids = "random", nstart = 100)
# regional_atb_norm_long = regional_atb_norm %>%
#   mutate(cluster = clusters$clustering) %>%
#   pivot_longer(matches("year_"), names_to = "year", values_to = "consumption") %>%
#   mutate(year = gsub("year_", "", year))
# 
# centers_norm_long = data.frame(clusters$medoids) %>%
#   rownames_to_column(var = "cluster") %>%
#   pivot_longer(cols = -cluster, names_to = "year", values_to = "consumption") %>%
#   mutate(year = gsub("year_", "", year))
# 
# p3 = ggplot() +
#   geom_line(data = regional_atb_norm_long,
#             aes(y = consumption, x = year,
#                 group = interaction(region, atb_class),
#                 col = region, linetype = atb_class)) +
#   facet_wrap(facets = vars(cluster), nrow = 1) +
#   geom_line(data = centers_norm_long, aes(y = consumption, x = year,
#                                           group = cluster),
#             col = "black", linewidth = 2) +
#   theme_bw() +
#   labs(x = "", y = "Normalized antibiotic consumption",
#        col = "Region",
#        linetype = "Antibiotic class")
# p3
# 
# p4 = regional_atb %>%
#   pivot_longer(matches("year_"), names_to = "year", values_to = "consumption") %>%
#   mutate(year = gsub("year_", "", year)) %>%
#   left_join(., regional_atb_norm_long %>% select(-consumption), by = c("region", "atb_class", "year")) %>%
#   ggplot(.,aes(y = consumption, x = year, group = interaction(region, atb_class),
#                col = region, linetype = atb_class)) +
#   geom_line() +
#   facet_wrap(facets = vars(cluster), nrow = 1) +
#   theme_bw() +
#   labs(x = "", y = "Antibiotic consumption\n(DDD/1,000 bed-days)",
#        col = "Region",
#        linetype = "Antibiotic class")
# p4
# 
# # Cluster composition
# p5 = regional_atb_norm %>%
#   mutate(cluster = clusters$clustering) %>%
#   group_by(cluster, region) %>%
#   summarise(n = n(), .groups = "drop") %>%
#   group_by(cluster) %>%
#   mutate(p = n/sum(n)) %>%
#   ungroup() %>%
#   ggplot(., aes(x = cluster, y = p, fill = region, label = n)) +
#   geom_bar(stat = "identity") +
#   geom_text(position = position_stack(vjust = 0.5)) +
#   theme_bw() +
#   labs(x = "Clusters", y = "Repartition", fill = "Region")
# 
# p6 = regional_atb_norm %>%
#   mutate(cluster = clusters$clustering) %>%
#   group_by(cluster, atb_class) %>%
#   summarise(n = n(), .groups = "drop") %>%
#   group_by(cluster) %>%
#   mutate(p = n/sum(n)) %>%
#   ungroup() %>%
#   ggplot(., aes(x = cluster, y = p, fill = atb_class, label = n)) +
#   geom_bar(stat = "identity") +
#   geom_text(position = position_stack(vjust = 0.5)) +
#   theme_bw() +
#   labs(x = "Clusters", y = "Repartition", fill = "Antibiotic class")
# 
# # Plot all data at the same time
# p_clusters = grid.arrange(
#   arrangeGrob(
#     arrangeGrob(
#       arrangeGrob(p1+theme(legend.position = "hidden"), NULL, NULL, ncol = 3),
#       arrangeGrob(p21, p22, p23, ncol = 3),
#       nrow = 2
#     ),
#     p3+theme(legend.position="hidden"),
#     p4+theme(legend.position="hidden"),
#     nrow = 3, heights = c(1.5,1, 1)
#   ),
#   get_legend(p3), ncol = 2, widths = c(3,1)
# )
# ggsave("plots/cons_atb/medoids_regional.png", p_clusters,
#        height = 10, width = 10)
# 
# p_comp = ggarrange(p5, p6, nrow = 2)
# ggsave("plots/cons_atb/medoids_regional_composition.png", p_comp,
#        height = 8, width = 8)
# 
# # Hierarchical clustering
# df = regional_atb_norm %>%
#   mutate(new_names = paste0(region, "_", atb_class)) %>%
#   select(-c(region, atb_class)) %>%
#   column_to_rownames(var = "new_names")
# regional_dist = dist(df, "euclidean")
# clust_fit = hclust(regional_dist, method = "ward.D")
# 
# png("plots/cons_atb/hierarchical_regional.png", height = 20, width = 20,
#     units = "cm", res = 360)
# par(mfrow = c(1,1))
# plot(clust_fit, xlab = "", main = "Cluster Dendogram - Ward's method", cex = 0.6)
# rect.hclust(clust_fit, k=4)
# dev.off()
# 
# p1 = regional_atb_norm %>%
#   mutate(cluster = cutree(clust_fit, k= 4)) %>%
#   pivot_longer(matches("year_"), values_to = "consumption", names_to = "year") %>%
#   mutate(year = gsub("year_", "", year)) %>%
#   ggplot(., aes(x = year, y = consumption, linetype = atb_class, col = region, group = interaction(region, atb_class))) +
#   geom_line() +
#   facet_grid(cols = vars(cluster)) +
#   theme_bw() +
#   labs(x = "", y = "Normalized annual antibiotic consumption\n(DDD per 1,000 bed-days)",
#        linetype = "Antibiotic class", col = "Geographical region")
# 
# p2 = regional_atb %>%
#   mutate(cluster = cutree(clust_fit, k= 4)) %>%
#   pivot_longer(matches("year_"), values_to = "consumption", names_to = "year") %>%
#   mutate(year = gsub("year_", "", year)) %>%
#   ggplot(., aes(x = year, y = consumption, linetype = atb_class, col = region, group = interaction(region, atb_class))) +
#   geom_line() +
#   facet_grid(cols = vars(cluster)) +
#   theme_bw() +
#   labs(x = "", y = "Annual antibiotic consumption\n(DDD per 1,000 bed-days)",
#        linetype = "Antibiotic class", col = "Geographical region")
# 
# p3 = regional_atb_norm %>%
#   mutate(cluster = cutree(clust_fit, k= 4)) %>%
#   group_by(cluster, atb_class) %>%
#   summarise(n = n(), .groups = "drop") %>%
#   group_by(cluster) %>%
#   mutate(p = n/sum(n)) %>%
#   ungroup() %>%
#   ggplot(., aes(x = cluster, y = p, fill = atb_class, label = n)) +
#   geom_bar(stat = "identity") +
#   geom_text(position = position_stack(vjust = 0.5)) +
#   theme_bw() +
#   labs(x = "Clusters", y = "Repartition", fill = "Antibiotic class")
# 
# 
# p4 = regional_atb_norm %>%
#   mutate(cluster = cutree(clust_fit, k= 4)) %>%
#   group_by(cluster, region) %>%
#   summarise(n = n(), .groups = "drop") %>%
#   group_by(cluster) %>%
#   mutate(p = n/sum(n)) %>%
#   ungroup() %>%
#   ggplot(., aes(x = cluster, y = p, fill = region, label = n)) +
#   geom_bar(stat = "identity") +
#   geom_text(position = position_stack(vjust = 0.5)) +
#   theme_bw() +
#   labs(x = "Clusters", y = "Repartition", fill = "Region")
# 
# pdata = ggarrange(p1,p2, nrow = 2, common.legend = T, legend = "right")
# ggsave("plots/cons_atb/hierarchical_clusters_regional.png", pdata,
#        height = 7, width = 10)
# pcomp = ggarrange(p3,p4, nrow = 2)
# ggsave("plots/cons_atb/hierarchical_comp_regional.png", pcomp,
#        height = 7, width = 10)

# K-medoids clustering on normalized data aggregated at the type level----------
# # Normalized data
# type_atb_norm = atb %>%
#   filter(secteur == "Total établissement") %>%
#   group_by(code, Date_year, atb_class, type) %>%
#   summarise(
#     consumption = sum(molDDD),
#     Nbhosp = unique(Nbhosp),
#     .groups = "drop"
#   ) %>%
#   group_by(Date_year, type, atb_class) %>%
#   summarise(consumption = sum(consumption)/sum(Nbhosp)*1000, .groups = "drop") %>%
#   pivot_wider(names_from = Date_year, values_from = consumption, names_prefix = "year_") %>%
#   pivot_longer(matches("year_"), names_to = "year", values_to = "consumption") %>%
#   group_by(type, atb_class) %>%
#   mutate(consumption = (consumption - mean(consumption)) / sd(consumption)) %>%
#   ungroup() %>%
#   pivot_wider(names_from = year, values_from = consumption)
# 
# # Plot of all data - regional consumption by antibiotic class
# p1 = type_atb_norm %>%
#   pivot_longer(matches("year_"), names_to = "year", values_to = "consumption") %>%
#   mutate(year= gsub("year_", "", year)) %>%
#   ggplot(., aes(x = year, y = consumption, group = interaction(type, atb_class),
#                 col = type, linetype = atb_class)) +
#   geom_line() +
#   theme_bw() +
#   labs(x = "", y = "Antibiotic consumption\n(DDD/1,000 bed-days)",
#        col = "Region",
#        linetype = "Antibiotic class", title = "All data")
# 
# # Define number of clusters
# df = select(type_atb_norm, -c(type, atb_class))
# summary(df)
# p21 = fviz_nbclust(df, pam, method = "wss", medoids = "random",
#                    nstart = 100, metric = "euclidean", k.max = 13) +
#   labs(title = "")
# p22 = fviz_nbclust(df, pam, method = "silhouette", medoids = "random",
#                    nstart = 100, metric = "euclidean", k.max = 13) +
#   labs(title = "")
# p23 = fviz_gap_stat(cluster::clusGap(df, FUN = pam, nstart = 100, K.max = 13,
#                                      B = 100)) +
#   labs(title = "")
# ggarrange(p21, p22, p23, ncol = 3)
# 
# # Add cluster assignation
# clusters = pam(df, 4, medoids = "random", nstart = 100)
# type_atb_norm_long = type_atb_norm %>%
#   mutate(cluster = clusters$clustering) %>%
#   pivot_longer(matches("year_"), names_to = "year", values_to = "consumption") %>%
#   mutate(year = gsub("year_", "", year))
# 
# centers_norm_long = data.frame(clusters$medoids) %>%
#   rownames_to_column(var = "cluster") %>%
#   pivot_longer(cols = -cluster, names_to = "year", values_to = "consumption") %>%
#   mutate(year = gsub("year_", "", year))
# 
# p3 = ggplot() +
#   geom_line(data = type_atb_norm_long,
#             aes(y = consumption, x = year,
#                 group = interaction(type, atb_class),
#                 col = type, linetype = atb_class)) +
#   facet_wrap(facets = vars(cluster), nrow = 1) +
#   geom_line(data = centers_norm_long, aes(y = consumption, x = year,
#                                           group = cluster),
#             col = "black", linewidth = 2) +
#   theme_bw() +
#   labs(x = "", y = "Antibiotic consumption\n(DDD/1,000 bed-days)",
#        col = "Region",
#        linetype = "Antibiotic class")
# p3
# 
# # Cluster composition
# p4 = type_atb_norm %>%
#   mutate(cluster = clusters$clustering) %>%
#   group_by(cluster, type) %>%
#   summarise(n = n(), .groups = "drop") %>%
#   group_by(cluster) %>%
#   mutate(p = n/sum(n)) %>%
#   ungroup() %>%
#   ggplot(., aes(x = cluster, y = p, fill = type, label = n)) +
#   geom_bar(stat = "identity") +
#   geom_text(position = position_stack(vjust = 0.5)) +
#   theme_bw() +
#   labs(x = "Clusters", y = "Repartition", fill = "Region")
# 
# p5 = type_atb_norm %>%
#   mutate(cluster = clusters$clustering) %>%
#   group_by(cluster, atb_class) %>%
#   summarise(n = n(), .groups = "drop") %>%
#   group_by(cluster) %>%
#   mutate(p = n/sum(n)) %>%
#   ungroup() %>%
#   ggplot(., aes(x = cluster, y = p, fill = atb_class, label = n)) +
#   geom_bar(stat = "identity") +
#   geom_text(position = position_stack(vjust = 0.5)) +
#   theme_bw() +
#   labs(x = "Clusters", y = "Repartition", fill = "Antibiotic class")
# 
# # Plot all data at the same time
# p_clusters = grid.arrange(
#   arrangeGrob(
#     arrangeGrob(
#       arrangeGrob(p1+theme(legend.position = "hidden"), NULL, NULL, ncol = 3),
#       arrangeGrob(p21, p22, p23, ncol = 3),
#       nrow = 2
#     ),
#     p3+theme(legend.position="hidden"), nrow = 2, heights = c(1.5,1)
#   ),
#   get_legend(p3), ncol = 2, widths = c(3,1)
# )
# ggsave("plots/cons_atb/medoids_type.png", p_clusters,
#        height = 10, width = 10)
# 
# p_comp = ggarrange(p4, p5, align = "v", nrow = 2)
# ggsave("plots/cons_atb/medoids_type_composition.png", p_comp,
#        height = 8, width = 8)
# 
# # Hierarchical clustering
# df = type_atb_norm %>%
#   mutate(new_names = paste0(type, "_", atb_class)) %>%
#   select(-c(type, atb_class)) %>%
#   column_to_rownames(var = "new_names")
# regional_dist = dist(df, "euclidean")
# clust_fit = hclust(regional_dist, method = "ward.D")
# 
# png("plots/cons_atb/hierarchical_type.png", height = 20, width = 20,
#     units = "cm", res = 360)
# par(mfrow = c(1,1))
# plot(clust_fit, xlab = "", main = "Cluster Dendogram - Ward's method", cex = 0.6)
# rect.hclust(clust_fit, k=4)
# dev.off()
# 
# p1 = type_atb_norm %>%
#   mutate(cluster = cutree(clust_fit, k= 4)) %>%
#   pivot_longer(matches("year_"), values_to = "consumption", names_to = "year") %>%
#   mutate(year = gsub("year_", "", year)) %>%
#   ggplot(., aes(x = year, y = consumption, linetype = atb_class, col = type,
#                 group = interaction(type, atb_class))) +
#   geom_line() +
#   facet_grid(cols = vars(cluster)) +
#   theme_bw() +
#   labs(x = "", y = "Normalized annual antibiotic consumption\n(DDD per 1,000 bed-days)",
#        linetype = "Antibiotic class", col = "Hospital type")
# 
# p2 = type_atb_norm %>%
#   mutate(cluster = cutree(clust_fit, k= 4)) %>%
#   pivot_longer(matches("year_"), values_to = "consumption", names_to = "year") %>%
#   mutate(year = gsub("year_", "", year)) %>%
#   ggplot(., aes(x = year, y = consumption, linetype = atb_class, col = type,
#                 group = interaction(type, atb_class))) +
#   geom_line() +
#   facet_grid(cols = vars(cluster)) +
#   theme_bw() +
#   labs(x = "", y = "Annual antibiotic consumption\n(DDD per 1,000 bed-days)",
#        linetype = "Antibiotic class", col = "Hospital type")
# 
# p3 = type_atb_norm %>%
#   mutate(cluster = cutree(clust_fit, k= 4)) %>%
#   group_by(cluster, atb_class) %>%
#   summarise(n = n(), .groups = "drop") %>%
#   group_by(cluster) %>%
#   mutate(p = n/sum(n)) %>%
#   ungroup() %>%
#   ggplot(., aes(x = cluster, y = p, fill = atb_class, label = n)) +
#   geom_bar(stat = "identity") +
#   geom_text(position = position_stack(vjust = 0.5)) +
#   theme_bw() +
#   labs(x = "Clusters", y = "Repartition", fill = "Antibiotic class")
# 
# 
# p4 = type_atb_norm %>%
#   mutate(cluster = cutree(clust_fit, k= 4)) %>%
#   group_by(cluster, type) %>%
#   summarise(n = n(), .groups = "drop") %>%
#   group_by(cluster) %>%
#   mutate(p = n/sum(n)) %>%
#   ungroup() %>%
#   ggplot(., aes(x = cluster, y = p, fill = type, label = n)) +
#   geom_bar(stat = "identity") +
#   geom_text(position = position_stack(vjust = 0.5)) +
#   theme_bw() +
#   labs(x = "Clusters", y = "Repartition", fill = "Hospital type")
# 
# pdata = ggarrange(p1,p2, nrow = 2, common.legend = T, legend = "right")
# ggsave("plots/cons_atb/hierarchical_clusters_type.png", pdata,
#        height = 7, width = 10)
# pcomp = ggarrange(p3,p4, nrow = 2, align = "v")
# ggsave("plots/cons_atb/hierarchical_comp_type.png", pcomp,
#        height = 7, width = 10)

# # K-means clustering on normalized data-------------------------
# # Normalized data
# hospital_atb_norm = atb %>%
#   filter(secteur == "Total établissement") %>%
#   group_by(code, Date_year, atb_class, region, type) %>%
#   summarise(consumption = sum(molDDD) / unique(Nbhosp) * 1000, .groups = "drop") %>%
#   group_by(code, type, atb_class, region) %>%
#   mutate(consumption = (consumption-mean(consumption))/sd(consumption)) %>%
#   ungroup() %>%
#   pivot_wider(names_from = Date_year, values_from = consumption, names_prefix = "year_") %>%
#   filter(!if_any(matches("year_"), is.na))
# 
# hospital_atb_norm %>%
#   filter(if_any(matches("year_"), is.na)) %>%
#   nrow(.)
# 
# # Plot of all data - regional consumption by antibiotic class
# p1 = hospital_atb_norm %>%
#   pivot_longer(matches("year_"), names_to = "year", values_to = "consumption") %>%
#   mutate(year= gsub("year_", "", year)) %>%
#   ggplot(., aes(x = year, y = consumption, group = interaction(code, atb_class),
#                 col = region, linetype = atb_class)) +
#   geom_line() +
#   theme_bw() +
#   labs(x = "", y = "Antibiotic consumption\n(DDD/1,000 bed-days)", 
#        col = "Region",
#        linetype = "Antibiotic class", title = "All data")
# 
# # Define number of clusters
# df = select(hospital_atb_norm, -c(region, type, atb_class, code))
# p21 = fviz_nbclust(df, kmeans, method = "wss", k.max = 13) +
#   labs(title = "")
# p22 = fviz_nbclust(df, kmeans, method = "silhouette", k.max = 13) +
#   labs(title = "")
# p23 = fviz_gap_stat(cluster::clusGap(df, FUN = kmeans, nstart = 20, K.max = 13, B = 100)) +
#   labs(title = "")
# ggarrange(p21, p22, p23, nrow = 1)
# 
# NbClust(data = df, method = "kmeans")
# 
# # Add cluster assignation  
# clusters = pam(df, k = 3, medoids = "random", nstart = 2)
# hospital_atb_norm_long = hospital_atb_norm %>% 
#   mutate(cluster = clusters$clustering) %>% 
#   pivot_longer(matches("year_"), names_to = "year", values_to = "consumption") %>%
#   mutate(year = gsub("year_", "", year))
# 
# centers_hosp_long = data.frame(clusters$medoids) %>%
#   rownames_to_column(var = "cluster") %>%
#   pivot_longer(cols = -cluster, names_to = "year", values_to = "consumption") %>%  
#   mutate(year = gsub("year_", "", year))
# 
# p3 = ggplot() +
#   geom_line(data = hospital_atb_norm_long, 
#             aes(y = consumption, x = year, group = interaction(code, atb_class),
#                 col = region, linetype = atb_class)) +
#   facet_wrap(facets = vars(cluster), nrow = 1) + 
#   geom_line(data = centers_hosp_long, aes(y = consumption, x = year, group = cluster), 
#             col = "black", linewidth = 2) +
#   theme_bw() +
#   labs(x = "", y = "Antibiotic consumption\n(DDD/1,000 bed-days)", 
#        col = "Region",
#        linetype = "Antibiotic class", title = "K-means clustering")
# p3
# 
# # Cluster composition
# p4 = hospital_atb_norm %>%
#   mutate(cluster = clusters$clustering) %>%
#   group_by(cluster, region) %>%
#   summarise(n = n(), .groups = "drop") %>%
#   group_by(cluster) %>%
#   mutate(p = n/sum(n)) %>%
#   ungroup() %>%
#   ggplot(., aes(x = cluster, y = p, fill = region, label = n)) +
#   geom_bar(stat = "identity") +
#   geom_text(position = position_stack(vjust = 0.5)) +
#   theme_bw() +
#   labs(x = "Clusters", y = "Repartition", fill = "Region")
# 
# p5 = hospital_atb_norm %>% 
#   mutate(cluster = clusters$clustering) %>%
#   group_by(cluster, atb_class) %>%
#   summarise(n = n(), .groups = "drop") %>%
#   group_by(cluster) %>%
#   mutate(p = n/sum(n)) %>%
#   ungroup() %>%
#   ggplot(., aes(x = cluster, y = p, fill = atb_class, label = n)) +
#   geom_bar(stat = "identity") +
#   geom_text(position = position_stack(vjust = 0.5)) +
#   theme_bw() +
#   labs(x = "Clusters", y = "Repartition", fill = "Antibiotic class")
# 
# p6 = hospital_atb_norm %>% 
#   mutate(cluster = clusters$clustering) %>%
#   group_by(cluster, type) %>%
#   summarise(n = n(), .groups = "drop") %>%
#   group_by(cluster) %>%
#   mutate(p = n/sum(n)) %>%
#   ungroup() %>%
#   ggplot(., aes(x = cluster, y = p, fill = type, label = n)) +
#   geom_bar(stat = "identity") +
#   geom_text(position = position_stack(vjust = 0.5)) +
#   theme_bw() +
#   labs(x = "Clusters", y = "Repartition", fill = "Hospital type")
# 
# ggarrange(p4, p5, p6, align = "v", nrow = 3)
# 
# # Plot all data at the same time
# p_clusters = grid.arrange(
#   arrangeGrob(
#     arrangeGrob(
#       arrangeGrob(p1+theme(legend.position = "hidden"), NULL, NULL, ncol = 3),
#       arrangeGrob(p21, p22, p23, ncol = 3),
#       nrow = 2
#     ),
#     p3+theme(legend.position="hidden"),nrow = 2,  heights = c(1.5,1)
#   ),
#   get_legend(p3), ncol = 2, widths = c(3,1)
# )
# ggsave("plots/cons_atb/medoids_hospital.png", p_clusters, 
#        height = 10, width = 10)
# 
# p_comp = ggarrange(p4, p5, p6, align = "v", nrow = 3)
# ggsave("plots/cons_atb/medoids_hospital_composition.png", p_comp, 
#        height = 11, width = 8)

