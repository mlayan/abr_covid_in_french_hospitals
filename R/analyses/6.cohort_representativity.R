##################################################
## EVALUATION OF THE REPRESENTATIVITY OF THE 
## COHORT
##################################################
library(tidyverse)
library(rstatix)
library(readxl)
library(haven)
library(DescTools)
library(RColorBrewer)
library(ggpubr)
library(gtsummary)
rm(list = ls())

load("data/cohort_final.rda")
load("data/metadata_beds.rda")
load("data/metadata_admin_espic.rda")
# load("data/metadata_all.rda")
source("R/helper/dictionaries.R")
source("R/helper/helper_functions.R")

# Load data on finess hospital by region in 2020
# Data found on ATIH website: https://www.scansante.fr/applications/analyse-activite-regionale
region_finess_all = data.frame() 
for (f in list.files("data-raw/atih/Hospital_activity_2020/", full.names = T)) {
  # Load and concat data 
  region_mco = read_excel(f, sheet = "MCO_finess", skip = 2) %>%
    dplyr::select(`Finess PMSI`, `Raison sociale`, `Finess de transmission`) %>%
    mutate(type = "MCO")
  region_ssr = read_excel(f, sheet = "SSR_finess", skip = 2) %>%
    dplyr::select(`Finess PMSI`, `Raison sociale`, `Finess de transmission`) %>%
    mutate(type = "SSR")
  region_finess = bind_rows(region_mco, region_ssr)
  
  # Add region name
  region_finess$region = case_when(
    grepl("_84", f) ~ "Auvergne-Rhône-Alpes",
    grepl("_27", f) ~ "Bourgogne-Franche-Comté",
    grepl("_53", f) ~ "Bretagne",
    grepl("_24", f) ~ "Centre-Val de Loire",
    grepl("_44", f) ~ "Grand-Est",
    grepl("_32", f) ~ "Hauts-de-France",
    grepl("_11", f) ~ "Île-de-France",
    grepl("_28", f) ~ "Normandie",
    grepl("_75", f) ~ "Nouvelle-Aquitaine",
    grepl("_76", f) ~ "Occitanie",
    grepl("_52", f) ~ "Pays de la Loire",
    grepl("_93", f) ~ "Provence-Alpes-Côte d'Azur"
  )
  
  # Duplicate rows if multiple transmission finess 
  rows_to_duplicate = region_finess[grepl("; ", region_finess$`Finess de transmission`),]
  if (nrow(rows_to_duplicate) > 0) {
    for (r in split(rows_to_duplicate, seq(nrow(rows_to_duplicate)))) {
      aggr_finess = r$`Finess de transmission`
      
      # Create as many rows as finess numbers 
      finess_l = strsplit(aggr_finess, ";")[[1]]
      finess_l = gsub(" ", "", finess_l)
      rows_duplicated = bind_rows(replicate(length(finess_l), r, simplify = F))
      rows_duplicated$`Finess de transmission` = finess_l
      
      # Add duplicated rows
      region_finess = region_finess %>%
        filter(`Finess de transmission` != aggr_finess) %>%
        bind_rows(., rows_duplicated)
    }
  }
  
  # Get final database
  region_finess_all = bind_rows(region_finess_all, region_finess)
  rm(region_finess)
}

################################################################################
# Hospital activity representativity in 2020 - MCO and SSR
################################################################################
# Load data on the number of hospitalisation days - MCO
region_mco= data.frame() 
for (f in list.files("data-raw/atih/Hospital_activity_2020/", full.names = T)) {
  
  if (grepl("_84", f)) region_bd_tmp = read_excel(f, range = "MCO_region!B77:M88") %>%
      mutate(region = "Auvergne-Rhône-Alpes")
  
  if (grepl("_27", f)) region_bd_tmp = read_excel(f, range = "MCO_region!B72:M83") %>%
      mutate(region = "Bourgogne-Franche-Comté")
  
  if (grepl("_53", f)) region_bd_tmp = read_excel(f, range = "MCO_region!B69:M80") %>%
      mutate(region = "Bretagne")
  
  if (grepl("_24", f)) region_bd_tmp = read_excel(f, range = "MCO_region!B67:M78") %>%
      mutate(region = "Centre-Val de Loire")
  
  if (grepl("_44", f)) region_bd_tmp = read_excel(f, range = "MCO_region!B74:M85") %>%
      mutate(region = "Grand-Est")
  
  if (grepl("_32", f)) region_bd_tmp = read_excel(f, range = "MCO_region!B85:M96") %>%
      mutate(region = "Hauts-de-France")
  
  if (grepl("_11", f)) region_bd_tmp = read_excel(f, range = "MCO_region!B84:M95") %>%
      mutate(region = "Île-de-France")
  
  if (grepl("_28", f)) region_bd_tmp = read_excel(f, range = "MCO_region!B69:M80") %>%
      mutate(region = "Normandie")
  
  if (grepl("_75", f)) region_bd_tmp = read_excel(f, range = "MCO_region!B73:M84") %>%
      mutate(region = "Nouvelle-Aquitaine")
  
  if (grepl("_76", f)) region_bd_tmp = read_excel(f, range = "MCO_region!B116:M127") %>%
      mutate(region = "Occitanie")
  
  if (grepl("_52", f)) region_bd_tmp = read_excel(f, range = "MCO_region!B66:M77") %>%
      mutate(region = "Pays de la Loire")
  
  if (grepl("_93", f)) region_bd_tmp = read_excel(f, range = "MCO_region!B67:M78") %>%
      mutate(region = "Provence-Alpes-Côte d'Azur")
  
  # Get final database
  region_mco = bind_rows(region_mco, region_bd_tmp)
  rm(region_bd_tmp)
}

colnames(region_mco) = c(
  "sector",
  "bd_2020",
  "bd_prop_2020",
  "bd_change_2019_2020",
  "bd_contribution_to_change",
  "bd_change_france",
  "stays_2020",
  "stays_prop_2020",
  "stays_change_2019_2020",
  "stays_contribution_to_change",
  "stays_change_france",
  "activity_rate",
  "region"
)

atih_mco = region_mco %>% 
  group_by(region) %>%
  summarise(bed_mco_atih = sum(bd_2020), .groups = "drop")

# Load SPARES data - MCO
spares_mco = metadata_beds %>%
  filter(year == 2020, code %in% cohort_final, !secteur %in% c("Total", "Rehabilitation care")) %>%
  left_join(., metadata_admin_espic %>% dplyr::select(code, region) %>% distinct(), by = "code") %>%
  group_by(region) %>%
  summarise(bed_mco_spares = sum(nbjh), .groups = "drop") %>%
  left_join(., atih_mco, by = "region")

# Load data on the number of hospitalisation days - SSR
atih_ssr = data.frame() 
for (f in list.files("data-raw/atih/Hospital_activity_2020/", full.names = T)) {
  
  if (grepl("_84", f)) region_bd_tmp = read_excel(f, range = "SSR_region!C61:C62") %>%
      mutate(region = "Auvergne-Rhône-Alpes")
  
  if (grepl("_27", f)) region_bd_tmp = read_excel(f, range = "SSR_region!C56:C57") %>%
      mutate(region = "Bourgogne-Franche-Comté")
  
  if (grepl("_53", f)) region_bd_tmp = read_excel(f, range = "SSR_region!C53:C54") %>%
      mutate(region = "Bretagne")
  
  if (grepl("_24", f)) region_bd_tmp = read_excel(f, range = "SSR_region!C51:C52") %>%
      mutate(region = "Centre-Val de Loire")
  
  if (grepl("_44", f)) region_bd_tmp = read_excel(f, range = "SSR_region!C58:C59") %>%
      mutate(region = "Grand-Est")
  
  if (grepl("_32", f)) region_bd_tmp = read_excel(f, range = "SSR_region!C69:C70") %>%
      mutate(region = "Hauts-de-France")
  
  if (grepl("_11", f)) region_bd_tmp = read_excel(f, range = "SSR_region!C68:C69") %>%
      mutate(region = "Île-de-France")
  
  if (grepl("_28", f)) region_bd_tmp = read_excel(f, range = "SSR_region!C53:C54") %>%
      mutate(region = "Normandie")
  
  if (grepl("_75", f)) region_bd_tmp = read_excel(f, range = "SSR_region!C57:C58") %>%
      mutate(region = "Nouvelle-Aquitaine")
  
  if (grepl("_76", f)) region_bd_tmp = read_excel(f, range = "SSR_region!C100:C101") %>%
      mutate(region = "Occitanie")
  
  if (grepl("_52", f)) region_bd_tmp = read_excel(f, range = "SSR_region!C50:C51") %>%
      mutate(region = "Pays de la Loire")
  
  if (grepl("_93", f)) region_bd_tmp = read_excel(f, range = "SSR_region!C51:C52") %>%
      mutate(region = "Provence-Alpes-Côte d'Azur")
  
  # Get final database
  atih_ssr = bind_rows(atih_ssr, region_bd_tmp)
  rm(region_bd_tmp)
}
colnames(atih_ssr) = c("bd_2020", "region")

# Load SPARES data - SSR
spares_ssr = metadata_beds %>%
  filter(year == 2020, code %in% cohort_final, secteur == "Rehabilitation care") %>%
  left_join(., metadata_admin_espic %>% dplyr::select(code, region) %>% distinct(), by = "code") %>%
  group_by(region) %>%
  summarise(bed_ssr_spares = sum(nbjh), .groups = "drop") %>%
  left_join(., atih_ssr %>% rename(bed_ssr_atih = bd_2020), by = "region")

# Merge databases 
beds = left_join(spares_ssr, spares_mco, by = c("region")) %>% 
  mutate(Counts_Cohort = bed_ssr_spares + bed_mco_spares, Counts_ATIH = bed_mco_atih + bed_ssr_atih) %>%
  mutate(Distribution_Cohort = Counts_Cohort / sum(Counts_Cohort), Distribution_ATIH = Counts_ATIH / sum(Counts_ATIH))

chisq.test(x=beds$Counts_Cohort, p=beds$Distribution_ATIH)

# Plot counts and distributions
p1=beds %>%
  dplyr::select(region, Counts_Cohort, Counts_ATIH, Distribution_Cohort, Distribution_ATIH) %>%
  pivot_longer(-region, names_to = "param", values_to = "value") %>%
  separate(param, into = c("param", "db"), sep = "_") %>%
  mutate(region = recode(region, !!!dict_regions)) %>%
  ggplot(., aes(x = region, y = value, fill = db)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.5), width = 0.5) +
  facet_grid(rows = vars(param), scales = "free_y") +
  scale_fill_manual(values =  c("Cohort" = "coral", "ATIH" = "darkorchid")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(fill = "", x = "", y = "Bed-days in 2020")

################################################################################
# Hospital counts representativity
################################################################################
# Get all finess available in France and verify that ATIH hospital finess 
# are in the national database (coherence between PMSI and national finess)
all_finess = read.csv2("data-raw/spares/hospital_location/etalab-cs1100502-stock-20240110-0338.csv", 
                       header = F)
colnames(all_finess) = c("section", "finess", "finess_juridique", "rs", "rslongue",
                         "complrs", "compldistrib", "numvoie", "typvoie", "voie",
                         "compvoie", "lieuditbp", "commune", "departement", 
                         "libdepartement", "ligneacheminement", "telephone",
                         "telecopie", "categetab", "libcategetab", "categagretab",
                         "libcategagretab", "siret", "codeape", "codemft", "libmft",
                         "codesph", "libsph", "dateouv", "dateautor", "datemaj", "numuai")

# Compare national and ATIH finess
length(unique(region_finess_all$`Finess de transmission`))
sum(unique(region_finess_all$`Finess de transmission`) %in% c(all_finess$finess, all_finess$finess_juridique))
region_finess_all[!region_finess_all$`Finess de transmission` %in% c(all_finess$finess, all_finess$finess_juridique),]

# Get number of hospitals used in ATIH in 2020
n_hospitals = data.frame()
for (f in list.files("data-raw/atih/Hospital_activity_2020/", full.names = T)) {
  
  n_tmp = bind_rows(
    read_excel(f, sheet = "SSR_finess", skip = 2)[,1:3] %>% mutate(sector = "SSR"),
    read_excel(f, sheet = "MCO_finess", skip = 2)[,1:3] %>% mutate(sector = "MCO")
  ) %>%
    dplyr::select(`Finess PMSI`, `Raison sociale`, `Finess de transmission`, sector) %>%
    rename(finess_pmsi = `Finess PMSI`, rs = `Raison sociale`, finess = `Finess de transmission`) %>%
    mutate(region = case_when(
      grepl("_84", f) ~ "Auvergne-Rhône-Alpes",
      grepl("_27", f) ~ "Bourgogne-Franche-Comté",
      grepl("_53", f) ~ "Bretagne",
      grepl("_24", f) ~ "Centre-Val de Loire",
      grepl("_44", f) ~ "Grand-Est",
      grepl("_32", f) ~ "Hauts-de-France",
      grepl("_11", f) ~ "Île-de-France", 
      grepl("_28", f) ~ "Normandie",
      grepl("_75", f) ~ "Nouvelle-Aquitaine", 
      grepl("_76", f) ~ "Occitanie",
      grepl("_52", f) ~ "Pays de la Loire", 
      grepl("_93", f) ~ "Provence-Alpes-Côte d'Azur"
      )
    )

  # Get final database
  n_hospitals = bind_rows(n_hospitals, n_tmp)
  rm(n_tmp)
}

# Create new rows for finess that have multiple transmission finess
n_hospitals_test = n_hospitals %>%
  separate(., finess, into = c("finess1", "finess2"), sep = "; ", fill = "right") %>%
  pivot_longer(c(finess1, finess2), names_to = "to_discard", values_to = "finess") %>%
  filter(!is.na(finess)) %>%
  dplyr::select(-c(to_discard, sector)) %>%
  distinct()

# Check which finess is available in the official data 
nrow(n_hospitals_test)
length(unique(n_hospitals_test$finess))
length(unique(n_hospitals_test$finess_pmsi))
sum(n_hospitals_test$finess_pmsi %in% chu_finess_juridique$finess_juridique)

sum(unique(n_hospitals_test$finess) %in% c(all_finess$finess, all_finess$finess_juridique))
sum(unique(n_hospitals_test$finess_pmsi) %in% c(all_finess$finess, all_finess$finess_juridique))

n_hospitals_test %>%
  mutate(
    finess_transmission_official = ifelse(finess %in% c(all_finess$finess, all_finess$finess_juridique), 1, 0),
    finess_pmsi_official = ifelse(finess_pmsi %in% c(all_finess$finess, all_finess$finess_juridique), 1, 0)
    ) %>%
  filter(finess_pmsi_official != finess_transmission_official)

# Consider finess PMSI -> one finess = one geographic hospital
# Except for CHU
atih_nb = n_hospitals_test %>%
  dplyr::select(finess_pmsi, region) %>%
  distinct() %>%
  group_by(region) %>%
  summarise(Counts_ATIH = n(), .groups = "drop")

# Number of hospitals in SPARES 
spares_nb = metadata_admin_espic %>%
  filter(code %in% cohort_final) %>%
  dplyr::select(code, finess_juridique, region) %>%
  distinct() %>%
  group_by(region) %>%
  summarise(Counts_Cohort = n(), .groups = "drop")

# Plot distributions
p2 = left_join(spares_nb, atih_nb, by = "region") %>%
  mutate(Distribution_Cohort = Counts_Cohort / sum(Counts_Cohort), 
         Distribution_ATIH = Counts_ATIH / sum(Counts_ATIH),
         region = recode(region, !!!dict_regions)) %>%
  pivot_longer(-region, names_to = "param", values_to = "value") %>%
  separate(param, into = c("param", "db"), sep = "_") %>%
  ggplot(., aes(x = region, y = value, fill = db)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.5), width = 0.5) +
  facet_grid(rows = vars(param), scales = "free_y") +
  theme_bw() +
  scale_fill_manual(values = c("Cohort" = "coral", "ATIH" = "darkorchid")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = "", y = "Legal entities in 2020", fill = "")

# Chi2 test of adequation
chisq.test(spares_nb$Counts_Cohort, 
  p = atih_nb$Counts_ATIH / sum(atih_nb$Counts_ATIH))

################################################################################
# Infection prevalence representativity 
################################################################################
# Save finess number of final cohort to retrieve denominator from PMSI
finess_point_prevalence_denominator = metadata_admin_espic %>%
  filter(code %in% cohort_final) %>%
  dplyr::select(finess) %>%
  distinct() %>%
  .$finess
writeLines(text = paste0('"', paste0(finess_point_prevalence_denominator, collapse = '","'), '"'), 
           con = "data-raw/atih/finess_point_prevalence.txt")

# Load denominator from ATIH
# It corresponds to the number of patients that were hospitalized in 
# 2019, 2020, 2021 or 2022 between May 15 and June 30 
denominator = read_sas("data-raw/atih/point_prevalence_denominateur.sas7bdat") %>%
  left_join(., metadata_admin_espic %>% dplyr::select(finess, region) %>% distinct(), by = c("FinessGeo"="finess")) %>%
  filter(!is.na(region), Date_year == 2022) %>%
  group_by(region) %>%
  summarise(npatients = sum(npatients), .groups = "drop")

metadata_admin_espic %>%
  filter(code %in% cohort_final, !finess %in% read_sas("data-raw/atih/point_prevalence_denominateur.sas7bdat")$FinessGeo) %>%
  dplyr::select(finess_juridique, name) %>%
  distinct()

metadata_admin_espic %>%
  filter(code %in% cohort_final, finess %in% read_sas("data-raw/atih/point_prevalence_denominateur.sas7bdat")$FinessGeo,
         type == "University hospital") %>%
  dplyr::select(finess_juridique, name) %>%
  distinct()

# Load final sample data 
# Get patients with incident HAI at the day level 
# If a patient has several incident HAI on the same day, we count as 1 patient  
all_bacterial_samples = read.table("data-raw/spares/combined/resistance_cohortfinal_representativity_2022.txt", header = T, sep = "\t") %>%
  dplyr::select(-c(molecule, atb_class)) %>%
  distinct() %>%
  dplyr::select(code, Date_day, Num_Patient) %>%
  distinct()
  
# Get average prevalence by region over the study period
spares_p = all_bacterial_samples %>%
  inner_join(., metadata_admin_espic %>% dplyr::select(code, region) %>% distinct(), by = "code") %>%
  group_by(region) %>%
  summarise(n = n(), .groups = "drop") %>%
  left_join(., denominator, by = "region") %>%
  mutate(p = n/npatients*100, database = "Cohort") %>%
  arrange(region)
spares_p$pmin = mapply(function(x,n) 100*prop.test(x=x,n=n,conf.level = 0.95, alternative = "two.sided")$conf.int[1],
                        spares_p$n, spares_p$npatients)
spares_p$pmax = mapply(function(x,n) 100*prop.test(x=x,n=n,conf.level = 0.95, alternative = "two.sided")$conf.int[2],
                        spares_p$n, spares_p$npatients)

# Get estimates from national prevalence survey
# These estimates do not account for nosocomial Covid-19 cases
enp = data.frame(
  region = c("Hauts-de-France", "Normandie", "Bretagne", "Île-de-France", "Grand-Est", "Pays de la Loire", "Centre-Val de Loire", "Bourgogne-Franche-Comté", "Nouvelle-Aquitaine", "Auvergne-Rhône-Alpes", "Occitanie", "Provence-Alpes-Côte d'Azur"), 
  p = c(5.27, 5.89, 3.96, 5.28, 6.32, 4.73, 4.05, 5.6, 5.06, 5.91, 5.56, 5.26),
  pmin = c(4.38, 5.03, 2.98, 4.78, 5.8, 4.23, 2.74, 4.97, 4.36, 5.45, 4.68, 4.07),
  pmax = c(6.34, 6.89, 5.25, 5.83, 6.88, 5.27, 5.95, 6.3, 5.87, 6.4, 6.58, 6.78),
  database = "SpF"
  ) %>%
  arrange(region)

# Comparison with SPARES data from final cohort
p3 = bind_rows(spares_p, enp) %>%
  mutate(region = recode(region, !!!dict_regions),
         database = factor(database, c("SpF", "Cohort"))) %>%
  ggplot(., aes(x = region, y = p, ymin = pmin, ymax = pmax, col = database)) +
  geom_pointrange(position = position_dodge(width = 0.2), fatten = 2) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_color_manual(values = c("Cohort"="coral", "SpF"="steelblue3")) +
  labs(x = "", y = "Prevalence of HAI (SpF) or\nbacterial episodes (Cohort), 95% CI", col = "") +
  expand_limits(y = 0)

p4 = bind_rows(spares_p, enp) %>%
  group_by(database) %>%
  mutate(dist = p/sum(p), 
         region = recode(region, !!!dict_regions),
         database = factor(database, c("SpF", "Cohort"))
         ) %>%
  ungroup() %>%
  ggplot(., aes(x = region, y = dist, fill = database)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.5), width = 0.5) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_manual(values = c("Cohort" = "coral", "SpF" = "steelblue3")) +
  labs(x = "", y = "Distribution of prevalence", fill = "") 

################################################################################
# Final supplementary figures 
################################################################################
# Final figure on the representativity of the cohort
p = ggarrange(
  ggarrange(p1, p2, ncol = 2, align = "hv", common.legend = T, legend = "right",
            labels = c("A", "B")),
  ggarrange(p3, p4, ncol = 2, align = "hv", common.legend = T, legend = "right",
            labels = c("C", "")),
  heights = c(1, 0.7),
  nrow = 2
)
p
ggsave("../Paper/Supplementary/hospital_cohort_representativity.png", 
       p, height = 7, width = 10)
ggsave("plots/cohort_description/hospital_cohort_representativity.png", 
       p, height = 7, width = 10)

################################################################################
# Comparison of included and excluded hospitals in terms of AMR
################################################################################
# All hospitals in hexagonal France--------------------------------------------- 
all_codes = bind_rows(
  read_excel("data-raw/spares/2019/SPARES_2019_BMRCovid.xlsx", sheet = "ADM", col_types = c("numeric", "text", "skip", "skip", "text", "text", "text")) %>% rename(code = IdEtablissement, type = groupe, Region = `Nouvelle-Region`, `finess géographique` = finess),
  read_excel("data-raw/spares/2020/atb2020_2021.xlsx", sheet = "Participants ES_2020", col_types = c("numeric", "text", "text", "skip", "skip", "text", "text")), 
  read_excel("data-raw/spares/2020/atb2020_2021.xlsx", sheet = "Participants ES_2021", col_types = c("numeric", "text", "text", "skip", "skip", "text", "text")),
  read_excel("data-raw/spares/2022/ATB2022.xlsx", sheet = "Participants ES_2022", col_types = c("numeric", "text", "text", "skip", "skip", "text", "text"))
) %>% 
  mutate(
    Region = case_when(
      Region == "Auvergne - Rhône Alpes" ~ "Auvergne-Rhône-Alpes",
      Region == "Bourgogne - Franche Comté" ~ "Bourgogne-Franche-Comté",
      Region == "Grand Est" ~ "Grand-Est",                 
      Region == "Hauts de France" ~ "Hauts-de-France",
      Region == "Ile de France" ~ "Île-de-France",
      Region == "Nouvelle Aquitaine" ~ "Nouvelle-Aquitaine",
      Region == "Pays de Loire" ~ "Pays de la Loire",
      Region == "Provence Alpes Côte d'Azur" ~ "Provence-Alpes-Côte d'Azur",
      Region == "Reunion - Mayotte" ~ "La Réunion-Mayotte",
      .default = Region
    ),
    type = case_when(code == 11014 ~ "CLCC",
                     code == 2013 ~ "MCO",
                     code == 1913 ~ "MCO", # physical Finess code: 210011847
                     code == 11180 ~ "MCO", # physical Finess code: 420000192
                     .default = type)
  ) %>%
  filter(
    !Region %in% c("Corse", "Guadeloupe", "Guyane", "Martinique", "Nouvelle Calédonie", "Reunion - Mayotte"),
    !type %in% c("CLCC", "ESLD", "HIA", "PSY")
  ) %>%
  distinct(code) %>%
  pull(code)

length(all_codes)

# Antibiotic use----------------------------------------------------------------
# Antibiotics included in the final analysis
load("data/atb_use_hospital_level.rda")
load("data/included_molecules.rda")
load("data/icu_cohort_final.rda")

# Antibiotic use all hospitals 
antibiotic_use = bind_rows(
  read_excel("data-raw/spares/2019/ATB2019_NAT.xlsx", sheet = "ATB2019") %>% mutate(Date_year = 2019), 
  read_excel("data-raw/spares/2020/atb2020_2021.xlsx", sheet = "ATB2020") %>% mutate(Date_year = 2020), 
  read_excel("data-raw/spares/2020/atb2020_2021.xlsx", sheet = "ATB2021") %>% mutate(Date_year = 2021), 
  read_excel("data-raw/spares/2022/ATB2022.xlsx", sheet = "ATB2022") %>% mutate(Date_year = 2022)
  ) %>%
  left_join(., included_molecules %>% mutate(included = 1), by = c("ATC", "DCI")) %>%
  filter(
    secteur %in% c("Chirurgie", "Gynécologie-Obstétrique", "Médecine", "Réanimation", "SSR"),
    !is.na(included)
    ) %>%
  mutate(
    cohort = code %in% cohort_final,
    excluded = code %in% all_codes
  )
  
# Compare annual antibiotic use in hospitals
bind_rows(
  antibiotic_use %>%
    group_by(Date_year, code, secteur) %>%
    summarise(molDDD = sum(molDDD), Nbhosp = unique(Nbhosp), .groups = "drop") %>%
    mutate(hospital = "all hospitals", atb_class = "Total"),
  antibiotic_use %>%
    filter(DCI %in% c("Imipenem","Meropenem", "Vancomycine", "Azithromycine") | ATC %in% vec_3gc) %>%
    mutate(atb_class = case_when(
      DCI %in% c("Meropenem", "Imipenem") ~ "Imipenem + Meropenem",
      ATC %in% vec_3gc ~ "3GC",
      .default = gsub("e$", "", DCI))
      ) %>%
    group_by(Date_year, atb_class, code, secteur) %>%
    summarise(molDDD = sum(molDDD), Nbhosp = unique(Nbhosp), .groups = "drop") %>%
    mutate(hospital = "all hospitals"),
  antibiotic_use %>%
    group_by(Date_year, code, secteur, atb_class) %>%
    summarise(molDDD = sum(molDDD), Nbhosp = unique(Nbhosp), .groups = "drop") %>%
    mutate(hospital = "all hospitals"),
  atb_use_hospital_level %>%
    mutate(hospital = "cohort", atb_class = ifelse(atb_class == "Third generation Cephalosporins", "3GC", atb_class))
  ) %>%
  group_by(Date_year, hospital, atb_class) %>%
  summarise(atb_use = sum(molDDD) / sum(Nbhosp) * 1000, .groups = "drop") %>%
  mutate(atb_class = factor(atb_class,
                            c("Aminoglycosides", "Carbapenems", "Imipenem + Meropenem", "Cephalosporins",
                              "3GC", "Fosfomycin", "Glycopeptides", "Vancomycin",
                              "Lipopeptides", "Macrolides", "Azithromycin", "Monobactams", "Oxazolidinones",
                              "Penicillins", "Polymyxins", "Quinolones", "Tetracyclines", "Trimethoprim",
                              "Total")
                            )) %>%
  ggplot(., aes(x = Date_year, y = atb_use, fill = hospital)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.5), width = 0.5) +
  facet_wrap(facets = vars(atb_class), ncol = 5, scales = "free_y") +
  theme_bw() +
  scale_fill_manual(values = c("all hospitals" = "grey80", "cohort" = "cornflowerblue")) +
  theme(axis.title.x = element_blank()) +
  labs(y = "Annual antibiotic use\n(in DDD for 1,000 occupied bed-days)")
ggsave("plots/antibiotic_use/total_use_sensitivity_cohort.png", height = 7, width = 11)

# Compare annual total antibiotic use
bind_rows(
  antibiotic_use %>%
    filter(secteur == "Réanimation") %>%
    group_by(Date_year, code) %>%
    summarise(molDDD = sum(molDDD), Nbhosp = unique(Nbhosp), .groups = "drop") %>%
    mutate(hospital = "all hospitals", atb_class = "Total"),
  antibiotic_use %>%
    filter(secteur == "Réanimation", DCI %in% c("Imipenem","Meropenem", "Vancomycine", "Azithromycine") | ATC %in% vec_3gc) %>%
    mutate(atb_class = case_when(
      DCI %in% c("Meropenem", "Imipenem") ~ "Imipenem + Meropenem",
      ATC %in% vec_3gc ~ "3GC",
      .default = gsub("e$", "", DCI))
    ) %>%
    group_by(Date_year, atb_class, code) %>%
    summarise(molDDD = sum(molDDD), Nbhosp = unique(Nbhosp), .groups = "drop") %>%
    mutate(hospital = "all hospitals"),
  antibiotic_use %>%
    filter(secteur == "Réanimation") %>%
    group_by(Date_year, code, atb_class) %>%
    summarise(molDDD = sum(molDDD), Nbhosp = unique(Nbhosp), .groups = "drop") %>%
    mutate(hospital = "all hospitals"),
  atb_use_hospital_level %>%
    filter(secteur == "ICU", code %in% icu_cohort_final) %>%
    mutate(hospital = "cohort", atb_class = ifelse(atb_class == "Third generation Cephalosporins", "3GC", atb_class))
) %>%
  group_by(Date_year, hospital, atb_class) %>%
  summarise(atb_use = sum(molDDD) / sum(Nbhosp) * 1000, .groups = "drop") %>%
  mutate(atb_class = factor(atb_class,
                            c("Aminoglycosides", "Carbapenems", "Imipenem + Meropenem", "Cephalosporins",
                              "3GC", "Fosfomycin", "Glycopeptides", "Vancomycin",
                              "Lipopeptides", "Macrolides", "Azithromycin", "Monobactams", "Oxazolidinones",
                              "Penicillins", "Polymyxins", "Quinolones", "Tetracyclines", "Trimethoprim",
                              "Total")
  )) %>%
  ggplot(., aes(x = Date_year, y = atb_use, fill = hospital)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.5), width = 0.5) +
  facet_wrap(facets = vars(atb_class), ncol = 5, scales = "free_y") +
  theme_bw() +
  scale_fill_manual(values = c("all hospitals" = "grey80", "cohort" = "cornflowerblue")) +
  theme(axis.title.x = element_blank()) +
  labs(y = "Annual antibiotic use\n(in DDD for 1,000 occupied bed-days)")
ggsave("plots/antibiotic_use/total_use_icu_sensitivity_cohort.png", height = 7, width = 11)

# Resistant isolates------------------------------------------------------------
# Cohort data
load("data/res_hospital.rda")
load("data/res_icu.rda")

# All hospitals
bact_names = c("Pseudomonas aeruginosa" = "CR P. aeruginosa", "Enterobacter cloacae complex" = "ESBL E. cloacae", 
               "Escherichia coli" = "ESBL E. coli", "Klebsiella pneumoniae" = "ESBL K. pneumoniae", 
               "Staphylococcus aureus" = "MRSA")

isolates_all = bind_rows(
  read.table("data-raw/spares/combined/resistance_allhospitals_2019.txt", sep = "\t", header = T),
  read.table("data-raw/spares/combined/resistance_allhospitals_2020.txt", sep = "\t", header = T),
  read.table("data-raw/spares/combined/resistance_allhospitals_2021.txt", sep = "\t", header = T),
  read.table("data-raw/spares/combined/resistance_allhospitals_2022.txt", sep = "\t", header = T)
) %>%
  group_by(code, site, Date_day, ligne, Num_Patient, age_cat, bacterie, secteur, Date_year) %>%
  summarise(Resultat = ifelse(any(Resultat %in% c("O", "R")), "R", "S"), .groups = "drop") %>%
  group_by(Date_year, bacterie) %>%
  summarise(
    n_res_icu = sum(Resultat %in% "R" & secteur == "ICU"), 
    n_tot_icu = sum(secteur == "ICU"), 
    n_res_hospital = sum(Resultat %in% "R"),
    n_tot_hospital = n(),
    .groups = "drop") %>%
  pivot_longer(matches("^n_.*_"), cols_vary = "slowest", names_to = c(".value", "level"), names_pattern = "(n_.*)_(.*)") %>%
  mutate(hospital = "all hospitals", bacterie = recode(bacterie, !!!bact_names))

# Compare prevalence
bind_rows(
  isolates_all, 
  res_hospital %>%
    mutate(Date_year = year(Date_week), level = "hospital", hospital = "cohort") %>%
    group_by(Date_year, bacterie, level, hospital) %>%
    summarise(n_res = sum(n_res), n_tot = sum(n_tot), .groups = "drop"),
  res_icu %>%
    mutate(Date_year = year(Date_week), level = "icu", hospital = "cohort") %>%
    group_by(Date_year, bacterie, level, hospital) %>%
    summarise(n_res = sum(n_res), n_tot = sum(n_tot), .groups = "drop")
  ) %>%
  nest(.by = c(Date_year, bacterie, level, hospital)) %>%
  mutate(cfint = map(data, function(.data) getBinomCI(.data, sides = "two.sided", method = "wilson"))) %>%
  unnest(cols = c(data, cfint)) %>%
  ggplot(., aes(x = Date_year, y = res_rate, ymin = res_rate_lwr, ymax = res_rate_upr, fill = hospital)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.5), width = 0.5) +
  geom_errorbar(position = position_dodge(width = 0.5), width = 0.3) +
  facet_grid(cols = vars(level), rows = vars(bacterie)) +
  scale_fill_manual(values = c("all hospitals" = "grey80", "cohort" = "cornflowerblue")) +
  theme_bw() +
  theme(axis.title.x = element_blank()) +
  labs(y = "Resistant isolate prevalence (95% Wilson CIs)")
ggsave("plots/bacterial_isolates/resistance_prevalence_sensitivity_cohort.png", height = 8, width = 6)

# # Variations in annual consumption (Tables in the main manuscript)
# hospital_tab = antibiotic_use %>%
#   group_by(Date_year, code, secteur, atb_class) %>%
#   summarise(molDDD = sum(molDDD), Nbhosp = unique(Nbhosp), .groups = "drop") %>%
#   group_by(Date_year, code, atb_class) %>%
#   summarise(molDDD = sum(molDDD), Nbhosp = sum(Nbhosp), .groups = "drop") %>%
#   group_by(atb_class) %>%
#   nest() %>%
#   mutate(mul_comp = map(data, atbmulcomp, level = "national", cohort = F)) %>%
#   select(-data) %>%
#   unnest(mul_comp) %>%
#   ungroup() %>%
#   mutate(
#     comparison = paste0(group1, "vs", group2),
#     conf_int = ifelse( -conf.low > -conf.high,
#                        paste0(-round(estimate, 1), " (", -round(conf.high, 1), ", ", -round(conf.low,1), ")"),
#                        paste0(-round(estimate, 1), " (", -round(conf.low, 1), ", ", -round(conf.high,1), ")")
#     ),
#     p.adj.bis = round(p.adj, 3),
#     p.adj = ifelse(p.adj < 0.001, "<0.001", round(p.adj, 3)),
#     # p_f = ifelse(p_f < 0.001, "<0.001", round(p_f, 3))
#   ) %>%
#   select(atb_class, comparison, conf_int, p.adj, p.adj.bis) %>%
#   pivot_wider(names_from = comparison, values_from = c(conf_int, p.adj, p.adj.bis)) %>%
#   arrange(atb_class) %>%
#   mutate(atb_class = as.character(atb_class)) %>%
#   gt(.) %>%
#   cols_hide(
#     columns = c(p.adj.bis_2019vs2020, p.adj.bis_2019vs2021, p.adj.bis_2019vs2022)
#   ) %>%
#   tab_spanner(
#     label = "2019 vs 2020",
#     columns = c(conf_int_2019vs2020, `p.adj_2019vs2020`)
#   ) %>%
#   tab_spanner(
#     label = "2019 vs 2021",
#     columns = c(conf_int_2019vs2021, `p.adj_2019vs2021`)
#   ) %>%
#   tab_spanner(
#     label = "2019 vs 2022",
#     columns = c(conf_int_2019vs2022, `p.adj_2019vs2022`)
#   ) %>%
#   cols_label(
#     # p_f = "Friedman test p-value",
#     conf_int_2019vs2020 = "Estimate (95% CI)",
#     `p.adj_2019vs2020` = "p-value",
#     conf_int_2019vs2021 = "Estimate (95% CI)",
#     `p.adj_2019vs2021` = "p-value",
#     conf_int_2019vs2022 = "Estimate (95% CI)",
#     `p.adj_2019vs2022` = "p-value",
#     atb_class = ""
#   ) %>%
#   tab_style(
#     style = list(
#       cell_text(weight = "bold")
#     ),
#     locations = list(
#       cells_body(
#         rows = p.adj.bis_2019vs2020 <= 0.05,
#         columns = conf_int_2019vs2020
#       ),
#       cells_body(
#         rows = p.adj.bis_2019vs2021 <= 0.05,
#         columns = conf_int_2019vs2021
#       ),
#       cells_body(
#         rows = p.adj.bis_2019vs2022 <= 0.05,
#         columns = conf_int_2019vs2022
#       ),
#       cells_column_spanners(),
#       cells_column_labels()
#     )
#   ) %>%
#   tab_options(table.font.size = 11) %>%
#   tab_footnote(
#     footnote = "Estimates and their 95% CI in bold have a p-value≤0.05",
#     locations = cells_column_labels(columns = c(conf_int_2019vs2020, conf_int_2019vs2021, conf_int_2019vs2022))
#   ) %>%
#   tab_footnote(
#     footnote = "Corresponds to trimethoprim and combinations of sulfanomides",
#     locations = cells_body(rows = atb_class == "Trimethoprim", columns = atb_class)
#   ) %>%
#   tab_footnote(
#     footnote = "P-value adjusted with the Benjamini-Hochberg procedure",
#     locations = cells_column_labels(columns = c(p.adj_2019vs2020, p.adj_2019vs2021, p.adj_2019vs2022))
#   ) #%>%
#   # sub_values(
#   #   columns = atb_class,
#   #   rows = DCI != "",
#   #   replacement = "",
#   #   pattern = "Imipenem \\+ Meropenem|Azithromycin|Vancomycin"
#   # )
# hospital_tab
# gtsave(hospital_tab, filename = "../Paper/Tables/national_hospital_bh_all-hospitals.docx")
# gtsave(hospital_tab, filename = "tables/Table2_bh_all-hospitals.docx")


# ################################################################################
# # Comparison of included and excluded hospitals
# ################################################################################
# # Finess of espic hospitals
# espic_list = read_excel("data-raw/spares/hospital_location/espic.xlsx", col_names = T)
# 
# # Get all finess numbers available in SPARES data
# all_finess = bind_rows(
#   read_excel("data-raw/spares/2019/SPARES_2019_BMRCovid.xlsx", sheet = "ADM", col_types = c("numeric", "text", "skip", "skip", "text", "text", "text")) %>% rename(code = IdEtablissement, type = groupe, Region = `Nouvelle-Region`, `finess géographique` = finess),
#   read_excel("data-raw/spares/2020/atb2020_2021.xlsx", sheet = "Participants ES_2020", col_types = c("numeric", "text", "text", "skip", "skip", "text", "text")), 
#   read_excel("data-raw/spares/2020/atb2020_2021.xlsx", sheet = "Participants ES_2021", col_types = c("numeric", "text", "text", "skip", "skip", "text", "text")),
#   read_excel("data-raw/spares/2022/ATB2022.xlsx", sheet = "Participants ES_2022", col_types = c("numeric", "text", "text", "skip", "skip", "text", "text"))
# ) %>% 
#   rename(finess_jur = finess_juridique, finess = `finess géographique`, region = Region) %>%
#   mutate(
#     region = case_when(
#       region == "Auvergne - Rhône Alpes" ~ "Auvergne-Rhône-Alpes",
#       region == "Bourgogne - Franche Comté" ~ "Bourgogne-Franche-Comté",
#       region == "Grand Est" ~ "Grand-Est",                 
#       region == "Hauts de France" ~ "Hauts-de-France",
#       region == "Ile de France" ~ "Île-de-France",
#       region == "Nouvelle Aquitaine" ~ "Nouvelle-Aquitaine",
#       region == "Pays de Loire" ~ "Pays de la Loire",
#       region == "Provence Alpes Côte d'Azur" ~ "Provence-Alpes-Côte d'Azur",
#       region == "Reunion - Mayotte" ~ "La Réunion-Mayotte",
#       .default = region
#     ),
#     type = case_when(code == 11014 ~ "CLCC",
#                      code == 2013 ~ "MCO",
#                      code == 1913 ~ "MCO", # physical Finess code: 210011847
#                      code == 11180 ~ "MCO", # physical Finess code: 420000192
#                      .default = type),
#     finess = ifelse(nchar(finess) < 9, paste0("0", finess), finess),
#     finess_jur = ifelse(nchar(finess_jur) < 9, paste0("0", finess_jur), finess_jur),
#   ) %>%
#   mutate(
#     type = recode(type, !!!dict_hospital_type),
#     finess = case_when(
#       code == 2406 ~ "850000647",
#       code == 12672 ~ "940110042",
#       code == 2453 ~ "440059319",
#       code == 2009 ~ "570026252",
#       code == 7115 ~ "630000404",
#       code == 9657 ~ "690781810", 
#       code == 9709 ~ "310025010",
#       .default = finess
#     ),
#     finess_jur = case_when(
#       code == 2478 ~ "440041895", 
#       code == 10597 ~ "920029527",
#       code == 8984 ~ "420784878",
#       code == 8836 ~ "350001137",
#       code == 9984 ~ "950042994",
#       code == 12672 ~ "940110042",
#       code == 2453 ~ "440059301",
#       code == 2009 ~ "570023630",
#       code == 9709 ~ "310025010",
#       .default = as.character(finess_jur)
#     ),
#     ) %>%
#   mutate(
#     finess = ifelse(type == "University hospital", finess_jur, finess),
#     type = ifelse(finess %in% espic_list$finess, "Private not-for-profit hospital", type)
#   ) %>%
#   distinct()
# 
# # Hospitals reporting antibiotics
# reporting_atb = bind_rows(
#   read_excel("data-raw/spares/2019/ATB2019_NAT.xlsx", sheet = "ATB2019") %>% distinct(code, secteur, Nbhosp, Nblits) %>% mutate(Date_year = 2019),
#   read_excel("data-raw/spares/2020/atb2020_2021.xlsx", sheet = "ATB2020") %>% distinct(code, secteur, Nbhosp, Nblits) %>% mutate(Date_year = 2020), 
#   read_excel("data-raw/spares/2020/atb2020_2021.xlsx", sheet = "ATB2021") %>% distinct(code, secteur, Nbhosp, Nblits) %>% mutate(Date_year = 2021),
#   read_excel("data-raw/spares/2022/ATB2022.xlsx", sheet = "ATB2022") %>% distinct(code,secteur, Nbhosp, Nblits) %>% mutate(Date_year = 2022)
# ) %>% 
#   distinct(code, secteur, Nbhosp, Nblits, Date_year) %>%
#   filter(secteur == "Total établissement")
# 
# # Hospitals reporting isolates
# reporting_isolates = bind_rows(
#   read.csv("data-raw/spares/2019/Souches_BMRCovid_2019_d1 copie.csv", sep = "|", encoding="latin1") %>% distinct(IdEtablissement),
#   read.csv("data-raw/spares/2020/BMR_Covid_souches_2020.csv", sep = "|", encoding = "latin1") %>% distinct(IdEtablissement),
#   read.csv("data-raw/spares/2021/BMR_Covid_souches_2021.csv", sep = "|", encoding = "latin1") %>% distinct(IdEtablissement),
#   read.csv("data-raw/spares/2022/BMR_Covid_souches_2022.csv", sep = "|", encoding = "latin1") %>% distinct(IdEtablissement)
# ) %>%
#   distinct(IdEtablissement) %>%
#   rename(code = IdEtablissement)
# 
# # Hospitals in SPARES for which we don't have finess numbers
# nrow(reporting_isolates)
# sum(reporting_isolates$code %in% all_finess$code)
# 
# length(unique(reporting_atb$code))
# sum(unique(reporting_atb$code) %in% all_finess$code)
# 
# all_codes = sort(unique(c(reporting_atb$code, reporting_isolates$code)))
# all_codes[!all_codes %in% all_finess$code] 
# 
# # Metadata
# metadata_all = bind_rows(
#   read_excel("data-raw/spares/2019/SPARES_2019_BMRCovid.xlsx", sheet = "ADM") %>%
#     rename(code = IdEtablissement, type = groupe, region = `Nouvelle-Region`) %>%
#     select(code, type, region) %>%
#     mutate(Date_year = 2019) %>%
#     left_join(., read_excel("data-raw/spares/2019/SPARES_2019_BMRCovid.xlsx", sheet = "JH_total") %>% rename(code = idetablissement, nbjh = nbJH), 
#               by = c("code")), 
#   read_excel("data-raw/spares/2020/BMR_Covid_2020.xlsx", sheet = "ADM") %>%
#     rename(code = idetablissement, type = groupe, region = Nouvelle_Region) %>%
#     select(code, type, region) %>%
#     mutate(Date_year = 2020) %>%
#     left_join(., read_excel("data-raw/spares/2020/BMR_Covid_2020.xlsx", sheet = "JH_total") %>% rename(code = idetablissement), 
#               by = c("code")),
#   read_excel("data-raw/spares/2021/BMR_Covid_2021.xlsx", sheet = "ADM") %>%
#     rename(code = idetablissement, type = groupe, region = Nouvelle_Region) %>%
#     select(code, type, region) %>%
#     mutate(Date_year = 2021) %>%
#     left_join(., read_excel("data-raw/spares/2021/BMR_Covid_2021.xlsx", sheet = "JH_total") %>% rename(code = idetablissement), 
#               by = c("code")),
#   read_excel("data-raw/spares/2022/BMR_Covid_2022.xlsx", sheet = "ADM") %>%
#     rename(code = idetablissement, type = groupe, region = Nouvelle_Region) %>%
#     select(code, type, region) %>%
#     mutate(Date_year = 2022) %>%
#     left_join(., read_excel("data-raw/spares/2022/BMR_Covid_2022.xlsx", sheet = "JH_total") %>% rename(code = idetablissement), 
#               by = c("code"))
# ) %>%
#   select(-c(type, region, secteur)) %>%
#   distinct() %>%
#   group_by(code) %>%
#   summarise(nbjh = round(mean(nbjh)), nb_lits = round(mean(nb_lits)), .groups = 'drop')
# 
# # Table with all metadata available
# k = all_finess %>%
#   filter(
#     !type %in% c("CLCC", "HIA", "PSY", "Long-term care facility"),
#     !region %in% c("Corse", "Guadeloupe", "Guyane", "La Réunion-Mayotte", "Martinique", "Nouvelle Calédonie")
#     ) %>%
#   left_join(., metadata_all, by = "code") %>%
#   mutate(
#     cohort = ifelse(code %in% cohort_final, "Selected", "Excluded")
#     ) 
# 
# chisq.test(k %>% filter(cohort == "Selected") %>% count(type) %>% pull(n),
#            p=k %>% count(type) %>% mutate(n = n/sum(n)) %>% pull(n))
# 
# chisq.test(k %>% filter(cohort == "Selected") %>% count(region) %>% pull(n),
#            p=k %>% count(region) %>% mutate(n = n/sum(n)) %>% pull(n))
# 
# 
# # Comparison in bed-days and number of hospitals per region
# # Plot distributions
# left_join(k %>% distinct(region, finess_jur) %>% count(region) %>% rename(Counts_SPARES = n), atih_nb, by = "region") %>%
#   mutate(Distribution_SPARES = Counts_SPARES / sum(Counts_SPARES), 
#          Distribution_ATIH = Counts_ATIH / sum(Counts_ATIH),
#          region = recode(region, !!!dict_regions)) %>%
#   pivot_longer(-region, names_to = "param", values_to = "value") %>%
#   separate(param, into = c("param", "db"), sep = "_") %>%
#   ggplot(., aes(x = region, y = value, fill = db)) +
#   geom_bar(stat = "identity", position = position_dodge(width = 0.5), width = 0.5) +
#   facet_grid(rows = vars(param), scales = "free_y") +
#   theme_bw() +
#   scale_fill_manual(values = c("SPARES" = "coral", "ATIH" = "darkorchid")) +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#   labs(x = "", y = "Legal entities", fill = "")
# 
# # Chi2 test of adequation
# chisq.test(k %>% distinct(region, finess_jur) %>% count(region) %>% pull(n), 
#            p = atih_nb$Counts_ATIH / sum(atih_nb$Counts_ATIH))
# 
# # Comparison in bed-days bed-days of hospitals per region
# k %>% filter(is.na(nbjh)) %>% count(region)
# k %>% filter(is.na(nbjh)) %>% count(type)
# 
# beds %>%
#   dplyr::select(region, Counts_ATIH, Distribution_ATIH) %>%
#   left_join(., 
#             k %>% group_by(region) %>% summarise(Counts_SPARES = sum(nbjh, na.rm = T), .groups = "drop") %>% mutate(Distribution_SPARES = Counts_SPARES / sum(Counts_SPARES)), 
#             by = "region") %>%
#   pivot_longer(-region, names_to = "param", values_to = "value") %>%
#   separate(param, into = c("param", "db"), sep = "_") %>%
#   mutate(region = recode(region, !!!dict_regions)) %>%
#   ggplot(., aes(x = region, y = value, fill = db)) +
#   geom_bar(stat = "identity", position = position_dodge(width = 0.5), width = 0.5) +
#   facet_grid(rows = vars(param), scales = "free_y") +
#   scale_fill_manual(values =  c("SPARES" = "coral", "ATIH" = "darkorchid")) +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#   labs(fill = "", x = "", y = "Bed-days in 2020")
# 
# # Chi2 test of adequation
# chisq.test(x=k %>% group_by(region) %>% summarise(nbjh = sum(nbjh, na.rm = T), .groups = "drop") %>% pull(nbjh), 
#            p=beds$Distribution_ATIH)
# 
