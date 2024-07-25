##################################################
## EVALUATION OF THE REPRESENTATIVITY OF THE 
## COHORT
##################################################
library(tidyverse) 
library(readxl)
library(haven)
library(RColorBrewer)
library(ggpubr)
rm(list = ls())

load("data/cohort_final.rda")
load("data/metadata_beds.rda")
load("data/metadata_admin_espic.rda")
source("R/helper/dictionaries.R")

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
