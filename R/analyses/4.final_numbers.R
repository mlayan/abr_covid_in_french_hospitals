##################################################
## CONSTITUTE FINAL COHORT AND FINAL PLOTS
## OF THE COHORT
##################################################
rm(list = ls())
library(tidyverse)
library(readxl)
library(sf)
library(RColorBrewer)
library(ggpubr)
library(ggnewscale)
library(geofacet)

load("data/metadata_admin.rda")
load("data/metadata_beds.rda")
load("data/cohort19202122.rda")
load("data/icu_cohort19202122.rda")

source("R/helper/helper_functions.R")
source("R/helper/dictionaries.R")

##################################################
# Data
##################################################
# Load data
res = bind_rows( 
  read.table("data-raw/spares/combined/resistance_cohort19202122_2019.txt", 
             sep = "\t", header = T),
  read.table("data-raw/spares/combined/resistance_cohort19202122_2020.txt",
             sep = "\t", header = T),
  read.table("data-raw/spares/combined/resistance_cohort19202122_2021.txt",
             sep = "\t", header = T),
  read.table("data-raw/spares/combined/resistance_cohort19202122_2022.txt",
             sep = "\t", header = T)
)

# Number of isolates that are in the final database
res %>%
  dplyr::select(-c(molecule, atb_class, Resultat)) %>% 
  distinct() %>%
  count(Date_year)

# Antibiotic data 
atb_sub = bind_rows( 
  read_excel("data-raw/spares/2019/ATB2019_NAT.xlsx", sheet = "ATB2019") %>% mutate(Date_year = 2019),
  read_excel("data-raw/spares/2020/atb2020_2021.xlsx", sheet = "ATB2020") %>% mutate(Date_year = 2020),
  read_excel("data-raw/spares/2020/atb2020_2021.xlsx", sheet = "ATB2021") %>% mutate(Date_year = 2021),
  read_excel("data-raw/spares/2022/ATB2022.xlsx", sheet = "ATB2022") %>% mutate(Date_year = 2022)
) %>%
  dplyr::select(-c(Type, type, region, Region, molG)) %>%
  filter(code %in% cohort19202122, !secteur %in% c("Hématologie", "Maladie inf", "SLD", "Psychiatrie", "Pédiatrie", "Total établissement")) %>%
  mutate(
    atb_class = recode(ATC, !!!dict_antibiotic_class),
    secteur = recode(secteur, !!!dict_secteur_spares)
    )

atb = bind_rows(
  atb_sub,
  atb_sub %>%
    group_by(code, DCI, ATC, DDJ, Date_year, atb_class) %>%
    summarise(Nbhosp = sum(Nbhosp), Nblits = sum(Nblits), molDDD = sum(molDDD), .groups = "drop") %>%
    mutate(secteur = "Total")
)

# PMSI hospitalisation days
jh_pmsi = read.table("data-raw/atih/pmsi_sejours.txt", header = T, sep = "\t") %>%
  mutate(secteur = recode(secteur, !!!dict_secteur_spares))

# PMSI COVID-19 hospitalisation days
covid_pmsi = read.table("data-raw/atih/pmsi_sejours_covid.txt", header = T, sep = "\t") %>%
  mutate(secteur = recode(secteur, !!!dict_secteur_spares))

##################################################
# Get final cohort - hospital level
##################################################
# No infection reported by species of interest
metadata_admin %>%
  filter(code %in% cohort19202122, !code %in% res$code) %>%
  count(type)

to_exclude = metadata_admin %>%
  filter(!code %in% cohort19202122 | !code %in% res$code) %>%
  dplyr::select(code) %>%
  distinct() %>%
  .$code
length(to_exclude)

# Hospitals not available in PMSI
to_exclude_pmsi = cohort19202122[!cohort19202122 %in% jh_pmsi$code]
length(unique(to_exclude_pmsi))
metadata_admin %>%
  filter(code %in% to_exclude_pmsi) %>%
  dplyr::select(code, type) %>%
  distinct() %>%
  count(type)

# Number of hospitals and ICUs
length(unique(cohort19202122))

# Get and save final cohort
cohort_final = cohort19202122[!cohort19202122 %in% c(to_exclude, to_exclude_pmsi)]
length(unique(cohort_final))
save(cohort_final, file = "data/cohort_final.rda")

length(cohort_final)/3008*100 # https://drees.solidarites-sante.gouv.fr/sites/default/files/2021-07/Vue%20d%27ensemble.pdf

# Number hospitals by type
metadata_admin %>%
  filter(code %in% cohort_final) %>%
  dplyr::select(code, type) %>%
  distinct() %>%
  count(type)

##################################################
# Get final cohort - ICU level
##################################################
# Number of ICUs in 2019-2022 cohort 
res_atb = bind_rows(
  read.csv("data-raw/spares/2019/Souches_BMRCovid_2019_d1 copie.csv", sep = "|", encoding="latin1") %>% dplyr::select(IdEtablissement, secteur) %>% rename(code = IdEtablissement) %>% mutate(db = "res", Date_year = 2019),
  read.csv("data-raw/spares/2020/BMR_Covid_souches_2020.csv", sep = "|", encoding = "latin1") %>% dplyr::select(IdEtablissement, secteur) %>% rename(code = IdEtablissement) %>% mutate(db = "res", Date_year = 2020),
  read.csv("data-raw/spares/2021/BMR_Covid_souches_2021.csv", sep = "|", encoding = "latin1") %>% dplyr::select(IdEtablissement, secteur) %>% rename(code = IdEtablissement) %>% mutate(db = "res", Date_year = 2021),
  read.csv("data-raw/spares/2022/BMR_Covid_souches_2022.csv", sep = "|", encoding = "latin1") %>% dplyr::select(IdEtablissement, secteur) %>% rename(code = IdEtablissement) %>% mutate(db = "res", Date_year = 2022)
) %>%
  bind_rows(.,
            read_excel("data-raw/spares/2019/ATB2019_NAT.xlsx", sheet = "ATB2019") %>% dplyr::select(code, secteur) %>% mutate(db = "atb", Date_year = 2019),
            read_excel("data-raw/spares/2020/atb2020_2021.xlsx", sheet = "ATB2020") %>% dplyr::select(code, secteur) %>% mutate(db = "atb", Date_year = 2020),
            read_excel("data-raw/spares/2020/atb2020_2021.xlsx", sheet = "ATB2021") %>% dplyr::select(code, secteur) %>% mutate(db = "atb", Date_year = 2021),
            read_excel("data-raw/spares/2022/ATB2022.xlsx", sheet = "ATB2022") %>% dplyr::select(code, secteur) %>% mutate(db = "atb", Date_year = 2022)
            ) %>%
  filter(secteur == "Réanimation") %>%
  distinct()

all_icus = res_atb %>% 
  filter(code %in% cohort19202122) %>%
  dplyr::select(code) %>%
  distinct() %>%
  nrow()

# Number that did not report their antibiotic consumption or resistance data
# for a given year 
all_icus - length(icu_cohort19202122)

# Number of ICUs from selected hospitals 
icu_cohort_final = cohort_final[cohort_final %in% icu_cohort19202122]
length(icu_cohort_final)

 # Verify that they reported infections in SPARES 
to_exclude_icu_res = metadata_admin %>%
  filter(code %in% icu_cohort_final) %>%
  filter(!code %in% res$code[res$secteur == "ICU"]) %>%
  dplyr::select(code) %>%
  distinct() %>%
  .$code
length(to_exclude_icu_res)

# Verify that they reported hospitalisation days in PMSI
to_exclude_icu_pmsi = metadata_admin %>%
  filter(code %in% icu_cohort_final) %>%
  filter(!code %in% jh_pmsi$code[jh_pmsi$secteur == "ICU"]) %>%
  dplyr::select(code) %>%
  distinct() %>%
  .$code
length(to_exclude_icu_pmsi)

# Save ICU cohort
icu_cohort_final = icu_cohort_final[!icu_cohort_final %in% c(to_exclude_icu_pmsi, to_exclude_icu_res)]
length(unique(icu_cohort_final))
save(icu_cohort_final, file = "data/icu_cohort_final.rda")

# Types of hospitals where ICUs are followed-up
metadata_admin %>%
  filter(code %in% icu_cohort_final) %>%
  dplyr::select(code, type) %>%
  distinct() %>%
  count(type)

##################################################
# Get final number of bacterial samples included
##################################################
# Total number of bacterial isolates - hospital level
res %>%
  filter(code %in% cohort_final) %>%
  dplyr::select(-c(Resultat, atb_class, molecule)) %>%
  distinct() %>%
  nrow(.)

# Total number of resistance bacterial isolates - hospital level
res %>%
  filter(code %in% cohort_final, Resultat %in% c("O", "R")) %>%
  dplyr::select(-c(Resultat, atb_class, molecule)) %>%
  distinct() %>%
  nrow(.)

# Total number of isolates - ICU level
res %>%
  filter(code %in% icu_cohort_final, secteur == "ICU") %>%
  dplyr::select(-c(Resultat, atb_class, molecule)) %>%
  distinct() %>%
  nrow(.)

# Total number of resistant isolates - ICU level
res %>%
  filter(code %in% icu_cohort_final, secteur == "ICU", Resultat %in% c("O", "R")) %>%
  dplyr::select(-c(Resultat, atb_class, molecule)) %>%
  distinct() %>%
  nrow(.)

##################################################
# Get information on ESPIC hospitals
##################################################
# ESPIC included in the final cohort
espic_list = read_excel("data-raw/spares/hospital_location/espic.xlsx", col_names = T)

# ESPIC final hospital cohort
metadata_admin %>%
  filter(code %in% cohort_final, finess %in% espic_list$finess | finess_juridique %in% espic_list$finess) %>%
  count(type)

# ESPIC final ICU cohort
metadata_admin %>%
  filter(code %in% icu_cohort_final, finess %in% espic_list$finess | finess_juridique %in% espic_list$finess) %>%
  count(type)

# ESPIC names in final hospital cohorts
metadata_admin %>%
  filter(
    code %in% cohort_final, 
    finess %in% espic_list$finess | finess_juridique %in% espic_list$finess,
    type == "Private for profit hospital"
  ) %>%
  .$name

metadata_admin %>%
  filter(
    code %in% cohort_final, 
    finess %in% espic_list$finess | finess_juridique %in% espic_list$finess,
    type == "Rehabilitation hospital"
  ) %>%
  .$name

# Codes corresponding to ESPIC hospitals
espic = metadata_admin %>%
  filter(
    code %in% cohort_final, 
    finess %in% espic_list$finess | finess_juridique %in% espic_list$finess
  ) %>%
  dplyr::select(code) %>%
  distinct() %>%
  .$code

# Save metadata admin with ESPIC info
metadata_admin_espic = metadata_admin %>%
  filter(code %in% cohort_final) %>%
  mutate(type = ifelse(code %in% espic, "Private not-for-profit hospital", type),
         icu = ifelse(code %in% icu_cohort_final, 1, 0))
save(metadata_admin_espic, file = "data/metadata_admin_espic.rda")

# Distribution by hospital type
metadata_admin_espic %>% 
  group_by(code) %>% 
  summarise(type = unique(type), .groups = "drop") %>% 
  count(type)

metadata_admin_espic %>% 
  filter(icu == 1) %>% 
  group_by(code) %>% 
  summarise(type = unique(type), .groups = "drop") %>% 
  count(type)

##################################################
# Map of France with abbreviated region names
##################################################
france_map = france %>%
  mutate(region_abbrev = recode(region, !!!dict_regions)) %>%
  ggplot(.) +
  geom_sf(fill = "white") +
  geom_sf_text(aes(label = region_abbrev), col = "black", size = 4) +
  theme_bw() +
  theme(axis.title = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm"))
ggsave("plots/cohort_description/map_regions_abbreviated.png", france_map,
       height = 5, width = 5)

##################################################
# Plot hospitals on map
##################################################
# Get GPS coordinates
data_gouv = readLines("data-raw/spares/hospital_location/etalab-cs1100507-stock-20230502-0337.csv")
gps_coor = data_gouv[sapply(data_gouv, grepl, pattern = "geolocalisation")]
gps_coor = data.frame(X = gps_coor) %>%
  tidyr::separate(X, c("section", "finess", "coordX", "coordY", "sourcecoordet", "datemaj"), ";") %>%
  dplyr::select(-c(section, sourcecoordet, datemaj))
rm(data_gouv)

# Map of final cohort
finess_hosp = read_excel("data-raw/spares/2020/BMR_Covid_2020.xlsx", sheet = "ADM") %>% 
  dplyr::select(idetablissement, finess, finess_juridique, etablissement, ville, groupe) %>%
  rename(code = idetablissement, type = groupe) %>%
  filter(code %in% cohort_final) %>%
  mutate(finess = ifelse(nchar(finess) < 9, paste0("0", finess), finess)) %>%
  left_join(., gps_coor, by = "finess") %>%
  mutate(
    coordX = case_when(
      code == 2079 ~ "774591.0", # 510000029 - 510002447 - 774591.0 6903471.4
      code == 2370 ~ "432382.3", # 490000031 - 490000049 - 432382.3 6714951.6
      code == 7115 ~ "707341.7", # 630780989 - 630000404 - 707341.7 6517496.7
      code == 8984 ~ "808974.5", # 420784878 - 420782559 - 808974.5 6480266.0
      code == 9243 ~ "350633.9", # 350005179 - 350000741 - 350633.9 6789911.9
      code == 10357 ~ "562764.9", # 870000015 - 870000064 - 562764.9 6525390.1
      code == 10726 ~ "806105.2", # 300780038 - 300782117 - 806105.2 6303500.4
      code == 12584 ~ "619539.2", # 450000088 - 450002613 - 619539.2 6749160.4
      .default = coordX
    ),
    coordY = case_when(
      code == 2079 ~ "6903471.4", # 510000029 - 510002447 - 774591.0 6903471.4
      code == 2370 ~ "6714951.6", # 490000031 - 490000049 - 432382.3 6714951.6
      code == 7115 ~ "6517496.7", # 630780989 - 630000404 - 707341.7 6517496.7
      code == 8984 ~ "6480266.0", # 420784878 - 420782559 - 808974.5 6480266.0
      code == 9243 ~ "6789911.9", # 350005179 - 350000741 - 350633.9 6789911.9
      code == 10357 ~ "6525390.1", # 870000015 - 870000064 - 562764.9 6525390.1
      code == 10726 ~ "6303500.4", # 300780038 - 300782117 - 806105.2 6303500.4
      code == 12584 ~ "6749160.4", # 450000088 - 450002613 - 619539.2 6749160.4
      .default = coordY
    ),
    type = case_when(
      type == "LOC" & !code %in% espic ~ "CH", 
      code %in% espic ~ "ESPIC",
      .default = type
      )
  ) %>%
  mutate(type = recode(type, !!!dict_hospital_type)) %>%
  st_as_sf(., coords = c("coordX", "coordY"), crs = 2154)

map_hosp = ggplot() +
  geom_sf(data = france, fill = "grey98") +
  new_scale_color() +
  geom_sf(data = finess_hosp, aes(col = type), size = 1.2) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  guides(colour = guide_legend(override.aes = list(size=3))) +
  labs(col = "", title = "Hospitals") +
  scale_color_manual(values = c("General public hospital" = "#8ecae6", 
                                "Private for profit hospital" = "#219ebc",
                                "Private not-for-profit hospital" = "#023047", 
                                "Rehabilitation hospital" = "#ffb703", 
                                "University hospital" = "#fb8500"))

map_icu = ggplot() +
  geom_sf(data = france, fill = "grey98") +
  new_scale_color() +
  geom_sf(data = finess_hosp %>% filter(code %in% icu_cohort_final), aes(col = type), size = 1.2) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  guides(colour = guide_legend(override.aes = list(size=3))) +
  labs(col = "", title = "ICUs") +
  scale_color_manual(values = c("General public hospital" = "#8ecae6", 
                                "Private for profit hospital" = "#219ebc",
                                "Private not-for-profit hospital" = "#023047", 
                                "Rehabilitation hospital" = "#ffb703", 
                                "University hospital" = "#fb8500"))

maps = ggarrange(map_hosp, map_icu, ncol = 2, common.legend = T, legend = "right") 
maps
ggsave("plots/cohort_description/map_cohortfinal.png", maps, height = 4, width = 10)
ggsave("plots/final_figures/Figure1D.pdf", maps, height = 4, width = 10)

##################################################
# Plot numbers of resistant isolates
# by bacterial species 
##################################################
# Dataframe of numbers
df_numbers = res %>%
  filter(code %in% cohort_final) %>%
  mutate(molecule = ifelse(molecule %in% c("Imipénème", "Méropénème"), "Carbapenems", molecule)) %>%
  group_by(across(-c(Resultat, atb_class))) %>%
  summarise(resistance = ifelse(any(Resultat %in% c("O", "R")), 1, 0), .groups = "drop") %>%
  group_by(Date_year, molecule, bacterie) %>%
  summarise(
    n_res = sum(resistance %in% 1),
    n_tot = n(),
    .groups = "drop"
  ) %>%
  mutate(
    bacterie = case_when(
      bacterie == "Staphylococcus aureus" ~ "MRSA",
      bacterie == "Escherichia coli" ~ "ESBL E. coli",
      bacterie == "Klebsiella pneumoniae" ~ "ESBL K. pneumoniae",
      bacterie == "Enterobacter cloacae complex" ~ "ESBL E. cloacae",
      bacterie == "Pseudomonas aeruginosa" ~ "CR P. aeruginosa"
    )
  )
sum(df_numbers$n_res)
sum(df_numbers$n_tot)

# Plot of numbers
p1 = ggplot(df_numbers, aes(x = Date_year, y = n_res, fill = bacterie)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.6), width = 0.4) +
  theme_bw() +
  scale_fill_manual(values = brewer.pal(8, "Set2")) +
  labs(x = "Year", y = "Number of resistant isolates", fill = "")

# Add annual bed days to get incidence
df_nbjh = metadata_beds %>%
  filter(secteur != "Total") %>%
  rename(Date_year = year) %>%
  group_by(Date_year) %>%
  summarise(nbjh = sum(nbjh), .groups = "drop")

# Plot incidence
p2 = left_join(df_numbers, df_nbjh, by = "Date_year") %>%
  mutate(incidence = n_res / nbjh * 1000) %>%
  ggplot(., aes(x = Date_year, y = incidence, fill = bacterie)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.6), width = 0.4) +
  theme_bw() +
  scale_fill_manual(values = brewer.pal(8, "Set2")) +
  labs(x = "Year", y = "Incidence of resistant infections\nfor 1,000 bed-days", fill = "")

# Merge the two plots
p = ggarrange(p1, p2, nrow = 2, common.legend = T, legend = "right", align = "hv")
p 
ggsave("plots/cohort_description/selected_samples_timeseries.png", p, 
       height = 5, width = 8)

##################################################
# Hospital type distribution by regions in hospital 
# cohort
##################################################
# List of University hospitals 
all_chu = read_excel("data-raw/communes/commune_chu.xlsx") %>%
  filter(!main_city_name %in% c("Basse-Terre", "Fort-de-France", "Saint-Denis")) %>%
  mutate(main_city_name = ifelse(main_city_name %in% c("Paris", "Marseille", "Lyon"), paste0(main_city_name, " 01"), main_city_name))

# GPS coordinates of French cities
gps_cities = read.csv("data-raw/communes/communes-departement-region.csv") %>%
  dplyr::select(latitude, longitude, nom_commune) %>%
  inner_join(., all_chu, by = c("nom_commune" = "main_city_name")) %>%
  distinct() %>%
  arrange(nom_commune) %>%
  mutate(nom_commune = gsub(" 01", "", nom_commune))

# CHUs in SPARES
chu_spares = bind_rows(
  read_excel("data-raw/spares/2019/SPARES_2019_BMRCovid.xlsx", sheet = "ADM") %>%
    rename(code = IdEtablissement, type = groupe, name = etablissement, city = ville, region = `Nouvelle-Region`), 
  read_excel("data-raw/spares/2020/BMR_Covid_2020.xlsx", sheet = "ADM") %>%
    rename(code = idetablissement, name = etablissement, city = ville, type = groupe, region = `Nouvelle_Region`),
  read_excel("data-raw/spares/2021/BMR_Covid_2021.xlsx", sheet = "ADM") %>%
    rename(code = idetablissement, name = etablissement, city = ville, type = groupe, region = `Nouvelle_Region`),
  read_excel("data-raw/spares/2022/BMR_Covid_2022.xlsx", sheet = "ADM") %>%
    rename(code = idetablissement, name = etablissement, city = ville, type = groupe, region = `Nouvelle_Region`)
) %>%
  filter(type == "CHU", code %in% cohort_final) %>%
  mutate(finess_juridique = ifelse(finess == "420784878", "420784878", finess_juridique),
         finess = case_when(finess_juridique == "630780989" ~ "630000404", 
                            finess_juridique == "690781810" ~ "690781810",
                            .default = finess),
         name = case_when(code == 7115 ~ "CHU G. MONTPIED", 
                          code == 9263 ~ "CHU GRENOBLE-HOPITAL NORD",
                          .default = name),
         reporting = ifelse(finess != finess_juridique, "Included in cohort - Physical reporting", "Included in cohort - Legal reporting"),
         city = case_when(
           city == "CLERMONT FERRAND" ~ "Clermont-Ferrand",
           city == "ST ETIENNE" ~ "Saint-Étienne",
           city == "NIMES" ~ "Nîmes",
           city == "ORLEANS" ~ "Orléans",
           .default = stringr::str_to_title(city)
         )) %>%
  distinct()

chu_inclusion = chu_spares %>%
  group_by(finess_juridique, city, reporting) %>%
  summarise(n = n(), .groups = "drop") %>%
  full_join(., gps_cities, by = c("city" = "nom_commune", "finess_juridique" = "finess_jur")) %>%
  mutate(n = ifelse(is.na(n), 0, n), 
         reporting = factor(
           ifelse(is.na(reporting), "Not included", reporting), 
           c("Included in cohort - Legal reporting", 'Included in cohort - Physical reporting', "Not included"))
         ) %>%
  st_as_sf(., coords = c("longitude", "latitude"), crs = 4326) 

chu_inclusion_n = chu_spares %>%
  filter(reporting == "Included in cohort - Physical reporting") %>%
  group_by(city, finess_juridique, reporting) %>%
  summarise(n = n(), .groups = "drop") %>%
  left_join(., gps_cities, by = c("city" = "nom_commune", "finess_juridique" = "finess_jur")) %>%
  mutate(latitude = latitude - 0.40) %>%
  st_as_sf(., coords = c("longitude", "latitude"), crs = 4326) 

# Distribution of university hospitals 
p1 = ggplot(data = chu_inclusion) +
  geom_sf(data = france, fill = "#FFFFFF00") +
  geom_sf(aes(col = reporting, shape = reporting), size = 2) +
  scale_color_manual(name = "University hospitals", 
                     labels = c("Included in cohort - Legal reporting",
                                "Included in cohort - Physical reporting",
                                "Not included"),
                     values = c("orchid1", "orchid1", "grey50")) +
  scale_shape_manual(name = "University hospitals", 
                     labels = c("Included in cohort - Legal reporting",
                                "Included in cohort - Physical reporting",
                                "Not included"),
                     values = c(19, 17, 15)) +
  geom_sf_label(data = chu_inclusion_n, aes(label = n), size = 3) +
  theme_minimal() +
  theme(legend.box = 'vertical', plot.title = element_text(hjust = 0.5)) +
  labs(x = "Longitude", y = "Latitude")

# Distribution by region and type
p2=metadata_admin_espic %>% 
  filter(code %in% cohort_final) %>%
  mutate(type = ifelse(type == "ESPIC", "Private not-for-profit hospital", type),
         code = ifelse(type == "University hospital" & city == "MARSEILLE", 11048, code),
         region = recode(region, !!!dict_regions)) %>%
  dplyr::select(code, type, region) %>%
  distinct() %>%
  group_by(region, type) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(region) %>%
  mutate(p = n/sum(n)) %>%
  ungroup() %>%
  ggplot(., aes(y = p, fill = type, x = type)) +
  geom_bar(stat = "identity") +
  geom_text(aes(y = p+0.05, label=n)) +
  facet_wrap(facets = vars(region)) +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  scale_fill_manual(values = brewer.pal(n = 6, name = "Dark2")) +
  labs(x = "", y = "Proportion (in %)", fill = "")

# Final figure 
p = ggarrange(p1, p2, labels = c("A", "B"), nrow = 2, heights = c(0.8, 1))
p
ggsave("plots/cohort_description/cohort_description.png", p, height = 8, width = 7)
ggsave("../Paper/Supplementary/cohort_description.png", p, height = 8, width = 7)
