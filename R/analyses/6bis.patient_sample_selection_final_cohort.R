rm(list = ls())
library(tidyverse)
library(foreach)
library(doParallel)

load("data/cohort_final.rda")
load("data/dict_id_site.rda")
load("data/dict_id_age.rda")
load("data/dict_molecule_class.rda")
load("data/enzymes.rda")
source("R/helper_functions.R")

##################################################
# Load resistance data 
##################################################
registerDoParallel(cores = 3)

foreach(year = 2019:2022, .combine = "list", .packages = c("tidyverse")) %dopar% {
  
  # Basic checks----------------------------------------------------------------
  if (year == 2019) res = read.table("data-raw/spares/2019/Souches_BMRCovid_2019_d1 copie.csv", 
                                     sep = "|", header = T, encoding = "latin1") %>%
    rename(code= IdEtablissement)
  if (year == 2020) res = read.table("data-raw/spares/2020/BMR_Covid_souches_2020.csv", 
                                     sep = "|", header = T, encoding = "latin1") %>%
    rename(code= IdEtablissement)
  if (year == 2021) res = read.table("data-raw/spares/2021/BMR_Covid_souches_2021.csv", 
                                     sep = "|", header = T, encoding = "latin1") %>%
    rename(code= IdEtablissement)
  if (year == 2022) res = read.table("data-raw/spares/2022/BMR_Covid_souches_2022.csv", 
                                     sep = "|", header = T, encoding = "latin1") %>%
      rename(code= IdEtablissement)
  
  # Filtering based on selected hospitals/departments, patient age, site, and 
  # correct coding of results variable------------------------------------------
  res = res %>%
    # Filter out excluded facilities
    filter(code %in% cohort_final) %>%
    # Filter out SLD because it is not in PMSI
    filter(secteur != c("SLD")) %>%
    # Filter out samples from patients above 5 but isolated in newborns
    filter(!(idtrancheage > 1 & IdSite == 13)) %>%
    # Filter out test results/phenotypes that are incorrectly coded 
    filter(
      !(Resultat %in% c("O", "N") & !molecule %in% enzymes) | 
        (Resultat %in% c("S", "I", "R") & molecule %in% enzymes)
    ) %>%
    # Filter out samples that were isolated before the 14th of April and after the 31st of July
    mutate(prelevement = as.Date(prelevement, "%d/%m/%Y")) %>%
    filter(prelevement >= as.Date(paste0(year, "/04/14")), prelevement <= as.Date(paste0(year, "/07/31"))) %>%
    # Remove unnecessary columns
    # that were previously checked for potential coding issues
    dplyr::select(-c(IdArchiveLaboratoire, IdPeriode, IdBacterie, IdMolecule, nosocomial, suppr_dedoublonnage2, 
              iduf, idactivite, code_de, code_ta, IdetablissementLaboratoire)) %>%
    mutate(
      age_cat = recode(idtrancheage, !!!dict_id_age),
      molecule_class = recode(molecule, !!!dict_molecule_class),
      site = recode(IdSite, !!!dict_id_site),
      year = year,
      bacterie = case_when(
        bacterie == "Enterobacter cloacae co" ~ "Enterobacter cloacae complex",
        bacterie == "Enterobacter cloacae" ~ "Enterobacter cloacae complex",
        bacterie == "Enterobacter asburiae" ~ "Enterobacter cloacae complex",
        bacterie == "Enterobacter hormaechei" ~ "Enterobacter cloacae complex",
        bacterie == "Enterobacter kobei" ~ "Enterobacter cloacae complex",
        bacterie == "Enterobacter ludwigii" ~ "Enterobacter cloacae complex",
        bacterie == "Enterobacter nimipressuralis" ~ "Enterobacter cloacae complex",
        bacterie == "Streptococcus pneumonia" ~ "Streptococcus pneumoniae",
        .default = bacterie)
    ) %>%
    dplyr::select(code, site, prelevement, ligne, Resultat, Num_Patient, 
           age_cat, bacterie, molecule, molecule_class, secteur, year) %>%
    distinct()
  
  # Get numbers-----------------------------------------------------------------
  # 1. Samples that are isolated in the same individual, the same day, at the same 
  # site for the same bacteria but with different sample IDs--------------------
  ids_duplicates = res %>% 
    dplyr::select(code, site, prelevement, Num_Patient, age_cat, bacterie, secteur, ligne) %>%
    distinct() %>%
    group_by(code, site, prelevement, Num_Patient, age_cat, bacterie, secteur) %>%
    mutate(n = n()) %>%
    filter(n>1) %>%
    ungroup() %>%
    dplyr::select(-n)
  
  exact_duplicates = res %>%
    inner_join(., ids_duplicates, by = c("code", "site", "prelevement", "ligne", 
                                         "Num_Patient", "age_cat", "bacterie", "secteur")) %>%
    dplyr::select(-c(molecule_class, year)) %>%
    group_by(code, site, prelevement, Num_Patient, age_cat, bacterie, secteur) %>%
    mutate(n_samples = length(unique(ligne))) %>%
    ungroup() %>%
    group_by(code, site, prelevement, Num_Patient, age_cat, bacterie, secteur, molecule) %>%
    mutate(n = n(), n_diff = length(unique(Resultat))) %>%
    ungroup() %>%
    group_by(code, site, prelevement, Num_Patient, age_cat, bacterie, secteur) %>%
    mutate(uniqueness = case_when(
      all(n > 1) & length(unique(n)) == 1 & all(n_diff == 1) ~ 1, # Exact duplicates
      any(n < n_samples) & all(n_diff == 1) ~ 2, # More antibiotics tested, otherwise identical phenotypes
      .default = 3 # At least one phenotype that is different
    )) %>%
    ungroup()
  
  duplicates_to_remove = exact_duplicates %>%
    group_by(code, site, prelevement, Num_Patient, age_cat, bacterie, secteur) %>%
    mutate(m = min(ligne)) %>%
    filter(ligne != m) %>%
    ungroup() %>%
    dplyr::select(code, site, prelevement, Num_Patient, age_cat, bacterie, secteur, ligne) %>%
    distinct() 
  
  rm(ids_duplicates)
  rm(exact_duplicates)
  
  # 2. Samples that are part of a 30-day sequence-------------------------------
  samples_30days = res %>%
    anti_join(., duplicates_to_remove, 
              by = c("code", "site", "prelevement", "Num_Patient", "age_cat", "bacterie", "secteur", "ligne")) %>%
    dplyr::select(code, site, Num_Patient, age_cat, bacterie, secteur, ligne, prelevement) %>%
    distinct() %>%
    mutate(prelevement = as.Date(prelevement, "%d/%m/%Y")) %>%
    arrange(code, site, Num_Patient, age_cat, bacterie, secteur, prelevement) %>%
    group_by(code, site, Num_Patient, age_cat, bacterie, secteur) %>%
    mutate(n = n(), prelevement_lag = lag(prelevement)) %>%
    ungroup() %>%
    filter(n > 1) %>%
    mutate(time_diff = ifelse(
      is.na(prelevement_lag), 
      0, 
      difftime(as.Date(prelevement, "%d/%m/%Y"), as.Date(prelevement_lag, "%d/%m/%Y"), units = "day")
    )) %>%
    filter(!is.na(prelevement), time_diff <= 30) %>%
    group_by(code, site, Num_Patient, age_cat, bacterie, secteur) %>%
    mutate(n2 = n()) %>%
    ungroup() %>%
    filter(n2 > 1)
  
  redundant_profiles = samples_30days %>%
    dplyr::select(-c(n, prelevement_lag, time_diff, n2)) %>%
    left_join(., res %>% mutate(prelevement = as.Date(prelevement, "%d/%m/%Y")), 
              by = c("code", "site", "Num_Patient", "age_cat", "bacterie", "secteur", "ligne", "prelevement")) %>%
    group_by(code, site, Num_Patient, age_cat, bacterie, secteur) %>%
    nest() %>%
    mutate(selected = map(data, selection_30days_sequence)) %>%
    dplyr::select(-data) %>%
    unnest(cols = selected) %>%
    ungroup()

  rm(samples_30days)
  
  # Samples and phenotypes selection--------------------------------------------
  # 1. Keep samples with lowest "ligne" ID among duplicates---------------------
  res = res %>%
    anti_join(., duplicates_to_remove, 
              by = c("code", "site", "prelevement", "Num_Patient", "age_cat", "bacterie", "secteur", "ligne"))  

  # 2. Exclude samples with sensitive phenotype when a patient
  # is tested multiple times for the same bacteria at the same site
  # Remove phenotype that is sensitive or non producing enzymes 
  # when multiple tests for the same patient-ligne-prelevement-site
  to_exclude = res %>%
    group_by(code, site, prelevement, ligne, Num_Patient, age_cat, bacterie, secteur, year, molecule) %>%
    summarise(n = n(), .groups = "drop") %>%
    filter(n > 1) %>%
    dplyr::select(-n)
  
  res = anti_join(res, to_exclude, 
                  by = c("code", "site", "prelevement", "ligne", "Num_Patient", "age_cat", "bacterie", 
                         "secteur", "year", "molecule"))    

  # 3. Remove samples isolated in Chirurgie when already isolated in Médecine
  # on day d, at site s, in patient p, for bacteria b
  sector_duplicates = res %>% 
    dplyr::select(code, site, prelevement, Num_Patient, age_cat, bacterie, year, secteur) %>%
    distinct() %>%
    group_by(code, site, prelevement, Num_Patient, age_cat, bacterie, year) %>%
    summarise(n = n(), 
              n_chirurgie = sum(secteur %in% "Chirurgie"), 
              n_medecine = sum(secteur %in% "Médecine"), 
              n_reanimation = sum(secteur %in% "Réanimation"),
              n_gynecologie = sum(secteur %in% "Gynécologie-Obstétrique"),
              n_ssr = sum(secteur %in% "SSR"),
              .groups = "drop") %>%
    filter(n > 1)
  
  if (any(sector_duplicates$n_chirurgie > 1) | any(sector_duplicates$n_medecine > 1) | 
      any(sector_duplicates$n_reanimation > 1) | any(sector_duplicates$n_gynecologie > 1) |
      any(sector_duplicates$n_ssr > 1)) 
    stop("Problem with the number of duplicates")
  

  res = anti_join(
    res,
    sector_duplicates %>% 
      mutate(secteur = case_when(
        n_reanimation == 1 & n_chirurgie == 1     ~ "Réanimation",
        n_reanimation == 1 & n_medecine == 1      ~ "Réanimation",
        n_reanimation == 1 & n_gynecologie == 1   ~ "Réanimation",
        n_reanimation == 1 & n_ssr == 1           ~ "Réanimation",
        
        n_chirurgie == 1 & n_medecine == 1        ~ "Chirurgie",
        n_chirurgie == 1 & n_gynecologie == 1     ~ "Chirurgie",
        n_chirurgie == 1 & n_ssr == 1             ~ "Chirurgie",
        
        n_gynecologie == 1 & n_medecine == 1      ~ "Gynécologie-Obstétrique",
        n_gynecologie == 1 & n_ssr == 1           ~ "Gynécologie-Obstétrique",
        
        n_medecine == 1 & n_ssr == 1              ~ "Médecine",
        
        .default = "Médecine"
      )) %>% 
      dplyr::select(-c(n, n_chirurgie, n_medecine, n_reanimation, n_gynecologie)),
    by = c("code", "site", "prelevement", "Num_Patient", "age_cat", "bacterie", "year", "secteur")
  )
 
  
  # 4. Remove the most recent samples when there are less than 2 phenotypic 
  # variations and carried out in the same facility
  res = res %>%
    anti_join(., redundant_profiles, by = c("code", "site", "Num_Patient", "age_cat", "bacterie",
                                            "secteur", "ligne", "prelevement"))    
  
  # Save data-------------------------------------------------------------------
  res %>% 
    filter(prelevement >= as.Date(paste0(year, "/05/15")), prelevement <= as.Date(paste0(year, "/06/30"))) %>%
    mutate(
      Date_week = as.Date(cut(prelevement, "week")),
      Date_month = as.Date(cut(prelevement, "month"))
    ) %>%
    rename(Date_day = prelevement, atb_class = molecule_class, Date_year = year) %>%
    dplyr::select(code, Date_day, ligne, Num_Patient, Date_year, atb_class, bacterie, molecule) %>%
    distinct() %>%
    write.table(., paste0("data-raw/spares/combined/resistance_cohortfinal_representativity_", year, ".txt"),
                sep = "\t", quote = F, row.names = F)
}
