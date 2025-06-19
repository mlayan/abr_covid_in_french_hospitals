##################################################
## SELECTION OF BACTERIAL SAMPLES FOR THE
## FIVE SPECIES OF INTEREST
##################################################

rm(list = ls())
library(tidyverse)
library(foreach)
library(doParallel)
library(gt)

load("data/cohort19202122.rda")
source("R/helper/helper_functions.R")
source("R/helper/dictionaries.R")

##################################################
# Process resistance data
##################################################
registerDoParallel(cores = 2)
all_years = 2019:2022

foreach(y = all_years, .combine = "list", .packages = c("tidyverse")) %dopar% {

  out = ""
  out = paste0(out,
               "###############################################################\n",
               "Resistance data - ", y, "\n",
               "###############################################################\n"
               )

  # Basic checks----------------------------------------------------------------
  if (y == 2019) res = read.table("data-raw/spares/2019/Souches_BMRCovid_2019_d1 copie.csv",
                                     sep = "|", header = T, encoding = "latin1") %>%
    rename(code= IdEtablissement)
  if (y == 2020) res = read.table("data-raw/spares/2020/BMR_Covid_souches_2020.csv",
                                     sep = "|", header = T, encoding = "latin1") %>%
    rename(code= IdEtablissement)
  if (y == 2021) res = read.table("data-raw/spares/2021/BMR_Covid_souches_2021.csv",
                                     sep = "|", header = T, encoding = "latin1") %>%
    rename(code= IdEtablissement)
  if (y == 2022) res = read.table("data-raw/spares/2022/BMR_Covid_souches_2022.csv",
                                     sep = "|", header = T, encoding = "latin1") %>%
    rename(code= IdEtablissement)

  out = paste0(out,
               "Number of isolates: ",
               res %>%
                 select(-c(IdMolecule, molecule, nosocomial,Resultat, IdArchiveLaboratoire)) %>%
                 distinct() %>%
                 nrow(.),
             "\nNumber of isolates when removing unwanted columns: ",
             res %>%
               select(-c(IdMolecule, molecule, nosocomial,Resultat, IdArchiveLaboratoire,
                         IdPeriode, IdBacterie, suppr_dedoublonnage2, iduf, idactivite, code_de, code_ta, IdetablissementLaboratoire)) %>%
               distinct() %>%
               nrow(.),
             "\nNumber of duplicated isolated (IdArchiveLaboratoire): ",
             nrow(res %>%
                    select(-c(IdMolecule, IdBacterie, IdetablissementLaboratoire, nosocomial, code_ta, code_de, idactivite, iduf, suppr_dedoublonnage2, IdPeriode)) %>%
                    distinct() %>%
                    group_by(across(-c(IdArchiveLaboratoire, Resultat))) %>%
                    mutate(n = n()) %>%
                    filter(n > 1)),
             "\nNumber of samples with multiple bacterial species for one isolate: ",
             res %>%
               select(code, IdSite, secteur, prelevement, ligne, Num_Patient, bacterie) %>%
               distinct() %>%
               group_by(across(-bacterie)) %>%
               mutate(n = n()) %>%
               filter(n > 1) %>%
               nrow(.),
             "\nNumber of samples in ICUs: ",
             res %>%
               filter(secteur == "Réanimation") %>%
               select(-c(IdMolecule, molecule, nosocomial, Resultat, IdArchiveLaboratoire)) %>%
               distinct() %>%
               nrow(.)
  )

  # Filtering based on selected hospitals/departments, patient age, site, and
  # correct coding of results variable------------------------------------------
  res = res %>%
    # Filter out Pediatry, Psychiatry and SSL
    filter(secteur %in% c("Chirurgie", "Gynécologie-Obstétrique", "Médecine", "Réanimation", "SSR")) %>%
    # Filter out excluded facilities
    filter(code %in% cohort19202122) %>%
    # Filter out samples isolated in new borns
    filter(IdSite != 13) %>%
    # Filter out samples from patients under 15 years old
    filter(idtrancheage >= 3) %>%
    # Filter out test results/phenotypes that are incorrectly coded
    filter(
      !(Resultat %in% c("O", "N") & !molecule %in% enzymes) |
      (Resultat %in% c("S", "I", "R") & molecule %in% enzymes)
      ) %>%
    # Remove unnecessary columns
    # that were previously checked for potential coding issues
    select(-c(IdArchiveLaboratoire, IdPeriode, IdBacterie, IdMolecule, nosocomial, suppr_dedoublonnage2,
              iduf, idactivite, code_de, code_ta, IdetablissementLaboratoire)) %>%
    mutate(
      age_cat = recode(idtrancheage, !!!dict_id_age),
      molecule_class = recode(molecule, !!!dict_molecule_class),
      secteur = recode(secteur, !!!dict_secteur_spares),
      site = recode(IdSite, !!!dict_id_site),
      Date_year = y,
      bacterie = case_when(
        bacterie == "Enterobacter cloacae co" ~ "Enterobacter cloacae complex",
        bacterie == "Enterobacter cloacae" ~ "Enterobacter cloacae complex",
        bacterie == "Enterobacter asburiae" ~ "Enterobacter cloacae complex",
        bacterie == "Enterobacter hormaechei" ~ "Enterobacter cloacae complex",
        bacterie == "Enterobacter kobei" ~ "Enterobacter cloacae complex",
        bacterie == "Enterobacter ludwigii" ~ "Enterobacter cloacae complex",
        bacterie == "Enterobacter nimipressuralis" ~ "Enterobacter cloacae complex",
        .default = bacterie)
      ) %>%
    select(code, site, prelevement, ligne, Resultat, Num_Patient,
           age_cat, bacterie, molecule, molecule_class, secteur, Date_year) %>%
    distinct()

  # Get basic information-------------------------------------------------------
  molecules = sort(unique(res$molecule))
  sites = sort(unique(res$site))
  bacteria = sort(unique(res$bacterie))
  nlines = res %>% select(code, site, secteur, prelevement, ligne, Num_Patient) %>% distinct() %>% nrow(.)
  npatients = res %>% select(code, Num_Patient) %>% distinct() %>% nrow(.)

  # Filter out samples that are not from bacteria of interest-------------------
  res = res %>%
    filter(bacterie %in% bacteria_of_interest)

  # Get number of samples per bacteria species stratified by whether they were
  # tested for the antibiotic of interest---------------------------------------
  numbers_tested = left_join(
    res %>%
      distinct(code, site, prelevement, ligne, Num_Patient, age_cat, bacterie, Date_year) %>%
      count(Date_year, bacterie) %>%
      rename(all_samples = n),
    res %>%
      filter(
        (bacterie %in% c("Escherichia coli", "Klebsiella pneumoniae", "Enterobacter cloacae complex") & molecule == "BLSE") |
          (bacterie == "Pseudomonas aeruginosa" & molecule %in% c("Imipénème", "Méropénème")) |
          (bacterie == "Staphylococcus aureus" & molecule == "Oxacilline")
      ) %>%
      distinct(code, site, prelevement, ligne, Num_Patient, age_cat, bacterie, Date_year) %>%
      count(Date_year, bacterie) %>%
      rename(tested_samples = n),
    by = c("Date_year", "bacterie")
  )

  # Excluded samples
  res %>%
    group_by(code, site, prelevement, ligne, Num_Patient, age_cat, bacterie, secteur, Date_year) %>%
    mutate(to_keep = case_when(
      all(bacterie %in% c("Escherichia coli", "Klebsiella pneumoniae", "Enterobacter cloacae complex")) & any(molecule == "BLSE") ~ 0,
      all(bacterie == "Pseudomonas aeruginosa") & any(molecule %in% c("Imipénème", "Méropénème")) ~ 0,
      all(bacterie == "Staphylococcus aureus") & any(molecule == "Oxacilline") ~ 0,
      .default = 1
    )) %>%
    ungroup() %>%
    filter(to_keep == 1) %>%
    select(-to_keep) %>%
    rename(Date_day = prelevement, atb_class = molecule_class) %>%
    mutate(
      Date_day = as.Date(Date_day, "%d/%m/%Y"),
      Date_week = as.Date(cut(Date_day, "week")),
      Date_month = as.Date(cut(Date_day, "month"))
    ) %>%
    write.table(., paste0("data-raw/spares/combined/not_tested_cohort19202122_detailed_", y, ".txt"),
                sep = "\t", quote = F, row.names = F)

  # All samples
  res %>%
    rename(Date_day = prelevement, atb_class = molecule_class) %>%
    mutate(
      Date_day = as.Date(Date_day, "%d/%m/%Y"),
      Date_week = as.Date(cut(Date_day, "week")),
      Date_month = as.Date(cut(Date_day, "month"))
    ) %>%
    write.table(., paste0("data-raw/spares/combined/cohort19202122_detailed_", y, ".txt"),
                sep = "\t", quote = F, row.names = F)

  ##############################################################################
  # Get numbers-----------------------------------------------------------------
  ##############################################################################
  # 1. Samples that are not tested for the antibiotics used to define resistances
  n_samples_bacteria_of_interest = res %>%
    select(code, site, prelevement, Num_Patient, bacterie, secteur, ligne) %>%
    distinct() %>%
    nrow()

  n_samples_bacteria_of_interest_included = res %>%
    filter(
      (bacterie %in% bacteria_of_interest[1:3] & molecule %in% "BLSE") |
        (bacterie == bacteria_of_interest[4] & molecule %in% "Oxacilline") |
        (bacterie == bacteria_of_interest[5] & molecule %in% c("Imipénème", "Méropénème")) |
        (bacterie == bacteria_of_interest[6] & molecule %in% c("Imipénème", "Méropénème")) |
        (bacterie %in% bacteria_of_interest[7:8] & molecule %in% "Vancomycine")
    ) %>%
    select(code, site, prelevement, Num_Patient, bacterie, secteur, ligne) %>%
    distinct() %>%
    nrow()

  out = paste0(out,
               "\nNumber of samples of bacteria of interest: ",
               n_samples_bacteria_of_interest,
               "\nNumber of samples of bacteria of interest excluded because they were not tested for specific antibiotics: ",
               n_samples_bacteria_of_interest - n_samples_bacteria_of_interest_included
               )

  # 2. Phenotype of bacteria that are not tested for the antibiotics used to
  # define resistance-----------------------------------------------------------
  res %>%
    group_by(code, site, prelevement, Num_Patient, bacterie, secteur, ligne) %>%
    summarise(atb_test = sum(
      (bacterie %in% bacteria_of_interest[1:3] & molecule %in% "BLSE") |
        (bacterie == bacteria_of_interest[4] & molecule %in% "Oxacilline") |
        (bacterie == bacteria_of_interest[5] & molecule %in% c("Imipénème", "Méropénème")) |
        (bacterie == bacteria_of_interest[6] & molecule %in% c("Imipénème", "Méropénème")) |
        (bacterie %in% bacteria_of_interest[7:8] & molecule %in% "Vancomycine")
    )
    , .groups = "drop") %>%
    group_by(bacterie) %>%
    summarise(n = n(),
              n_not_tested = sum(atb_test == 0),
              p = sum(atb_test == 0) / n() *100,
              .groups = "drop") %>%
    mutate(Date_year = y) %>%
    write.table(., paste0("data-raw/spares/combined/not_tested_cohort19202122_", y, ".txt"),
                sep = "\t", quote = F, row.names = F)


  res %>%
    filter(bacterie %in% c("Escherichia coli", "Enterobacter cloacae complex")) %>%
    group_by(code, site, prelevement, Num_Patient, bacterie, secteur, ligne) %>%
    summarise(
      blse_tested = ifelse("BLSE" %in% molecule, "tested", "not tested"),
      blse_production = ifelse(any(molecule %in% "BLSE" & Resultat == "O"), 1, 0),
      n_beta_lactamins = sum(molecule_class %in% c("Broad spectrum Penicillins", "Carbapenems", "Cephalosporins", "Narrow spectrum Penicillins")),
      n_beta_res = sum(molecule_class %in% c("Broad spectrum Penicillins", "Carbapenems", "Cephalosporins", "Narrow spectrum Penicillins") & Resultat == "R"),
      n_beta_int = sum(molecule_class %in% c("Broad spectrum Penicillins", "Carbapenems", "Cephalosporins", "Narrow spectrum Penicillins") & Resultat == "I"),
      n_beta_s = sum(molecule_class %in% c("Broad spectrum Penicillins", "Carbapenems", "Cephalosporins", "Narrow spectrum Penicillins") & Resultat == "S"),
      .groups = "drop"
      ) %>%
    filter(blse_tested == "not tested") %>%
    group_by(bacterie) %>%
    summarise(
      n_r = sum(n_beta_lactamins == n_beta_res),
      n_s = sum(n_beta_lactamins == n_beta_s),
      n_int = sum(n_beta_lactamins != n_beta_s & n_beta_lactamins != n_beta_res),
      .groups = "drop"
    ) %>%
    mutate(Date_year = y) %>%
    write.table(., paste0("data-raw/spares/combined/enterobacter_c3g_cohort19202122_", y, ".txt"),
                sep = "\t", quote = F, row.names = F)

  # 3. Samples that are I among the bacteria of interest and the antibiotics used
  # to define resistance--------------------------------------------------------
  res %>%
    filter(
      (bacterie %in% bacteria_of_interest[1:3] & molecule %in% "BLSE") |
        (bacterie == bacteria_of_interest[4] & molecule %in% "Oxacilline") |
        (bacterie == bacteria_of_interest[5] & molecule %in% c("Imipénème", "Méropénème")) |
        (bacterie == bacteria_of_interest[6] & molecule %in% c("Imipénème", "Méropénème")) |
        (bacterie %in% bacteria_of_interest[7:8] & molecule %in% "Vancomycine")
    ) %>%
    group_by(code, site, prelevement, Num_Patient, bacterie, secteur, ligne, molecule) %>%
    summarise(phenotype = case_when(
      "R" %in% Resultat || "O" %in% "Resultat" ~ "R",
      "I" %in% Resultat ~ "I",
      .default = "S"
    ), .groups = "drop") %>%
    group_by(bacterie, molecule) %>%
    summarise(
      n_tot = n(),
      n_res = sum(phenotype == "R"),
      n_int = sum(phenotype %in% "I"),
      p_int_by_tot = sum(phenotype %in% "I") / n() * 100,
      p_int_by_nots = sum(phenotype %in% "I") / sum(!phenotype %in% "S") * 100,
      n_s = sum(phenotype %in% "S"),
      .groups = "drop"
      ) %>%
    mutate(Date_year = y) %>%
    write.table(., paste0("data-raw/spares/combined/phenotype_int_cohort19202122_", y, ".txt"),
                sep = "\t", quote = F, row.names = F)

  # 4. Verify that there is only one phenotype per sample defined at the "ligne"
  # level ----------------------------------------------------------------------
  ligne_level = res %>%
    group_by(code, site, prelevement, ligne, Num_Patient, age_cat, bacterie, secteur, Date_year, molecule) %>%
    summarise(n = n(), .groups = "drop") %>%
    filter(n > 1)

  if(nrow(ligne_level) > 0) stop("Some isolates defined at the ligne level have multiple phenotypes for an antibiotic")

  # 5. Samples that are isolated in the same individual, the same day, at the same
  # site for the same bacteria but with different sample IDs--------------------
  ids_duplicates = res %>%
    select(code, site, prelevement, Num_Patient, age_cat, bacterie, secteur, ligne) %>%
    distinct() %>%
    group_by(code, site, prelevement, Num_Patient, age_cat, bacterie, secteur) %>%
    mutate(n = n()) %>%
    filter(n>1) %>%
    ungroup() %>%
    select(-n)

  out = paste0(out,
               "\nNumber of potential duplicates: ",
               ids_duplicates %>% select(-ligne) %>% distinct() %>% nrow(.)
  )

  exact_duplicates = res %>%
    inner_join(., ids_duplicates, by = c("code", "site", "prelevement", "ligne",
                                         "Num_Patient", "age_cat", "bacterie", "secteur")) %>%
    group_by(code, site, prelevement, Num_Patient, age_cat, bacterie, secteur) %>%
    mutate(n_samples = length(unique(ligne))) %>%
    ungroup() %>%
    group_by(code, site, prelevement, Num_Patient, age_cat, bacterie, secteur, molecule) %>%
    mutate(n = n(), n_diff = length(unique(Resultat))) %>%
    ungroup() %>%
    group_by(code, site, prelevement, Num_Patient, age_cat, bacterie, secteur) %>%
    mutate(
      uniqueness = case_when(
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
    select(code, site, prelevement, Num_Patient, age_cat, bacterie, secteur, ligne) %>%
    distinct()

  rm(ids_duplicates)

  # 6. Samples that are part of a 30-day sequence-------------------------------
  samples_30days = res %>%
    anti_join(., duplicates_to_remove,
              by = c("code", "site", "prelevement", "Num_Patient", "age_cat", "bacterie", "secteur", "ligne")) %>%
    select(code, site, Num_Patient, age_cat, bacterie, secteur, ligne, prelevement) %>%
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

  out = paste0(out,
               "\nNumber of sequences of samples separated by less than 30 days: ",
               samples_30days %>%
                 select(code, site, Num_Patient, age_cat, bacterie, secteur) %>%
                 distinct() %>%
                 nrow(.)
  )

  redundant_profiles = samples_30days %>%
    select(-c(n, prelevement_lag, time_diff, n2)) %>%
    left_join(., res %>% mutate(prelevement = as.Date(prelevement, "%d/%m/%Y")),
              by = c("code", "site", "Num_Patient", "age_cat", "bacterie", "secteur", "ligne", "prelevement")) %>%
    group_by(code, site, Num_Patient, age_cat, bacterie, secteur) %>%
    nest() %>%
    mutate(selected = map(data, selection_30days_sequence)) %>%
    select(-data) %>%
    unnest(cols = selected) %>%
    ungroup()

  rm(samples_30days)

  ##############################################################################
  # Samples and phenotypes selection--------------------------------------------
  ##############################################################################
  # 1. Keep samples with lowest "ligne" ID among exact duplicates---------------
  exact_duplicates_to_remove = exact_duplicates %>%
    filter(uniqueness == 1) %>%
    group_by(code, site, prelevement, Num_Patient, age_cat, bacterie, secteur) %>%
    mutate(m = min(ligne)) %>%
    filter(ligne != m) %>%
    ungroup() %>%
    select(code, site, prelevement, Num_Patient, age_cat, bacterie, secteur, ligne) %>%
    distinct()

  res = res %>%
    anti_join(., exact_duplicates_to_remove,
              by = c("code", "site", "prelevement", "Num_Patient", "age_cat", "bacterie", "secteur", "ligne"))

  out = paste0(out,
               "\nNumber of exact duplicates: ",
               exact_duplicates %>%
                 filter(uniqueness == 1) %>%
                 select(code, site, prelevement, Num_Patient, age_cat, bacterie, secteur, ligne) %>%
                 distinct() %>%
                 nrow(),
               "\nNumber of samples from exact duplicates that were removed: ",
               nrow(exact_duplicates_to_remove)
  )

  # 2. Merge samples with same phenotype but different antibiotics tested-------
  uncomplete_duplicates_to_remove = exact_duplicates %>%
    filter(uniqueness == 2) %>%
    select(code, site, prelevement, Num_Patient, age_cat, bacterie, secteur, ligne) %>%
    distinct()

  uncomplete_duplicates_merged = exact_duplicates %>%
    filter(uniqueness == 2) %>%
    group_by(code, site, prelevement, Num_Patient, age_cat, bacterie, secteur) %>%
    mutate(ligne = min(ligne)) %>%
    ungroup() %>%
    group_by(code, site, prelevement, Num_Patient, age_cat, bacterie, molecule, molecule_class,
             ligne, secteur, Date_year) %>%
    summarise(Resultat = unique(Resultat), .groups = "drop") %>%
    distinct()

  res = res %>%
    anti_join(., uncomplete_duplicates_to_remove,
              by = c("code", "site", "prelevement", "Num_Patient", "age_cat", "bacterie", "secteur", "ligne")) %>%
    bind_rows(., uncomplete_duplicates_merged)

  out = paste0(out,
               "\nNumber of duplicated samples with same phenotypes but different combination of antibiotics: ",
               nrow(uncomplete_duplicates_to_remove),
               "\nNumber of samples with same phenotypes but different combination of antibiotics tested that were removed: ",
               nrow(uncomplete_duplicates_to_remove) -
                uncomplete_duplicates_merged %>%
                select(code, site, prelevement, Num_Patient, age_cat, bacterie, secteur, ligne) %>%
                distinct() %>%
                nrow(.)
  )

  # 3. Merge samples isolated the same day that have different phenotypes-------
  different_phenotypes = exact_duplicates %>%
    filter(uniqueness == 3) %>%
    group_by(code, site, prelevement, Num_Patient, age_cat, bacterie, secteur) %>%
    mutate(ligne = min(ligne)) %>%
    ungroup() %>%
    group_by(code, site, prelevement, Num_Patient, age_cat, bacterie, molecule, molecule_class, ligne, secteur, Date_year) %>%
    mutate( n = length(unique(Resultat))) %>%
    filter(n > 1, molecule == "BLSE" & bacterie %in% bacteria_of_interest[1:3] | molecule %in% c("Imipénème", "Méropénème") & bacterie %in% bacteria_of_interest[5:6] | molecule == "Vancomycine" & bacterie %in% bacteria_of_interest[7:8])

  different_duplicates_to_remove = exact_duplicates %>%
    filter(uniqueness == 3) %>%
    select(code, site, prelevement, Num_Patient, age_cat, bacterie, secteur, ligne) %>%
    distinct()

  different_duplicates_merged = exact_duplicates %>%
    filter(uniqueness == 3) %>%
    group_by(code, site, prelevement, Num_Patient, age_cat, bacterie, secteur) %>%
    mutate(ligne = min(ligne)) %>%
    ungroup() %>%
    group_by(code, site, prelevement, Num_Patient, age_cat, bacterie, molecule, molecule_class,
             ligne, secteur, Date_year) %>%
    summarise(Resultat = case_when(
      any(Resultat %in% "O") ~ "O",
      any(Resultat %in% "R") ~ "R",
      all(Resultat %in% "N") ~ "N",
      all(Resultat %in% c("I", "S")) ~ "S",
      all(Resultat %in% "S") ~ "S"
    ),
    .groups = "drop") %>%
    distinct()

  res = res %>%
    anti_join(., different_duplicates_to_remove,
              by = c("code", "site", "prelevement", "Num_Patient", "age_cat", "bacterie", "secteur", "ligne")) %>%
    bind_rows(., different_duplicates_merged)

  out = paste0(out,
          "\nNumber of duplicated samples with same antibiotics tested but different phenotypes: ",
          nrow(different_duplicates_to_remove),
          "\nNumber of samples with same antibiotics tested but different phenotypes that were removed: ",
          nrow(different_duplicates_to_remove) -
            different_duplicates_merged %>%
            select(code, site, prelevement, Num_Patient, age_cat, bacterie, secteur, ligne) %>%
            distinct() %>%
            nrow(.),
          "\nNumber of different phenotypes for same antibiotic in duplicates: ",
           different_phenotypes %>%
            nrow(.),
          "\nNumber of phenotypes with same antibiotic tested that were removed: ",
          different_phenotypes %>%
            anti_join(., different_duplicates_merged, by = c("code", "site", "prelevement", "ligne", "Resultat", "Num_Patient", "age_cat", "bacterie", "molecule", "molecule_class", "secteur", "Date_year")) %>%
            nrow(.)
  )

  # 4. Remove samples isolated in Chirurgie when already isolated in Médecine
  # on day d, at site s, in patient p, for bacteria b---------------------------
  sector_duplicates = res %>%
    select(code, site, prelevement, Num_Patient, age_cat, bacterie, Date_year, secteur) %>%
    distinct() %>%
    group_by(code, site, prelevement, Num_Patient, age_cat, bacterie, Date_year) %>%
    summarise(n = n(),
              n_chirurgie = sum(secteur %in% "Surgery"),
              n_medecine = sum(secteur %in% "Medicine"),
              n_reanimation = sum(secteur %in% "ICU"),
              n_gynecologie = sum(secteur %in% "Obstetrics and gynaecology"),
              n_ssr = sum(secteur %in% "Rehabilitation care"),
              .groups = "drop") %>%
    filter(n > 1)

  if (any(sector_duplicates$n_chirurgie > 1) |
      any(sector_duplicates$n_medecine > 1) |
      any(sector_duplicates$n_reanimation > 1) |
      any(sector_duplicates$n_gynecologie > 1) |
      any(sector_duplicates$n_ssr > 1))
    stop("Deduplication steps did not work properly")

  out = paste0(out,
               "\nNumber of isolates on the same day in two different sectors: ",
               sum(sector_duplicates$n)
               )

  res = anti_join(
    res,
    sector_duplicates %>%
      mutate(secteur = case_when(
        n_reanimation > 0                                         ~ "ICU",

        n_chirurgie > 0 & n_reanimation == 0                      ~ "Surgery",

        n_gynecologie > 0 & n_reanimation == 0 & n_chirurgie == 0 ~ "Obstetrics and gynaecology",

        .default = "Medicine"
      )) %>%
      select(-c(n, n_chirurgie, n_medecine, n_reanimation, n_gynecologie)),
    by = c("code", "site", "prelevement", "Num_Patient", "age_cat", "bacterie", "Date_year", "secteur")
    )

  # 5. Remove the most recent samples when there are less than 2 phenotypic
  # variations and carried out in the same facility
  out = paste0(out,
               "\nNumber of isolates that are removed due to similar phenotypic profile: ",
               nrow(redundant_profiles)
               )

  res = res %>%
    mutate(prelevement = as.Date(prelevement, "%d/%m/%Y")) %>%
    anti_join(., redundant_profiles, by = c("code", "site", "Num_Patient", "age_cat", "bacterie",
                                            "secteur", "ligne", "prelevement"))

  # 6. Remove samples from the first week of 2019 and the last week from 2021---
  res = res %>%
    mutate(prelevement = as.Date(prelevement, "%d/%m/%Y")) %>%
    mutate(
      Date_week = as.Date(cut(prelevement, "week")),
      Date_month = as.Date(cut(prelevement, "month"))
    )

  out = paste0(out,
               "\nNumber of isolates in the 1st week and last week of the study period: ",
               res %>%
                 filter(Date_week %in% c("2018-12-31", "2022-12-26")) %>%
                 select(code, site, prelevement, ligne, Num_Patient, age_cat, bacterie, secteur) %>%
                 distinct() %>%
                 nrow()
  )

  res = res %>%
    filter(!Date_week %in% c("2018-12-31", "2022-12-26"))

  ##############################################################################
  # Save data-------------------------------------------------------------------
  ##############################################################################
  # Save data-------------------------------------------------------------------
  res = res %>%
    rename(Date_day = prelevement, atb_class = molecule_class) %>%
    filter(
      (bacterie %in% bacteria_of_interest[1:3] & molecule %in% "BLSE") |
        (bacterie == bacteria_of_interest[4] & molecule %in% "Oxacilline") |
        (bacterie == bacteria_of_interest[5] & molecule %in% c("Imipénème", "Méropénème")) |
        (bacterie == bacteria_of_interest[6] & molecule %in% c("Imipénème", "Méropénème")) |
        (bacterie %in% bacteria_of_interest[7:8] & molecule %in% "Vancomycine")
    )

  res %>%
    write.table(., paste0("data-raw/spares/combined/resistance_cohort19202122_", y, ".txt"),
                sep = "\t", quote = F, row.names = F)

  # Save final number of samples------------------------------------------------
  numbers_tested = left_join(
    numbers_tested,
    res %>%
      distinct(code, site, Date_day, ligne, Num_Patient, age_cat, bacterie, Date_year) %>%
      count(Date_year, bacterie) %>%
      rename(final_samples = n),
    by = c("Date_year", "bacterie")
    )
  write.table(numbers_tested, paste0("data-raw/spares/combined/numbers_samples_", y, ".txt"),
              sep = "\t", quote = F, row.names = F)

  # Get updated basic information-----------------------------------------------
  new_nlines = res %>% select(code, site, secteur, Date_day, ligne, Num_Patient) %>% distinct() %>% nrow(.)
  new_npatients = res %>% select(code, Num_Patient) %>% distinct() %>% nrow(.)
  out = paste0(out,
          "\nFinal number of isolates: ",
          res %>%
            select(-c(molecule, atb_class, Resultat)) %>%
            distinct() %>%
            nrow(.),
          "\nFinal number of isolates in ICUs : ",
          res %>%
            select(-c(molecule, atb_class, Resultat)) %>%
            filter(secteur == "Réanimation") %>%
            distinct() %>%
            nrow(.)
          )
  writeLines(out, paste0("data-raw/spares/combined/sample_selection_", y, ".txt"))

  list(molecules, sites, bacteria, nlines, npatients, new_nlines, new_npatients)
}


##################################################
# Table of selected samples
##################################################
n_samples = data.frame()
for (f in list.files("data-raw/spares/combined/", "^numbers_samples_", full.names = T)) {
  n_samples = bind_rows(
    n_samples,
    read.table(f, header = T, sep = "\t")
  )
}

# Table with raw numbers
n_samples %>%
  pivot_longer(c(all_samples, tested_samples, final_samples), names_to = "metric", values_to = "n") %>%
  pivot_wider(names_from = Date_year, values_from = "n") %>%
  mutate(
    metric = case_when(
      metric == "all_samples" ~ "All",
      metric == "tested_samples" ~ "Tested",
      metric == "final_samples" ~ "Included"
    )
  ) %>%
  gt(groupname_col = "bacterie") %>%
  cols_label(metric = "Sample number")

# Table with raw numbers and percentages 
n_samples %>%
  mutate(p = paste0(tested_samples, " / ", all_samples, " (", round(tested_samples/all_samples*100), ")")) %>%
  select(-matches("samples")) %>%
  pivot_wider(names_from = Date_year, values_from = "p") %>%
  gt(rowname_col = "bacterie") 

##################################################
# Carbapenem phenotypes - EUCAST guidelines
##################################################
# Final cohort 
load("data/cohort_final.rda")

# Load resistance data for pseudomonas aeruginosa
pa = bind_rows( 
  read.table("data-raw/spares/combined/resistance_cohort19202122_2019.txt", 
             sep = "\t", header = T),
  read.table("data-raw/spares/combined/resistance_cohort19202122_2020.txt",
             sep = "\t", header = T),
  read.table("data-raw/spares/combined/resistance_cohort19202122_2021.txt",
             sep = "\t", header = T),
  read.table("data-raw/spares/combined/resistance_cohort19202122_2022.txt",
             sep = "\t", header = T)
) %>% 
  filter(bacterie == "Pseudomonas aeruginosa", code %in% cohort_final)
  
# Antibiotics used to define phenotypes 
pa_atb = pa %>%
  pivot_wider(names_from = molecule, values_from = Resultat) %>%
  mutate(
    atb = case_when(
      !is.na(Imipénème) & is.na(Méropénème) ~ "imi",
      is.na(Imipénème) & !is.na(Méropénème) ~ "mero",
      !is.na(Imipénème) & !is.na(Méropénème) ~ "imi+mero"
    ),
    phenotype = case_when(
      Imipénème == Méropénème & !is.na(Imipénème) & Imipénème == "R" ~ "R",
      Imipénème != Méropénème & !is.na(Imipénème) & !is.na(Méropénème) ~ "discordant",
      is.na(Imipénème) & Méropénème == "R" ~ "R",
      is.na(Méropénème) & Imipénème == "R" ~ "R",
      .default = "S"
    )
  )

# Check antibiotic-resistant phenotypes 
pa_atb %>% count(atb)
pa_atb %>% filter(Date_year == 2022) %>% filter(Méropénème == "R" | Imipénème == "R") %>% nrow()
pa_atb %>% filter(Date_year == 2022) %>% count(atb, phenotype)
pa_atb %>% mutate(atb_tested = ifelse(atb == "mero", atb, "other")) %>% count(atb_tested)
pa_atb %>% filter(phenotype == "discordant") %>% count(Imipénème, Méropénème)
pa_atb %>% filter(atb == "imi+mero") %>% group_by(Date_year) %>% 
  summarise(n_mero_r = sum(Imipénème != "R" & Méropénème == "R"), n_tot = n(), .groups = "drop") %>%
  mutate(p = n_mero_r / n_tot)

# Meningitis samples
pa %>% 
  distinct(code, Num_Patient, site, Date_day, ligne) %>%
  mutate(site = ifelse(site == "Liquide céphalorachidien", "Meningitis", "Other")) %>%
  count(site)

##################################################
# Stability of the proportion of excluded samples 
##################################################
# Final cohort
load("data/cohort_final.rda")

# All isolates
all_samples = bind_rows( 
  read.table("data-raw/spares/combined/cohort19202122_detailed_2019.txt", 
             sep = "\t", header = T),
  read.table("data-raw/spares/combined/cohort19202122_detailed_2020.txt",
             sep = "\t", header = T),
  read.table("data-raw/spares/combined/cohort19202122_detailed_2021.txt",
             sep = "\t", header = T),
  read.table("data-raw/spares/combined/cohort19202122_detailed_2022.txt",
             sep = "\t", header = T)
) %>%
  filter(code %in% cohort_final, code != 9657) %>%
  distinct(code, site, Num_Patient, ligne, bacterie, Date_day, Date_week, Date_month) %>%
  group_by(bacterie, Date_day, Date_week, Date_month) %>%
  summarise(n_all = n(), .groups = "drop")

# Load excluded isolates
excluded = bind_rows(
  read.table("data-raw/spares/combined/not_tested_cohort19202122_detailed_2019.txt",
             sep = "\t", header = T),
  read.table("data-raw/spares/combined/not_tested_cohort19202122_detailed_2020.txt",
             sep = "\t", header = T),
  read.table("data-raw/spares/combined/not_tested_cohort19202122_detailed_2021.txt",
             sep = "\t", header = T),
  read.table("data-raw/spares/combined/not_tested_cohort19202122_detailed_2022.txt",
             sep = "\t", header = T)
) %>% 
  filter(code %in% cohort_final, code != 9657) %>%
  distinct(code, site, Num_Patient, ligne, bacterie, Date_day, Date_week, Date_month) %>%
  group_by(bacterie, Date_day, Date_week, Date_month) %>%
  summarise(n_excluded = n(), .groups = "drop")
  
# Plot proportion by week
full_join(excluded, all_samples, by = c("bacterie", "Date_day", "Date_week", "Date_month")) %>%
  filter(!Date_week %in% c("2018-12-31", "2022-12-26")) %>%
  replace_na(list(n_all = 0, n_excluded = 0)) %>%
  group_by(Date_week, bacterie) %>%
  summarise(n_excluded = sum(n_excluded), n_all = sum(n_all), .groups = "drop") %>%
  pivot_longer(cols = c(n_all, n_excluded), names_to = "sample_type", values_to = "n") %>%
  mutate(Date_week = as_date(Date_week),
         sample_type = ifelse(sample_type == "n_all", "All isolats", "Excluded isolates")) %>%
  ggplot(., aes(x = Date_week, y = n, col = bacterie, linetype = sample_type)) +
  facet_wrap(facets = vars(bacterie), scales = "free_y", ncol = 2) +
  geom_line() +
  theme_bw() +
  theme(axis.title.x = element_blank(), legend.position = "bottom", legend.title = element_blank()) +
  scale_color_manual(values = c("Enterobacter cloacae complex" = "aquamarine", 
                                "Escherichia coli" = "aquamarine4",
                                "Klebsiella pneumoniae" = "bisque2", 
                                "Pseudomonas aeruginosa" = "indianred1", 
                                "Staphylococcus aureus" = "indianred4"),
                    ) +
  guides(color="none") +
  labs(y = "Weekly number of isolates") +
  expand_limits(ymin = 0)
ggsave("plots/bacterial_isolates/sample_numbers_exclusion.png", height = 7, width = 10)

# Plot ratio by week
full_join(excluded, all_samples, by = c("bacterie", "Date_day", "Date_week", "Date_month")) %>%
  filter(!Date_week %in% c("2018-12-31", "2022-12-26")) %>%
  replace_na(list(n_all = 0, n_excluded = 0)) %>%
  group_by(Date_week, bacterie) %>%
  summarise(n_excluded = sum(n_excluded), n_all = sum(n_all), .groups = "drop") %>%
  mutate(Date_week = as_date(Date_week)) %>%
  mutate(r = n_excluded / n_all ) %>%
  ggplot(., aes(x = Date_week, y = r, col = bacterie)) +
  facet_wrap(facets = vars(bacterie), scales = "free_y", ncol = 2) +
  geom_line() +
  theme_bw() +
  theme(legend.position = "none", axis.title.x = element_blank()) +
  scale_color_manual(values = c("Enterobacter cloacae complex" = "aquamarine", 
                                "Escherichia coli" = "aquamarine4",
                                "Klebsiella pneumoniae" = "bisque2", 
                                "Pseudomonas aeruginosa" = "indianred1", 
                                "Staphylococcus aureus" = "indianred4")) +
  labs(y = "Weekly proportion of excluded isolates") +
  expand_limits(ymin = 0)
ggsave("plots/bacterial_isolates/sample_proportion_exclusion.png", height = 7, width = 10)

# Investigation of dynamics of ESBL testing
bind_rows(
  read.table("data-raw/spares/combined/not_tested_cohort19202122_detailed_2019.txt",
             sep = "\t", header = T),
  read.table("data-raw/spares/combined/not_tested_cohort19202122_detailed_2020.txt",
             sep = "\t", header = T),
  read.table("data-raw/spares/combined/not_tested_cohort19202122_detailed_2021.txt",
             sep = "\t", header = T),
  read.table("data-raw/spares/combined/not_tested_cohort19202122_detailed_2022.txt",
             sep = "\t", header = T)
) %>% 
  filter(code %in% cohort_final, bacterie %in% c("Enterobacter cloacae complex", "Escherichia coli", "Klebsiella pneumoniae")) %>%
  distinct(code, site, Num_Patient, ligne, age_cat, bacterie, Date_day, Date_week, Date_month) %>%
  left_join(., metadata_admin %>% distinct(code, region), by = "code") %>%
  filter(region == "Auvergne-Rhône-Alpes") %>%
  group_by(bacterie, code, Date_week, Date_month) %>%
  summarise(n_excluded = n(), .groups = "drop") %>%
  ggplot(., aes(x = as_date(Date_week), y = n_excluded)) +
  geom_line() +
  facet_grid(rows = vars(code), cols = vars(bacterie), scales = "free_y") +
  theme_bw()

# Bump of samples not tested for BLSE in 2021
data21 = read.table("data-raw/spares/2021/BMR_Covid_souches_2021.csv", 
           sep = "|", header = T, encoding = "latin1") %>%
  rename(code= IdEtablissement) %>%
  filter(secteur %in% c("Chirurgie", "Gynécologie-Obstétrique", "Médecine", "Réanimation", "SSR")) %>%
  # Filter out excluded facilities
  filter(code %in% cohort19202122) %>%
  # Filter out samples isolated in new borns
  filter(IdSite != 13) %>%
  # Filter out samples from patients under 15 years old
  filter(idtrancheage >= 3) %>%
  # Filter out test results/phenotypes that are incorrectly coded
  filter(
    !(Resultat %in% c("O", "N") & !molecule %in% enzymes) |
      (Resultat %in% c("S", "I", "R") & molecule %in% enzymes)
  ) %>%
  # Remove unnecessary columns
  # that were previously checked for potential coding issues
  select(-c(IdArchiveLaboratoire, IdPeriode, IdBacterie, IdMolecule, nosocomial, suppr_dedoublonnage2,
            iduf, idactivite, code_de, code_ta, IdetablissementLaboratoire)) %>%
  mutate(
    age_cat = recode(idtrancheage, !!!dict_id_age),
    molecule_class = recode(molecule, !!!dict_molecule_class),
    secteur = recode(secteur, !!!dict_secteur_spares),
    site = recode(IdSite, !!!dict_id_site),
    Date_year = 2021,
    bacterie = case_when(
      bacterie == "Enterobacter cloacae co" ~ "Enterobacter cloacae complex",
      bacterie == "Enterobacter cloacae" ~ "Enterobacter cloacae complex",
      bacterie == "Enterobacter asburiae" ~ "Enterobacter cloacae complex",
      bacterie == "Enterobacter hormaechei" ~ "Enterobacter cloacae complex",
      bacterie == "Enterobacter kobei" ~ "Enterobacter cloacae complex",
      bacterie == "Enterobacter ludwigii" ~ "Enterobacter cloacae complex",
      bacterie == "Enterobacter nimipressuralis" ~ "Enterobacter cloacae complex",
      .default = bacterie)
  ) %>%
  select(code, site, prelevement, ligne, Resultat, Num_Patient,
         age_cat, bacterie, molecule, molecule_class, secteur, Date_year) %>%
  distinct() %>%
  filter(bacterie %in% bacteria_of_interest)

##################################################
# Susceptibility test results for macrolides 
# other than azithromycin
##################################################
# All samples
selected_samples = bind_rows( 
  read.table("data-raw/spares/combined/cohort19202122_detailed_2019.txt", 
             sep = "\t", header = T),
  read.table("data-raw/spares/combined/cohort19202122_detailed_2020.txt",
             sep = "\t", header = T),
  read.table("data-raw/spares/combined/cohort19202122_detailed_2021.txt",
             sep = "\t", header = T),
  read.table("data-raw/spares/combined/cohort19202122_detailed_2022.txt",
             sep = "\t", header = T)
) %>%
  filter(code %in% cohort_final) %>%
  mutate(
    Date_week = as_date(Date_week),
    molecule = case_when(molecule == "Pristinamycine" ~ "Pristinamycin", 
                         molecule == "Erythromycine" ~ "Erythromycin",
                         .default = molecule)
    )

# Which samples are tested for macrolides other than azithromycin
selected_samples %>% 
  filter(atb_class == "Macrolides") %>%
  count(bacterie, molecule)

# Plot resistance phenotype counts 
selected_samples %>% 
  filter(atb_class == "Macrolides") %>%
  group_by(code, site, Date_day, Date_week, ligne, Num_Patient, age_cat, bacterie, molecule) %>%
  summarise(phenotype = ifelse(any(Resultat == "R"), "R", "S"), .groups = "drop") %>%
  group_by(Date_week, bacterie, molecule) %>%
  summarise(Resistant = sum(phenotype == "R"), Total = n(), .groups = "drop") %>%
  pivot_longer(c(Resistant, Total), names_to = "Phenotype", values_to = "Count") %>%
  ggplot(., aes(x = Date_week, y = Count, col = Phenotype)) +
  geom_line() +
  facet_grid(rows = vars(molecule)) +
  scale_color_manual(values = c("Resistant" = "red3", "Total" = "navyblue")) +
  theme_bw()

# Plot resistance incidence rate
load("data/bd_pmsi_hospital.rda")
selected_samples %>% 
  filter(atb_class == "Macrolides") %>%
  group_by(code, site, Date_day, Date_week, ligne, Num_Patient, age_cat, bacterie, molecule) %>%
  summarise(phenotype = ifelse(any(Resultat == "R"), "R", "S"), .groups = "drop") %>%
  group_by(Date_week, bacterie, molecule) %>%
  summarise(Resistant = sum(phenotype == "R"), Total = n(), .groups = "drop") %>%
  right_join(., bd_pmsi_hospital, by = "Date_week") %>%
  mutate(i = ifelse(is.na(Resistant), 0, Resistant / nbjh * 1000)) %>%
  ggplot(., aes(x = as_date(Date_week), y = i, col = molecule)) +
  geom_line() +
  scale_color_manual(values = c("Erythromycin" = "paleturquoise3", "Pristinamycin" = "saddlebrown")) +
  theme_bw() +
  theme(axis.title.x = element_blank()) +
  labs(y = "Incidence of resistant isolates\n(for 1,000 bed-days)", col = "Resistance")
ggsave("plots/antibiotic_resistance/macrolide_saureus.png", height = 2.5, width = 6)

