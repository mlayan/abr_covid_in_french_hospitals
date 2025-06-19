# Impacts of the COVID-19 pandemic on antibiotic use and resistance in hospitals : a retrospective ecological analysis of French national surveillance data over 2019-2022

<!--a href="https://doi.org/10.1101/2024.12.04.24317990"><img src="" alt="DOI"></a-->

Full text is available [here](https://doi.org/10.1101/2024.12.04.24317990), and Supplementary Materials [here](https://www.medrxiv.org/content/10.1101/2024.12.04.24317990v1.supplementary-material).

## Usage

### Dependencies
The code has been tested using R version 4.4.1. It depends on the following R libraries: 
- Data management: `readxl`, `tidyverse`, `gt`
- Plots : `ggpubr`, `cowplot`, `sf`
- Descriptive statistical analyses: `rstatix`
- Statistical modelling : `MASS`, `DescTools`, `performance`

### Bacterial sample and hospital selection procedures
Codes used to analyze the raw data from the national surveillance system on antibiotic resistance in hospitals (SPARES) and the national hospital discharge database (PMSI) are not made available to ensure hospital anonymity. 

Briefly, these codes performed the following procedures:
- Selection of hospitals according to our inclusion and exclusion criteria
- Selection of bacterial samples according to our inclusion and exclusion criteria
- Assessment of the regional representativeness of our cohort (in terms of hospital numbers, hospital activity in bed-days, and prevalence of healthcare-associated infections)
- Assessment of the impact of the hospital selection procedure on annual antibiotic use and resistance proportion
- Aggregation of the antibiotic resistance, antibiotic use, and bed-days data at the national and regional level

Several figures were generated with these codes and are available in the `plots/cohort_description` directory. 

### Descriptive analyses
One can find in the `R/analyses` directory the R codes to reproduce the figures of the paper as well as the statistical analyses. All descriptive analyses can be reproduced using the following scripts:

- `8.bed_days_time_series.R`: plot total and COVID-19 bed-days at the hospital and ICU levels ;
- `9.resistance_time_series.R`: plot incidence and resistance proportion of resistant bacteria at the hospital and ICU levels ;
- `10.antibiotic_consumption.R`: perform the descriptive analyses of antibiotic use. This code uses antibiotic consumption at the hospital or ICU level. For anonymity reasons, we did not make the data available.

### Statistical modelling
Aggregated data at the national or regional level are available in the `data` directory. They were analyzed with the following codes:
- `11.temporal_analyses_national.R`: all statistical analyses at the national level for both hospitals and ICUs, including sensitivity analyses for carbapenem-resistant *Pseudomonas aeruginosa* (CR-PA) where we tested the association between CR-PA incidence and COVID-19 prevalence at weeks *w+1*, *w+2*, and *w+3* or where we investigated the association between CR-PA isolates collected in wards different from the ICU and COVID-19 variables ;
- `12.temporal_analyses_national_sensitivity.R`: sensitivity analysis at the national level where annual antibiotic use has not been taken into account for both hospitals and ICUs ;
- `13.temporal_analyses_national_bloodstream.R`: statistical analysis of the bloodstream infections at the national level for both hospitals and ICUs ;
- `14.temporal_analyses_regions.R`: statistical analyses of the association between CR-PA, MRSA, or ESBL-EC incidence in hospitals and COVID-19 variables across 12 French regions ;
- `15.temporal_analyses_regions_sensitivity.R`: sensitivity analysis at the regional level where antibiotic use has not been taken into account.   

