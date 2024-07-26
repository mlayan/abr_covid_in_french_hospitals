# Impact of the COVID-19 pandemic on antibiotic resistance in French hospitals : a retrospective ecological study using national surveillance data

<a href=""><img src="" alt="DOI"></a>

Full text and supplementary materials are available [here]().

## Usage

### Dependencies
The code has been tested using R version 4.3.0. It depends on the following R libraries: 
- Data management: `readxl`, `tidyverse`, `gt`
- Plots : `RColorBrewer`, `ggpubr`, `sf`, `geofacet`
- Statistical modelling : `MASS`

### Bacterial sample and hospital selection procedures
Codes used to analyze the raw data from the national surveillance system on antibiotic resistance in hospitals, SPARES, are not made available to ensure hospital anonymity. 

Briefly, these codes performed the following procedures:
- Selection of hospitals according to our inclusion and exclusion criteria
- Selection of bacterial samples according to our inclusion and exclusion criteria
- Assessment of the regional representativity of our cohort (in terms of hospital numbers, hospital activity in bed-days, and prevalence of healthcare-associated infections)
- Aggregation of the antibiotic resistance, antibiotic use, and bed-days data at the national and regional level

Several figures were generated with these codes and are available in the `plots/cohort_description` directory. 

### Descriptive analyses
Briefly, these codes performed the following procedures:
- Selection of hospitals according to our inclusion and exclusion criteria
- Selection of bacterial samples according to our inclusion and exclusion criteria
- Assessment of the regional representativity of our cohort (in terms of hospital numbers, hospital activity in bed-days, and prevalence of healthcare-associated infections)
- Aggregation of the antibiotic resistance, antibiotic use, and bed-days data at the national and regional level

### Statistical modelling
