# Impact of the COVID-19 pandemic on antibiotic resistance in French hospitals : a retrospective ecological study using national surveillance data

<!--a href="https://doi.org/10.1101/2024.12.04.24317990"><img src="" alt="DOI"></a-->

Full text is available [here](https://doi.org/10.1101/2024.12.04.24317990), and Supplementary Materials [here](https://www.medrxiv.org/content/10.1101/2024.12.04.24317990v1.supplementary-material).

## Usage

### Dependencies
The code has been tested using R version 4.3.0. It depends on the following R libraries: 
- Data management: `readxl`, `tidyverse`, `gt`
- Plots : `RColorBrewer`, `ggpubr`, `cowplot`, `sf`, `geofacet`
- Descriptive statistical analyses: `rstatix`
- Statistical modelling : `MASS`, `DescTools`, `performance`

### Bacterial sample and hospital selection procedures
Codes used to analyze the raw data from the national surveillance system on antibiotic resistance in hospitals, SPARES, are not made available to ensure hospital anonymity. 

Briefly, these codes performed the following procedures:
- Selection of hospitals according to our inclusion and exclusion criteria
- Selection of bacterial samples according to our inclusion and exclusion criteria
- Assessment of the regional representativity of our cohort (in terms of hospital numbers, hospital activity in bed-days, and prevalence of healthcare-associated infections)
- Aggregation of the antibiotic resistance, antibiotic use, and bed-days data at the national and regional level

Several figures were generated with these codes and are available in the `plots/cohort_description` directory. 

### Descriptive analyses
One can find in the `R/analyses` directory the R codes to reproduce the figures of the paper as well as the statistical analyses. All descriptive analyses can be reproduced using the following scripts:

- `8.bed_days_time_series.R`: plot total and intubated COVID-19 bed-days at the hospital and ICU levels ;
- `9.resistance_time_series.R`: plot incidence and resistance proportion of resistant bacteria at the hospital and ICU levels ;
- `10.antibiotic_consumption.R`: perform the descriptive analyses of antibiotic use. This code uses antibiotic consumption at the hospital or ICU level. For anonymity reasons, we did not make the data available ;
- `11.temporal_analyses_correlations.R`: perform univariate analyses of the association between COVID-19 variables and incidence of resistant bacteria at the ICU and hospital levels.

### Statistical modelling
Aggregated data at the national or regional level are available in the `data` directory. They were analyzed with the following codes:
- `12.temporal_analyses_national.R`: all statistical analyses at the national level for both hospitals and ICUs, including sensitivity analyses for carbapenem-resistant *Pseudomonas aeruginosa* (CR-PA) ;
- `13.temporal_analyses_regions.R`: statistical analyses of the association between CR-PA incidence in hospitals and the prevalence of intubated COVID-19 patients at week *w-2* across 12 French regions.   

