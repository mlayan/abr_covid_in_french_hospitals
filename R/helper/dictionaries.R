# THIS DICTIONARY FILE USES THE FOLLOWING DEPENDENCIES
  # dplyr
  # geofacet
  # sf
  # readxl
  # RColorBrewer

# Bacteria of interest 
bacteria_of_interest = c("Enterobacter cloacae complex", 
                         "Escherichia coli", 
                         "Klebsiella pneumoniae",
                         "Staphylococcus aureus", 
                         "Pseudomonas aeruginosa")

# All week dates of the study period
all_dates = seq(as.Date("2019-01-07"), as.Date("2022-12-19"), 7)

# Enzymes
enzymes = c("BLSE", "Carbapénèmase", "Céphalosporinase déréprimée")

# Beta-lactamins
beta_lactams = c("BLSE", "Carbapénèmase", "Céphalosporinase déréprimée")

# Last line antibiotics
last_line = c("Cefiderocol", 
              "Ceftazidime avibactam", 
              "Ceftolozane tazobactam", 
              "Imipenem+relebactam", 
              "Meropenem+vaborbactam")

# IdSite
dict_id_site = c(
  "1"	= "Hémoculture",
  "2" =	"Urine sauf matériel de sondage",
  "3"	= "Dispositif intravasculaire",
  "4"	= "Liquide céphalorachidien",
  "5"	= "Liquide pleural",
  "6"	= "Liquide articulaire",
  "7"	= "Liquide d ascite",
  "8"	= "Prélèvement profond",
  "9"	= "Prélèvement respiratoire protégé ou distal",
  "10" = "Prélèvement respiratoire non protégé", 
  "11" = "Coproculture",
  "12" = "Prélèvement génital urétral cervico vaginal endocol",
  "13" = "Prélèvement nouveau né liquide gastrique multi sites",
  "14" = "Pus superficiel",
  "15" = "Autre prélèvement"
)


# site plot
dict_site_plot = c(
  "Hémoculture" = "Hémoculture",
  "Urine : à l'exclusion du matériel de sondage" = "Urine",
  "Dispositif intravasculaire" = "Dispositif IV",
  "Liquide céphalorachidien" = "LCR",
  "Liquide pleural" = "Liquide pleural",
  "Liquide articulaire" = "Liquide articulaire",
  "Liquide d'ascite" = "Liquide d'ascite",
  "Prélèvement profond" = "Pvt profond",
  "Prélèvement respiratoire protégé ou distal" = "Pvt respi distal",
  "Prélèvement respiratoire non protégé" = "Pvt respi", 
  "Coproculture" = "Coproculture",
  "Prélèvement génital: urétral, cervico-vaginal, endocol" = "Pvt génital",
  "Prélèvement nouveau-né : liquide gastrique, multi-sites" = "Pvt nouveau-né",
  "Pus superficiel" = "Pus",
  "Autre prélèvement" = "Autre"
)

# idactivite
dict_id_activite = c(
  "1" = "Médecine",
  "2" = "Réanimation",
  "3" = "Psychatrie",
  "4" = "SSR", 
  "5" = "SLD",
  "6" = "Gynéco-Obstétrique", 
  "7" = "Pédiatrie",
  "8" = "Chirurgie",
  "14" = "Urgences"
)

# idtrancheage
dict_id_age = c(
  "1" = "[0,5[",
  "2" = "[5,15[",
  "3" = "[15,25[",
  "4" = "[25,45[",
  "5" = "[45,65[",
  "6"	= "65+"
)

# nosocomial
dict_nosocomial = c(
  "0" = "community-onset", 
  "1" = "hospital-onset", 
  "9" = "undetermined"
)

# Hospital type
dict_hospital_type = c(
  "CH" = "General public hospital",
  "CHU" = "University hospital",
  "ESLD" = "Long-term care facility",
  "ESPIC" = "Private not-for-profit hospital",
  "ESSR" = "Rehabilitation hospital",
  "LOC" = "General public hospital",
  "MCO" = "Private for profit hospital"
)

# Antibiotic class full (https://www.whocc.no/atc_ddd_index/)
dict_antibiotic_class = c(
  "J01GB01" = "Aminoglycosides",
  "J01GB06" = "Aminoglycosides",
  "J01GB03" = "Aminoglycosides",
  
  "J01DH03" = "Carbapenems",
  "J01DH51" = "Carbapenems",
  "J01DH02" = "Carbapenems",
  
  "J01DB05" = "Cephalosporins",
  "J01DB01" = "Cephalosporins",
  "J01DB04" = "Cephalosporins",
  "J01DC04" = "Cephalosporins",
  "J01DC03" = "Cephalosporins",
  "J01DC07" = "Cephalosporins",
  "J01DC01" = "Cephalosporins",
  "J01DC02" = "Cephalosporins",
  "J01DD08" = "Cephalosporins",
  "J01DD01" = "Cephalosporins",
  "J01DD13" = "Cephalosporins",
  "J01DD02" = "Cephalosporins",
  "J01DD52" = "Cephalosporins",
  "J01DD04" = "Cephalosporins", 
  "J01DE01" = "Cephalosporins",
  "J01DI02" = "Cephalosporins",
  "J01DI01" = "Cephalosporins",
  "J01DI54" = "Cephalosporins", 
  
  "J01XA04" = "Glycopeptides",
  "J01XA01" = "Glycopeptides",
  "J01XA02" = "Glycopeptides",

  "J01XX09" = "Lipopeptides",
  
  "J01FA01" = "Macrolides",
  "J01FA07" = "Macrolides", 
  "J01FA10" = "Macrolides",
  "J01FA09" = "Macrolides", 
  "J01FA03" = "Macrolides",
  "J01FA06" = "Macrolides",
  "J01FA02" = "Macrolides", 
  "J01FA15" = "Macrolides",
  "J01FF02" = "Macrolides", 
  "J01FF01" = "Macrolides",
  "J01FG01" = "Macrolides", 
  
  "J01DF01" = "Monobactams",
  
  "J01XX08" = "Oxazolidinones", 
  
  "J01CA04" = "Penicillins",
  "J01CA01" = "Penicillins",
  "J01CA08" = "Penicillins",
  "J01CA17" = "Penicillins",
  "J01CA13" = "Penicillins",
  "J01CA12" = "Penicillins",
  "J01CR02" = "Penicillins",
  "J01CR01" = "Penicillins",
  "J01CR03" = "Penicillins", 
  "J01CR05" = "Penicillins",
  "J01CE08" = "Penicillins",
  "J01CE01" = "Penicillins",
  "J01CE02" = "Penicillins",
  "J01CF02" = "Penicillins",
  "J01CF04" = "Penicillins",
  
  "J01XX01" = "Fosfomycin",
  
  "J01XB01" = "Polymyxins",
  
  "J01AA01" = "Tetracyclines",
  "J01AA02" = "Tetracyclines",
  "J01AA12" = "Tetracyclines", 
  "J01AA04" = "Tetracyclines", 
  "J01AA05" = "Tetracyclines",
  "J01AA08" = "Tetracyclines",
  
  "J01MA12" = "Quinolones",
  "J01MA23" = "Quinolones",
  "J01MA07" = "Quinolones",
  "J01MA02" = "Quinolones",
  "J01MA06" = "Quinolones",
  "J01MA01" = "Quinolones", 
  "J01MA14" = "Quinolones",
  
  "J01EA01" = "Trimethoprim",
  "J01EC02" = "Trimethoprim",
  "J01EE01" = "Trimethoprim",
  
  "J04AB02" = "Other", 
  "J01BA02" = "Other",
  "J01EB02" = "Other", 
  "J01GA01" = "Other", 
  "J01RA04" = "Other",
  "J01RA02" = "Other",
  "J01XC01" = "Other",
  "J01XD01" = "Other",
  "J01XD03" = "Other",
  "J01XE01" = "Other",
  "J01XX11" = "Other",
  "A07AA12" = "Other",
  "P01AB01" = "Other", 
  "P01AB03" = "Other", 
  "P01AB02" = "Other",
  
  "J01DH56" = "Anti-MDR GNB",
  "J01DH52" = "Anti-MDR GNB",
  "J01DI04" = "Anti-MDR GNB"
)

# 3GC
vec_3gc = c("J01DD08", "J01DD01", "J01DD13", "J01DD02", "J01DD52", "J01DD04")

# Dictionary of antibiotic class and antibiotic name (for microbiological data)
dict_molecule_class = c(
  "Acide fusidique" = "Other", 
  "Amikacine" = "Aminoglycosides", 
  "Amoxicilline - acide clavulanique" = "Broad spectrum Penicillins",
  "Amoxicilline" = "Broad spectrum Penicillins", 
  "Ampicilline" = "Broad spectrum Penicillins",
  "Ampicilline Sulbactam" = "Broad spectrum Penicillins",
  "Azithromycine" = "Macrolides",
  "Aztréonam" = "Monobactams",
  "Benzathine" = "Narrow spectrum Penicillins",
  "Benzylpénicilline - pénicilline G" = "Narrow spectrum Penicillins",
  "Céfaclor" = "Cephalosporins",
  "Céfadroxil" = "Cephalosporins",
  "Céfalexine" = "Cephalosporins",
  "Céfamandole" = "Cephalosporins",
  "Céfazoline" = "Cephalosporins",
  "Céfépime" = "Cephalosporins",
  "Céfixime" = "Cephalosporins",
  "Céfotaxime" = "Cephalosporins",
  "Céfotiam" = "Cephalosporins",
  "Céfoxitine" = "Cephalosporins",
  "Cefpodoxime" = "Cephalosporins",
  "Ceftaroline" = "Cephalosporins",
  "Ceftazidime" = "Cephalosporins",
  "Ceftazidime avibactam" = "Cephalosporins",
  "Ceftobiprole" = "Cephalosporins",
  "Ceftolozane tazobactam" = "Cephalosporins",
  "Ceftriaxone" = "Cephalosporins",
  "Céfuroxime" = "Cephalosporins",
  "Ciprofloxacine" = "Quinolones",
  "Clarithromycine" = "Macrolides",
  "Clindamycine" = "Macrolides",
  "Cloxacilline" = "Narrow spectrum Penicillins",
  "Colistine" = "Polymyxins",
  "Colistine/colimycine" = "Polymyxins",
  "Cotrimoxazole" = "Trimethoprim",
  "Dalbavancine" = "Glycopeptides",
  "Daptomycine" = "Lipopeptides",
  "Demeclocycline" = "Tetracyclines",
  "Doxycycline" = "Tetracyclines",
  "Ertapénème" = "Carbapenems",
  "Erythromycine" = "Macrolides",
  "Erythromycine+sulfafurazole" = "Other",
  "Fidaxomicine" = "Other",
  "Fosfomycine" = "Phosphonics",
  "Gentamicine" = "Aminoglycosides",
  "Imipénème" = "Carbapenems",
  "Josamycine" = "Macrolides",
  "Kanamycine" = "Aminoglycosides",
  "Lévofloxacine" = "Quinolones",
  "Lincomycine" = "Macrolides",
  "Linézolide" = "Oxazolidinones",
  "Loméfloxacine" = "Quinolones",
  "Lymécycline" = "Tetracyclines",
  "Méropénème" = "Carbapenems",
  "Meropenem+vaborbactam" = "Carbapenems",
  "Metacycline" = "Tetracyclines",
  "Métronidazole" = "Other",
  "Midécamycine" = "Macrolides",
  "Minocycline" = "Tetracyclines",
  "Moxifloxacine" = "Quinolones",
  "Nétilmicine" = "Aminoglycosides",
  "Nitrofurantoïne" = "Other",
  "Norfloxacine" = "Quinolones",
  "Ofloxacine" = "Quinolones",
  "Ornidazole" = "Other", 
  "Oxacilline" = 'Narrow spectrum Penicillins',
  "Pénicilline V" = "Narrow spectrum Penicillins",
  "Pipéracilline" = "Broad spectrum Penicillins", 
  "Pipéracilline - tazobactam" = "Broad spectrum Penicillins",
  "Pivmécillinam" = "Broad spectrum Penicillins",
  "Pristinamycine" = "Macrolides",
  "Rifampicine" = 'Other',
  "Roxithromycine" = "Macrolides",
  "Spiramycine" = "Macrolides",
  "Spiramycine+métronidazole" = "Other",
  "Streptomycine" = "Other", 
  "Sulfadiazine" = "Trimethoprim",
  "Sulfaméthoxazole - triméthoprime" = "Trimethoprim",
  "Sulfamethizol" = "Other",
  "Tédizolide" = "Other",
  "Teicoplanine" = "Glycopeptides",
  "Télithromycine" = "Macrolides",
  "Tétracycline" = "Tetracyclines",
  "Témocilline" = "Broad spectrum Penicillins",
  "Thiamphénicol" = "Other",
  "Ticarcilline" = "Broad spectrum Penicillins",
  "Ticarcilline ac clavulanique" = "Broad spectrum Penicillins",
  "Tigecycline" = "Tetracyclines",
  "Tinidazole" = "Other", 
  "Tobramycine" = "Aminoglycosides",
  "Triméthoprime" = "Trimethoprim",
  "Vancomycine" = "Glycopeptides"
)

# Bacterial species names in EARS-NET
dict_bacterie_earsnet = c(
  "ACISPP" = "CR A spp", 
  "ENCFAE" = "VR E. faecalis", 
  "ENCFAI" = "VR E. faecium", 
  "ESCCOL" = "ESBL E. coli", 
  "KLEPNE" = "ESBL K. pneumoniae", 
  "PSEAER" = "CR P. aeruginosa", 
  "STAAUR" = "MRSA"
)

# Secteurs hôpitaux 
dict_secteur_sivic = c(
  "Hospitalisation conventionnelle" = "Hospitalisation conventionnelle", 
  "Hospitalisation en SSR" = "SSR", 
  "Hospitalisation psychiatrique" = "Psychiatrie", 
  "Hospitalisation réanimatoire (réa ou SI)" = "Réanimation", 
  "Soins aux urgences" = "Urgences",
  "Total" = "Total", 
  "USLD"="SLD"
)

# Secteurs hôpitaux 
dict_secteur_spares = c(
  "Chirurgie" = "Surgery", 
  "Gynécologie-Obstétrique" = "Obstetrics and gynaecology", 
  "Hématologie" = "Hematology",
  "Maladie inf" = "Infectious diseases", 
  "Médecine" = "Medicine", 
  "Pédiatrie" = "Pediatry", 
  "Psychiatrie" = "Psychiatry",
  "Réanimation" = "ICU",
  "SLD"="LTC",
  "SSR" = "Rehabilitation care", 
  "Total établissement" = "Total"
)

# Regions
dict_regions = c(
  "Auvergne-Rhône-Alpes" = "ARA",
  "Bourgogne-Franche-Comté" = "BFC",
  "Bretagne" = "BRE",
  "Centre-Val de Loire" = "CVL",
  "Grand-Est" = "GES",
  "Hauts-de-France" = "HDF",
  "Île-de-France" = "IDF",
  "Normandie" = "NOR",
  "Nouvelle-Aquitaine" = "NAQ",
  "Occitanie" = "OCC",
  "Pays de la Loire" = "PDL",
  "Provence-Alpes-Côte d'Azur" = "PAC"
)

dict_regions_inv = c(
  "ARA" = "Auvergne-Rhône-Alpes",
  "BFC" = "Bourgogne-Franche-Comté",
  "BRE" = "Bretagne",
  "CVL" = "Centre-Val de Loire",
  "GES" = "Grand-Est",
  "HDF" = "Hauts-de-France",
  "IDF" = "Île-de-France",
  "NOR" = "Normandie",
  "NAQ" = "Nouvelle-Aquitaine",
  "OCC" = "Occitanie",
  "PDL" = "Pays de la Loire",
  "PAC" = "Provence-Alpes-Côte d'Azur"
)

dict_regions_nb = c(
  "84" = "Auvergne-Rhône-Alpes",
  "27" = "Bourgogne-Franche-Comté",
  "53" = "Bretagne",
  "24" = "Centre-Val de Loire",
  "94" = "Corse",
  "44" = "Grand-Est",
  "1" = "Guadeloupe",
  "3" = "Guyane",
  "32" = "Hauts-de-France",
  "11" = "Île-de-France",
  "2" = "Martinique",
  "6" = "Mayotte",
  "28" = "Normandie",
  "75" = "Nouvelle-Aquitaine",
  "76" = "Occitanie",
  "52" = "Pays de la Loire",
  "93" = "Provence-Alpes-Côte d'Azur",
  "4" = "La Réunion"
)

# Correspondences departments/regions France
corres = read.csv("data-raw/geography/departement_2022.csv") %>%
  dplyr::select(DEP, REG) %>%
  rename(dep = DEP, reg = REG) %>%
  mutate(region = recode(reg, !!!dict_regions_nb)) %>%
  dplyr::select(-reg) %>%
  filter(!dep %in% c("2A", "2B", "971", "972", "973", "974", "976"))

# Regional grid for geofacet
my_regional_grid = geofacet::fr_regions_grid1
my_regional_grid = subset(my_regional_grid, !name %in% c("Guadeloupe", "Martinique", "Guyane", "Réunion", "Mayotte", "Corse"))
my_regional_grid$name[my_regional_grid$name == "Grand Est"] = "Grand-Est"

# Database on CHU in France
chu_france = readxl::read_excel("data-raw/spares/finess_issues/all_chu_finess.xlsx")

# Hexagonal France map
france = sf::st_read("data-raw/geography/data_gouv/regions-20161121.shp") %>%
  rename(region = name) %>% 
  filter(!region %in% c("Guadeloupe", "Martinique", "Guyane", "La Réunion", "Mayotte", "Corse")) %>% 
  dplyr::select(region, geometry) %>%
  mutate(region = ifelse(region == "Pays-de-la-Loire", "Pays de la Loire", region))

# finess juridique for CHU
chu_finess_juridique = data.frame(
  finess_juridique = c("800000044", "590780193", "290000017", "350005179", "140000100", "760780239", "750712184", "210780581", "250000015", "540023264", "510000029", "670780055", "490000031", "440000289", "370000481", "450000088", "690781810", "380780080", "630780989", "420784878", "330781196", "870000015", "860010560", "310781406", "300780038", "340780477", "060785011", "130786049", "570005165"),
  city = c("AMIENS", "LILLE", "BREST", "RENNES", "CAEN", "ROUEN", "PARIS", "DIJON", "BESANCON", "NANCY", "REIMS", "STRASBOURG", "ANGERS", "NANTES", "TOURS", "ORLEANS", "LYON", "GRENOBLE", "CLERMONT FERRAND", "ST ETIENNE", "BORDEAUX", "LIMOGES", "POITIERS", "TOULOUSE", "NIMES", "MONTPELLIER", "NICE", "MARSEILLE", "METZ"),
  region = c("HDF", "HDF", "BRE", "BRE", "NOR", "NOR", "IDF", "BFC", "BFC", "GES", "GES", "GES", "PDL", "PDL", "CVL", "CVL", "ARA", "ARA", "ARA", "ARA", "NAQ", "NAQ", "NAQ", "OCC", "OCC", "OCC", "PAC", "PAC", "GES")
)

# Dictionary of departments
dict_departments = read.csv("data-raw/geography/departement_2022.csv")$LIBELLE
names(dict_departments) = read.csv("data-raw/geography/departement_2022.csv")$DEP

# Count regression model names
model_names = c(
  "model0" = "No Covid-19 variable",
  "model1" = "Pandemic periods w",
  "model2" = "Pandemic periods w-1",
  "model3" = "Pandemic periods w-2",
  "model4" = "Covid-19 prevalence w", 
  "model5" = "Covid-19 prevalence w-1", 
  "model6" = "Covid-19 prevalence w-2"
)

# Hospital types
hospital_types_dict = c(
"Centre Hospitalier (C.H.)",
"Centre Hospitalier Régional (C.H.R.)",
"Centre de Médecine Sportive",
"Centre de Médecine Universitaire",
"Centre de Médecine collective",
"Centre de Santé",
"Centre hospitalier, ex Hôpital local",
"Etablissement Soins Obstétriques Chirurgico-Gynécologiques",
"Etablissement de Soins Chirurgicaux",
"Etablissement de Soins Médicaux",
"Etablissement de Soins Pluridisciplinaire",
"Etablissement de santé privé autorisé en SSR",
"Groupement de coopération sanitaire - Etablissement de santé",
"Groupement de coopération sanitaire de moyens",
"Groupement de coopération sanitaire de moyens - Exploitant",
"Installation autonome de chirurgie esthétique"
)

# COVID-19 parameters names
covid_var_names = c(
  'periodsfirst wave' = "First wave w",
  'periodsstrong res' = "Strong w",
  'periodsmild res' = "Mild w",
  'periodslow to no res' = "Low to none w",
  
  'lag1_periodsfirst wave' = "First wave w-1",
  'lag1_periodsstrong res' = "Strong w-1",
  'lag1_periodsmild res' = "Mild w-1",
  'lag1_periodslow to no res' = "Low to none w-1",
  
  "lag2_periodsfirst wave" = "First wave w-2",
  "lag2_periodsstrong res" = "Strong w-2",
  "lag2_periodsmild res" = "Mild w-2",
  "lag2_periodslow to no res" = "Low to none w-2",
  
  "lag2_covid_prev" = "Covid-19 prev. w-2"
)

# Model parameters names
var_names = c(
  "(Intercept)" = "Intercept",
  "lag1_i_res" = "Incidence w-1",
  
  "periodsfirst wave" = "Pandemic periods",
  "periodsstrong res" = "Pandemic periods",
  "periodsmild res" = "Pandemic periods",
  "periodslow to no res" = "Pandemic periods",
  
  "lag1_periodsfirst wave" = "Pandemic periods w-1",
  "lag1_periodsstrong res" = "Pandemic periods w-1",
  "lag1_periodsmild res" = "Pandemic periods w-1",
  "lag1_periodslow to no res" = "Pandemic periods w-1",
  
  "lag2_covid_prev" = "Covid-19 prev. w-2",
  
  "Penicillins" = "Penicillins",
  "TGC" = "3rd generation Cephalosporins",
  "Carbapenems" = "Imipenem + Meropenem"
)

