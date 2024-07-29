##################################################
## PLOTS OF BED-DAYS DATA
##################################################
rm(list = ls())
library(tidyverse)
library(ggpubr)
source("R/helper/dictionaries.R")
source("R/helper/helper_plots.R")

load("data/int_national_start_end.rda")

load("data/bd_pmsi_hospital.rda")
load("data/bd_pmsi_icu.rda")

load("data/covid_intub_icu.rda")
load("data/covid_intub_hospital.rda")

##################################################
# Paper figure 3 - Bed days - Total and Covid-19
##################################################
# Total bed-days
plot_bd = bind_rows(
  bd_pmsi_hospital %>% mutate(setting = "Hospital"),
  bd_pmsi_icu %>% mutate(setting = "ICU")
) %>%
  ggplot(., aes(x = Date_week, y = nbjh)) +
  geom_rect(data = int_national_start_end, 
            aes(NULL, NULL, xmin=start, xmax=end, fill=restrictions, ymin=-Inf, ymax=Inf)) +
  geom_line() +
  facet_wrap(facets = vars(setting), scales = "free_y") +
  theme_bw() +
  theme(legend.position = "bottom",  
        legend.key = element_rect(colour = "black"),
        legend.title = element_text(hjust = 0.5)) +
  guides(fill = guide_legend(title.position = "top")) +
  scale_fill_manual(
    name = "Anti-Covid-19 interventions",
    labels = c("First wave", "Strong", "Intermediary", "Light to none"),
    breaks = c("p_first_wave", "p_strong_res", "p_mild_res", "p_no_res"),
    values = c("p_first_wave" = col_interventions(1), "p_strong_res" = col_interventions(2), "p_mild_res" = col_interventions(3), "p_no_res" = col_interventions(4))
  ) +
  scale_y_continuous(labels = scales::label_comma()) +
  expand_limits(y = 0) +
  labs(x = "", y = "Weekly bed-days")

# Intubated Covid-19 bed-days
all_intubated_covid = bind_rows(
  covid_intub_hospital %>% mutate(setting = "Hospital"),
  covid_intub_icu %>% mutate(setting = "ICU")
) 

zero_intub_covid = data.frame(
  expand.grid(
    Date_week = all_dates[!all_dates %in% all_intubated_covid$Date_week],
    setting = c("ICU", "Hospital"),
    covid_intub = 0 
  )
)
  

plot_covid = bind_rows(zero_intub_covid, all_intubated_covid) %>%
  ggplot(., aes(x = Date_week, y = covid_intub)) +
  geom_rect(data = int_national_start_end, 
            aes(NULL, NULL, xmin=start, xmax=end, fill=restrictions, ymin=-Inf, ymax=Inf)) +
  geom_line() +
  facet_wrap(facets = vars(setting), scales = "free_y") +
  theme_bw() +
  theme(legend.position = "bottom",  
        legend.key = element_rect(colour = "black"),
        legend.title.align = 0.5) +
  guides(fill = guide_legend(title.position = "top")) +
  scale_fill_manual(
    name = "Anti-Covid-19 interventions",
    labels = c("First wave", "Strong", "Intermediary", "Light to none"),
    breaks = c("p_first_wave", "p_strong_res", "p_mild_res", "p_no_res"),
    values = c("p_first_wave" = col_interventions(1), "p_strong_res" = col_interventions(2), "p_mild_res" = col_interventions(3), "p_no_res" = col_interventions(4))
  ) +
  expand_limits(y = 0) +
  scale_y_continuous(labels = scales::label_comma()) +
  labs(x = "", y = "Weekly intubated Covid-19 bed-days")

# Final figure
figure2 = ggarrange(plot_bd, plot_covid, nrow = 2, 
                    labels = c("A", "B"), hjust = -1.5,
                    common.legend = T, legend = "bottom")
figure2
ggsave("../Paper/Figures/Figure2.png", figure2, height = 5, width = 8)
ggsave("plots/Figure2.png", figure2, height = 5, width = 8)
