# Function to plot antibiotic consumption per year and region
atbplot = function(df, atbclass, level = "regional", facet_type = "geographic") {
  
  if (level == "regional") {
    # Consumption per antibiotic class
    df_final = df %>%
      group_by(code, Date_year, region) %>%
      summarise(consumption = sum(molDDD)/unique(Nbhosp)*1000, .groups = "drop")
    m = max(df_final$consumption)
    
    # P value of Friedman test
    pvals = df_final %>%
      group_by(region) %>%
      friedman_test(consumption ~ Date_year | code) %>%
      mutate(p_sign = factor(case_when(
        p > 0.05 ~ "n.s.",
        p <= 0.05 & p > 0.01 ~ "*",
        p <= 0.01 & p > 0.001 ~ "**",
        p <= 0.001 ~ "***"
      ), 
      levels = c("n.s.", "*", "**", "***"),
      labels = c("n.s.", "*", "**", "***")
    ))
    
    # Mean consumption per region and year
    df_mean = df_final %>%
      group_by(Date_year, region) %>%
      summarise(consumption = mean(consumption), .groups = "drop") %>%
      left_join(., pvals %>% select(region, p_sign), by = "region")
    
    # Significance levels of multiple comparison tests
    df_comp = df_final %>%
      arrange(code, Date_year) %>%
      group_by(region) %>%
      wilcox_test(
        consumption ~ Date_year, 
        paired = T, 
        p.adjust.method = "bonferroni", 
        alternative = "two.sided",
        conf.level = 1-0.05/6,
        detailed = T
      ) %>%
      mutate(group1 = factor(group1), group2 = factor(group2)) %>%
      filter(p.adj <= 0.05, group1 == "2019") %>%
      mutate(y = case_when(group2 == "2020" ~ m/2, group2 == "2021" ~ m*2/3, group2 == "2022" ~ m*5/6)) %>%
      mutate(p.adj = ifelse(p.adj<0.001, "<0.001", round(p.adj, 3)))
    
    # Plot
    colfunc <- colorRampPalette(c("floralwhite", "deepskyblue4"))
    out = ggplot(df_final, aes(x = as.factor(Date_year))) +
      geom_line(alpha = 0.1, aes(y = consumption, group = code)) +
      geom_point(data = df_mean, aes(x = as.factor(Date_year), y = consumption, fill = p_sign), 
                 col = "black", pch=21, size = 3, stroke = 0.2) +
      geom_signif(data = df_comp, aes(xmin = group1, xmax = group2, annotations = p.adj,
                                      y_position = y),
                  manual = T, textsize = 3, vjust = -0.2) +
      theme_bw() +
      scale_fill_manual(
        limits = c("n.s.", "*", "**", "***"),
        values = colfunc(4),
        name = "Friedman test p-value") +
      theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom", legend.title = element_text(hjust=0.5)) +
      guides(fill = guide_legend(title.position = "top")) +
      labs(x = "", y = "Annual antibiotic consumption\n(DDD/1,000 bed days)")
    
    if (facet_type == "geographic") out = out + facet_geo(~ region, grid = my_regional_grid)
    if (facet_type == "alphabetical") out = out + facet_wrap(facets = vars(region), ncol = 4)
    # ggsave(paste0("plots/cons_atb/", tolower(gsub(" ", "", atbclass)), ".png"), 
    #        out, 
    #        height = 7, width = 10)
    return(out)
  }
  
  if (level == "type") {
    # Consumption per antibiotic class
    df = df %>%
      mutate(consumption = molDDD/Nbhosp*1000)
    df_max = df %>%
      group_by(type) %>%
      summarise(m = max(consumption), .groups = "drop")
    
    # P value of Friedman test
    pvals = df %>%
      group_by(type) %>%
      friedman_test(consumption ~ Date_year | code) %>%
      mutate(p_sign = factor(case_when(
        p > 0.05 ~ "n.s.",
        p <= 0.05 & p > 0.01 ~ "*",
        p <= 0.01 & p > 0.001 ~ "**",
        p <= 0.001 ~ "***"
      ), 
      levels = c("n.s.", "*", "**", "***"),
      labels = c("n.s.", "*", "**", "***")
      ))
    
    # Mean consumption per type and year
    df_mean = df %>%
      group_by(Date_year, type) %>%
      summarise(consumption = mean(consumption), .groups = "drop") %>%
      left_join(., pvals %>% select(type, p_sign), by = "type")
    
    # Significance levels of multiple comparison tests
    df_comp = df %>%
      arrange(code, Date_year) %>%
      group_by(type) %>%
      wilcox_test(
        consumption ~ Date_year, 
        paired = T, 
        p.adjust.method = "bonferroni", 
        alternative = "two.sided",
        conf.level = 1-0.05/6,
        detailed = T
      ) %>%
      mutate(group1 = factor(group1), group2 = factor(group2)) %>%
      filter(p.adj <= 0.05, group1 == "2019") %>%
      left_join(., df_max, by = "type") %>%
      mutate(y = case_when(group2 == "2020" ~ m/2, group2 == "2021" ~ m*2/3, group2 == "2022" ~ m*5/6))
    
    # Plot
    colfunc <- colorRampPalette(c("floralwhite", "deepskyblue4"))
    out = ggplot(df, aes(x = as.factor(Date_year))) +
      geom_line(alpha = 0.1, aes(y = consumption, group = code)) +
      geom_point(data = df_mean, aes(x = as.factor(Date_year), y = consumption, fill = p_sign), 
                 col = "black", pch=21, size = 2, stroke = 0.2) +
      geom_signif(data = df_comp, aes(xmin = group1, xmax = group2, annotations = round(p.adj, 4),
                                      y_position = y),
                  manual = T, textsize = 3, vjust = -0.2) +
      facet_wrap(facets = vars(type), ncol = 3, scales = "free_y") +
      theme_bw() +
      scale_fill_manual(
        limits = c("n.s.", "*", "**", "***"),
        values = colfunc(4),
        name = "Friedman test p-value") +
      theme(plot.title = element_text(hjust = 0.5), legend.title = element_text(hjust = 0.5), legend.position = "bottom") +
      guides(col = guide_legend(title.position = "top")) +
      labs(x = "", y = "Annual antibiotic consumption\n(DDD/1,000 bed days)", 
           title = atbclass)

    # ggsave(paste0("plots/cons_atb/type_", tolower(gsub(" ", "", atbclass)), ".png"), 
    #        out, 
    #        height = 7, width = 10)
    return(out)
  }
 
  if (level == "national") {
    df_final = df %>%
      mutate(consumption = molDDD / Nbhosp * 1000) %>%
      arrange(code, Date_year) %>%
      group_by(code) %>%
      mutate(n = n()) %>%
      ungroup() %>%
      filter(n == 4)
    m = max(df$consumption)
    
    df_mean = df %>%
      group_by(Date_year) %>%
      summarise(consumption = mean(consumption), .groups = "drop")
    
    df_comp = df %>%
      arrange(code, Date_year) %>%
      wilcox_test(
        consumption ~ Date_year, 
        paired = T, 
        p.adjust.method = "bonferroni", 
        alternative = "two.sided",
        conf.level = 1-0.05/6,
        detailed = T
      ) %>%
      mutate(group1 = factor(group1), group2 = factor(group2)) %>%
      filter(p.adj <= 0.05, group1 == "2019") %>%
      mutate(y = case_when(group2 == "2020" ~ m/2, group2 == "2021" ~ m*2/3, group2 == "2022" ~ m*5/6))
    
    out = ggplot() +
      geom_line(data = df_final, aes(x = factor(Date_year), y = consumption, group = code), alpha = 0.1) +
      geom_point(data = df_mean, aes(x = factor(Date_year), y = consumption), col = "red") +
      geom_signif(data = df_comp, aes(xmin = group1, xmax = group2, annotations = p.adj, y_position = y), 
                  manual = T) +
      theme_bw() +
      labs(x = "", y = "Annual antibiotic consumption\n(DDD/1,000 hospitalisation days)", 
           title = atb_class)
    # ggsave(paste0("plots/cons_atb/", tolower(gsub(" ", "", atbclass)), ".png"), 
    #        out, height = 4, width = 5)
    return(out)
  }
}


col_interventions <- function(i) { 
  paste0(colorRampPalette(c("steelblue4", "ivory"))(4)[i], "99")
}
