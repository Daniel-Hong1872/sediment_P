setwd("C:/Users/dan91/OneDrive/桌面/Field work")
library(tidyverse)
library(readxl)
library(lme4)
library(lmerTest)
library(emmeans)
library(patchwork)
library(performance)

sediment <- read_excel("CNP_data.xlsx", sheet = "sediment", na = "NA")
P <- read_excel("CNP_data.xlsx", sheet = "P", na = "NA")
IO <- read_excel("CNP_data.xlsx", sheet = "Inorganic vs. Organic", na = "NA")

sediment <- sediment %>% 
  rename(
    sed_dish = `sample A(g) sediment+dish`,
    dish = `dish(g)`
  ) %>%
  mutate(
    sed_wg = sed_dish - dish
  )

# mg/g = (mg/L * mL / 1000) / (mg / 1000) = mg/L * mL / mg
V_ml <- 5
P <- P %>%
  rename(
    A_cmgL = `sample A_P concentration(mg/L)`,
    B_cmgL = `sample B_P concentration(mg/L)`,
    A_wmg = `A subsample_P weight(mg)`,
    B_wmg = `B subsample_P weight(mg)`
  ) %>%
  mutate(
    bulk_P_cmgg = (A_cmgL * V_ml / 1000) / A_wmg * 1000,
    algae_P_cmgg = (B_cmgL * V_ml / 1000)/ B_wmg * 1000
  )

IO <- IO %>%
  rename(
    tube = `tube(g)`,
    subsample = `subsample(g)`,
    wo_salt = `tube+subsample w/o salt`,
    tube_inorg  = `tube+inorganic weight(g)`
  ) %>%
  mutate(
    salt_g = tube + subsample - wo_salt,
    inorg_g = tube_inorg - tube,
    org_g = wo_salt - tube_inorg,
    salt_pct = salt_g  / subsample * 100,
    inorg_pct = inorg_g / subsample * 100,
    org_pct = org_g   / subsample * 100
  )

region_color <- c(
  "GI" = "#4CAF50",
  "NE" = "#3288BD",
  "XLQ" = "#D53E4F"
)

sed_site_mean <- sediment %>% 
  group_by(Region, site) %>%
  summarise(
    depth = mean(depth, na.rm = T),
    temperature = mean(temperature, na.rm = T),
    sediment_weight = mean(sed_wg, na.rm = T),
    .groups = "drop"
  )

sed_depth <- 
  ggplot(sed_site_mean, aes(x = depth, y = sediment_weight)) + 
  geom_point() +
  geom_smooth(method = "lm", se = F) +
  facet_wrap(~Region) +
  theme_bw()
sed_depth

sed_temp <- 
  ggplot(sed_site_mean, aes(x = temperature, y = sediment_weight)) + 
  geom_point() +
  geom_smooth(method = "lm", se = F) +
  facet_wrap(~Region) +
  theme_bw()
sed_temp

sed_region <- 
  ggplot(sed_site_mean, aes(Region, sediment_weight, fill = Region)) +
  geom_boxplot() +
  scale_fill_manual(values = region_color) +
  guides(fill = "none") +
  geom_jitter() +
  theme_bw()
sed_region

#sea urchin----

XLQ_ur <- read_excel("urchin density.xlsx", sheet = "XLQ") %>% 
  mutate(Region = "XLQ")
GI_ur <- read_excel("urchin density.xlsx", sheet = "GI") %>%
  mutate(Region = "GI")
NE_ur <- read_excel("urchin density.xlsx", sheet = "NE") %>%
  mutate(Region = "NE")

urchin_all <- bind_rows(XLQ_ur, GI_ur, NE_ur)

herb_urchin <- urchin_all %>%
  filter(Genus %in% c("Diadema", "Echinothrix", "Stomopneustes"))

herb_urchin_region <- herb_urchin %>%
  group_by(Region, Genus) %>%
  summarise(Total = n(), .groups = "drop")

urchin_color <- c(
  "Diadema" = "#6BAED6",
  "Echinothrix" = "#8FCB9B",
  "Stomopneustes" = "#6C6F9A"
)

herb_urchin_sum_region <- 
  ggplot(herb_urchin_region, aes(Region, Total, fill = Genus)) +
  geom_col() +
  scale_fill_manual(values = urchin_color, drop = F) +
  theme_bw() +
  labs(x = "Region", y = "Total number of urchins", fill = "Genus")
herb_urchin_sum_region

# Relationship between P concentration & Region----

P_site_mean <- P %>% 
  group_by(Region, site) %>%
  summarise(
    bulk_P = mean(bulk_P_cmgg, na.rm = T),
    algae_P = mean(algae_P_cmgg, na.rm = T),
    .groups = "drop"
  )

bulk_P_region <- 
  ggplot(P_site_mean, aes(Region, bulk_P, fill = Region)) +
  geom_boxplot() +
  scale_fill_manual(values = region_color) +
  guides(fill = "none") +
  geom_jitter() +
  theme_bw()
bulk_P_region

algae_P_region <- 
  ggplot(P_site_mean, aes(Region, algae_P, fill = Region)) +
  geom_boxplot() +
  scale_fill_manual(values = region_color) +
  guides(fill = "none") +
  geom_jitter() +
  theme_bw()
algae_P_region

kruskal.test(bulk_P ~ Region, data = P_site_mean)
pairwise.wilcox.test(
  P_site_mean$bulk_P, 
  P_site_mean$Region,
  p.adjust.method = "BH"
)

kruskal.test(algae_P ~ Region, data = P_site_mean)
pairwise.wilcox.test(
  P_site_mean$algae_P, 
  P_site_mean$Region,
  p.adjust.method = "BH"
)

cor.test(P_site_mean$bulk_P,
         P_site_mean$algae_P,
         method = "spearman")

# Fish bite----

bite <- read_excel("Fish bite_data.xlsx")

area <- 1
duration <- 35
bite <- bite %>%
  mutate(bite_rate = bites / area * 60 / duration, 
         `trophic group` = fct_relevel(
           `trophic group`,
           "Herbivore", "Benthic invertivore", "Corallivore", "Omnivore"
         ))

bite_quad_trophic <- bite %>%
  group_by(Region, site, Camera, `trophic group`) %>%
  summarise(
    total_bites = sum(bites, na.rm = T),
    .groups = "drop"
  ) %>%
  mutate(
    bite_rate = total_bites * 60 / duration) %>%
  filter(!is.na(`trophic group`))

trophic_color <- c(
  "Herbivore" = "#4DB6AC",
  "Benthic invertivore" = "#FDB863",
  "Corallivore" = "#80B1D3",
  "Omnivore" = "#BEAED4"
)

trophic_region <- 
  ggplot(bite_quad_trophic, 
         aes(`trophic group`, bite_rate, fill = `trophic group`)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.6) +
  stat_summary(fun = mean, geom = "point", color = "red") +
  facet_wrap(~ Region, nrow = 1) +
  scale_fill_manual(values = trophic_color, drop = F) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  labs(
    x = "Trophic group",
    y = expression(Bite~rate~(bites~m^{-2}~hr^{-1}))
  )
trophic_region

bite_quad_func <- bite %>%
  filter(`trophic group` == "Herbivore") %>%
  group_by(Region, site, Camera, `herbivore functional group`) %>%
  summarise(
    total_bites = sum(bites, na.rm = T),
    .groups = "drop"
  ) %>%
  mutate(
    bite_rate = total_bites * 60 / duration
  )

func_group_color <- c(
  "browsers" = "#BC9C6D",
  "farmers (territorial croppers)" = "#D8B365",
  "grazers (croppers)" = "#7FBF7B",
  "scrapers" = "#80B1D3",
  "small excavators" = "#BEAED4"
)

herb_func_region <- 
  ggplot(bite_quad_func, aes(`herbivore functional group`, bite_rate,
      fill = `herbivore functional group`)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.6) +
  stat_summary(fun = mean, geom = "point", color = "red") +
  facet_wrap(~ Region, nrow = 1) +
  scale_fill_manual(values = func_group_color, drop = F) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none") +
  labs(
    x = "Herbivore functional group",
    y = expression(Bite~rate~(bites~m^{-2}~hr^{-1}))
  )
herb_func_region

# Bite rate by event----

all_camera <- expand.grid(
  Region = unique(bite$Region),
  site = unique(bite$site),
  Camera = c("C1", "C2", "C3")
)

bite_by_quadrat <- bite %>%
  filter(`trophic group` == "Herbivore") %>%
  group_by(Region, site, Camera) %>%
  summarise(
    total_bites = sum(bites, na.rm = T),
    n_events = n(),
    mean_bites_per_event = mean(bites, na.rm = T),
    .groups = "drop"
  ) 

bite_by_quadrat <- all_camera %>%
  left_join(bite_by_quadrat, by = c("Region", "site", "Camera"))

bite_by_quadrat <- bite_by_quadrat %>%
  mutate(
    total_bites = replace_na(total_bites, 0),
    n_events = replace_na(n_events, 0),
    total_bite_rate = total_bites * 60 / duration,
    event_rate = n_events * 60 / duration
  )

# occurrence of herbivore feeding among regions----

camera_zero_summary <- bite_by_quadrat %>%
  mutate(has_herbivore = n_events > 0) %>%
  group_by(Region) %>%
  summarise(
    n_camera = n(),
    n_with_herbivore = sum(has_herbivore),
    prop_with_herbivore = mean(has_herbivore),
    prop_zero = mean(!has_herbivore),
    .groups = "drop"
  )

ggplot(camera_zero_summary, aes(x = Region, y = prop_with_herbivore, fill = Region)) +
  geom_col() +
  scale_fill_manual(values = region_color) +
  scale_y_continuous(labels = scales::percent_format()) +
  guides(fill = "none") +
  theme_bw() +
  labs(
    x = "Region",
    y = "Proportion of cameras with herbivore feeding",
    title = "Occurrence of herbivore feeding among regions"
  )
# herbivore grazing pressure----

herb_bite_positive <- bite_by_quadrat %>%
  filter(n_events > 0)

bite_by_quadrat_region <- herb_bite_positive %>%
  select(Region, total_bite_rate, event_rate, mean_bites_per_event) %>%
  pivot_longer(
    cols = c(total_bite_rate, event_rate, mean_bites_per_event),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(metric = factor(
    metric, levels = c("total_bite_rate","event_rate","mean_bites_per_event"), 
    labels = c("Total bite rate\n(bites/hr)",
               "Feeding event rate\n(events/hr)",
               "Bites per event\n(bites/event)")
  ))

herb_grazing_pressure <- 
  ggplot(bite_by_quadrat_region, aes(Region, value, fill = Region)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.6) +
  stat_summary(fun = mean, geom = "point", color = "red") +
  scale_fill_manual(values = region_color) +
  facet_wrap(~metric, scales = "free_y", nrow = 1) +
  guides(fill = "none") +
  theme_bw() +
  labs(
    title = "Herbivore Grazing Pressure",
    x = "Region",
    y = NULL
  )
herb_grazing_pressure

kruskal.test(total_bite_rate ~ Region, data = bite_by_quadrat)
kruskal.test(event_rate ~ Region, data = bite_by_quadrat)
kruskal.test(mean_bites_per_event ~ Region, data = herb_bite_positive)

pairwise.wilcox.test(
  bite_by_quadrat$total_bite_rate,
  bite_by_quadrat$Region,
  p.adjust.method = "BH"
)

pairwise.wilcox.test(
  bite_by_quadrat$event_rate,
  bite_by_quadrat$Region,
  p.adjust.method = "BH"
)

# Bite rate considered of organic/inorganic percentage----

# 整理 IO data

IO <- IO %>%
  rename(sediment_id = sediment)

IO_quadrat <- IO %>%
  mutate(Camera = gsub("S", "C", sediment_id)) %>%
  select(Region, site, Camera, org_pct, inorg_pct, salt_pct)

bite_by_quadrat <- bite_by_quadrat %>% 
  left_join(IO_quadrat, by = c("Region", "site", "Camera"))

herb_bite_positive <- herb_bite_positive %>%
  left_join(IO_quadrat, by = c("Region", "site", "Camera"))

#check the relationship between org and total bite rate/event rate/bite per event
cor.test(bite_by_quadrat$total_bite_rate,
         bite_by_quadrat$org_pct,
         method = "spearman")
cor.test(bite_by_quadrat$event_rate,
         bite_by_quadrat$org_pct,
         method = "spearman")
cor.test(herb_bite_positive$mean_bites_per_event,
         herb_bite_positive$org_pct,
         method = "spearman")

IO_site <- IO %>%
  group_by(Region, site) %>%
  summarise(
    org_pct = mean(org_pct, na.rm = T),
    inorg_pct = mean(inorg_pct, na.rm = T),
    salt_pct = mean(salt_pct, na.rm = T),
    .groups = "drop"
  )

IO_site_long <- IO_site %>%
  select(Region, site, org_pct, inorg_pct, salt_pct) %>%
  pivot_longer(
    cols = c(org_pct, inorg_pct, salt_pct), 
    names_to = "component",
    values_to = "percentage"
  )

sediment_color <- c(
  "org_pct" = "#FFBB39",
  "inorg_pct" = "#6FA8DC",
  "salt_pct" = "#999999"
)

IO_pct <- 
  ggplot(IO_site_long, aes(site, percentage, fill = component)) + 
  geom_bar(
    stat = "identity",
    position = "stack",
    color = "black",
    width = 0.8
  ) +
  facet_wrap(~ Region, scales = "free_x") +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  scale_fill_manual(
    values = sediment_color,
    labels = c(
      "org_pct" = "Organic matter",
      "inorg_pct" = "Inorganic matter",
      "salt_pct" = "Salt"
    )) +
  theme_bw()+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    x = "Site",
    y = "Percentage (%)",
    fill = "Sediment component"
  )
IO_pct  

#Bite rate consider bulk_P and algae_P----

P <- P %>%
  rename(sediment_id = sediment)

P_data <- P %>%
  mutate(Camera = gsub("S", "C", sediment_id)) %>%
  select(Region, site, Camera, bulk_P_cmgg, algae_P_cmgg)

bite_by_quadrat <- bite_by_quadrat %>%
  left_join(P_data, by = c("Region", "site", "Camera"))

herb_bite_positive <- herb_bite_positive %>%
  left_join(P_data, by = c("Region", "site", "Camera"))

cor.test(bite_by_quadrat$total_bite_rate,
         bite_by_quadrat$algae_P_cmgg,
         method = "spearman")
cor.test(bite_by_quadrat$event_rate,
         bite_by_quadrat$algae_P_cmgg,
         method = "spearman")
cor.test(herb_bite_positive$mean_bites_per_event,
         herb_bite_positive$algae_P_cmgg,
         method = "spearman")

#Modeling----
total_bite_rate_org_lmer <- 
  lmer(total_bite_rate ~ org_pct + (1|site), data = bite_by_quadrat)
event_rate_org_lmer <- 
  lmer(event_rate ~ org_pct + (1|site), data = bite_by_quadrat)
bite_per_event_org_lmer <- 
  lmer(mean_bites_per_event ~ org_pct + (1|site), data = herb_bite_positive)

total_bite_rate_P_lmer <- 
  lmer(total_bite_rate ~ algae_P_cmgg + (1|site), data = bite_by_quadrat)
event_rate_P_lmer <- 
  lmer(event_rate ~ algae_P_cmgg + (1|site), data = bite_by_quadrat)
bite_per_event_P_lmer <- 
  lmer(mean_bites_per_event ~ algae_P_cmgg + (1|site), data = herb_bite_positive)

#plot
bite_by_quadrat_env <- bite_by_quadrat %>%
  pivot_longer(
    cols = c(total_bite_rate, event_rate, mean_bites_per_event),
    names_to = "response",
    values_to = "response_value"
  ) %>%
  pivot_longer(
    cols = c(org_pct, algae_P_cmgg),
    names_to = "predictor",
    values_to = "predictor_value"
  ) %>%
  mutate(
    response = factor(
      response,
      levels = c("total_bite_rate", "event_rate", "mean_bites_per_event"),
      labels = c("Total bite rate", "Event rate", "Bites per event")
    ),
    predictor = factor(
      predictor,
      levels = c("algae_P_cmgg", "org_pct"),
      labels = c("Algal phosphorus (mg/g)", "Organic matter (%)")
    )
  )

grazing_org_P_label <- bite_by_quadrat_env %>%
  group_by(response, predictor) %>%
  summarise(
    rho = cor(predictor_value, response_value,
              method = "spearman", use = "complete.obs"),
    p = cor.test(predictor_value, response_value,
                 method = "spearman")$p.value,
    .groups = "drop"
  ) %>%
  mutate(label = paste0("ρ = ", round(rho, 2),
                        "\np = ", signif(p, 2)))
  
grazing_org_P <- 
  ggplot(bite_by_quadrat_env, aes(predictor_value, response_value)) +
  geom_point(alpha = 0.8) +
  geom_smooth(method = "lm", linewidth = 0.8) +
  geom_text(
    data = grazing_org_P_label,
    aes(label = label),
    x = Inf,
    y = Inf,
    hjust = 1.1,
    vjust = 1.1
  ) +
  facet_grid(response ~ predictor, scales = "free") +
  labs(
    x = NULL,
    y = NULL,
    title = "Relationship between environmental variables and herbivore grazing"
  )
grazing_org_P

#bite rate per P/org----

#check P and organic correlation
cor.test(
  ~ algae_P_cmgg + org_pct,
  data = bite_by_quadrat,
  method = "spearman"
)

algal_P_org <- 
  ggplot(bite_by_quadrat, aes(algae_P_cmgg, org_pct, colour = Region)) +
  geom_point() +
  scale_colour_manual(values = region_color) +
  geom_smooth(method = "lm", se = F) +
  labs(
    x = "Algal phosphorous (mg/g)",
    y = "Organic matter (%)",
    title = "Relationship between algal P and organic matter"
  ) +
  theme_bw()
algal_P_org

bite_by_quadrat <- bite_by_quadrat %>%
  mutate(
    Region = factor(Region, levels = c("GI", "NE", "XLQ")),
    site = factor(site)
  )

# Region-only models
m_total_region <- lmer(total_bite_rate ~ Region + (1 | site), data = bite_by_quadrat)
m_event_region <- lmer(event_rate ~ Region + (1 | site), data = bite_by_quadrat)
m_bite_region  <- lmer(mean_bites_per_event ~ Region + (1 | site), data = herb_bite_positive)

summary(m_total_region)
anova(m_total_region)

summary(m_event_region)
anova(m_event_region)

summary(m_bite_region)
anova(m_bite_region)

# Region + algal P models
m_total_Padj <- lmer(total_bite_rate ~ Region + algae_P_cmgg + (1 | site), 
                     data = bite_by_quadrat)

m_event_Padj <- lmer(event_rate ~ Region + algae_P_cmgg + (1 | site), 
                     data = bite_by_quadrat)

m_bite_Padj <- lmer(mean_bites_per_event ~ Region + algae_P_cmgg + (1 | site), 
                    data = herb_bite_positive)

summary(m_total_Padj)
anova(m_total_Padj)

summary(m_event_Padj)
anova(m_event_Padj)

summary(m_bite_Padj)
anova(m_bite_Padj)

emm_total_Padj <- emmeans(m_total_Padj, ~ Region)
emm_event_Padj <- emmeans(m_event_Padj, ~ Region)
emm_bite_Padj  <- emmeans(m_bite_Padj, ~ Region)

emm_total_Padj
pairs(emm_total_Padj)

emm_event_Padj
pairs(emm_event_Padj)

emm_bite_Padj
pairs(emm_bite_Padj)

emm_total_df <- as.data.frame(emm_total_Padj) %>%
  mutate(metric = "Total bite rate\n(bites/hr)")

emm_event_df <- as.data.frame(emm_event_Padj) %>%
  mutate(metric = "Feeding event rate\n(events/hr)")

emm_bite_df <- as.data.frame(emm_bite_Padj) %>%
  mutate(metric = "Bites per event\n(bites/event)")

emm_all <- bind_rows(emm_total_df, emm_event_df, emm_bite_df)

grazing_P_adjusted <- 
  ggplot(emm_all, aes(x = Region, y = emmean, colour = Region)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.15, linewidth = 0.7) +
  facet_wrap(~ metric, scales = "free_y", nrow = 1) +
  scale_colour_manual(values = region_color) +
  theme_bw() +
  guides(colour = "none") +
  labs(
    title = "Herbivore Grazing Pressure adjusted for algal P",
    x = "Region",
    y = "Estimated marginal mean"
  )
grazing_P_adjusted

#testing P influence on bite regardless of region----

m_total_P <- lmer(total_bite_rate ~ algae_P_cmgg + (1 | site), 
                     data = bite_by_quadrat)

m_event_P <- lmer(event_rate ~ algae_P_cmgg + (1 | site), 
                     data = bite_by_quadrat)

m_bite_P <- lmer(mean_bites_per_event ~ algae_P_cmgg + (1 | site), 
                    data = herb_bite_positive)

summary(m_total_P)
anova(m_total_P)

summary(m_event_P)
anova(m_event_P)

summary(m_bite_P)
anova(m_bite_P)
