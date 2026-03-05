setwd("C:/Users/dan91/OneDrive/桌面/Field work")
library(readxl)
library(dplyr)
library(forcats)
library(lme4)
library(ggplot2)
library(patchwork)
library(tidyr)

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

# Relationship between sediment weight & Region/environmental parameters----

sed_model1 <- lmer(sed_wg ~ Region + depth + temperature + (1 | site), 
                data = sediment)
summary(sed_model1)

sed_model2 <- lmer(sed_wg ~ Region * depth + Region * temperature + 
                  (1 | site), data = sediment)
summary(sed_model2)

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

bite_by_quadrat <- bite %>%
  filter(`trophic group` == "Herbivore") %>%
  group_by(Region, site, Camera) %>%
  summarise(
    total_bites = sum(bites, na.rm = T),
    n_events = n(),
    mean_bites_per_event = mean(bites, na.rm = T),
    .groups = "drop"
  ) %>%
  mutate(
    total_bite_rate = total_bites * 60 / duration,
    event_rate = n_events * 60 / duration
  )

bite_per_event <- 
  ggplot(bite_by_quadrat, aes(Region, mean_bites_per_event, fill = Region)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = region_color) +
  guides(fill = "none") +
  geom_jitter(width = 0.15, size = 2, alpha = 0.6) +
  stat_summary(fun = mean, geom = "point", color = "red") +
  theme_bw() +
  labs(
    y = expression(Bites~per~event~(bites~events^{-1}))
  )
bite_per_event


total_bite_rate <- 
  ggplot(bite_by_quadrat, aes(Region, total_bite_rate, fill = Region)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = region_color) +
  guides(fill = "none") +
  geom_jitter(width = 0.15, size = 2, alpha = 0.6) +
  stat_summary(fun = mean, geom = "point", color = "red") +
  theme_bw() +
  labs(
    y = expression(Total~bite~rate~(bites~m^{-2}~hr^{-1}))
  )
total_bite_rate


event_rate <- 
  ggplot(bite_by_quadrat, aes(Region, event_rate, fill = Region)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = region_color) +
  guides(fill = "none") +
  geom_jitter(width = 0.15, size = 2, alpha = 0.6) +
  stat_summary(fun = mean, geom = "point", color = "red") +
  theme_bw() +
  labs(
    y = expression(Feeding~events~rate~(events~m^{-2}~hr^{-1}))
  )
event_rate

(total_bite_rate|event_rate|bite_per_event) +
  plot_annotation(
    title = "Herbivore Grazing Pressure"
  )

kruskal.test(total_bite_rate ~ Region, data = bite_by_quadrat)
kruskal.test(event_rate ~ Region, data = bite_by_quadrat)
kruskal.test(mean_bites_per_event ~ Region, data = bite_by_quadrat)

pairwise.wilcox.test(
  bite_by_quadrat$event_rate,
  bite_by_quadrat$Region,
  p.adjust.method = "BH"
)

# 檢查 total bite rate 是由哪個 component 驅動
cor.test(bite_by_quadrat$total_bite_rate,
         bite_by_quadrat$event_rate,
         method = "spearman")

cor.test(bite_by_quadrat$total_bite_rate,
         bite_by_quadrat$mean_bites_per_event,
         method = "spearman")

# Bite rate considered of organic/inorganic percentage----

# 整理 IO data 到 site level
IO_site <- IO %>%
  group_by(Region, site) %>%
  summarise(
    org_pct = mean(org_pct, na.rm = T),
    inorg_pct = mean(inorg_pct, na.rm = T),
    salt_pct = mean(salt_pct, na.rm = T),
    .groups = "drop"
  )

bite_by_quadrat <- bite_by_quadrat %>% 
  left_join(IO_site, by = c("Region", "site"))

#check the relationship between org and total bite rate/event rate/bite per event
cor.test(bite_by_quadrat$total_bite_rate,
         bite_by_quadrat$org_pct,
         method = "spearman")
cor.test(bite_by_quadrat$event_rate,
         bite_by_quadrat$org_pct,
         method = "spearman")
cor.test(bite_by_quadrat$mean_bites_per_event,
         bite_by_quadrat$org_pct,
         method = "spearman")

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

#Bite rate consider bulk_P and algae_P
P_data <- P %>%
  mutate(Camera = gsub("S", "C", sediment)) %>%
  select(Region, site, Camera, bulk_P_cmgg, algae_P_cmgg)

bite_by_quadrat <- bite_by_quadrat %>%
  left_join(P_data, by = c("Region", "site", "Camera"))

cor.test(P_data$bulk_P_cmgg,
         P_data$algae_P_cmgg,
         method = "spearman")

cor.test(bite_by_quadrat$total_bite_rate,
         bite_by_quadrat$algae_P_cmgg,
         method = "spearman")
cor.test(bite_by_quadrat$event_rate,
         bite_by_quadrat$algae_P_cmgg,
         method = "spearman")
cor.test(bite_by_quadrat$mean_bites_per_event,
         bite_by_quadrat$algae_P_cmgg,
         method = "spearman")

event_rate_P <- 
  ggplot(bite_by_quadrat, aes(algae_P_cmgg, event_rate)) +
  geom_point() +
  geom_smooth(method = "lm")
event_rate_P

lmer(event_rate ~ algae_P_cmgg + (1|site), data=bite_by_quadrat)

#event rate ~ P (each region)
event_rate_P_region <- 
  ggplot(bite_by_quadrat, aes(algae_P_cmgg, event_rate, color = Region)) +
  geom_point(size = 2) +
  scale_color_manual(values = region_color) +
  geom_smooth(method = "lm", se = F)
event_rate_P_region
