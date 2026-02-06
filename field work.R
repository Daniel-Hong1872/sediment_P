setwd("C:/Users/dan91/OneDrive/桌面/Field work")
library(readxl)
library(dplyr)
library(forcats)
library(lme4)
library(ggplot2)
library(patchwork)

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

sed_site_mean <- sediment %>% 
  group_by(Region, site) %>%
  summarise(
    depth = mean(depth, na.rm = T),
    temperature = mean(temperature, na.rm = T),
    sediment_weight = mean(sed_wg, na.rm = T),
    .groups = "drop"
  )

sed_depth <- ggplot(sed_site_mean, aes(x = depth, y = sediment_weight)) + 
  geom_point() +
  geom_smooth(method = "lm", se = F) +
  facet_wrap(~Region)

sed_temp <- ggplot(sed_site_mean, aes(x = temperature, y = sediment_weight)) + 
  geom_point() +
  geom_smooth(method = "lm", se = F) +
  facet_wrap(~Region)

sed_region <- ggplot(sed_site_mean, aes(x = Region, y = sediment_weight)) +
  geom_boxplot() +
  geom_jitter()

sed_depth
sed_temp
sed_region

# Relationship between P concentration & Region----

bulk_P_model1 <- lmer(bulk_P_cmgg ~ Region + (1 | site), data = P)
summary(bulk_P_model1)

algae_P_model1 <- lmer(algae_P_cmgg ~ Region + (1 | site), data = P)
summary(algae_P_model1)

P_site_mean <- P %>% 
  group_by(Region, site) %>%
  summarise(
    bulk_P = mean(bulk_P_cmgg, na.rm = T),
    algae_P = mean(algae_P_cmgg, na.rm = T),
    .groups = "drop"
  )

bulk_P_region <- ggplot(P_site_mean, aes(x = Region, y = bulk_P)) +
  geom_boxplot() +
  geom_jitter()

algae_P_region <- ggplot(P_site_mean, aes(x = Region, y = algae_P)) +
  geom_boxplot() +
  geom_jitter()

bulk_P_region
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
  group_by(Region, Site, Camera, `trophic group`) %>%
  summarise(
    total_bites = sum(bites, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    bite_rate = total_bites * 60 / duration   # m-2 hr-1
  )

trophic_color <- c(
  "Herbivore" = "#7FDEA0",
  "Benthic invertivore" = "#47A0BD",
  "Corallivore" = "#FF6D69",
  "Omnivore" = "#AF84F3",
  "NA" = "#D0D2D4"
)

trophic_region_p1 <- ggplot(
  bite_quad_trophic, 
  aes(x = `trophic group`, y = bite_rate, fill = `trophic group`)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.6) +
  stat_summary(fun = mean, geom = "point", color = "red", size = 2.5) +
  facet_wrap(~ Region, nrow = 1) +
  scale_fill_manual(values = trophic_color, drop = FALSE) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  labs(
    x = "Trophic group",
    y = expression(Bite~rate~(bites~m^{-2}~hr^{-1}))
  )
trophic_region_p1

bite_quad_func <- bite %>%
  filter(`trophic group` == "Herbivore") %>%
  group_by(Region, Site, Camera, `herbivore functional group`) %>%
  summarise(
    total_bites = sum(bites, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    bite_rate = total_bites * 60 / duration
  )

bite_herb <- bite %>%
  filter(`trophic group` == "Herbivore") %>%
  rename(site = Site)

region_color <- c(
  "GI" = "lightgreen",
  "NE" = "skyblue1",
  "XLQ" = "indianred1"
)

bite_rate_indiv_p3 <- ggplot(bite_herb, aes(Region, bite_rate)) +
  geom_boxplot(outlier.shape = NA, fill = region_color) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.6) +
  stat_summary(fun = mean, geom = "point", color = "red", size = 3) +
  theme_bw()
bite_rate_indiv_p3

func_group_color <- c(
  "browsers" = "#8A6646",
  "farmers (territorial cropper)" = "#E499C4",
  "grazers (croppers)" = "#7CCE7C",
  "scrapers" = "#6CA6CC",
  "small excavators" = "#D89D65"
)

herb_func_region_p4 <- ggplot(
  bite_quad_func,
  aes(x = `herbivore functional group`, y = bite_rate,
      fill = `herbivore functional group`)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.6) +
  stat_summary(fun = mean, geom = "point", color = "red") +
  facet_wrap(~ Region, nrow = 1) +
  scale_fill_manual(values = func_group_color, drop = FALSE) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none") +
  labs(
    x = "Herbivore functional group",
    y = expression(Bite~rate~(bites~m^{-2}~hr^{-1}))
  )
herb_func_region_p4

bite_herb_site <- bite_herb %>%
  group_by(Region, site) %>%
  summarise(
    bite_rate = mean(bite_rate, na.rm = T),
    .groups = "drop"
  )
kruskal.test(bite_rate ~ Region, data = bite_herb_site)
pairwise.wilcox.test(
  bite_herb_site$bite_rate,
  bite_herb_site$Region,
  p.adjust.method = "BH"
)

# Bite rate considered of organic/inorganic percentage----

IO_site <- IO %>%
  group_by(Region, site) %>%
  summarise(
    org_pct = mean(org_pct, na.rm = T),
    inorg_pct = mean(inorg_pct, na.rm = T),
    .groups = "drop"
  )

bite_IO <- bite_herb_site %>%
  left_join(IO_site, by = c("Region", "site"))

bite_region_IO <- lm(
  bite_rate ~ Region + org_pct,
  data = bite_IO
)
summary(bite_region_IO)

bite_IO$adj_bite <- resid(bite_region_IO)

p5 <- ggplot(bite_IO, aes(x = Region, y = adj_bite)) +
  geom_boxplot() +
  geom_jitter(width = 0.1) +
  theme_bw() +
  labs(
    y = "Herbivore bite rate (residuals)"
  )
p5

cor.test(
  bite_IO$bite_rate,
  bite_IO$org_pct,
  method = "spearman"
)
cor.test(
  bite_IO$bite_rate,
  bite_IO$inorg_pct,
  method = "spearman"
)

p6 <- ggplot(bite_IO, aes(x = org_pct, y = bite_rate)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = F) +
  facet_wrap(~Region)+
  theme_bw()+
  labs(
    x = "Organic matter (%)",
    y = expression(Herbivore~bite~rate~(m^{-2}~hr^{-1}))
  )
p6

p7 <- ggplot(bite_IO, aes(x = inorg_pct, y = bite_rate)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = F) +
  facet_wrap(~Region)+
  theme_bw()+
  labs(
    x = "Inorganic matter (%)",
    y = expression(Herbivore~bite~rate~(m^{-2}~hr^{-1}))
  )
p7

# Bite rate by individual----

bite_by_quadrat <- bite %>%
  filter(`trophic group` == "Herbivore") %>%
  group_by(Region, Site, Camera) %>%
  summarise(
    total_bites = sum(bites, na.rm = TRUE),
    n_events = n(),
    mean_bites_per_event = mean(bites, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    total_bite_rate = total_bites * 60 / duration,
    event_rate = n_events * 60 / duration
  )

bite_per_event_p8 <- ggplot(bite_by_quadrat, aes(Region, mean_bites_per_event)) +
  geom_boxplot(outlier.shape = NA, fill = region_color) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.6) +
  stat_summary(fun = mean, geom = "point", color = "red", size = 3) +
  theme_bw() +
  labs(
    y = expression(Bites~per~event~(bites~events^{-1}))
  )
bite_per_event_p8


total_bite_rate_p9 <- ggplot(bite_by_quadrat, aes(Region, total_bite_rate)) +
  geom_boxplot(outlier.shape = NA, fill = region_color) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.6) +
  stat_summary(fun = mean, geom = "point", color = "red", size = 3) +
  theme_bw() +
  labs(
    y = expression(Total~bite~rate~(bites~m^{-2}~hr^{-1}))
  )
total_bite_rate_p9


event_rate_p10 <- ggplot(bite_by_quadrat, aes(Region, event_rate)) +
  geom_boxplot(outlier.shape = NA, fill = region_color) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.6) +
  stat_summary(fun = mean, geom = "point", color = "red", size = 3) +
  theme_bw() +
  labs(
    y = expression(Feeding~events~rate~(events~m^{-2}~hr^{-1}))
  )
event_rate_p10

(total_bite_rate_p9|event_rate_p10|bite_per_event_p8) +
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
