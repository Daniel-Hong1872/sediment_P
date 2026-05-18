# Field work analysis

setwd("C:/Users/dan91/OneDrive/桌面/Field work")
library(tidyverse)
library(readxl)
library(lme4)
library(lmerTest)
library(emmeans)
library(patchwork)
library(performance)
library(glmmTMB)
library(DHARMa)
library(coin)
library(sf)
library(terra)
library(geodata)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggspatial)
library(ggrepel)
library(ragg)
library(cowplot)

# Import data----
sed <- read_excel("CNP_data.xlsx", sheet = "sediment", na = "NA")
P <- read_excel("CNP_data.xlsx", sheet = "P", na = "NA")
IO <- read_excel("CNP_data.xlsx", sheet = "Inorganic vs. Organic", na = "NA")

# Sediment weight across regions----

sed <- sed %>% 
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
    sediment_id = sediment,
    bulk_P_mgL = `sample A_P concentration(mg/L)`,
    algal_P_mgL = `sample B_P concentration(mg/L)`,
    bulk_P_mg = `A subsample_P weight(mg)`,
    algal_P_mg = `B subsample_P weight(mg)`
  ) %>%
  mutate(
    bulk_P_mgg = (bulk_P_mgL * V_ml / 1000) / bulk_P_mg * 1000,
    algal_P_mgg = (algal_P_mgL * V_ml / 1000)/ algal_P_mg * 1000
  )

IO <- IO %>%
  rename(
    sediment_id = sediment,
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

# setup region color
region_color <- c(
  "GI" = "#4CAF50",
  "NE" = "#3288BD",
  "XLQ" = "#D53E4F"
)

sed_site_mean <- sed %>% 
  group_by(Region, site) %>%
  summarise(
    depth = mean(depth, na.rm = T),
    temperature = mean(temperature, na.rm = T),
    sediment_weight = mean(sed_wg, na.rm = T),
    .groups = "drop"
  )

## sediment weight vs. depth----
sed_depth <- 
  ggplot(sed_site_mean, aes(x = depth, y = sediment_weight)) + 
  geom_point() +
  geom_smooth(method = "lm", se = F) +
  facet_wrap(~Region) +
  theme_bw()
sed_depth

## sediment weight vs. temperature----
sed_temp <- 
  ggplot(sed_site_mean, aes(x = temperature, y = sediment_weight)) + 
  geom_point() +
  geom_smooth(method = "lm", se = F) +
  facet_wrap(~Region) +
  theme_bw()
sed_temp

## sediment weight vs. region----
sed_region <- 
  ggplot(sed_site_mean, aes(Region, sediment_weight, fill = Region)) +
  geom_boxplot() +
  scale_fill_manual(values = region_color) +
  guides(fill = "none") +
  geom_jitter() +
  theme_bw()
sed_region


# Sea urchin density and composition----

XLQ_ur <- read_excel("urchin density.xlsx", sheet = "XLQ") %>% 
  mutate(Region = "XLQ")

GI_ur <- read_excel("urchin density.xlsx", sheet = "GI") %>%
  mutate(Region = "GI")

NE_ur <- read_excel("urchin density.xlsx", sheet = "NE") %>%
  mutate(Region = "NE")

urchin_all <- bind_rows(XLQ_ur, GI_ur, NE_ur)

# only select herbivore sea urchin that can move and graze
herb_urchin <- urchin_all %>%
  filter(Genus %in% c("Diadema", "Echinothrix", "Stomopneustes"))

herb_urchin_region <- herb_urchin %>%
  group_by(Region, Genus) %>%
  summarise(Total = n(), .groups = "drop")

urchin_color <- c(
  "Diadema" = "#F2A541",
  "Echinothrix" = "#C97B84",
  "Stomopneustes" = "#6C6F9A"
)

herb_urchin_sum_region <- 
  ggplot(herb_urchin_region, aes(Region, Total, fill = Genus)) +
  geom_col() +
  scale_fill_manual(
    values = urchin_color,
    breaks = c("Diadema", "Echinothrix", "Stomopneustes"),
    labels = c(
      expression(italic(Diadema)),
      expression(italic(Echinothrix)),
      expression(italic(Stomopneustes))
    ),
    drop = FALSE
  ) +
  theme_bw() +
  labs(x = "Region", y = "Total number of urchins", fill = "Genus")

herb_urchin_sum_region

all_samples <- urchin_all %>%
  distinct(Region, Site)

urchin_density <- herb_urchin %>%
  count(Region, Site, name = "count") %>%
  right_join(all_samples, by = c("Region", "Site")) %>%
  mutate(
    count = replace_na(count, 0),
    density = count / 10
  ) %>%
  arrange(Region, Site)
urchin_density

urchin_summary <- urchin_density %>%
  group_by(Region) %>%
  summarise(
    mean_density = mean(density, na.rm = TRUE),
    sd_density = sd(density, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )
urchin_summary

urchin_density <- urchin_density %>%
  mutate(
    Region = as.factor(Region),
    density = as.numeric(density)
  )

kruskal_test(
  density ~ Region, 
  data = urchin_density,
  distribution = coin::approximate(nresample = 9999)
)

pairwise.wilcox.test(
  x = urchin_density$density,
  g = urchin_density$Region,
  p.adjust.method = "BH"
)

urchin_plot <- 
  ggplot(urchin_density, aes(x = Region, y = density, color = Region)) +
  stat_summary(
    fun.data = mean_sdl,
    fun.args = list(mult = 1),
    geom = "errorbar",
    width = 0.12
  ) +
  stat_summary(
    fun = mean,
    geom = "point",
    size = 4
  ) +
  scale_color_manual(values = region_color) +
  guides(color = "none") +
  theme_bw() +
  labs(
    x = NULL,
    y = expression("Density (ind m"^{-2}*")")
  )
urchin_plot

urchin_plot | herb_urchin_sum_region

# Nutrient concentration and stoichiometric ratios----

P_site_mean <- P %>% 
  group_by(Region, site) %>%
  summarise(
    bulk_P = mean(bulk_P_mgg, na.rm = TRUE),
    algal_P = mean(algal_P_mgg, na.rm = TRUE),
    .groups = "drop"
  )

CN <- read_excel("CNP_data.xlsx", sheet = "CN", na = "NA")

CN <- CN %>%
  rename(
    sediment_id = sediment,
    bulk_mg = `sample_a_weight_mg`,
    bulk_N_pct = `sample_a_N%`,
    bulk_C_pct = `sample_a_C%`,
    algal_mg = `sample_b_weight_mg`,
    algal_N_pct = `sample_b_N%`,
    algal_C_pct = `sample_b_C%`
  ) %>%
  mutate(
    # concentration (mg/g)
    bulk_N_mgg  = bulk_N_pct * 10,
    bulk_C_mgg  = bulk_C_pct * 10,
    algal_N_mgg = algal_N_pct * 10,
    algal_C_mgg = algal_C_pct * 10,
    
    # molar concentration (mmol/g)
    bulk_C_mmolg  = bulk_C_mgg / 12,
    bulk_N_mmolg  = bulk_N_mgg / 14,
    algal_C_mmolg = algal_C_mgg / 12,
    algal_N_mmolg = algal_N_mgg / 14
  )

## CN adjusted back to bulk basis ----

P_IO_adj <- P %>%
  left_join(
    IO %>% select(Region, site, sediment_id, inorg_pct, org_pct),
    by = c("Region", "site", "sediment_id")
  ) %>%
  mutate(
    bulk_P_mmolg = bulk_P_mgg / 31,
    algal_P_mmolg = algal_P_mgg / 31
  )

CNP_adj_2 <- CN %>%
  left_join(
    P_IO_adj %>% select(
      Region, site, sediment_id, org_pct,
      algal_P_mgg, algal_P_mmolg,
      bulk_P_mgg, bulk_P_mmolg
    ),
    by = c("Region", "site", "sediment_id")
  ) %>%
  mutate(
    Camera = gsub("S", "C", sediment_id),
    bulk_C_mgg_adj = bulk_C_mgg * org_pct / 100,
    bulk_N_mgg_adj = bulk_N_mgg * org_pct / 100,
    bulk_C_mmolg_adj = bulk_C_mmolg * org_pct / 100,
    bulk_N_mmolg_adj = bulk_N_mmolg * org_pct / 100
  )

CNP_site_mean_adj_2 <- CNP_adj_2 %>%
  group_by(Region, site) %>%
  summarise(
    across(
      c(
        bulk_C_mgg_adj, bulk_N_mgg_adj, bulk_P_mgg,
        algal_C_mgg, algal_N_mgg, algal_P_mgg,
        bulk_C_mmolg_adj, bulk_N_mmolg_adj, bulk_P_mmolg,
        algal_C_mmolg, algal_N_mmolg, algal_P_mmolg
      ),
      ~ mean(.x, na.rm = TRUE),
      .names = "{.col}_mean"
    ),
    .groups = "drop"
  ) %>%
  mutate(
    bulk_CN_ratio_adj = bulk_C_mmolg_adj_mean / bulk_N_mmolg_adj_mean,
    algal_CN_ratio    = algal_C_mmolg_mean / algal_N_mmolg_mean,
    bulk_CP_ratio_adj = bulk_C_mmolg_adj_mean / bulk_P_mmolg_mean,
    algal_CP_ratio    = algal_C_mmolg_mean / algal_P_mmolg_mean,
    bulk_NP_ratio_adj = bulk_N_mmolg_adj_mean / bulk_P_mmolg_mean,
    algal_NP_ratio    = algal_N_mmolg_mean / algal_P_mmolg_mean
  )

CNP_conc_long_adj_2 <- CNP_site_mean_adj_2 %>%
  pivot_longer(
    cols = c(
      bulk_C_mgg_adj_mean, bulk_N_mgg_adj_mean, bulk_P_mgg_mean,
      algal_C_mgg_mean, algal_N_mgg_mean, algal_P_mgg_mean
    ),
    names_to = "type",
    values_to = "value"
  )

bulk_CNP_ratio_table_adj_2 <- CNP_site_mean_adj_2 %>%
  group_by(Region) %>%
  summarise(
    bulk_C = mean(bulk_C_mmolg_adj_mean, na.rm = TRUE),
    bulk_N = mean(bulk_N_mmolg_adj_mean, na.rm = TRUE),
    bulk_P = mean(bulk_P_mmolg_mean, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    C_rel = bulk_C / bulk_P,
    N_rel = bulk_N / bulk_P,
    P_rel = 1,
    C_int = round(C_rel),
    N_int = round(N_rel),
    P_int = 1
  ) %>%
  transmute(
    Region,
    `C:N:P` = paste0(C_int, " : ", N_int, " : ", P_int)
  )
bulk_CNP_ratio_table_adj_2

CNP_site_mean_adj_2 <- CNP_site_mean_adj_2 %>%
  mutate(
    Region = as.factor(Region),
    bulk_C_mgg_adj_mean = as.numeric(bulk_C_mgg_adj_mean),
    bulk_N_mgg_adj_mean = as.numeric(bulk_N_mgg_adj_mean),
    bulk_P_mgg_mean = as.numeric(bulk_P_mgg_mean),
    bulk_CN_ratio_adj = as.numeric(bulk_CN_ratio_adj),
    bulk_CP_ratio_adj = as.numeric(bulk_CP_ratio_adj),
    bulk_NP_ratio_adj = as.numeric(bulk_NP_ratio_adj),
  )

kruskal_test(bulk_C_mgg_adj_mean ~ Region, 
             data = CNP_site_mean_adj_2,
             distribution = coin::approximate(nresample = 9999)
)

pairwise.wilcox.test(
  CNP_site_mean_adj_2$bulk_C_mgg_adj_mean,
  CNP_site_mean_adj_2$Region,
  p.adjust.method = "BH"
)

kruskal_test(bulk_N_mgg_adj_mean ~ Region, 
             data = CNP_site_mean_adj_2,
             distribution = coin::approximate(nresample = 9999)
)

pairwise.wilcox.test(
  CNP_site_mean_adj_2$bulk_N_mgg_adj_mean,
  CNP_site_mean_adj_2$Region,
  p.adjust.method = "BH"
)

kruskal_test(bulk_P_mgg_mean ~ Region, 
             data = CNP_site_mean_adj_2,
             distribution = coin::approximate(nresample = 9999)
)

pairwise.wilcox.test(
  CNP_site_mean_adj_2$bulk_P_mgg_mean,
  CNP_site_mean_adj_2$Region,
  p.adjust.method = "BH"
)

kruskal_test(bulk_CN_ratio_adj ~ Region, 
             data = CNP_site_mean_adj_2,
             distribution = coin::approximate(nresample = 9999)
)

kruskal_test(bulk_CP_ratio_adj ~ Region, 
             data = CNP_site_mean_adj_2,
             distribution = coin::approximate(nresample = 9999)
)

kruskal_test(bulk_NP_ratio_adj ~ Region, 
             data = CNP_site_mean_adj_2,
             distribution = coin::approximate(nresample = 9999)
)

## bulk only plot----

p1_data <- CNP_site_mean_adj_2 %>%
  select(
    Region,
    bulk_C_mgg_adj_mean,
    bulk_N_mgg_adj_mean,
    bulk_P_mgg_mean
  ) %>%
  pivot_longer(
    cols = -Region,
    names_to = "variable",
    values_to = "value"
  )

p1 <- ggplot(p1_data, aes(x = Region, y = value, fill = Region)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.8) +
  scale_fill_manual(values = region_color) +
  guides(fill = "none") +
  facet_wrap(~ variable, scales = "free_y", ncol = 3,
             labeller = as_labeller(c(
               "bulk_C_mgg_adj_mean" = "Carbon (mg g⁻¹)",
               "bulk_N_mgg_adj_mean" = "Nitrogen (mg g⁻¹)",
               "bulk_P_mgg_mean"     = "Phosphorus (mg g⁻¹)"
             ))) +
  theme_bw() +
  labs(
    x = NULL,
    y = "Concentration (mg/g)"
  )

p2_data <- CNP_site_mean_adj_2 %>%
  select(
    Region,
    bulk_CN_ratio_adj,
    bulk_CP_ratio_adj,
    bulk_NP_ratio_adj
  ) %>%
  pivot_longer(
    cols = -Region,
    names_to = "variable",
    values_to = "value"
  )

p2 <- ggplot(p2_data, aes(x = Region, y = value, fill = Region)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.8) +
  scale_fill_manual(values = region_color) +
  guides(fill = "none") +
  facet_wrap(~ variable, scales = "free_y", ncol = 3,
             labeller = as_labeller(c(
               "bulk_CN_ratio_adj" = "C:N ratio",
               "bulk_CP_ratio_adj" = "C:P ratio",
               "bulk_NP_ratio_adj" = "N:P ratio"
             ))) +
  theme_bw() +
  labs(
    x = "Region",
    y = "Molar ratio"
  )
p1/p2

# Organic, inorganic, and salt fractions across regions----

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

IO_site <- IO_site %>%
  mutate(
    Region = as.factor(Region),
    org_pct = as.numeric(org_pct)
  )

kruskal_test(org_pct ~ Region, 
             data = IO_site,
             distribution = coin::approximate(nresample = 9999)
)

pairwise.wilcox.test(
  IO_site$org_pct,
  IO_site$Region,
  p.adjust.method = "BH"
)

# Fish bite data----

bite <- read_excel("Fish bite_data.xlsx")

area <- 1 # m^2
duration <- 35 # 35 mins
bite <- bite %>%
  mutate(bite_rate = bites / area * 60 / duration, # bites/hr*m^2
         `trophic group` = fct_relevel(
           `trophic group`,
           "Herbivore", "Benthic invertivore", "Corallivore", "Omnivore"
         ))

# Bite rate by trophic group----
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
  geom_jitter(width = 0.15, alpha = 0.6) +
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

# Bite rate by herbivore functional group----

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

# Quadrat-level grazing pressure----
# Retained for total bite rate and event rate only

bite_by_quadrat <- bite %>%
  filter(`trophic group` == "Herbivore") %>%
  group_by(Region, site, Camera) %>%
  summarise(
    total_bites = sum(bites, na.rm = T),
    n_events = sum(bites > 0, na.rm = T),
    mean_bites_per_event = ifelse(n_events > 0, total_bites / n_events, NA_real_),
    .groups = "drop"
  ) 

all_camera <- bite %>%
  distinct(Region, site) %>%
  tidyr::crossing(Camera = c("C1", "C2", "C3"))

bite_by_quadrat <- all_camera %>%
  left_join(bite_by_quadrat, by = c("Region", "site", "Camera"))

bite_by_quadrat <- bite_by_quadrat %>%
  mutate(
    total_bites = replace_na(total_bites, 0),
    n_events = replace_na(n_events, 0),
    total_bite_rate = total_bites * 60 / duration,
    event_rate = n_events * 60 / duration
  )

# Occurrence of herbivore feeding----

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

# Herbivore grazing pressure: total bite rate and event rate

bite_herb_all <- bite_by_quadrat %>%
  select(Region, total_bite_rate, event_rate) %>%
  pivot_longer(
    cols = c(total_bite_rate, event_rate),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(metric = factor(
    metric, levels = c("total_bite_rate", "event_rate"), 
    labels = c("Total bite rate\n(bites/hr)",
               "Feeding event rate\n(events/hr)")
  ))

herb_grazing_pressure <- 
  ggplot(bite_herb_all, aes(Region, value, fill = Region)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.6) +
  stat_summary(fun = mean, geom = "point", color = "red") +
  scale_fill_manual(values = region_color) +
  facet_wrap(~metric, scales = "free_y", nrow = 1) +
  guides(fill = "none") +
  theme_bw() +
  labs(
    x = "Region",
    y = NULL
  )
herb_grazing_pressure

bite_by_quadrat <- bite_by_quadrat %>%
  mutate(
    Region = as.factor(Region),
    total_bite_rate = as.numeric(total_bite_rate),
    event_rate = as.numeric(event_rate)
  )

kruskal_test(total_bite_rate ~ Region, 
             data = bite_by_quadrat,
             distribution = coin::approximate(nresample = 9999)
)

pairwise.wilcox.test(
  bite_by_quadrat$total_bite_rate,
  bite_by_quadrat$Region,
  p.adjust.method = "BH"
)

kruskal_test(event_rate ~ Region, 
             data = bite_by_quadrat,
             distribution = coin::approximate(nresample = 9999)
)

pairwise.wilcox.test(
  bite_by_quadrat$event_rate,
  bite_by_quadrat$Region,
  p.adjust.method = "BH"
)

# exclude zero bite data for bite per event
herb_bite_positive <- bite_by_quadrat %>%
  filter(n_events > 0)

bites_per_event <- herb_bite_positive %>%
  select(Region, mean_bites_per_event) %>%
  pivot_longer(
    cols = c(mean_bites_per_event),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(metric = factor(
    metric, levels = c("mean_bites_per_event"), 
    labels = c("Bites per event\n(bites/event)")
  ))

bites_per_event_plot <- 
  ggplot(bites_per_event, aes(Region, value, color = Region)) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
  stat_summary(fun = mean, geom = "point", color = "black") +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width = 0.12) +
  scale_color_manual(values = region_color) +
  guides(fill = "none") +
  theme_bw() +
  labs(
    title = "Bites per event",
    x = "Region",
    y = NULL
  )
bites_per_event_plot

herb_bite_positive <- herb_bite_positive %>%
  mutate(
    Region = as.factor(Region),
    mean_bites_per_event = as.numeric(mean_bites_per_event)
  )

kruskal_test(mean_bites_per_event ~ Region, 
             data = herb_bite_positive,
             distribution = coin::approximate(nresample = 9999)
)


# Site-level organic matter data for BPE analysis----

IO_site_mean <- IO %>%
  group_by(Region, site) %>%
  summarise(
    org_pct_mean   = mean(org_pct, na.rm = TRUE),
    inorg_pct_mean = mean(inorg_pct, na.rm = TRUE),
    salt_pct_mean  = mean(salt_pct, na.rm = TRUE),
    org_g_mean     = mean(org_g, na.rm = TRUE),
    inorg_g_mean   = mean(inorg_g, na.rm = TRUE),
    salt_g_mean    = mean(salt_g, na.rm = TRUE),
    .groups = "drop"
  )

# Site-level bites per event
# BPE is retained only at site level in this cleaned version

bpe_site <- bite_by_quadrat %>%
  filter(n_events > 0) %>%
  group_by(Region, site) %>%
  summarise(bpe_site_mean = mean(mean_bites_per_event, na.rm = T), .groups = "drop") %>%
  left_join(
    IO_site_mean, 
    by = c("Region", "site")) %>%
  left_join(
    CNP_site_mean_adj_2, 
    by = c("Region", "site")) %>%
  mutate(
    bites_per_event_org_adj = bpe_site_mean / org_g_mean,
    bites_per_event_CN_adj = bpe_site_mean / bulk_CN_ratio_adj,
    bites_per_event_CP_adj = bpe_site_mean / bulk_CP_ratio_adj,
    bites_per_event_NP_adj = bpe_site_mean / bulk_NP_ratio_adj
  )

# Statistical test for site-level BPE among regions
bpe_site <- bpe_site %>%
  mutate(
    Region = as.factor(Region),
    bpe_site_mean = as.numeric(bpe_site_mean)
  )

kruskal_test(
  bpe_site_mean ~ Region,
  data = bpe_site,
  distribution = coin::approximate(nresample = 9999)
)

pairwise.wilcox.test(
  bpe_site$bpe_site_mean,
  bpe_site$Region,
  p.adjust.method = "BH"
)


bpe_site_adj <- bpe_site %>%
  select(
    Region,
    bites_per_event_org_adj,
    bites_per_event_CN_adj,
    bites_per_event_CP_adj,
    bites_per_event_NP_adj
  ) %>%
  pivot_longer(
    cols = c(
      bites_per_event_org_adj,
      bites_per_event_CN_adj,
      bites_per_event_CP_adj,
      bites_per_event_NP_adj
    ),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = factor(
      metric,
      levels = c(
        "bites_per_event_org_adj",
        "bites_per_event_CN_adj",
        "bites_per_event_CP_adj",
        "bites_per_event_NP_adj"
      ),
      labels = c(
        "adjusted by organic matter",
        "adjusted by C:N",
        "adjusted by C:P",
        "adjusted by N:P"
      )
    )
  )

NE_bpe_site <- bpe_site %>%
  filter(Region == "NE")

cor.test(NE_bpe_site$bpe_site_mean,
         NE_bpe_site$org_pct_mean,
         method = "spearman")

cor.test(NE_bpe_site$bpe_site_mean,
         NE_bpe_site$bulk_CN_ratio_adj,
         method = "spearman")

cor.test(NE_bpe_site$bpe_site_mean,
         NE_bpe_site$bulk_CP_ratio_adj,
         method = "spearman")

cor.test(NE_bpe_site$bpe_site_mean,
         NE_bpe_site$bulk_NP_ratio_adj,
         method = "spearman")

# Site-level BPE adjusted by nutrient quality----

bpe_site_adj_plot <- 
  ggplot(bpe_site_adj, aes(Region, value, color = Region)) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
  stat_summary(fun = mean, geom = "point", color = "black") +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width = 0.12) +
  scale_color_manual(values = region_color) +
  facet_wrap(~metric, scales = "free_y", nrow = 1) +
  guides(fill = "none") +
  theme_bw() +
  labs(
    title = "Bites per event site mean adjusted by nutrient quality",
    x = "Region",
    y = NULL
  )
bpe_site_adj_plot

bpe_site$log_bpe_site_mean <- log(bpe_site$bpe_site_mean)

log_bpe_site_org <- bpe_site %>%
  select(log_bpe_site_mean, org_pct_mean, Region)

# Site-level log BPE vs stoichiometric ratios----

log_bpe_site_CNP_ratio <- bpe_site %>%
  select(
    log_bpe_site_mean,
    org_pct_mean,
    bulk_CN_ratio_adj,
    bulk_CP_ratio_adj,
    bulk_NP_ratio_adj,
    Region
  ) %>%
  pivot_longer(
    cols = c(
      org_pct_mean,
      bulk_CN_ratio_adj,
      bulk_CP_ratio_adj,
      bulk_NP_ratio_adj
    ),
    names_to = "Ratio",
    values_to = "Value"
  ) %>%
  mutate(
    Ratio = recode(
      Ratio,
      "org_pct_mean" = "organic matter",
      "bulk_CN_ratio_adj" = "C:N ratio",
      "bulk_CP_ratio_adj" = "C:P ratio",
      "bulk_NP_ratio_adj" = "N:P ratio"
    ),
    Ratio = factor(
      Ratio,
      levels = c(
        "organic matter",
        "C:N ratio",
        "C:P ratio",
        "N:P ratio"
      )
    )
  )


log_bpe_site_CNP_ratio_plot <- 
  ggplot(
    log_bpe_site_CNP_ratio,
    aes(x = Value, y = log_bpe_site_mean, color = Region)
  ) +
  geom_point(
    aes(size = Region)
  ) +
  facet_wrap(~ Ratio, scales = "free_x", ncol = 2) +
  scale_color_manual(values = region_color) +
  scale_size_manual(
    values = c("GI" = 1,
               "XLQ" = 1,
               "NE" = 3.5)
  ) +
  labs(
    y = "site mean bite per event log-transformed",
    color = "Region",
    size = "Region"
  ) +
  theme_bw()
log_bpe_site_CNP_ratio_plot

# Sampling site map----

sites <- read_excel("sites.xlsx")

sites_sf <- st_as_sf(
  sites,
  coords = c("lon", "lat"),
  crs = 4326,
  remove = F
)

tw_county <- st_read("TW map/COUNTY_MOI_1140318.shp")
tw_county <- st_transform(tw_county, 4326)

taiwan <- tw_county %>%
  summarise(geometry = st_union(geometry))

extent_TW <- list(
  xlim = c(119.8, 122.1),
  ylim = c(21.7, 25.5)
)

extent_NE <- list(
  xlim = c(121.77, 122.02),
  ylim = c(24.99, 25.17)
)

extent_GI <- list(
  xlim = c(121.45, 121.53),
  ylim = c(22.62, 22.69)
)

extent_XLQ <- list(
  xlim = c(120.34, 120.40),
  ylim = c(22.31, 22.37)
)

boxes <- tibble(
  region = c("NE", "GI", "XLQ"),
  xmin = c(extent_NE$xlim[1], extent_GI$xlim[1], extent_XLQ$xlim[1]),
  xmax = c(extent_NE$xlim[2], extent_GI$xlim[2], extent_XLQ$xlim[2]),
  ymin = c(extent_NE$ylim[1], extent_GI$ylim[1], extent_XLQ$ylim[1]),
  ymax = c(extent_NE$ylim[2], extent_GI$ylim[2], extent_XLQ$ylim[2])
)

p_main <- ggplot() +
  geom_sf(
    data = taiwan,
    fill = "olivedrab3",
    color = NA
  ) +
  geom_rect(
    data = boxes,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = NA,
    color = "black",
    linewidth = 0.6
  ) +
  annotate(
    "text",
    x = 121.57, y = 25.25,
    label = "Northeast\n(NE)",
    size = 5,
    fontface = "bold"
  ) +
  annotate(
    "text",
    x = 121.5, y = 22.45,
    label = "Green Island\n(GI)",
    size = 5,
    fontface = "bold"
  ) +
  annotate(
    "text",
    x = 120.30, y = 22.15,
    label = "Xiaoliuqiu\n(XLQ)",
    size = 5,
    fontface = "bold"
  ) +
  annotation_scale(
    location = "bl",
    width_hint = 0.35,
    text_cex = 0.8
  ) +
  annotation_north_arrow(
    location = "tl",
    which_north = "true",
    style = north_arrow_fancy_orienteering,
    height = unit(1.3, "cm"),
    width = unit(1.3, "cm"),
    pad_x = unit(0.5, "cm"),
    pad_y = unit(0.5, "cm")
  ) +
  coord_sf(
    xlim = extent_TW$xlim,
    ylim = extent_TW$ylim,
    expand = F
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.major = element_line(color = "gray90", linewidth = 0.4),
    panel.background = element_rect(fill = "white", color = NA),
    axis.title = element_text(size = 15, color = "black"),
    axis.text = element_text(size = 11, color = "black")
  )

make_inset_map <- function(region_name, xlim, ylim, scale_width = 0.35) {
  
  site_sub <- sites %>%
    filter(region == region_name)
  
  site_sub_sf <- sites_sf %>%
    filter(region == region_name)
  
  ggplot() +
    geom_sf(
      data = taiwan,
      fill = "olivedrab3",
      color = NA
    ) +
    geom_sf(
      data = site_sub_sf,
      color = "steelblue",
      size = 2.5
    ) +
    geom_text_repel(
      data = site_sub,
      aes(x = lon, y = lat, label = site),
      size = 3.7,
      color = "black",
      box.padding = 0.4,
      point.padding = 0.3,
      min.segment.length = 0,
      segment.color = "grey40",
      max.overlaps = Inf
    ) +
    annotation_scale(
      location = "bl",
      width_hint = scale_width,
      text_cex = 0.7
    ) +
    coord_sf(
      xlim = xlim,
      ylim = ylim,
      expand = F
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position = "none",
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid.major = element_line(color = "gray90", linewidth = 0.4),
      panel.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(2, 2, 2, 2)
    )
}

p_NE <- make_inset_map(
  region_name = "NE",
  xlim = extent_NE$xlim,
  ylim = extent_NE$ylim,
  scale_width = 0.5
)

p_GI <- make_inset_map(
  region_name = "GI",
  xlim = extent_GI$xlim,
  ylim = extent_GI$ylim,
  scale_width = 0.4
)

p_XLQ <- make_inset_map(
  region_name = "XLQ",
  xlim = extent_XLQ$xlim,
  ylim = extent_XLQ$ylim,
  scale_width = 0.35
)

right_col <- p_NE / p_GI / p_XLQ

final_map <- ggdraw() +
  draw_plot(
    p_main,
    x = 0.00, y = 0.00,
    width = 0.60, height = 1.00
  ) +
  draw_plot(
    right_col,
    x = 0.54, y = 0.03,
    width = 0.43, height = 0.97
  )

final_map

