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

# import data
sed <- read_excel("CNP_data.xlsx", sheet = "sediment", na = "NA")
P <- read_excel("CNP_data.xlsx", sheet = "P", na = "NA")
IO <- read_excel("CNP_data.xlsx", sheet = "Inorganic vs. Organic", na = "NA")

## Sediment across region----

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
    A_cmgL = `sample A_P concentration(mg/L)`,
    B_cmgL = `sample B_P concentration(mg/L)`,
    A_wmg = `A subsample_P weight(mg)`,
    B_wmg = `B subsample_P weight(mg)`
  ) %>%
  mutate(
    bulk_P_mgg = (A_cmgL * V_ml / 1000) / A_wmg * 1000,
    algal_P_mgg = (B_cmgL * V_ml / 1000)/ B_wmg * 1000
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

# sediment weight vs. depth
sed_depth <- 
  ggplot(sed_site_mean, aes(x = depth, y = sediment_weight)) + 
  geom_point() +
  geom_smooth(method = "lm", se = F) +
  facet_wrap(~Region) +
  theme_bw()
sed_depth

# sediment weight vs. temperature
sed_temp <- 
  ggplot(sed_site_mean, aes(x = temperature, y = sediment_weight)) + 
  geom_point() +
  geom_smooth(method = "lm", se = F) +
  facet_wrap(~Region) +
  theme_bw()
sed_temp

# sediment weight vs. region
sed_region <- 
  ggplot(sed_site_mean, aes(Region, sediment_weight, fill = Region)) +
  geom_boxplot() +
  scale_fill_manual(values = region_color) +
  guides(fill = "none") +
  geom_jitter() +
  theme_bw()
sed_region

## Sea urchins----

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
  scale_fill_manual(values = urchin_color, drop = F) +
  theme_bw() +
  labs(x = "Region", y = "Total number of urchins", fill = "Genus")
herb_urchin_sum_region


## Nutrient concentration across region----


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


P_join <- P %>%
  select(Region, site, sediment, bulk_P_mgg, algal_P_mgg)

CNP_sample <- CN %>%
  left_join(P_join, by = c("Region", "site", "sediment")) %>%
  mutate(
    Camera = gsub("S", "C", sediment),
    
    # P molar
    bulk_P_mmolg  = bulk_P_mgg / 31,
    algal_P_mmolg = algal_P_mgg / 31,
    
    # molar ratios
    bulk_CN_ratio  = bulk_C_mmolg / bulk_N_mmolg,
    algal_CN_ratio = algal_C_mmolg / algal_N_mmolg,
    bulk_CP_ratio  = bulk_C_mmolg / bulk_P_mmolg,
    algal_CP_ratio = algal_C_mmolg / algal_P_mmolg,
    bulk_NP_ratio  = bulk_N_mmolg / bulk_P_mmolg,
    algal_NP_ratio = algal_N_mmolg / algal_P_mmolg
  )


CNP_site_mean <- CNP_sample %>%
  group_by(Region, site) %>%
  summarise(
    # concentration
    bulk_C_mean = mean(bulk_C_mgg, na.rm = TRUE),
    bulk_N_mean = mean(bulk_N_mgg, na.rm = TRUE),
    bulk_P_mean = mean(bulk_P_mgg, na.rm = TRUE),
    algal_C_mean = mean(algal_C_mgg, na.rm = TRUE),
    algal_N_mean = mean(algal_N_mgg, na.rm = TRUE),
    algal_P_mean = mean(algal_P_mgg, na.rm = TRUE),
    
    # molar ratios
    bulk_CN_ratio = mean(bulk_CN_ratio, na.rm = TRUE),
    algal_CN_ratio = mean(algal_CN_ratio, na.rm = TRUE),
    bulk_CP_ratio = mean(bulk_CP_ratio, na.rm = TRUE),
    algal_CP_ratio = mean(algal_CP_ratio, na.rm = TRUE),
    bulk_NP_ratio = mean(bulk_NP_ratio, na.rm = TRUE),
    algal_NP_ratio = mean(algal_NP_ratio, na.rm = TRUE),
    .groups = "drop"
  )


CNP_conc_long <- CNP_site_mean %>%
  pivot_longer(
    cols = c(
      bulk_C_mean, bulk_N_mean, bulk_P_mean,
      algal_C_mean, algal_N_mean, algal_P_mean
    ),
    names_to = "type",
    values_to = "value"
  )

ggplot(CNP_conc_long, aes(x = Region, y = value, fill = Region)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.8) +
  scale_fill_manual(values = region_color) +
  guides(fill = "none") +
  facet_wrap(~ type, scales = "free_y") +
  theme_bw() +
  labs(
    x = "Region",
    y = "Concentration (mg/g)"
  )


CNP_ratio_long <- CNP_site_mean %>%
  pivot_longer(
    cols = c(
      bulk_CN_ratio, algal_CN_ratio,
      bulk_CP_ratio, algal_CP_ratio,
      bulk_NP_ratio, algal_NP_ratio
    ),
    names_to = "type",
    values_to = "value"
  )

ggplot(CNP_ratio_long, aes(x = Region, y = value, fill = Region)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.8) +
  scale_fill_manual(values = region_color) +
  guides(fill = "none") +
  facet_wrap(~ type, scales = "free_y") +
  theme_bw() +
  labs(
    x = "Region",
    y = "Molar ratio"
  )

## IO across region----

IO <- IO %>%
  rename(sediment_id = sediment)

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

## Fish bites----

### Trophic group/Herbivore functional group bite rate

bite <- read_excel("Fish bite_data.xlsx")

area <- 1 # m^2
duration <- 35 # 35 mins
bite <- bite %>%
  mutate(bite_rate = bites / area * 60 / duration, # bites/hr*m^2
         `trophic group` = fct_relevel(
           `trophic group`,
           "Herbivore", "Benthic invertivore", "Corallivore", "Omnivore"
         ))

# bite rate of different trophic groups
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

# bite rate of different herbivore functional groups
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

# Bite rate by quadrat

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

# 1. occurrence of herbivore feeding

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


# 2. Herbivore grazing pressure (total bite rate/event rate/bite per event)

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
    title = "Herbivore Grazing Pressure",
    x = "Region",
    y = NULL
  )
herb_grazing_pressure

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
  ggplot(bites_per_event, aes(Region, value, fill = Region)) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.6) +
  stat_summary(fun = mean, geom = "point", color = "red") +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width = 0.12) +
  scale_fill_manual(values = region_color) +
  guides(fill = "none") +
  theme_bw() +
  labs(
    title = "Bites per event",
    x = "Region",
    y = NULL
  )
bites_per_event_plot


# 3. Grazing pressure vs. region model, site as random effect

# total bite/number of event -> GLMMs
# bite per event -> LMMs

total_bite_region <- glmmTMB(
  total_bites ~ Region + (1 | site),
  family = nbinom2,
  data = bite_by_quadrat
)

event_region <- glmmTMB(
  n_events ~ Region + (1 | site),
  family = nbinom2,
  data = bite_by_quadrat
)

check_overdispersion(total_bite_region)
check_zeroinflation(total_bite_region)
res_total <- simulateResiduals(total_bite_region)
plot(res_total)
testDispersion(res_total)
testZeroInflation(res_total)
testOutliers(res_total)

check_overdispersion(event_region)
check_zeroinflation(event_region)
res_event <- simulateResiduals(event_region)
plot(res_event)
testDispersion(res_event)
testZeroInflation(res_event)
testOutliers(res_event)

summary(total_bite_region)
summary(event_region)

emmeans(total_bite_region, pairwise ~ Region, type = "response")
emmeans(event_region, pairwise ~ Region, type = "response")

bite_per_event_region <- glmmTMB(
  mean_bites_per_event ~ Region + (1 | site),
  family = Gamma(link = "log"),
  data = herb_bite_positive
)

sim_res <- simulate_residuals(bite_per_event_region)
plot(sim_res)
testDispersion(sim_res)

summary(bite_per_event_region)

emmeans(bite_per_event_region, pairwise ~ Region, type = "response")

## Bite per event adjust for nutrient----

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

herb_bite_positive <- herb_bite_positive %>%
  left_join(IO_site_mean, by = c("Region", "site")) %>%
  left_join(CNP_site_mean, by = c("Region", "site")) %>%
  mutate(
    bites_per_event_org_adj = mean_bites_per_event / org_g_mean,
    bites_per_event_C_adj   = mean_bites_per_event / bulk_C_mean,
    bites_per_event_N_adj   = mean_bites_per_event / bulk_N_mean,
    bites_per_event_P_adj   = mean_bites_per_event / algal_P_mean
  )

# Long format for plotting

bite_per_event_adj <- herb_bite_positive %>%
  select(
    Region,
    bites_per_event_org_adj,
    bites_per_event_C_adj,
    bites_per_event_N_adj,
    bites_per_event_P_adj
  ) %>%
  pivot_longer(
    cols = c(
      bites_per_event_org_adj,
      bites_per_event_C_adj,
      bites_per_event_N_adj,
      bites_per_event_P_adj
    ),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = factor(
      metric,
      levels = c(
        "bites_per_event_org_adj",
        "bites_per_event_C_adj",
        "bites_per_event_N_adj",
        "bites_per_event_P_adj"
      ),
      labels = c(
        "adjusted for organic matter",
        "adjusted for bulk C",
        "adjusted for bulk N",
        "adjusted for algal P"
      )
    )
  )

# Plot

bite_per_event_ad_plot <- 
  ggplot(bite_per_event_adj, aes(Region, value, fill = Region)) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.6) +
  stat_summary(fun = mean, geom = "point", color = "red") +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width = 0.12) +
  scale_fill_manual(values = region_color) +
  facet_wrap(~metric, scales = "free_y", nrow = 1) +
  guides(fill = "none") +
  theme_bw() +
  labs(
    title = "Bites per event adjusted for nutrient content",
    x = "Region",
    y = NULL
  )

bite_per_event_ad_plot