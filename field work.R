# Result

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

# import data----
sed <- read_excel("CNP_data.xlsx", sheet = "sediment", na = "NA")
P <- read_excel("CNP_data.xlsx", sheet = "P", na = "NA")
IO <- read_excel("CNP_data.xlsx", sheet = "Inorganic vs. Organic", na = "NA")

# Sediment across region----

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


# Sea urchins----

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
  scale_fill_manual(values = urchin_color, drop = FALSE) +
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

kruskal.test(density ~ Region, data = urchin_density)

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

urchin_plot / herb_urchin_sum_region

# Nutrient concentration across region----

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

## P adjusted to OM basis ----

P_org <- P %>%
  left_join(
    IO %>% select(Region, site, sediment_id, inorg_pct, org_pct),
    by = c("Region", "site", "sediment_id")
  ) %>%
  mutate(
    bulk_P_org_mgg = if_else(org_pct > 0, bulk_P_mgg / org_pct * 100, NA_real_),
    bulk_P_org_mmolg = bulk_P_org_mgg / 31,
    algal_P_mmolg = algal_P_mgg / 31
  )

CNP_site_mean_adj <- CN %>%
  left_join(
    P_org %>% select(
      Region, site, sediment_id,
      algal_P_mgg, algal_P_mmolg,
      bulk_P_org_mgg, bulk_P_org_mmolg
    ),
    by = c("Region", "site", "sediment_id")
  ) %>%
  mutate(
    Camera = gsub("S", "C", sediment_id)
  ) %>%
  group_by(Region, site) %>%
  summarise(
    across(
      c(
        bulk_C_mgg, bulk_N_mgg, bulk_P_org_mgg,
        algal_C_mgg, algal_N_mgg, algal_P_mgg,
        bulk_C_mmolg, bulk_N_mmolg, bulk_P_org_mmolg,
        algal_C_mmolg, algal_N_mmolg, algal_P_mmolg
      ),
      ~ mean(.x, na.rm = TRUE),
      .names = "{.col}_mean"
    ),
    .groups = "drop"
  ) %>%
  mutate(
    bulk_CN_ratio_adj = bulk_C_mmolg_mean / bulk_N_mmolg_mean,
    algal_CN_ratio    = algal_C_mmolg_mean / algal_N_mmolg_mean,
    bulk_CP_ratio_adj = bulk_C_mmolg_mean / bulk_P_org_mmolg_mean,
    algal_CP_ratio    = algal_C_mmolg_mean / algal_P_mmolg_mean,
    bulk_NP_ratio_adj = bulk_N_mmolg_mean / bulk_P_org_mmolg_mean,
    algal_NP_ratio    = algal_N_mmolg_mean / algal_P_mmolg_mean
  )

CNP_conc_long_adj <- CNP_site_mean_adj %>%
  pivot_longer(
    cols = c(
      bulk_C_mgg_mean, bulk_N_mgg_mean, bulk_P_org_mgg_mean,
      algal_C_mgg_mean, algal_N_mgg_mean, algal_P_mgg_mean
    ),
    names_to = "type",
    values_to = "value"
  )

CNP_conc_adj <- 
  ggplot(CNP_conc_long_adj, aes(x = Region, y = value, fill = Region)) +
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
CNP_conc_adj

CNP_ratio_long_adj <- CNP_site_mean_adj %>%
  pivot_longer(
    cols = c(
      bulk_CN_ratio_adj, algal_CN_ratio,
      bulk_CP_ratio_adj, algal_CP_ratio,
      bulk_NP_ratio_adj, algal_NP_ratio
    ),
    names_to = "type",
    values_to = "value"
  )

CNP_ratio_adj <- 
  ggplot(CNP_ratio_long_adj, aes(x = Region, y = value, fill = Region)) +
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
CNP_ratio_adj

bulk_CNP_ratio_table_adj <- CNP_site_mean_adj %>%
  group_by(Region) %>%
  summarise(
    bulk_C = mean(bulk_C_mmolg_mean, na.rm = TRUE),
    bulk_N = mean(bulk_N_mmolg_mean, na.rm = TRUE),
    bulk_P = mean(bulk_P_org_mmolg_mean, na.rm = TRUE),
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
bulk_CNP_ratio_table_adj

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

CNP_conc_adj_2 <- 
  ggplot(CNP_conc_long_adj_2, aes(x = Region, y = value, fill = Region)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.8) +
  scale_fill_manual(values = region_color) +
  guides(fill = "none") +
  facet_wrap(~ type, scales = "free_y",
             labeller = as_labeller(c(
               "bulk_C_mgg_adj_mean" = "Carbon (mg g⁻¹)",
               "bulk_N_mgg_adj_mean" = "Nitrogen (mg g⁻¹)",
               "bulk_P_mgg_mean"     = "Phosphorus (mg g⁻¹)"
             ))) +
  theme_bw() +
  labs(
    x = "Region",
    y = "Concentration (mg/g)"
  )
CNP_conc_adj_2

CNP_ratio_long_adj_2 <- CNP_site_mean_adj_2 %>%
  pivot_longer(
    cols = c(
      bulk_CN_ratio_adj, algal_CN_ratio,
      bulk_CP_ratio_adj, algal_CP_ratio,
      bulk_NP_ratio_adj, algal_NP_ratio
    ),
    names_to = "type",
    values_to = "value"
  )

CNP_ratio_adj_2 <- 
  ggplot(CNP_ratio_long_adj_2, aes(x = Region, y = value, fill = Region)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.8) +
  scale_fill_manual(values = region_color) +
  guides(fill = "none") +
  facet_wrap(~ type, scales = "free_y",
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
CNP_ratio_adj_2

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

kruskal.test(bulk_C_mgg_adj_mean ~ Region, data = CNP_site_mean_adj_2)

pairwise.wilcox.test(
  CNP_site_mean_adj_2$bulk_C_mgg_adj_mean,
  CNP_site_mean_adj_2$Region,
  p.adjust.method = "BH"
)

kruskal.test(bulk_N_mgg_adj_mean ~ Region, data = CNP_site_mean_adj_2)

pairwise.wilcox.test(
  CNP_site_mean_adj_2$bulk_N_mgg_adj_mean,
  CNP_site_mean_adj_2$Region,
  p.adjust.method = "BH"
)

kruskal.test(bulk_P_mgg_mean ~ Region, data = CNP_site_mean_adj_2)

pairwise.wilcox.test(
  CNP_site_mean_adj_2$bulk_P_mgg_mean,
  CNP_site_mean_adj_2$Region,
  p.adjust.method = "BH"
)

kruskal.test(bulk_CN_ratio_adj ~ Region, data = CNP_site_mean_adj_2)

kruskal.test(bulk_CP_ratio_adj ~ Region, data = CNP_site_mean_adj_2)

kruskal.test(bulk_NP_ratio_adj ~ Region, data = CNP_site_mean_adj_2)

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
p1

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
p2

p1/p2

## IO across region----

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

kruskal.test(org_pct ~ Region, data = IO_site)

pairwise.wilcox.test(
  IO_site$org_pct,
  IO_site$Region,
  p.adjust.method = "BH"
)

# Fish bites----

bite <- read_excel("Fish bite_data.xlsx")

area <- 1 # m^2
duration <- 35 # 35 mins
bite <- bite %>%
  mutate(bite_rate = bites / area * 60 / duration, # bites/hr*m^2
         `trophic group` = fct_relevel(
           `trophic group`,
           "Herbivore", "Benthic invertivore", "Corallivore", "Omnivore"
         ))

## bite rate of different trophic groups----
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

## bite rate of different herbivore functional groups----

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

## Bite rate by quadrat----

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

## 1.  occurrence of herbivore feeding----

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

## 2.  Herbivore grazing pressure (total bite rate/event rate/bite per event)----

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

kruskal.test(total_bite_rate ~ Region, data = bite_by_quadrat)

pairwise.wilcox.test(
  bite_by_quadrat$total_bite_rate,
  bite_by_quadrat$Region,
  p.adjust.method = "BH"
)

kruskal.test(event_rate ~ Region, data = bite_by_quadrat)

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

kruskal.test(mean_bites_per_event ~ Region, data = herb_bite_positive)

## Bite per event adjust by nutrient----

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
  left_join(
    IO_site_mean, 
    by = c("Region", "site")) %>%
  left_join(
    CNP_site_mean_adj_2, 
    by = c("Region", "site")) %>%
  mutate(
    bites_per_event_org_adj = mean_bites_per_event / org_g_mean,
    bites_per_event_CN_adj = mean_bites_per_event / bulk_CN_ratio_adj,
    bites_per_event_CP_adj = mean_bites_per_event / bulk_CP_ratio_adj,
    bites_per_event_NP_adj = mean_bites_per_event / bulk_NP_ratio_adj
  )

# Long format for plotting

bite_per_event_adj <- herb_bite_positive %>%
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

# Plot

bite_per_event_adj_plot <- 
  ggplot(bite_per_event_adj, aes(Region, value, fill = Region)) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.6) +
  stat_summary(fun = mean, geom = "point", color = "red") +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width = 0.12) +
  scale_fill_manual(values = region_color) +
  facet_wrap(~metric, scales = "free_y", nrow = 1) +
  guides(fill = "none") +
  theme_bw() +
  labs(
    title = "Bites per event adjusted by nutrient quality",
    x = "Region",
    y = NULL
  )

bite_per_event_adj_plot

## correlation of nutritional quality and bite per event----

cor.test(herb_bite_positive$mean_bites_per_event,
         herb_bite_positive$org_pct_mean,
         method = "spearman")

cor.test(herb_bite_positive$mean_bites_per_event,
         herb_bite_positive$bulk_CN_ratio_adj,
         method = "spearman")

cor.test(herb_bite_positive$mean_bites_per_event,
         herb_bite_positive$bulk_CP_ratio_adj,
         method = "spearman")

cor.test(herb_bite_positive$mean_bites_per_event,
         herb_bite_positive$bulk_NP_ratio_adj,
         method = "spearman")

bite_vs_ratio <- herb_bite_positive %>%
  select(mean_bites_per_event, bulk_CN_ratio_adj, bulk_CP_ratio_adj, bulk_NP_ratio_adj, Region) %>%
  pivot_longer(
    cols = c(bulk_CN_ratio_adj, bulk_CP_ratio_adj, bulk_NP_ratio_adj),
    names_to = "Ratio",
    values_to = "Value"
  )

bite_vs_ratio$Ratio <- recode(bite_vs_ratio$Ratio,
                          "bulk_CP_ratio_adj" = "C:P ratio",
                          "bulk_NP_ratio_adj" = "N:P ratio",
                          "bulk_CN_ratio_adj" = "C:N ratio")

bite_vs_ratio_plot <- 
  ggplot(bite_vs_ratio, aes(x = Value, y = mean_bites_per_event, color = Region)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~ Ratio, scales = "free_x") +
  scale_color_manual(values = region_color) +
  labs(
    x = "Stoichiometric ratio",
    y = "Bite per event",
    color = "Region"
  ) +
  theme_bw()
bite_vs_ratio_plot

NE_bite_per_event <- herb_bite_positive %>%
  filter(Region == "NE")

cor.test(NE_bite_per_event$mean_bites_per_event,
         NE_bite_per_event$bulk_CN_ratio_adj,
         method = "spearman")

cor.test(NE_bite_per_event$mean_bites_per_event,
         NE_bite_per_event$bulk_CP_ratio_adj,
         method = "spearman")

cor.test(NE_bite_per_event$mean_bites_per_event,
         NE_bite_per_event$bulk_NP_ratio_adj,
         method = "spearman")

GI_bite_per_event <- herb_bite_positive %>%
  filter(Region == "GI")

cor.test(GI_bite_per_event$mean_bites_per_event,
         GI_bite_per_event$bulk_CN_ratio_adj,
         method = "spearman")


# modeling to control nutrient for bite per event----

herb_bite_positive$log_bpe <- log(herb_bite_positive$mean_bites_per_event)


hist(herb_bite_positive$log_bpe)

m_CN <- lmer(
  log_bpe ~ bulk_CN_ratio_adj + Region + (1|site),
  data = herb_bite_positive
)
summary(m_CN)

plot(m_CN)

qqnorm(residuals(m_CN))
qqline(residuals(m_CN))

m_org <- lmer(
  log_bpe ~ org_pct_mean + Region + (1|site),
  data = herb_bite_positive
)
summary(m_org)

plot(m_org)

qqnorm(residuals(m_org))
qqline(residuals(m_org))

m_CP <- lmer(
  log_bpe ~ bulk_CP_ratio_adj + Region + (1|site),
  data = herb_bite_positive
)
summary(m_CP)

plot(m_CP)

qqnorm(residuals(m_CP))
qqline(residuals(m_CP))

m_NP <- lmer(
  log_bpe ~ bulk_NP_ratio_adj + Region + (1|site),
  data = herb_bite_positive
)
summary(m_NP)

plot(m_NP)

qqnorm(residuals(m_NP))
qqline(residuals(m_NP))

# log bite per event plot----

bite_log_vs_ratio <- herb_bite_positive %>%
  select(log_bpe, bulk_CN_ratio_adj, bulk_CP_ratio_adj, bulk_NP_ratio_adj, Region) %>%
  pivot_longer(
    cols = c(bulk_CN_ratio_adj, bulk_CP_ratio_adj, bulk_NP_ratio_adj),
    names_to = "Ratio",
    values_to = "Value"
  )

bite_log_vs_ratio$Ratio <- recode(bite_log_vs_ratio$Ratio,
                              "bulk_CP_ratio_adj" = "C:P ratio",
                              "bulk_NP_ratio_adj" = "N:P ratio",
                              "bulk_CN_ratio_adj" = "C:N ratio")

bite_log_vs_ratio_plot <- 
  ggplot(bite_log_vs_ratio, aes(x = Value, y = log_bpe, color = Region)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.8) +
  facet_wrap(~ Ratio, scales = "free_x") +
  scale_color_manual(values = region_color) +
  labs(
    x = "Stoichiometric ratio",
    y = "Bite per event log transformed",
    color = "Region"
  ) +
  theme_bw()
bite_log_vs_ratio_plot

# all region smooth line
bite_log_vs_ratio_all_plot <- 
  ggplot(bite_log_vs_ratio, aes(x = Value, y = log_bpe)) +
  geom_point(aes(color = Region), size = 2) +
  geom_smooth(aes(group = 1), method = "lm", se = FALSE, linewidth = 0.8, colour = "grey40") +
  facet_wrap(~ Ratio, scales = "free_x") +
  scale_color_manual(values = region_color) +
  labs(
    x = "Stoichiometric ratio",
    y = "Bite per event log transformed",
    color = "Region"
  ) +
  theme_bw()
bite_log_vs_ratio_all_plot

## log bpe vs org pct

bite_log_vs_org_pct <- herb_bite_positive %>%
  select(log_bpe, org_pct_mean, Region)

bite_log_vs_org_pct_plot <- 
  ggplot(bite_log_vs_org_pct, aes(x = org_pct_mean, y = log_bpe, color = Region)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.8) +
  scale_color_manual(values = region_color) +
  labs(
    x = "Organic matter percentage",
    y = "Bite per event log transformed",
    color = "Region"
  ) +
  theme_bw()
bite_log_vs_org_pct_plot

bite_log_vs_org_pct_all_plot <- 
  ggplot(bite_log_vs_org_pct, aes(x = org_pct_mean, y = log_bpe)) +
  geom_point(aes(color = Region), size = 2) +
  geom_smooth(aes(group = 1), method = "lm", se = FALSE, linewidth = 0.8, colour = "grey40") +
  scale_color_manual(values = region_color) +
  labs(
    x = "Organic matter percentage",
    y = "Bite per event log transformed",
    color = "Region"
  ) +
  theme_bw()
bite_log_vs_org_pct_all_plot


## Quadrat-level BPE vs. organic matter and CNP ratios ----

herb_bite_positive_clean <- bite_by_quadrat %>%
  filter(n_events > 0)

IO_quadrat <- IO %>%
  mutate(
    Camera = gsub("S", "C", sediment_id)
  ) %>%
  select(
    Region,
    site,
    Camera,
    org_pct
  )

CNP_quadrat <- CNP_adj_2 %>%
  mutate(
    bulk_CN_ratio_adj = bulk_C_mmolg_adj / bulk_N_mmolg_adj,
    bulk_CP_ratio_adj = bulk_C_mmolg_adj / bulk_P_mmolg,
    bulk_NP_ratio_adj = bulk_N_mmolg_adj / bulk_P_mmolg
  ) %>%
  select(
    Region,
    site,
    Camera,
    bulk_CN_ratio_adj,
    bulk_CP_ratio_adj,
    bulk_NP_ratio_adj
  )

herb_bite_quadrat <- herb_bite_positive_clean %>%
  left_join(
    IO_quadrat,
    by = c("Region", "site", "Camera")
  ) %>%
  left_join(
    CNP_quadrat,
    by = c("Region", "site", "Camera")
  ) %>%
  mutate(
    log_bpe = log(mean_bites_per_event)
  )

## Check missing values after joining ----

herb_bite_quadrat %>%
  summarise(
    n_total = n(),
    n_org_NA = sum(is.na(org_pct)),
    n_CN_NA = sum(is.na(bulk_CN_ratio_adj)),
    n_CP_NA = sum(is.na(bulk_CP_ratio_adj)),
    n_NP_NA = sum(is.na(bulk_NP_ratio_adj))
  )

## Spearman correlation ----

cor_org <- herb_bite_quadrat %>%
  select(log_bpe, org_pct) %>%
  drop_na()

cor.test(
  cor_org$log_bpe,
  cor_org$org_pct,
  method = "spearman"
)

cor_CN <- herb_bite_quadrat %>%
  select(log_bpe, bulk_CN_ratio_adj) %>%
  drop_na()

cor.test(
  cor_CN$log_bpe,
  cor_CN$bulk_CN_ratio_adj,
  method = "spearman"
)

cor_CP <- herb_bite_quadrat %>%
  select(log_bpe, bulk_CP_ratio_adj) %>%
  drop_na()

cor.test(
  cor_CP$log_bpe,
  cor_CP$bulk_CP_ratio_adj,
  method = "spearman"
)

cor_NP <- herb_bite_quadrat %>%
  select(log_bpe, bulk_NP_ratio_adj) %>%
  drop_na()

cor.test(
  cor_NP$log_bpe,
  cor_NP$bulk_NP_ratio_adj,
  method = "spearman"
)

## Plot: log BPE vs. C:N / C:P / N:P ----

bite_log_vs_ratio_quadrat <- herb_bite_quadrat %>%
  select(
    log_bpe,
    bulk_CN_ratio_adj,
    bulk_CP_ratio_adj,
    bulk_NP_ratio_adj,
    Region
  ) %>%
  pivot_longer(
    cols = c(
      bulk_CN_ratio_adj,
      bulk_CP_ratio_adj,
      bulk_NP_ratio_adj
    ),
    names_to = "Ratio",
    values_to = "Value"
  ) %>%
  drop_na(log_bpe, Value)

bite_log_vs_ratio_quadrat$Ratio <- recode(
  bite_log_vs_ratio_quadrat$Ratio,
  "bulk_CN_ratio_adj" = "C:N ratio",
  "bulk_CP_ratio_adj" = "C:P ratio",
  "bulk_NP_ratio_adj" = "N:P ratio"
)

bite_log_vs_ratio_quadrat_plot <- 
  ggplot(
    bite_log_vs_ratio_quadrat,
    aes(x = Value, y = log_bpe, color = Region)
  ) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.8) +
  facet_wrap(~ Ratio, scales = "free_x") +
  scale_color_manual(values = region_color) +
  labs(
    x = "Stoichiometric ratio",
    y = "Log bite per event",
    color = "Region"
  ) +
  theme_bw()
bite_log_vs_ratio_quadrat_plot

bite_log_vs_ratio_quadrat_all_plot <- 
  ggplot(
    bite_log_vs_ratio_quadrat,
    aes(x = Value, y = log_bpe)
  ) +
  geom_point(aes(color = Region), size = 2) +
  geom_smooth(
    aes(group = 1),
    method = "lm",
    se = FALSE,
    linewidth = 0.8,
    colour = "grey40"
  ) +
  facet_wrap(~ Ratio, scales = "free_x") +
  scale_color_manual(values = region_color) +
  labs(
    x = "Stoichiometric ratio",
    y = "Log bite per event",
    color = "Region"
  ) +
  theme_bw()

bite_log_vs_ratio_quadrat_all_plot

## Plot: log BPE vs. organic matter percentage ----

bite_log_vs_org_quadrat <- herb_bite_quadrat %>%
  select(log_bpe, org_pct, Region) %>%
  drop_na(log_bpe, org_pct)

bite_log_vs_org_quadrat_plot <- 
  ggplot(
    bite_log_vs_org_quadrat,
    aes(x = org_pct, y = log_bpe, color = Region)
  ) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.8) +
  scale_color_manual(values = region_color) +
  labs(
    x = "Organic matter percentage",
    y = "Log bite per event",
    color = "Region"
  ) +
  theme_bw()

bite_log_vs_org_quadrat_plot

bite_log_vs_org_quadrat_all_plot <- 
  ggplot(
    bite_log_vs_org_quadrat,
    aes(x = org_pct, y = log_bpe)
  ) +
  geom_point(aes(color = Region), size = 2) +
  geom_smooth(
    aes(group = 1),
    method = "lm",
    se = FALSE,
    linewidth = 0.8,
    colour = "grey40"
  ) +
  scale_color_manual(values = region_color) +
  labs(
    x = "Organic matter percentage",
    y = "Log bite per event",
    color = "Region"
  ) +
  theme_bw()

bite_log_vs_org_quadrat_all_plot

# change bite per event into site-level----

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

# bpe site mean adj plot

bpe_site_adj_plot <- 
  ggplot(bpe_site_adj, aes(Region, value, fill = Region)) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.6) +
  stat_summary(fun = mean, geom = "point", color = "red") +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width = 0.12) +
  scale_fill_manual(values = region_color) +
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

# log bpe site mean vs. org
  
log_bpe_site_org_all_plot <- 
  ggplot(
    bpe_site_org,
    aes(x = org_pct_mean, y = log_bpe_site_mean)
  ) +
  geom_point(aes(color = Region), size = 2) +
  geom_smooth(
    aes(group = 1),
    method = "lm",
    se = FALSE,
    linewidth = 0.8,
    colour = "grey40"
  ) +
  scale_color_manual(values = region_color) +
  labs(
    x = "Organic matter percentage",
    y = "Log bite per event site mean",
    color = "Region"
  ) +
  theme_bw()
log_bpe_site_org_all_plot

log_bpe_site_org_plot <- 
  ggplot(
    bpe_site_org,
    aes(x = org_pct_mean, y = log_bpe_site_mean, color = Region)
  ) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.8) +
  scale_color_manual(values = region_color) +
  labs(
    x = "Organic matter percentage",
    y = "Log bite per event site mean",
    color = "Region"
  ) +
  theme_bw()
log_bpe_site_org_plot

# log bpe site mean vs. CNP ratio

log_bpe_site_CNP_ratio <- bpe_site %>%
  select(
    log_bpe_site_mean,
    bulk_CN_ratio_adj,
    bulk_CP_ratio_adj,
    bulk_NP_ratio_adj,
    Region
  ) %>%
  pivot_longer(
    cols = c(
      bulk_CN_ratio_adj,
      bulk_CP_ratio_adj,
      bulk_NP_ratio_adj
    ),
    names_to = "Ratio",
    values_to = "Value"
  )

log_bpe_site_CNP_ratio$Ratio <- recode(
  log_bpe_site_CNP_ratio$Ratio,
  "bulk_CN_ratio_adj" = "C:N ratio",
  "bulk_CP_ratio_adj" = "C:P ratio",
  "bulk_NP_ratio_adj" = "N:P ratio"
)

log_bpe_site_CNP_ratio_all_plot <- 
  ggplot(
    log_bpe_site_CNP_ratio,
    aes(x = Value, y = log_bpe_site_mean)
  ) +
  geom_point(aes(color = Region), size = 2) +
  geom_smooth(
    aes(group = 1),
    method = "lm",
    se = FALSE,
    linewidth = 0.8,
    colour = "grey40"
  ) +
  facet_wrap(~ Ratio, scales = "free_x") +
  scale_color_manual(values = region_color) +
  labs(
    x = "Stoichiometric ratio",
    y = "Log bite per event site mean",
    color = "Region"
  ) +
  theme_bw()
log_bpe_site_CNP_ratio_all_plot

log_bpe_site_CNP_ratio_plot <- 
  ggplot(
    log_bpe_site_CNP_ratio,
    aes(x = Value, y = log_bpe_site_mean, color = Region)
  ) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.8) +
  facet_wrap(~ Ratio, scales = "free_x") +
  scale_color_manual(values = region_color) +
  labs(
    x = "Stoichiometric ratio",
    y = "Log bite per event site mean",
    color = "Region"
  ) +
  theme_bw()
log_bpe_site_CNP_ratio_plot
