setwd("C:/Users/dan91/OneDrive/桌面/Field work")
library(readxl)
library(dplyr)
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

P <- P %>%
  rename(
    A_cmgL = `sample A_P concentration(mg/L)`,
    B_cmgL = `sample B_P concentration(mg/L)`,
    A_wmg = `A subsample_P weight(mg)`,
    B_wmg = `B subsample_P weight(mg)`
  ) %>%
  mutate(
    bulk_P_cmgg = A_cmgL * 5 / A_wmg,
    algae_P_cmgg = B_cmgL * 5 / B_wmg
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

library(lme4)
library(ggplot2)

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

# Relationship between P concentration & Region/environmental parameters----

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

