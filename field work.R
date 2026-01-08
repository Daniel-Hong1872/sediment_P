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

# Relationship between sediment weight & environmental parameters----

library(lme4)
library(ggplot2)

model_1 <- lmer(sed_wg ~ Region + depth + temperature + (1 | Region/site), 
                data = sediment)
summary(model_1)

model_2 <- lmer(sed_wg ~ Region * depth + Region * temperature + 
                  (1 | Region/site), data = sediment)
summary(model_2)

site_mean <- sediment %>% 
  group_by(Region, site) %>%
  summarise(
    depth = mean(depth, na.rm = T),
    temperature = mean(temperature, na.rm = T),
    sed_wg = mean(sed_wg, na.rm = T),
    .groups = "drop"
  )

ggplot(site_mean, aes(x = depth, y = sed_wg)) + 
  geom_point() +
  geom_smooth(method = "lm", se = F) +
  facet_wrap(~Region)

ggplot(site_mean, aes(x = temperature, y = sed_wg)) + 
  geom_point() +
  geom_smooth(method = "lm", se = F) +
  facet_wrap(~Region)

ggplot(site_mean, aes(x = Region, y = sed_wg)) +
  geom_boxplot() +
  geom_jitter(width = 0.1, size = 2)





