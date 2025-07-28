library(readxl)
library(ggplot2)
df <- read_excel("C:/Users/dan91/OneDrive/桌面/urchin_sediment.xlsx", sheet = "sediment", na = "NA")

#Total sediment weight different among site----
# ANOVA
anova_tot_W <- aov(`sample A sediment(g)` ~ XLQ, data = df)
summary(anova_tot_W)

# Tukey post-hoc
TukeyHSD(anova_tot_W)

# boxplot
ggplot(df, aes(x = XLQ, y = `sample A sediment(g)`)) +
  geom_boxplot(fill = "lightblue") +
  labs(title = "Total sediment weight", x = "Site", y = "Total weight (g)")

#P concentration different among site----
# ANOVA
anova_p <- aov(`sample A_P (μmol/L)` ~ XLQ, data = df)
summary(anova_p)

# Tukey
TukeyHSD(anova_p)

# boxplot
ggplot(df, aes(x = XLQ, y = `sample A_P (μmol/L)`)) +
  geom_boxplot(fill = "lightgreen") +
  labs(title = "P concentration", x = "Site", y = "P concentration (μmol/L)")

#Relationship between P concentration and subsample weight----
library(ggpubr)

# correlation
cor.test(df$`A subsample_P weight(mg)`, df$`sample A_P (μmol/L)`, method = "spearman")

# plot
ggplot(df, aes(x = `A subsample_P weight(mg)`, y = `sample A_P (μmol/L)`)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, color = "red") +
  stat_cor(method = "spearman", label.x = min(df$`A subsample_P weight(mg)`, na.rm = TRUE),
           label.y = max(df$`sample A_P (μmol/L)`, na.rm = TRUE)) +
  labs(title = "Subsample Weight vs P Concentration",
       x = "Subsample Weight (mg)", y = "P Concentration (μmol/L)")

#Relationship between total weight and depth/temperature----
# correlation
cor.test(df$depth, df$`sample A sediment(g)`)
df_temp <- df[!is.na(df$temperature), ]
cor.test(df_temp$temperature, df_temp$`sample A sediment(g)`)

# plot
ggplot(df, aes(x = depth, y = `sample A sediment(g)`)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  stat_cor(method = "spearman", label.x = min(df$depth, na.rm = TRUE), 
           label.y = max(df$`sample A sediment(g)`, na.rm = TRUE)) +
  labs(title = "Total Sediment Weight vs Depth",
       x = "Depth (m)", y = "Total Sediment Weight (g)")

ggplot(df_temp, aes(x = temperature, y = `sample A sediment(g)`)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  stat_cor(method = "spearman", label.x = min(df_temp$temperature), 
           label.y = max(df_temp$`sample A sediment(g)`, na.rm = TRUE)) +
  labs(title = "Total Sediment Weight vs Temperature",
       x = "Temperature (°C)", y = "Total Sediment Weight (g)")

#Relationship between P concentration and depth/temperature----
cor.test(df$depth, df$`sample A_P (μmol/L)`)
cor.test(df_temp$temperature, df_temp$`sample A_P (μmol/L)`)

ggplot(df, aes(x = depth, y = `sample A_P (μmol/L)`)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  stat_cor(method = "spearman", label.x = min(df$depth, na.rm = TRUE), 
           label.y = max(df$`sample A_P (μmol/L)`, na.rm = TRUE)) +
  labs(title = "P Concentration vs Depth",
       x = "Depth (m)", y = "P Concentration (μmol/L)")

ggplot(df, aes(x = temperature, y = `sample A_P (μmol/L)`)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  stat_cor(method = "spearman", label.x = min(df_temp$temperature), 
           label.y = max(df$`sample A_P (μmol/L)`, na.rm = TRUE)) +
  labs(title = "P Concentration vs Temperature",
       x = "Temperature (°C)", y = "P Concentration (μmol/L)")
