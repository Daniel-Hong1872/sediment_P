setwd("C:/Users/dan91/OneDrive/桌面/Field work")
library(readxl)
library(dplyr)
sediment <- read_excel("CNP_data.xlsx", sheet = "sediment", na = "NA")

