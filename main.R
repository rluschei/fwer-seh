

rm(list=ls())

source("libraries.R")

tar_make()

files <- list.files(path = "_targets/objects", pattern = "^results.*", full.names = TRUE)
results <- lapply(files, readRDS) %>% bind_rows()

summary <- results %>%
  group_by(N, EHF, CHF, alloc, case, test, rate) %>%
  summarise(
    n = sum(!is.na(value)),
    mean = ifelse(n == 0, NA, round(mean(value, na.rm = TRUE), 4)),
    sd = ifelse(n == 0, NA, round(sd(value, na.rm = TRUE), 4)),
    min = ifelse(n == 0, NA, round(min(value, na.rm = TRUE), 4)),
    q1 = ifelse(n == 0, NA, round(quantile(value, 0.25, na.rm = TRUE), 4)),
    median = ifelse(n == 0, NA, round(median(value, na.rm = TRUE), 4)),
    q3 = ifelse(n == 0, NA, round(quantile(value, 0.75, na.rm = TRUE), 4)),
    max = ifelse(n == 0, NA, round(max(value, na.rm = TRUE), 4)),
    .groups = "drop"
  ) %>%
  dplyr::select(-n)
