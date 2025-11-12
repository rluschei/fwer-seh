

rm(list=ls())

source("libraries.R")

tar_make()

files <- list.files(path = "_targets/objects", pattern = "^results.*", full.names = TRUE)
results <- lapply(files, readRDS) %>% bind_rows()

means <- results %>%
  group_by(N, EHF, CHF, alloc, case, test, rate) %>%
  summarise(
    mean = round(mean(value), 4),
    .groups = "drop"
  )
