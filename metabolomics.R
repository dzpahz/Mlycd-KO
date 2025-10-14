library(data.table)
library(dplyr)
library(ggplot2)
library(stringr)


#load the transition results
tran_res <- fread("Transition results.csv")


#sample info
#fill in accordingly
WT <- paste("Organ", "sample ids", sep = "_")
KO <- paste("Organ", "sample ids", sep = "_")


data <- trans_res %>%
  filter(`Transition Result Is Quantitative`== 1) %>%
  select(Compound = "Peptide", 
         Replicate,
         Prec = `Precursor Mz`,
         Prod = `Product Mz`,
         Area,
         Background,
         RT = `Retention Time`) %>%
  group_by(Compound) %>%
  mutate(group = case_when(
    str_detect(Replicate, "QC") ~ "QC",
    Replicate %in% WT ~ "WT",
    Replicate %in% KO ~ "KO"),
  ) %>%
  group_by(Compound, Prod) %>%
 
  mutate(transition_area = median(Area)) %>%
  ungroup() %>%
  group_by(Compound) %>%
  mutate(transition_rank = frank(-transition_area, ties.method = "dense"),
         keep = transition_rank == 1) %>%
  filter(keep == 1) %>%
  select(Compound:RT, group)

samples <- data %>%
  filter(!str_detect(group, "QC")) %>%
  mutate(group = factor(group, levels = c("WT", "KO")))


save(samples, file = "sample_name.rdata")