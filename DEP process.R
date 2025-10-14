library(data.table)
library(dplyr)
library(tidyr)
library(DEP)
library(queryup)
library(readxl)
library(stringr)

res <- fread("Spectronaut report.tsv")

res$PG.Quantity[which(is.na(res$PG.Quantity))] <- NA

#change the name of the samples
name_data <- res %>%
  dplyr::select(R.Condition, R.FileName) %>%
  distinct() %>%
  group_by(R.Condition) %>%
  mutate(Sample = paste(R.Condition, row_number(), sep = "_"))

#add new names to the results
data <- left_join(res, name_data, by = "R.FileName") %>%
  dplyr::select(Sample, Protein = PG.ProteinAccessions,
                Quantity = PG.Quantity) %>%
  spread(Sample, value = Quantity)


#remove duplicated uniprot accessions
Reviewed_protein <- query_uniprot(list("organism_id" = c("10116"), "reviewed" = "true")) 

uni_protein_names <- data %>%
  dplyr::select(Protein) %>%
  mutate(sep_Protein = Protein) %>%
  separate_rows(sep_Protein) %>%
  left_join(Reviewed_protein, by = c("sep_Protein" = "Entry")) %>%
  group_by(Protein) %>%
  mutate(new_Protein = case_when(
    n() == 1 ~ sep_Protein[1],
    any(Reviewed == "reviewed") ~ sep_Protein[which(Reviewed == "reviewed")[1]],
    TRUE ~ sep_Protein[1]
  )) %>%
  dplyr::select(Protein, new_Protein) %>%
  distinct() 

format_data <- res %>%
  dplyr::select(-R.Condition, -PG.IBAQ, -PG.ProteinDescriptions) %>%
  right_join(uni_protein_names, by = c("PG.ProteinAccessions" = "Protein")) %>%
  spread(R.FileName, PG.Quantity)

data_unique <- make_unique(format_data, "Gene", "new_Protein", delim = ";") %>%
  as.data.frame()

# Are there any duplicated names?
data_unique$name %>% duplicated() %>% any()

# Generate a SummarizedExperiment object using an experimental design
# get quant column numbers;fill in according to the column names "09_m"
quant_columns <- grep("", colnames(data_unique)) 

experimental_design <- name_data %>%
  mutate(condition = R.Condition,
         replicate = str_sub(Sample, -1, -1),
         label = R.FileName)

#remove empty rows (proteins identified in 0 samples)
data_clean <- data_unique %>%
  mutate(empty = rowSums(.[quant_columns], na.rm = T)) %>%
  filter(empty > 0) %>%
  dplyr::select(-empty)

data_se <- make_se(data_clean, quant_columns, experimental_design)

# Plot a barplot of the protein identification overlap between samples
plot_frequency(data_se)


# Filter for proteins identified in most replicates within each condition
data_filt <- filter_missval(data_se, thr = 1)
plot_numbers(data_filt)

# Plot a barplot of the protein identification overlap between samples
plot_coverage(data_filt)

# Normalize the data
data_norm <- normalize_vsn(data_filt)
meanSdPlot(data_norm)
meanSdPlot(data_filt)
# Visualize normalization by boxplots for all samples before and after normalization
plot_normalization(data_filt, data_norm)

# Plot a heatmap of proteins with missing values
plot_missval(data_filt)

# Plot intensity distributions and cumulative fraction of proteins with and without missing values
plot_detect(data_filt)


# Impute missing data using random draws from a manually defined left-shifted Gaussian distribution 
data_imp <- impute(data_norm, fun = "man", shift = 1.8, scale = 0.3)
imp_data <- as.data.frame(assay(data_imp))

#Plot intensity distributions before and after imputation
plot_imputation(data_norm, data_imp)

# Test every sample versus control
data_diff <- test_diff(data_imp, type = "control", control = "WT")


# Denote significant proteins based on user defined cutoffs
dep <- add_rejections(data_diff, alpha = 0.05, lfc = log2(1.5))


# Plot the first and second principal components
plot_pca(dep)

# Plot the Pearson correlation matrix
plot_cor(dep, significant = TRUE, lower = 0, upper = 1)

# Plot a heatmap of all significant proteins with the data centered per protein
plot_heatmap(dep, type = "centered", kmeans = TRUE, 
             k = 6, show_row_names = T, indicate = c("condition"))


# Plot a volcano plot for the contrast 
plot_volcano(dep, contrast = "KO_vs_WT", add_names = T)

diff <- as.data.frame(rowData(data_diff)) %>%
  dplyr::select(new_Protein:Entrez, contains("KO"))

