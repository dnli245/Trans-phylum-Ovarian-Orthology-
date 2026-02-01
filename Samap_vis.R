# Load required libraries
library(dplyr)
library(networkD3)
library(ggsankeyfier)
library(ggplot2)
library(ggalluvial)

# Read the mapping table and set row names
# whole: 'hmsp_no_ribosomal_no_stromal_w_germ_8.14/SAMap_results_8.15_final/mappingtables/'
# germ: 'hm_germ_cell_experiment_8.14/SAMap_results_final/mappingtables/'
folder <- 'hmsp_no_ribosomal_no_stromal_w_germ_8.14/SAMap_results_8.15_final/mappingtables/'
mappingtable_whole <- read.csv(paste(folder, 'hmsp_whole_mappingtable_preannotated.csv', sep = ''))
mappingtable_7 <- read.csv(paste(folder, 'hm_7_sp_mappingtable_preannotated.csv', sep = ''))
mappingtable_9 <- read.csv(paste(folder, 'hm_9_sp_mappingtable_preannotated.csv', sep = ''))
mappingtable_10 <- read.csv(paste(folder, 'hm_10_sp_mappingtable_preannotated.csv', sep = ''))
mappingtable_13 <- read.csv(paste(folder, 'hm_13_sp_mappingtable_preannotated.csv', sep = ''))
mappingtable_16 <- read.csv(paste(folder, 'hm_16_sp_mappingtable_preannotated.csv', sep = ''))
mappingtable_adult <- read.csv(paste(folder, 'hm_adult_sp_mappingtable_preannotated.csv', sep = ''))
mappingtable_fetal <- read.csv(paste(folder, 'hm_fetal_sp_mappingtable_preannotated.csv', sep = ''))



mp_list <- list(
  hmsp_whole = mappingtable_whole,
  hm7_sp = mappingtable_7,
  hm9_sp = mappingtable_9,
  hm10_sp = mappingtable_10,
  hm13_sp = mappingtable_13,
  hm16_sp = mappingtable_16,
  hmadult_sp = mappingtable_adult,
  hmfetal_sp = mappingtable_fetal
)


# Define function to clean and prepare Sankey data
clean_sankey_data <- function(mappingtable, threshold = 0.5) {
  # Apply threshold to filter out low values
  mappingtable[mappingtable < threshold] <- 0
  
  # Remove empty rows and columns
  mappingtable <- mappingtable[rowSums(mappingtable) > 0, colSums(mappingtable) > 0]
  
  # Extract non-zero elements (source, target, value)
  non_zero_indices <- which(mappingtable > 0, arr.ind = TRUE)
  long_df <- data.frame(
    Source = rownames(mappingtable)[non_zero_indices[, 1]],
    Target = colnames(mappingtable)[non_zero_indices[, 2]],
    Value = mappingtable[non_zero_indices]
  )
  
  return(long_df)
}

rename_hm_to_hs <- function(df) {
  colnames(df) <- gsub("^hm_", "hs_", colnames(df))
  return(df)
}

for (mp_name in names(mp_list)) {
  mappingtable <- mp_list[[mp_name]]
  mappingtable <- rename_hm_to_hs(mappingtable)
  rownames(mappingtable) <- colnames(mappingtable)
  # Apply function and filter half the rows
  sankey_data <- clean_sankey_data(mappingtable, threshold = 0.1)
  sankey_data <- sankey_data[1:(nrow(sankey_data) / 2), ]
  
  # Rename columns and add category
  colnames(sankey_data) <- c("sp", "lv", "value")
  sankey_data$category <- rep("lvsp", nrow(sankey_data))
  
  # Plot Sankey diagram using ggalluvial-style ggplot
  san_plot <- ggplot(sankey_data, aes(axis1 = sp, axis2 = lv, y = value)) +
    geom_alluvium(aes(fill = value), alpha = 0.7) +
    geom_stratum() +
    geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
    scale_fill_gradient(low = "lightblue", high = "darkblue") +
    theme_minimal() +
    labs(title = "Sankey Diagram using ggalluvial")
  
  
  pdf(file = paste(folder, mp_name, "_sankey.pdf", sep = ""),width=10,height=24)
  print(san_plot)
  dev.off()
}