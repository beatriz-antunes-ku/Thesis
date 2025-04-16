#Load in required packages
library(tidyverse)        # Includes ggplot2, dplyr, purrr, stringr, etc.
library(ComplexHeatmap)   # For advanced heatmaps (includes circlize)
library(pROC)             # For ROC analysis
library(ggpubr)           # Publication-quality ggplot2 enhancements
library(ggrepel)          # Repel text labels in ggplot2
library(readxl)           # For reading Excel files
library(Biostrings)       # For biological string manipulation
library(ggExtra)
library(dplyr)

#Read in dataframes
mave_abundance <- read_excel(file.path("/Users/rnq676/Documents/thesis/KRAS/mavedb_4obe.xlsx"), 
                        sheet = 2) %>%
  select(aa_seq, fitness, sigma, WT, assay) %>%
  filter(assay == "AbundancePCA")

mave_interaction <- read_excel(file.path("/Users/rnq676/Documents/thesis/KRAS/mavedb_4obe.xlsx"), 
                        sheet = 2) %>%
  select(aa_seq, fitness, sigma, WT, assay) %>%
  filter(assay == "BindingPCA RAF1RBD")

data_2struct <- read_csv(file.path("/Users/rnq676/Documents/thesis/KRAS/0000_Sequence.csv")) %>% select(c(n,q3))

colnames(data_2struct)[1] <- "Pos"

SASA_data <- read_csv(file.path("/Users/rnq676/Documents/thesis/KRAS/SASA_4obe.csv"), col_names = TRUE)

colnames(SASA_data)[2] <- "Pos"
colnames(SASA_data)[8] <- "SASA"

chemical_order <- c(
  "A", "V", "I", "L", "M", "F", "Y", "W", #hydrophobic (non-polar)
  "R", "H", "K", "D", "E", "S", "T", "N", "Q", "C", "G", "P" #hydrophilic (polar)
)

wildType_abundance <- mave_abundance %>% filter(WT == "TRUE") %>% select(aa_seq)
wildType_interaction <- mave_interaction %>% filter(WT == "TRUE") %>% select(aa_seq)

mut_abundance <- mave_abundance %>% filter(WT!=TRUE)
mut_interaction <- mave_interaction %>% filter(WT!=TRUE)

changes <- function(mut, wildType) {
  # Ensure wildType is a character
  wildSplit <- strsplit(as.character(wildType), "")[[1]]  # Split the wild type sequence
  
  # Iterate through all rows and compare each mutated sequence
  mut <- mut %>%
    rowwise() %>%
    mutate(variant = {
      mutSeq <- as.character(aa_seq)  # Get mutated sequence for the current row
      mutSplit <- strsplit(mutSeq, "")[[1]]  # Split the mutated sequence into characters
      
      # Compare wild-type and mutated sequence positions
      changes <- map_chr(seq_along(wildSplit), function(i) {
        if (i <= length(mutSplit) && wildSplit[i] != mutSplit[i]) {
          paste0(wildSplit[i], i+1, mutSplit[i])  # Format as "WT Position Mut"
        } else {
          NA_character_  # Keep NA for unchanged positions
        }
      })
      
      changes <- na.omit(changes)  # Remove NAs (unchanged positions)
      
      # Combine multiple changes into a single string
      if (length(changes) > 0) {
        paste(changes, collapse = ";")  # Concatenate all mutations for a single row
      } else {
        NA_character_  # If no changes, set as NA
      }
    }) %>%
    ungroup()  # Restore dataframe structure after row-wise operation
  
  return(mut)
}

diff_abundance <- changes(mut_abundance, as.character(wildType_abundance[1,1]))
diff_interaction <- changes(mut_interaction, as.character(wildType_interaction[1,1]))

path <- "/Users/rnq676/Documents/thesis/KRAS"

data_esm_1 <- read_csv(file.path(path,"ESM/esm1v_1.csv")) %>% select(c(-1))

data_esm_1$Pos <- rownames(data_esm_1)

data_long_1 <- melt(data_esm_1, id.vars = "Pos")

data_esm_2 <- read_csv(file.path(path,"ESM/esm1v_2.csv")) %>% select(c(-1))

data_esm_2$Pos <- rownames(data_esm_2)

data_long_2 <- melt(data_esm_2, id.vars = "Pos")

data_esm_3 <- read_csv(file.path(path,"ESM/esm1v_3.csv")) %>% select(c(-1))

data_esm_3$Pos <- rownames(data_esm_3)

data_long_3 <- melt(data_esm_3, id.vars = "Pos")

data_esm_4 <- read_csv(file.path(path,"ESM/esm1v_4.csv")) %>% select(c(-1))

data_esm_4$Pos <- rownames(data_esm_4)

data_long_4 <- melt(data_esm_4, id.vars = "Pos")

data_esm_5 <- read_csv(file.path(path,"ESM/esm1v_5.csv")) %>% select(c(-1))

data_esm_5$Pos <- rownames(data_esm_5)

data_long_5 <- melt(data_esm_5, id.vars = "Pos")

merged_esm <- merge(data_long_1, data_long_2, by = c("Pos", "variable"), suffixes = c("_1", "_2"))

merged_esm <- merge(merged_esm, data_long_3, by = c("Pos", "variable"), suffixes = c("_3"))

merged_esm <- merge(merged_esm, data_long_4, by = c("Pos", "variable"), suffixes = c("_3", "_4"))

merged_esm <- merge(merged_esm, data_long_5, by = c("Pos", "variable"), suffixes = c("_4", "_5"))

colnames(merged_esm)[7] <- "value_5"

merged_esm <- merged_esm %>%
  mutate(mean_value = rowMeans(select(., starts_with("value_")), na.rm = TRUE)) %>%
  arrange(Pos)

ggplot(merged_esm, aes(x = variable, y = Pos, fill = mean_value)) +
  geom_tile(color = "gray80", linewidth = 0.3) +  # Add lines between tiles
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, name="Score Range") +
  labs(title = "Mutation Effect Scores by Position and Amino Acid",
       x = "Amino Acid", y = "Sequence Position") +
  theme_minimal() 

#Abundance MAVE vs ESM1v
mave_abundance_esm <- diff_abundance %>%
  filter(!str_detect(variant, "\\*|;")) %>%
  mutate(wild_type = substr(variant, 1, 1),  # Extract the wild-type amino acid
         Pos = as.integer(substr(variant, 2, nchar(variant) - 1))+1,  # Extract the position (integer)
         mutant = substr(variant, nchar(variant), nchar(variant))) %>%
  left_join(data_2struct, by = "Pos") %>%
  left_join(SASA_data %>% select(Pos, SASA), by = "Pos") %>%
  mutate(
    SASA = replace_na(SASA, "o"),  # Replace NA with "o"
    SASA = recode(SASA, "i" = "Core", "o" = "Surface"),  # Rename categories
    q3 = recode(q3, "C" = "Coil", "E" = "Strand", "H" = "Helix")
  )

merged_abundance_esm <- merge(mave_abundance_esm, merged_esm, by.x = c("Pos", "mutant"), by.y = c("Pos", "variable")) %>%
  mutate(prop = factor(mutant, levels = chemical_order))

ggplot(merged_abundance_esm, aes(x = fitness, y = mean_value)) +
  geom_point(color = "dodgerblue2", alpha = 0.5, size=1.5) +  # This adds points to the plot
  labs(
    title = "MAVE Abundance score comparison to ESM1v",
    x = "MAVE Score",
    y = "ESM1v Score"
  ) +
  geom_smooth(aes(group = 1), method = "lm", se = FALSE, color = "black", linetype = "dashed", size = 0.5) +
  stat_cor(aes(label = paste(..r.label..)), method = "pearson", size = 4) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

ggplot(merged_abundance_esm, aes(x = fitness, y = mean_value)) +
  geom_point(aes(fill = q3), size = 1.5, shape = 21, color = "grey20", alpha = 0.7) +  # Use `fill` for coloring by `In/Out`
  scale_fill_manual(values = c("Coil" = "darkorchid", "Helix" = "dodgerblue2", "Strand" = "coral"), name="Secondary Structure") + 
  labs(x = "MAVE Abundance Score", y = "ESM1v Score", title = "MAVE Abundance score comparison to ESM1v") +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed", size = 0.5) +
  stat_cor(aes(label = paste(..r.label..)), method = "pearson", size = 4) +  # Remove the `group` argument to avoid grouping issues
  theme_bw() + theme(plot.title = element_text(hjust = 0.5))

ggplot(merged_abundance_esm, aes(x = fitness, y = mean_value)) +
  geom_point(aes(fill = SASA), size = 1.5, shape = 21, color = "grey20", alpha = 0.7) +  # Use `fill` for coloring by `In/Out`
  scale_fill_manual(values = c("Core" = "darkorchid", "Surface" = "dodgerblue2"), name="Solvent Accessible\nSurface Area") + 
  labs(x = "MAVE Abundance Score", y = "ESM1v Score", title = "MAVE Abundance score comparison to ESM1v") +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed", size = 0.5) +
  stat_cor(aes(label = paste(..r.label..)), method = "pearson", size = 4) +  # Remove the `group` argument to avoid grouping issues
  theme_bw() + theme(plot.title = element_text(hjust = 0.5))

ggplot(merged_abundance_esm, aes(x = fitness, y = mean_value, fill = q3)) +
  geom_point(size = 1.5, shape = 21, color = "grey20", alpha = 0.7) +  # Fill controls point color
  scale_fill_manual(values = c("Coil" = "darkorchid", "Helix" = "dodgerblue2", "Strand" = "coral"), name="Secondary Structure") + 
  labs(x = "MAVE Abundance Score", y = "ESM1v Score", title = "MAVE Abundance Score comparison to ESM1v") +
  
  # Facet by amino acid chemistry, ensuring correct order
  facet_wrap(~prop, scales = "free_y", labeller = labeller(prop = merged_abundance_esm$prop)) +
  
  geom_smooth(aes(group = 1), method = "lm", se = FALSE, color = "black", linetype = "dashed", size = 0.5) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

ggplot(merged_abundance_esm, aes(x = fitness, y = mean_value, fill = SASA)) +
  geom_point(size = 1.5, shape = 21, color = "grey20", alpha = 0.7) +  # Fill controls point color
  scale_fill_manual(values = c("Core" = "darkorchid", "Surface" = "dodgerblue2"), name="Solvent Accessible\nSurface Area") + 
  labs(x = "MAVE Abundance Score", y = "ESM1v Score", title = "MAVE Abundance Score comparison to ESM1v") +
  
  # Facet by amino acid chemistry, ensuring correct order
  facet_wrap(~prop, scales = "free_y", labeller = labeller(prop = merged_abundance_esm$prop)) +
  
  geom_smooth(aes(group = 1), method = "lm", se = FALSE, color = "black", linetype = "dashed", size = 0.5) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

#Interaction MAVE vs ESM1v
mave_interaction_esm <- diff_interaction %>%
  filter(!str_detect(variant, "\\*|;")) %>%
  mutate(wild_type = substr(variant, 1, 1),  # Extract the wild-type amino acid
         Pos = as.integer(substr(variant, 2, nchar(variant) - 1))+1,  # Extract the position (integer)
         mutant = substr(variant, nchar(variant), nchar(variant))) %>%
  left_join(data_2struct, by = "Pos") %>%
  left_join(SASA_data %>% select(Pos, SASA), by = "Pos") %>%
  mutate(
    SASA = replace_na(SASA, "o"),  # Replace NA with "o"
    SASA = recode(SASA, "i" = "Core", "o" = "Surface"),  # Rename categories
    q3 = recode(q3, "C" = "Coil", "E" = "Strand", "H" = "Helix")
  )

merged_interaction_esm <- merge(mave_interaction_esm, merged_esm, by.x = c("Pos", "mutant"), by.y = c("Pos", "variable")) %>%
  mutate(prop = factor(mutant, levels = chemical_order))

ggplot(merged_interaction_esm, aes(x = fitness, y = mean_value)) +
  geom_point(color = "dodgerblue2", alpha = 0.5, size=1.5) +  # This adds points to the plot
  labs(
    title = "MAVE Interaction Score from RAF1 comparison to ESM1v",
    x = "MAVE Interaction Score",
    y = "ESM1v Score"
  ) +
  geom_smooth(aes(group = 1), method = "lm", se = FALSE, color = "black", linetype = "dashed", size = 0.5) +
  stat_cor(aes(label = paste(..r.label..)), method = "pearson", size = 4) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

ggplot(merged_interaction_esm, aes(x = fitness, y = mean_value)) +
  geom_point(aes(fill = q3), size = 1.5, shape = 21, color = "grey20", alpha = 0.7) +  # Use `fill` for coloring by `In/Out`
  scale_fill_manual(values = c("Coil" = "darkorchid", "Helix" = "dodgerblue2", "Strand" = "coral"), name="Secondary Structure") + 
  labs(x = "MAVE Interaction Score", y = "ESM1v Score", title = "MAVE Interaction Score from RAF1 comparison to ESM1v") +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed", size = 0.5) +
  stat_cor(aes(label = paste(..r.label..)), method = "pearson", size = 4) +  # Remove the `group` argument to avoid grouping issues
  theme_bw() + theme(plot.title = element_text(hjust = 0.5))

ggplot(merged_interaction_esm, aes(x = fitness, y = mean_value)) +
  geom_point(aes(fill = SASA), size = 1.5, shape = 21, color = "grey20", alpha = 0.7) +  # Use `fill` for coloring by `In/Out`
  scale_fill_manual(values = c("Core" = "darkorchid", "Surface" = "dodgerblue2"), name="Solvent Accessible\nSurface Area") + 
  labs(x = "MAVE Interaction Score", y = "ESM1v Score", title = "MAVE Interaction Score from RAF1 comparison to ESM1v") +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed", size = 0.5) +
  stat_cor(aes(label = paste(..r.label..)), method = "pearson", size = 4) +  # Remove the `group` argument to avoid grouping issues
  theme_bw() + theme(plot.title = element_text(hjust = 0.5))

ggplot(merged_interaction_esm, aes(x = fitness, y = mean_value, fill = q3)) +
  geom_point(size = 1.5, shape = 21, color = "grey20", alpha = 0.7) +  # Fill controls point color
  scale_fill_manual(values = c("Coil" = "darkorchid", "Helix" = "dodgerblue2", "Strand" = "coral"), name="Secondary Structure") + 
  labs(x = "MAVE Interaction Score", y = "ESM1v Score", title = "MAVE Interaction Score from RAF1 comparison to ESM1v") +
  
  # Facet by amino acid chemistry, ensuring correct order
  facet_wrap(~prop, scales = "free_y", labeller = labeller(prop = merged_interaction_esm$prop)) +
  
  geom_smooth(aes(group = 1), method = "lm", se = FALSE, color = "black", linetype = "dashed", size = 0.5) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

ggplot(merged_interaction_esm, aes(x = fitness, y = mean_value, fill = SASA)) +
  geom_point(size = 1.5, shape = 21, color = "grey20", alpha = 0.7) +  # Fill controls point color
  scale_fill_manual(values = c("Core" = "darkorchid", "Surface" = "dodgerblue2"), name="Solvent Accessible\nSurface Area") + 
  labs(x = "MAVE Interaction Score", y = "ESM1v Score", title = "MAVE Interaction Score from RAF1 comparison to ESM1v") +
  
  # Facet by amino acid chemistry, ensuring correct order
  facet_wrap(~prop, scales = "free_y", labeller = labeller(prop = merged_interaction_esm$prop)) +
  
  geom_smooth(aes(group = 1), method = "lm", se = FALSE, color = "black", linetype = "dashed", size = 0.5) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

