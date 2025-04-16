#Load in required packages
library(tidyverse)        # Includes ggplot2, dplyr, purrr, stringr, etc.
library(ComplexHeatmap)   # For advanced heatmaps (includes circlize)
library(pROC)             # For ROC analysis
library(ggpubr)           # Publication-quality ggplot2 enhancements
library(ggrepel)          # Repel text labels in ggplot2
library(readxl)           # For reading Excel files
library(Biostrings)       # For biological string manipulation
library(ggExtra)

#Change the path to your directory
path <- "/Users/rnq676/Documents/thesis/KRAS - PIK3CG"

#Read in dataframes
data_crystal <- read.table(file.path(path, "prism_rosetta_XXX_1he8_mut.txt"), header = TRUE) %>%
  mutate(Pos = as.numeric(str_extract(variant, "\\d+"))-749,
         variant = str_replace(variant, "\\d+", function(x) as.numeric(x) - 749))

data_mave <- read_excel(file.path("/Users/rnq676/Documents/thesis/KRAS/mavedb_4obe.xlsx"), 
                        sheet = 2) %>%
  select(aa_seq, fitness, sigma, WT, assay) %>%
  filter(assay == "BindingPCA PIK3CGRBD")

data_2struct <- read_csv(file.path("/Users/rnq676/Documents/thesis/KRAS/2struct_4obe.csv")) %>% select(c(n,q3))

colnames(data_2struct)[1] <- "Pos"

SASA_data <- read_csv(file.path("/Users/rnq676/Documents/thesis/KRAS/SASA_4obe.csv"), col_names = TRUE)

colnames(SASA_data)[2] <- "Pos"
colnames(SASA_data)[8] <- "SASA"

#Extract wildtype sequences
wildType <- data_mave %>% filter(WT == "TRUE") %>% select(aa_seq)

#All variants except wildtype
mut <- data_mave %>% filter(WT!=TRUE)

#function to convert sequence into variant notation (WildtypePositionMutation)
changes <- function(mut, wildType) {
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

diff <- changes(mut, as.character(wildType[1,1]))

data <- diff %>%
  merge(data_crystal, by = "variant") %>%
  mutate(cutoff_mave = ifelse(fitness > -0.5, 1, 0))%>%
  filter(!str_detect(variant, "\\*|;"))

# Distribution MAVE score
ggplot(data, aes(x=fitness)) +
  geom_histogram(binwidth=.05, colour="black", fill="dodgerblue2", alpha=0.7)+
  geom_vline(xintercept = -0.5, color = "darkgrey", linetype="dashed",size=0.5)+ 
  labs(title = "Mave scores distribution for KRAS - PIK3CG", x = "MAVE Score", y = "Frequency") +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5))

#ROC curve
ggroc(list(roc_crystal = roc(data$cutoff_mave, data$mean_ddG)), 
      legacy.axes = TRUE, 
      linewidth = 1,
      aes = c("color", "linetype")) +
  labs(title="ROC curve for KRAS - PIK3CG",x = "False Positive Rate", y = "True Positive Rate") +
  scale_linetype_manual("", 
                        labels = c(bquote(crystal*Delta*Delta*G (AUC == .(round(auc(data$cutoff_mave, data$mean_ddG),2))))),
                        values = 1) +  # Both lines solid (no dashed line)
  scale_color_manual("", 
                     labels = c(bquote(crystal*Delta*Delta*G (AUC == .(round(auc(data$cutoff_mave, data$mean_ddG),2))))),
                     values = "darkorchid") +
  geom_abline(linetype = "solid", linewidth = 1) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        text = element_text(size = 14), 
        legend.text.align = 0, 
        legend.position = "bottom", 
        legend.title = element_blank(), 
        aspect.ratio = 1,
        plot.margin = margin(10, 10, 10, 10))+ theme(plot.title = element_text(hjust = 0.5))

#MAVE vs ddG
data <- data %>%
  mutate(Pos = as.numeric(str_extract(variant, "\\d+"))) %>%
  left_join(data_2struct, by = "Pos") %>%
  left_join(SASA_data %>% select(Pos, SASA), by = "Pos") %>%
  mutate(
    SASA = replace_na(SASA, "o"),  # Replace NA with "o"
    SASA = recode(SASA, "i" = "Core", "o" = "Surface"),  # Rename categories
    SASA = ifelse(Pos %in% c(3, 21, 24, 25, 29, 33, 36, 37, 38, 39, 40, 41, 63, 64, 70, 73), "Interaction", SASA),
    q3 = recode(q3, "C" = "Coil", "E" = "Strand", "H" = "Helix")
  )

df_long <- data %>%
  mutate(mutant_aa = substr(variant, nchar(variant), nchar(variant)))

# Merge correlation values back into the dataset
chemical_order <- c(
  "A", "V", "I", "L", "M", "F", "Y", "W", #hydrophobic (non-polar)
  "R", "H", "K", "D", "E", "S", "T", "N", "Q", "C", "G", "P" #hydrophilic (polar)
)

# Merge correlations and add title including core/surface info
df_long <- df_long %>%
  mutate(
    # Ensure proper order of amino acids based on chemistry
    prop = factor(mutant_aa, levels = chemical_order),
    Pos = as.numeric(str_extract(variant, "\\d+")))

cor_values <- df_long %>%
  group_by(mutant_aa) %>%
  summarise(R = round(cor(fitness, mean_ddG, use = "complete.obs"), 2), .groups = "drop")  # Keeps separate R for each amino acid

# Add correlation as a single title (not per mutant_aa)
df_long <- df_long %>%
  left_join(cor_values, by = "mutant_aa") %>%
  mutate(title = paste0(mutant_aa, " (R = ", R, ")"))  # Each facet gets its own R

ggplot(df_long, aes(x = fitness, y = mean_ddG, fill = SASA)) +
  geom_point(size = 1.5, shape = 21, color = "grey20", alpha = 0.7) +  # Fill controls point color
  scale_fill_manual(values = c("Core" = "darkorchid", "Surface" = "dodgerblue2", "Interaction" = "coral"), name="Solvent Accessible\nSurface Area") +  # Custom colors
  labs(x = "MAVE Score", y = "ddG", title = "MAVE Scores vs Crystal ddG for KRAS - PIK3CG") +
  
  # Facet by amino acid chemistry, ensuring correct order
  facet_wrap(~prop, scales = "free_y", labeller = labeller(prop = setNames(df_long$title, df_long$prop))) +
  
  geom_smooth(aes(group = 1), method = "lm", se = FALSE, color = "black", linetype = "dashed", size = 0.5) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

ggplot(df_long, aes(x = fitness, y = mean_ddG, fill = q3)) +
  geom_point(color="grey20",size = 1.5, shape = 21, alpha = 0.7) +  # Fill controls point color
  scale_fill_manual(values = c("Coil" = "darkorchid", "Helix" = "dodgerblue2", "Strand" = "coral"), name="Secondary Structure") +  # Custom colors
  labs(x = "MAVE Score", y = "ddG", title = "MAVE Scores vs Crystal ddG for KRAS - PIK3CG") +
  
  # Facet by amino acid chemistry, ensuring correct order
  facet_wrap(~prop, scales = "free_y", labeller = labeller(prop = setNames(df_long$title, df_long$prop))) +
  
  geom_smooth(aes(group = 1), method = "lm", se = FALSE, color = "black", linetype = "dashed", size = 0.5) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

#MAVE vs ddG with colors
ggplot(df_long, aes(x = fitness, y = mean_ddG)) +
  geom_point(aes(fill = SASA), size = 1.5, shape = 21, color = "grey20", alpha = 0.7) +  # Use `fill` for coloring by SASA
  scale_fill_manual(values = c("Core" = "darkorchid", "Surface" = "dodgerblue2", "Interaction" = "coral"), name="Secondary Structure") +  # Custom colors for Core, Surface, and Interaction
  labs(x = "MAVE Score", y = "ddG", title = "MAVE Scores vs Crystal ddG for KRAS - PIK3CG") +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed", size = 0.5) +
  stat_cor(aes(label = paste(..r.label..)), method = "pearson", size = 4) +  # Remove the `group` argument to avoid grouping issues
  theme_bw() + theme(plot.title = element_text(hjust = 0.5))

ggplot(df_long, aes(x = fitness, y = mean_ddG)) +
  geom_point(aes(fill = q3), size = 1.5, shape = 21, color = "grey20", alpha = 0.7) +  # Use `fill` for coloring by SASA
  scale_fill_manual(values = c("Coil" = "darkorchid", "Helix" = "dodgerblue2", "Strand" = "coral"), name="Solvent Accessible\nSurface Area") + 
  labs(x = "MAVE Score", y = "ddG", title = "MAVE Scores vs Crystal ddG for KRAS - PIK3CG") +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed", size = 0.5) +
  stat_cor(aes(label = paste(..r.label..)), method = "pearson", size = 4) +  # Remove the `group` argument to avoid grouping issues
  theme_bw() + theme(plot.title = element_text(hjust = 0.5))

#MAVE abundance vs MAVE interaction

data_mave_merge <- read_excel(file.path("/Users/rnq676/Documents/thesis/KRAS/mavedb_4obe.xlsx"), 
                              sheet = 2) %>%
  select(aa_seq, fitness, sigma, WT, assay) %>%
  filter(assay == "AbundancePCA") %>%
  filter(WT!=TRUE) %>%
  distinct(aa_seq, .keep_all = TRUE) %>%
  inner_join(diff, by = "aa_seq")%>%
  mutate(Pos = as.numeric(str_extract(variant, "\\d+"))) %>%
  left_join(data_2struct, by = "Pos") %>%
  left_join(SASA_data %>% select(Pos, SASA), by = "Pos") %>%
  filter(Pos<170)%>%
  mutate(
    SASA = replace_na(SASA, "o"),  # Replace NA with "o"
    SASA = recode(SASA, "i" = "Core", "o" = "Surface"),  # Rename categories
    SASA = ifelse(Pos %in% c(3, 21, 24, 25, 29, 33, 36, 37, 38, 39, 40, 41, 63, 64, 70, 73), "Interaction", SASA),
    q3 = recode(q3, "C" = "Coil", "E" = "Strand", "H" = "Helix")
  ) %>%
  filter(!str_detect(variant, "\\*|;"))

p <- ggplot(data_mave_merge, aes(x = fitness.x, y = fitness.y)) +
  geom_point(size = 1.5, shape = 21, fill = "dodgerblue2", color = "grey20", alpha = 0.5) +  
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "black", linetype = "dashed", size = 0.5) +
  labs(x = "Abundace MAVE Score", y = "Interaction MAVE Score") +
  annotate("text", x = 0.75, y = 0.1, 
           label = paste("R =", round(cor(data_mave_merge$fitness.y, data_mave_merge$fitness.x), 2)), 
           size = 4) +
  theme_bw() +
  coord_fixed(ratio = 1) +
  xlim(-1.5, 1) +  # Set same limits for both axes
  ylim(-1.5, 0.5) 

p_marginal <- ggMarginal(p, type = "density", fill = "dodgerblue2", color= "dodgerblue2",alpha = 0.5, margins = "both", size = 5)

title <- ggplot() +
  ggtitle("MAVE Scores from Abundance vs Interaction for KRAS - PIK3CG") +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5, size = 14))

grid.arrange(title, p_marginal, ncol = 1, heights = c(0.1, 1))

#colored with secondary structure
base_plot <- ggplot(data_mave_merge, aes(fitness.x, fitness.y)) +
  geom_point(aes(color=q3), size=1.5, alpha=0.7) +
  scale_colour_manual(values = c("Coil" = "darkorchid", "Helix" = "dodgerblue2", "Strand" = "coral"), name="Secondary Structure") +
  labs(
    x = "Abundance MAVE Score",
    y = "Interaction MAVE Score") +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed", size = 0.5) +
  stat_cor(aes(label = paste(..r.label..)), method = "pearson", size = 4) + 
  coord_fixed(ratio = 1)+
  theme_bw()

p_marginal <- ggMarginal(base_plot, groupColour = TRUE, groupFill = TRUE)

title <- ggdraw() +
  draw_label(
    "MAVE Scores from Abundance vs Interaction for KRAS - PIK3CG", 
    hjust = 0.5
  )

plot_grid(title, p_marginal, ncol = 1, rel_heights = c(0.1, 1))

#colored with SASA

base_plot <- ggplot(data_mave_merge, aes(fitness.x, fitness.y)) +
  geom_point(aes(color=SASA), size=1.5, alpha=0.7) +
  scale_colour_manual(values = c("Core" = "darkorchid", "Surface" = "dodgerblue2", "Interaction" = "coral"), name="Solvent Accessible\nSurface Area" ) +
  labs(
    x = "Abundance MAVE Score",
    y = "Interaction MAVE Score") +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed", size = 0.5) +
  stat_cor(aes(label = paste(..r.label..)), method = "pearson", size = 4) + 
  coord_fixed(ratio = 1)+
  theme_bw()

p_marginal <- ggMarginal(base_plot, groupColour = TRUE, groupFill = TRUE)

title <- ggdraw() +
  draw_label(
    "MAVE Scores from Abundance vs Interaction for KRAS - PIK3CG", 
    hjust = 0.5
  )

plot_grid(title, p_marginal, ncol = 1, rel_heights = c(0.1, 1))

