#Load in libraries
library(tidyverse)        # Includes ggplot2, dplyr, purrr, stringr, etc.
library(ComplexHeatmap)   # For advanced heatmaps (includes circlize)
library(pROC)             # For ROC analysis
library(ggpubr)           # Publication-quality ggplot2 enhancements
library(ggrepel)          # Repel text labels in ggplot2
library(readxl)           # For reading Excel files
library(Biostrings)       # For biological string manipulation
library(ggExtra)          # Add marginal histograms or boxplots to ggplot2
library(gridExtra)        # Arrange multiple grid-based plots
library(cowplot)          # Streamlined plot themes and plot composition

#Path to data
path <- "/Users/rnq676/Documents/thesis/KRAS - RAL"

#Read in data
data_crystal <- read.table(file.path(path, "prism_rosetta_XXX_1lfd_mut.txt"), header = TRUE) %>%
  mutate(Pos = as.numeric(str_extract(variant, "\\d+"))-87,
         variant = str_replace(variant, "\\d+", function(x) as.numeric(x) - 87))

data_alphafold <- read.table(file.path(path, "prism_rosetta_XXX_fold_1lfd.txt"), header = TRUE) %>%
  mutate(Pos = as.numeric(str_extract(variant, "\\d+"))-87,
         variant = str_replace(variant, "\\d+", function(x) as.numeric(x) - 87))

data_mave <- read_excel(file.path("/Users/rnq676/Documents/thesis/KRAS/mavedb_4obe.xlsx"), 
                        sheet = 2) %>%
  select(aa_seq, fitness, sigma, WT, assay) %>%
  filter(assay == "BindingPCA RALGDSRBD")

data_2struct <- read_csv(file.path("/Users/rnq676/Documents/thesis/KRAS/2struct_4obe.csv")) %>% select(c(n,q3))

colnames(data_2struct)[1] <- "Pos"

SASA_data <- read_csv(file.path("/Users/rnq676/Documents/thesis/KRAS/SASA_4obe.csv"), col_names = TRUE)

colnames(SASA_data)[2] <- "Pos"
colnames(SASA_data)[8] <- "SASA"

#Merge crystal and alphafold data in the same df
data_ddG <- data_crystal %>%
  merge(data_alphafold, by = "variant", suffixes = c(".crystal", ".alphafold"))

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
          paste0(wildSplit[i], i+1, mutSplit[i])
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

#Define cutoff for ROC curve and exclude double mutations
data <- diff %>%
  merge(data_ddG, by = "variant") %>%
  mutate(cutoff_mave = ifelse(fitness > -0.5, 1, 0)) %>%
  filter(!str_detect(variant, "\\*|;"))

# Distribution MAVE score
ggplot(data, aes(x=fitness)) +
  geom_histogram(binwidth=.05, colour="black", fill="dodgerblue2", alpha=0.7)+
  geom_vline(xintercept = -0.5, color = "darkgrey", linetype="dashed",size=0.5)+ 
  labs(title = "Mave scores distribution for KRAS - RAL", x = "MAVE Score", y = "Frequency") +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5))

#ROC curve
ggroc(list(roc_crystal = roc(data$cutoff_mave, data$mean_ddG.crystal), roc_alphafold = roc(data$cutoff_mave, data$mean_ddG.alphafold)), 
      legacy.axes = TRUE, 
      linewidth = 1,
      aes = c("color", "linetype")) +
  labs(title="ROC curve for KRAS - RAL",x = "False Positive Rate", y = "True Positive Rate") +
  scale_linetype_manual("", 
                        labels = c(bquote(crystal*Delta*Delta*G (AUC == .(round(auc(data$cutoff_mave, data$mean_ddG.crystal),2)))),
                                   bquote(alphafold*Delta*Delta*G (AUC == .(round(auc(data$cutoff_mave, data$mean_ddG.alphafold),2))))),
                        values = c(1, 1)) +
  scale_color_manual("", 
                     labels = c(bquote(crystal*Delta*Delta*G (AUC == .(round(auc(data$cutoff_mave, data$mean_ddG.crystal),2)))),
                                bquote(alphafold*Delta*Delta*G (AUC == .(round(auc(data$cutoff_mave, data$mean_ddG.alphafold),2))))),
                     values = c("darkorchid", "dodgerblue2")) +
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

#ddG crystal vs ddG alphafold (with std errors and distributions)
p <- ggplot(data, aes(x = mean_ddG.crystal, y = mean_ddG.alphafold)) +
  geom_errorbarh(aes(xmin = mean_ddG.crystal - std_ddG.crystal, 
                     xmax = mean_ddG.crystal + std_ddG.crystal), 
                 height = 0, color = "grey20", alpha = 0.5) +
  geom_errorbar(aes(ymin = mean_ddG.alphafold - std_ddG.alphafold, 
                    ymax = mean_ddG.alphafold + std_ddG.alphafold), 
                width = 0, color = "grey20", alpha = 0.5) +
  geom_point(size = 1.5, shape = 21, fill = "dodgerblue2", color = "grey20", alpha = 0.5) +  
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "black", linetype = "dashed", size = 0.5) +
  labs(x = "Crystal ddG", y = "Alphafold ddG") +
  stat_cor(aes(label = paste(..r.label..)), method = "pearson", size = 4) +
  theme_bw() +
  coord_fixed(ratio = 1) 

p_marginal <- ggMarginal(p, type = "density", fill = "dodgerblue2",color="dodgerblue2", alpha = 0.5, margins = "both", size = 5)

title <- ggplot() +
  ggtitle("Crystal Structure ddG vs Alphafold Prediction ddG for KRAS - RAL") +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5, size = 14))

grid.arrange(title, p_marginal, ncol = 1, heights = c(0.1, 1))

#ddG crystal vs ddG alphafold (with std errors and distributions) label in outliers
lm_model <- lm(mean_ddG.alphafold ~ mean_ddG.crystal, data = data)

data$residuals <- residuals(lm_model)

residual_threshold <- 6 * sd(data$residuals)
data$outlier <- abs(data$residuals) > residual_threshold

outliers_data <- data %>% filter(outlier == TRUE)

p <- ggplot(data, aes(x = mean_ddG.crystal, y = mean_ddG.alphafold)) +
  geom_errorbarh(aes(xmin = mean_ddG.crystal - std_ddG.crystal, 
                     xmax = mean_ddG.crystal + std_ddG.crystal), 
                 height = 0, color = "grey20", alpha = 0.5) +
  geom_errorbar(aes(ymin = mean_ddG.alphafold - std_ddG.alphafold, 
                    ymax = mean_ddG.alphafold + std_ddG.alphafold), 
                width = 0, color = "grey20", alpha = 0.5) +
  geom_point(size = 1, shape = 21, fill = "dodgerblue2", color = "grey20", alpha = 0.5) +  
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "black", linetype = "dashed", size = 0.5) +
  labs(x = "Crystal ddG", y = "Alphafold ddG") +
  stat_cor(aes(label = paste(..r.label..)), method = "pearson", size = 4) +
  theme_bw() +
  coord_fixed(ratio = 1)

p_outliers <- p + geom_text_repel(data = outliers_data, aes(label = variant), 
                                  size = 3, color = "black", box.padding = 0.5, max.overlaps = 5)

p_marginal <- ggMarginal(p_outliers, type = "density", fill = "dodgerblue2",color="dodgerblue2", alpha = 0.5, margins = "both", size = 5)

title <- ggplot() +
  ggtitle("Crystal Structure ddG vs Alphafold Prediction ddG for KRAS - RAL") +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5, size = 14))

grid.arrange(title, p_marginal, ncol = 1, heights = c(0.1, 1))

data <- data %>%
  mutate(Pos = as.numeric(str_extract(variant, "\\d+"))) %>%
  left_join(data_2struct, by = "Pos") %>%
  left_join(SASA_data %>% select(Pos, SASA), by = "Pos") %>%
  mutate(
    SASA = replace_na(SASA, "o"),  # Replace NA with "o"
    SASA = recode(SASA, "i" = "Core", "o" = "Surface"),  # Rename categories
    SASA = ifelse(Pos %in% c(24, 25, 31, 33, 36, 37, 38, 39, 40, 41, 56, 64, 67), "Interaction", SASA),
    q3 = recode(q3, "C" = "Coil", "E" = "Strand", "H" = "Helix")
  )

#ddG crystal vs ddG alphafold (with distributions) colored with secondary structure
base_plot <- ggplot(data, aes(mean_ddG.crystal, mean_ddG.alphafold)) +
  geom_point(aes(color=q3), size=1.5,alpha=0.7) +
  scale_colour_manual(values = c("Helix" = "dodgerblue2", "Coil" = "darkorchid", "Strand" = "coral"), name="Secondary Structure") +
  labs(
    x = "Crystal ddG",
    y = "Alphafold ddG" ) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed", size = 0.5) +
  stat_cor(aes(label = paste(..r.label..)), method = "pearson", size = 4) + 
  theme(legend.title = element_blank())+
  theme_bw()+
  coord_fixed(ratio = 1)

p_marginal <- ggMarginal(base_plot, groupColour = TRUE, groupFill = TRUE)

title <- ggdraw() +
  draw_label(
    "Crystal Structure ddG vs Alphafold Prediction ddG for KRAS - RAL", 
    hjust = 0.5
  )

plot_grid(title, p_marginal, ncol = 1, rel_heights = c(0.1, 1))

#ddG crystal vs ddG alphafold (with distributions) colored with SASA
base_plot <- ggplot(data, aes(mean_ddG.crystal, mean_ddG.alphafold)) +
  geom_point(aes(color=SASA), size=1.5, alpha=0.7) +
  scale_colour_manual(values = c("Core" = "darkorchid", "Surface" = "dodgerblue2", "Interaction" = "coral"), name="Solvent Accessible\nSurface Area") +
  labs(
    x = "Crystal ddG",
    y = "Alphafold ddG") +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed", size = 0.5) +
  stat_cor(aes(label = paste(..r.label..)), method = "pearson", size = 4) + 
  theme(legend.title = element_blank())+
  theme_bw()+
  coord_fixed(ratio = 1)

p_marginal <- ggMarginal(base_plot, groupColour = TRUE, groupFill = TRUE)

title <- ggdraw() +
  draw_label(
    "Crystal Structure ddG vs Alphafold Prediction ddG for KRAS - RAL", 
    hjust = 0.5
  )

plot_grid(title, p_marginal, ncol = 1, rel_heights = c(0.1, 1))

df_long <- data %>%
  gather(key = "model", value = "ddg", mean_ddG.crystal, mean_ddG.alphafold)%>%
  mutate(mutant_aa = substr(variant, nchar(variant), nchar(variant)))

#MAVE vs ddG
ggplot(df_long, aes(x = fitness, y = ddg)) +
  geom_point(size = 1.5, shape = 21, fill = "dodgerblue2", color = "grey20", alpha = 0.5) +
  labs(x = "MAVE Score", y = "ddG", title = "MAVE Scores vs ddG for KRAS - RAL") +
  facet_wrap(~model, scales = "free_y", labeller = labeller(model = c("mean_ddG.crystal" = "Crystal Structure ddG", 
                                                                      "mean_ddG.alphafold" = "AlphaFold Prediction ddG"))) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed", size = 0.5) +
  stat_cor(aes(label = ..r.label..), method = "pearson", size = 4) +  # Add correlation for each facet
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5))

# Merge correlation values back into the dataset
chemical_order <- c(
  "A", "V", "I", "L", "M", "F", "Y", "W", #hydrophobic (non-polar)
  "R", "H", "K", "D", "E", "S", "T", "N", "Q", "C", "G", "P" #hydrophilic (polar)
)

# Merge correlations and add title including core/surface info
df_long <- df_long %>%
  filter(model == "mean_ddG.alphafold") %>%
  mutate(
    # Ensure proper order of amino acids based on chemistry
    prop = factor(mutant_aa, levels = chemical_order),
    Pos = as.numeric(str_extract(variant, "\\d+")))

cor_values <- df_long %>%
  group_by(mutant_aa) %>%
  summarise(R = round(cor(fitness, ddg, use = "complete.obs"), 2), .groups = "drop")  # Keeps separate R for each amino acid

# Add correlation as a single title (not per mutant_aa)
df_long <- df_long %>%
  left_join(cor_values, by = "mutant_aa") %>%
  mutate(title = paste0(mutant_aa, " (R = ", R, ")"))  # Each facet gets its own R

#MAVE vs ddG clustered by amino acid and colored by SASA
ggplot(df_long, aes(x = fitness, y = ddg, fill = SASA)) +
  geom_point(size = 1.5, shape = 21, color = "grey20", alpha = 0.7) +  # Fill controls point color
  scale_fill_manual(values = c("Core" = "darkorchid", "Surface" = "dodgerblue2", "Interaction" = "coral"), name="Solvent Accessible\nSurface Area") +  # Custom colors
  labs(x = "MAVE Score", y = "ddG", title = "MAVE Scores vs Alphafold ddG for KRAS - RAL") +
  
  # Facet by amino acid chemistry, ensuring correct order
  facet_wrap(~prop, scales = "free_y", labeller = labeller(prop = setNames(df_long$title, df_long$prop))) +
  
  geom_smooth(aes(group = 1), method = "lm", se = FALSE, color = "black", linetype = "dashed", size = 0.5) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

#MAVE vs ddG clustered by amino acid and colored by secondary structure
ggplot(df_long, aes(x = fitness, y = ddg, fill = q3)) +
  geom_point(color="grey20",size = 1.5, shape = 21, alpha = 0.7) +  # Fill controls point color
  scale_fill_manual(values = c("Coil" = "darkorchid", "Helix" = "dodgerblue2", "Strand" = "coral"), name="Secondary Structure") +  # Custom colors
  labs(x = "MAVE Score", y = "ddG", title = "MAVE Scores vs Alphafold ddG for KRAS - RAL") +
  
  # Facet by amino acid chemistry, ensuring correct order
  facet_wrap(~prop, scales = "free_y", labeller = labeller(prop = setNames(df_long$title, df_long$prop))) +
  
  geom_smooth(aes(group = 1), method = "lm", se = FALSE, color = "black", linetype = "dashed", size = 0.5) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

#MAVE vs ddG colored by SASA
ggplot(df_long, aes(x = fitness, y = ddg)) +
  geom_point(aes(fill = SASA), size = 1.5, shape = 21, color = "grey20", alpha = 0.7) +  # Use `fill` for coloring by SASA
  scale_fill_manual(values = c("Core" = "darkorchid", "Surface" = "dodgerblue2", "Interaction" = "coral"), name="Secondary Structure") +  # Custom colors for Core, Surface, and Interaction
  labs(x = "MAVE Score", y = "ddG", title = "MAVE Scores vs Alphafold ddG for KRAS - RAL") +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed", size = 0.5) +
  stat_cor(aes(label = paste(..r.label..)), method = "pearson", size = 4) +  # Remove the `group` argument to avoid grouping issues
  theme_bw() + theme(plot.title = element_text(hjust = 0.5))

#MAVE vs ddG colored by secondary structure
ggplot(df_long, aes(x = fitness, y = ddg)) +
  geom_point(aes(fill = q3), size = 1.5, shape = 21, color = "grey20", alpha = 0.7) +  # Use `fill` for coloring by SASA
  scale_fill_manual(values = c("Coil" = "darkorchid", "Helix" = "dodgerblue2", "Strand" = "coral"), name="Solvent Accessible\nSurface Area") + 
  labs(x = "MAVE Score", y = "ddG", title = "MAVE Scores vs Alphafold ddG for KRAS - RAL") +
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
    SASA = ifelse(Pos %in% c(24, 25, 31, 33, 36, 37, 38, 39, 40, 41, 56, 64, 67), "Interaction", SASA),
    q3 = recode(q3, "C" = "Coil", "E" = "Strand", "H" = "Helix")
  ) %>%
  filter(!str_detect(variant, "\\*|;"))

p <- ggplot(data_mave_merge, aes(x = fitness.x, y = fitness.y)) +
  geom_point(size = 1.5, shape = 21, fill = "dodgerblue2", color = "grey20", alpha = 0.5) +  
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "black", linetype = "dashed", size = 0.5) +
  labs(x = "Abundace MAVE Score", y = "Interaction MAVE Score") +
  annotate("text", x = 0.5, y = 0.4, 
           label = paste("R =", round(cor(data_mave_merge$fitness.y, data_mave_merge$fitness.x), 2)), 
           size = 4) +
  theme_bw() +
  coord_fixed(ratio = 1) +
  xlim(-1.5, 1) +  # Set same limits for both axes
  ylim(-1.5, 0.5) 

p_marginal <- ggMarginal(p, type = "density", fill = "dodgerblue2", color= "dodgerblue2",alpha = 0.5, margins = "both", size = 5)

title <- ggplot() +
  ggtitle("MAVE Scores from Abundance vs Interaction for KRAS - RAL") +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5, size = 14))

grid.arrange(title, p_marginal, ncol = 1, heights = c(0.1, 1))

#colored by secondary structure
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
    "MAVE Scores from Abundance vs Interaction for KRAS - RAL", 
    hjust = 0.5
  )

plot_grid(title, p_marginal, ncol = 1, rel_heights = c(0.1, 1))

#colored by SASA
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
    "MAVE Scores from Abundance vs Interaction for KRAS - RAL", 
    hjust = 0.5
  )

plot_grid(title, p_marginal, ncol = 1, rel_heights = c(0.1, 1))
