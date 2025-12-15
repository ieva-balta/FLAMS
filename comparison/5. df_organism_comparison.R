"""
@author: ieva-balta, majocava, naaattella
"""

library(readr)
library(dplyr)
library(stringr)
library(VennDiagram)
library(ggplot2)

"""
This script compares species compositions between the dbPTM/CPLM and UniProt versions of FLAMS.
It will generate three files:
1. Venn Digram plot showing unique, and shared species between UniProt and dbPTM/CPLM
2. CSV file showing the counts of unique accessions per species in both versions of FLAMS.
3. Dot Plot showing the accession distribution per species.
"""

light_blue = "#E4A788"
light_red = "#5361B1"

palette = c(light_blue, light_red)

# Read the organism comparison data, change the working directory to your working directoy
setwd("C:/Users/majoc/Downloads/OneDrive_1_11-8-2025")

# Read the datasets, change the path if needed. 
ptm_old <- read_csv("PTM_old_organism.csv", show_col_types = FALSE)
uniprot_df <- read_csv("uniprot_organism_data_test.csv", show_col_types = FALSE)

# Update for the first accession
uniprot_df[uniprot_df$UniProt_Accession == "A0A5S8WF70", "Organism"] <- "Leptolyngbya sp. JSC-1"

# Update for the second accession
uniprot_df[uniprot_df$UniProt_Accession == "Q845W8", "Organism"] <- "Pseudomonas dacunhae"


clean_species <- function(x) {
  if (is.na(x) || x == "") return(NA)
  
  x <- tolower(trimws(x))
  x <- strsplit(x, "\\(")[[1]][1]
  
  parts <- unlist(strsplit(x, "\\s+"))
  if (length(parts) == 0) return(NA)
  
  # Python behavior: keep 1-word species!!
  if (startsWith(parts[1], "recombinant")) {
    parts <- parts[-1]
  }
  
  # Keep up to first 2 words (Python does this)
  parts <- parts[1:min(2, length(parts))]
  
  paste(parts, collapse = " ")
}

uniprot_df$Species <- sapply(uniprot_df$Organism, clean_species)
ptm_old$Species <- sapply(ptm_old$Organism, clean_species)

# Now compare the organism distributions
old_acc <- unique(ptm_old$UniProt_Accession)
uniprot_acc <- unique(uniprot_df$UniProt_Accession)
shared_acc <- intersect(old_acc, uniprot_acc)
unique_old_acc <- setdiff(old_acc, uniprot_acc)
unique_uniprot_acc <- setdiff(uniprot_acc, old_acc)

species_old <- unique(ptm_old$Species)
species_uniprot <- unique(uniprot_df$Species)
organism_shared <- intersect(species_old, species_uniprot)
organism_old_only <- setdiff(species_old, species_uniprot)
organism_uniprot_only <- setdiff(species_uniprot, species_old)

print(paste("Organisms in both datasets:", length(organism_shared)))
print(paste("\nOrganisms unique to old dataset:", length(organism_old_only)))
print(paste("\nOrganisms unique to UniProt dataset:", length(organism_uniprot_only)))

venn.diagram(
  x = list(`UniProt` = species_uniprot, `dbPTM/CPLM` = species_old),
  filename = "organism_comparison_venn_new.png",
  height = 800, 
  width = 1000
)

x = list(`UniProt` = species_uniprot, `dbPTM/CPLM` = species_old)
venn_plot = venn.diagram(
  x,
  filename = NULL,
  fill = palette,
  alpha = 0.5,
  cex = 0.9,
  fontfamily = "Arial",
  fontface = "bold",
  cat.cex = 0,
  margin = 0.2,              # prevents clipping
  #  cat.pos = c(0, 0),      # dbPTM label (left), UniProt label (right)
  #  cat.dist = c(0.05, 0.05),
  #  cat.fontfamily = "Arial",
  lty= "blank",
  lwd = 0, 
  bg = NULL
)

svg("Final_pictures/organims_venn_diagram_set.svg")
grid.draw(venn_plot)
dev.off()

# Group by species and count unique accessions

# First I want to group by species and count unique accessions
uniprot_counts <- uniprot_df %>%
  group_by(Species) %>%
  summarise(UniProt = n_distinct(UniProt_Accession)) %>%
  ungroup()

old_counts <- ptm_old %>%
  group_by(Species) %>%
  summarise(`dbPTM/CPLM` = n_distinct(UniProt_Accession)) %>%
  ungroup()

merged_counts <- full_join(uniprot_counts, old_counts, by = "Species") %>%
  replace_na(list(UniProt = 0, `dbPTM/CPLM` = 0)) %>%
  mutate(Total = UniProt + `dbPTM/CPLM`) %>%
  arrange(desc(UniProt))

# Save the merged counts to a CSV for reference
write_csv(merged_counts, "organism_protein_counts_comparison_total.csv")


###################################################################
######################### Dot Plot ################################
###################################################################

org_df <- read.csv("organism_protein_counts_comparison_total.csv")

org_df <- org_df %>%
  rename(
    Species   = Species,           # PTM name
    dbPTM_CPLM = dbPTM.CPLM,  # dbPTM column
    UniProt    = UniProt      # UniProt column
  )

library(stringr)
org_df$Species <- str_to_sentence(org_df$Species)

org_df$Species[org_df$Species == "Human immunodeficiency"] <- "HIV"
org_df$Species[org_df$Species == "Influenza a"] <- "Influenza A"


org_df_long <- org_df %>%
  mutate(
    dbPTM_log = log10(dbPTM_CPLM + 1),
    UniProt_log = log10(UniProt + 1)
  ) %>%
  select(Species, dbPTM_log, UniProt_log) %>%
  pivot_longer(
    cols = ends_with("_log"),
    names_to = "database",
    values_to = "log_count"
  ) %>%
  mutate(
    database = recode(database,
                      "dbPTM_log" = "dbPTM/CPLM",
                      "UniProt_log" = "UniProt"))
# spcies order
organism_order <- org_df %>%
  arrange(desc(UniProt)) %>%      # highest -> lowest
  pull(Species)
org_df_long <- org_df_long %>%
  mutate(
    Species = factor(Species, levels = organism_order)
  )

# top 50
top50_species <- levels(org_df_long$Species)[1:50]

org_df_top50 <- org_df_long |>
  dplyr::filter(Species %in% top50_species) |>
  dplyr::mutate(Species = droplevels(Species))



p_org = ggplot(org_df_top50, aes(x = Species, y = log_count, color = database)) +
  geom_line(aes(group = Species), color="grey70", size = 1.2) +
  geom_point(size = 5, alpha = 0.85) +
  scale_color_manual(values = palette,
  ) + 
  labs(
    y = "Number of PTMs (log-scale)",
    color = NULL,
  ) +
  theme_minimal(base_size = 10) +
  theme(panel.grid = element_blank(),
        text = element_text(family = "Arial"),
        axis.text.x = element_text(angle=90, hjust=1, size = 16, face = "italic"),
        axis.text.y = element_text(size = 14),
        legend.position = c(0.95, 0.95),     # << move here,
        legend.key = element_rect(fill = "white", color = NA),
        legend.text = element_text(size = 12),
        legend.key.size = unit(0.8, "cm"),
        panel.background  = element_rect(fill = NA, colour = NA),   # key
        plot.background   = element_rect(fill = NA, colour = NA),   # key
        legend.background = element_rect(fill = NA, colour = NA),
        # keep axes
        axis.line         = element_line(color = "black", linewidth = 0.5),
        axis.ticks        = element_line(color = "black"),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank())
p_org

ggsave(
  "organismtop50_dotplot.svg",
  p_org,
  width = 14,
  height = 7,
  bg = "transparent"
)