################################################################################
############################### ONCOPRINT #######################################
# Tips on how to draw an oncoPrint using ComplexHeatmap.
#
# The oncoPrint coded here visualizes copy number amplifications across a
# selected gene panel and patient cohort.
#
# In this plot:
#  - Rows represent genes
#  - Columns represent patients
#  - Red rectangles indicate amplification events
#  - The top barplot shows the number of amplified genes per patient
#  - The right barplot shows the number of patients amplified per gene
#
# To replicate this plot you need:
#  - genePanel_proteinCodingGenes.RData
#    (containing amp_cnaTCGA_pc and cnaTCGA_pc_Panel objects)
#
# NOTE:
# This oncoPrint only visualizes amplification events.
################################################################################


############################### IMPORTS #########################################

library(ComplexHeatmap)
library(circlize)
library(grid)


############################### DATA LOADING ####################################

# Load the Example matrix
load("Example_Data/example_genePanel_amplifications.RData")
# CNA matrix with only genes selected for the panel
# Rows: genes
# Columns: patients


###################### FORMAT DATA FOR ONCOPRINT ################################

# Convert the numeric CNA matrix into a character matrix
# oncoPrint requires character strings for alterations
#
# Here:
#   1  -> "Amplification"
#   all other values -> "" (empty, no alteration shown)
cna_pcPanel_char <- apply(cna_pcPanel, 2, function(col) {
  ifelse(col == 1, "Amplification", "")
})


############################ COLOR DEFINITIONS ##################################

# Define colors for each alteration type
# Only amplifications are plotted
col <- c("Amplification" = "#C41E3A")


######################## DEFINE ALTERATION DRAWING ##############################

# alter_fun defines how each CNA event is drawn inside the heatmap cells
alter_fun <- list(
  # Background of each cell (light gray rectangle)
  background = function(x, y, w, h) {
    grid.rect(x, y, w, h * 0.7, gp = gpar(fill = "grey92", col = NA))},
  # Amplification event (red rectangle)
  Amplification = function(x, y, w, h) {
    grid.rect(x, y, w, h * 0.85, gp = gpar(fill = col["Amplification"], col = NA))}
)


########################## SUMMARY STATISTICS ###################################

# Count the number of amplification events per gene (row-wise)
gene_counts <- rowSums(cna_pcPanel_char == "Amplification")

# Count the number of amplified genes per patient (column-wise)
patient_counts <- colSums(cna_pcPanel_char == "Amplification")


######################## CREATE ANNOTATIONS #####################################

# Top annotation:
# Barplot showing number of amplified genes per patient
top_anno <- HeatmapAnnotation(
  top_bar = anno_barplot(
    patient_counts,
    border = FALSE,
    gp = gpar(fill = "#C41E3A", col = NA),
    height = unit(2.5, "cm"),
    labels = " "
  ),
  show_annotation_name = FALSE
)

# Right annotation:
# Barplot showing number of patients with amplification per gene
right_anno <- rowAnnotation(
  barplot = anno_barplot(
    gene_counts,
    border = FALSE,
    gp = gpar(fill = "#1F4F4F", col = NA)
  ),
  show_annotation_name = FALSE
)


############################ DRAW AND SAVE PLOT #################################

# Open a high-resolution PNG device suitable for publication
# Adapt size settings for your plot
png("oncoPrint.png", width = 2800, height = 1800, res = 300)

# Draw the oncoPrint
oncoPrint(
  cna_pcPanel_char,
  alter_fun = alter_fun,
  col = col,
  column_title = "Gene Panel Amplification Distribution in Patients",
  top_annotation = top_anno,
  right_annotation = right_anno,
  show_heatmap_legend = FALSE,
  pct_side = "right",
  row_names_side = "left",
  show_column_names = FALSE,
  row_names_gp = gpar(fontsize = 8, fontface = "bold"),
  column_title_gp = gpar(fontsize = 14, fontface = "bold"),
  
  # Preserve original gene order
  row_order = rownames(cna_pcPanel_char)
)

# Close the graphics device
dev.off()
