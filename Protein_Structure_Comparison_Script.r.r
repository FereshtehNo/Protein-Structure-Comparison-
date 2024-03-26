# Load necessary libraries
library(Biostrings)
library(msa)
library(Matrix)

# Step 1: Set the project directory
project_dir <- "D:/Phd-classes/Final-project"
input_dir <- file.path(project_dir, "ProcessedData")
output_dir <- file.path(project_dir, "AnalysisResults")

# Step 2: Create input and output directories if they don't exist
dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Step3: Define the protein sequences
sequence1 <- "PVIPLDPARRPVIKAQVDTQTSHPKTIEALLDTGADMTVIPIALFSSNTPLKNTSVLGAGGQTQDHFKLTSLPVLIRLPFRTTPIVLTSCLVDTKNNWAIIGRDALQQCQGVLYLP"
sequence2 <- "PQITLWQRPIVTIKIGGQLKEALLNTGADDTVLEEVNLPGRWKPKLIGGIGGFVKVRQYDQVPIEICGHKVIGTVLVGPTPANVIGRNLMTQIGCTLNF"

# Create amino acids String Set objects
protein_strings1 <- AAStringSet(sequence1)
protein_strings2 <- AAStringSet(sequence2)

# Perform pairwise alignment
alignment <- pairwiseAlignment(protein_strings1, protein_strings2)

# Print the alignment
print(alignment)

# Step4: Define the sequences for multiple sequence alignment
sequence3 <- "NPSPHQVYNVTWTITNLVTGTKANATSMLGTLTDAFPTMYFDLCDIIGNTWNPSDQEPFPGYGCDQPMRRWQQRNTPFYVCPGHANRKQCGGPQDGFCAVWGCETTGETYWRPTSSWDYITVKKGVTQGIYQCSGGGWCGPCYDKAVHSSTTGASEGGRCNPLILQFTQKGRQTSWDGPKSWGLRLYRSGYDPIALFSVSRQVMTITPPQA"

# Create amino acids String Set objects
sequences <- AAStringSet(c(sequence1, sequence2, sequence3))

# Perform multiple sequence alignment
alignment <- msa(sequences)

# Print the alignment
print(alignment)

# Step 5: TM-Score
matrix_file1 <- file.path(input_dir, "distance_matrix_HIV-DOMAIN-A-ALL-ATOMS.txt")
matrix_file2 <- file.path(input_dir, "distance_matrix_HTLV-1-DOMAIN-A-ALL-ATOMS.txt")

# Read the distance matrices from files (assuming they are tab-separated text files)
dist_matrix1 <- as.matrix(read.table(matrix_file1, header = FALSE, sep = "\t"))
dist_matrix2 <- as.matrix(read.table(matrix_file2, header = FALSE, sep = "\t"))

tm_score <- calculate_tm_score(dist_matrix1, dist_matrix2)

# Print the TM-score value
cat("TM-score between the two protein structures:", tm_score, "\n")

# Compare TM-score to the thresholds and provide interpretation
if (tm_score <= 0.17) {
  cat("TM-score indicates unrelated proteins (TM-score <= 0.17)\n")
} else if (tm_score > 0.17 && tm_score <= 0.5) {
  cat("TM-score suggests some structural similarity (0.17 < TM-score <= 0.5)\n")
} else {
  cat("TM-score indicates generally the same fold (TM-score > 0.5)\n")
}
# Calculate TM-score for the distance matrices
tm_score <- calculate_tm_score(dist_matrix1, dist_matrix2)

# Compare TM-score to the thresholds and provide interpretation
if (tm_score <= 0.17) {
  cat("TM-score indicates unrelated proteins (TM-score <= 0.17)\n")
} else if (tm_score > 0.17 && tm_score <= 0.5) {
  cat("TM-score suggests some structural similarity (0.17 < TM-score <= 0.5)\n")
} else {
  cat("TM-score indicates generally the same fold (TM-score > 0.5)\n")
}

# Step6: RMSD
dist_matrix1 <- as.matrix(read.table("HIV-RMSD.txt", header = FALSE))
dist_matrix2 <- as.matrix(read.table("HTLV-1-RMSD.txt", header = FALSE))

# calculate RMSD between two distances
calculate_rmsd <- function(matrix1, matrix2) {
  num_atoms <- min(nrow(matrix1), nrow(matrix2))
  
  # Calculate squared differences between distances
  squared_diffs <- (matrix1[1:num_atoms, ] - matrix2[1:num_atoms, ])^2
  
  # Sum of squared differences
  sum_squared_diffs <- sum(squared_diffs)
  
  # Calculate RMSD
  rmsd <- sqrt(sum_squared_diffs / num_atoms)
  return(rmsd)
}

# Calculate RMSD between the two protein structures
rmsd_value <- calculate_rmsd(dist_matrix1, dist_matrix2)

cat("RMSD between the two protein structures:", rmsd_value, "Ångström\n")

# Create scatter plot with RMSD value displayed
png("scatter_plot.png")
plot(dist_matrix1, dist_matrix2, main = "RMSD VALUE", xlab = "HIV", ylab = "HTLV-1", col = "yellow", pch = 16)
abline(0, 1, col = "green")

# Add text with RMSD value to the plot (further adjust X-coordinate for positioning)
text(25, 2.8, paste("RMSD =", round(rmsd_value, 2), "Ångström"), col = "black")

dev.off()

cat("Scatter plot saved as scatter_plot.png\n")

