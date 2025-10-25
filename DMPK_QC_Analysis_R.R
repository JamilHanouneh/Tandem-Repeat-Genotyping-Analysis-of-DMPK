# DMPK Repeat Length QC & Visualization — R Version
# Converted from your Python script; logic preserved.

# Packages
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2", repos = "https://cloud.r-project.org")
}
library(ggplot2)

# 1) Input path
qc_file <- 'C:\\Users\\jamel\\Desktop\\ALS\\data\\QC-Target_Genotype_Table.csv'

# 2) Load CSV
qc_df <- read.csv(qc_file, row.names = 1, check.names = FALSE)
cat("Raw CSV head:\n"); print(head(qc_df))

# 3) Transpose
qc_df <- as.data.frame(t(qc_df))
qc_df$Sample <- rownames(qc_df); rownames(qc_df) <- NULL
qc_df <- qc_df[, c("Sample", setdiff(names(qc_df), "Sample"))]
cat("\nTransposed DataFrame head:\n"); print(head(qc_df))

# 4) Extract DMPK consensus sizes
dmpk_df <- data.frame(Sample = character(), Allele = character(), Consensus.Size = integer(), stringsAsFactors = FALSE)
for (sample in qc_df$Sample) {
  sample_row <- qc_df[qc_df$Sample == sample, , drop = FALSE]
  if ("DMPK consensus size allele 0" %in% names(qc_df)) {
    value <- sample_row[["DMPK consensus size allele 0"]][1]
    if (!is.na(value)) dmpk_df <- rbind(dmpk_df, data.frame(Sample = sample, Allele = "Allele 0", Consensus.Size = as.integer(value)))
  }
  if ("DMPK consensus size allele 1" %in% names(qc_df)) {
    value <- sample_row[["DMPK consensus size allele 1"]][1]
    if (!is.na(value)) dmpk_df <- rbind(dmpk_df, data.frame(Sample = sample, Allele = "Allele 1", Consensus.Size = as.integer(value)))
  }
}
cat("\nDMPK DataFrame:\n"); print(dmpk_df)

# 5) Classify
if (nrow(dmpk_df) == 0) { stop("Error: No DMPK data extracted. Check 'DMPK consensus size' columns.") }
dmpk_df$Status <- ifelse(dmpk_df$Consensus.Size > 50, "Expanded (Potential Myotonic Dystrophy Risk)", "Normal")
cat("\nDMPK Repeat Lengths Across Samples:\n")
print(dmpk_df[order(-dmpk_df$Consensus.Size), c("Sample","Allele","Consensus.Size","Status")])

# 6) Plot
dmpk_df$Sample_Allele <- paste(dmpk_df$Sample, dmpk_df$Allele, sep = "_")
p <- ggplot(dmpk_df, aes(x = Sample_Allele, y = Consensus.Size, fill = Status)) +
  geom_col() + geom_hline(yintercept = 50, linetype = "dashed", color = "red") +
  geom_text(aes(label = Consensus.Size), vjust = -0.3, size = 3) +
  scale_fill_manual(values = c("Normal" = "skyblue", "Expanded (Potential Myotonic Dystrophy Risk)" = "salmon")) +
  labs(x = "Sample and Allele", y = "DMPK Repeat Length (bp)",
       title = "DMPK Repeat Expansions in Coriell Samples (Mimicking Variant Discovery)", fill = "Status") +
  theme_minimal(base_size = 12) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p)

# 7) Pathogenic subset
pathogenic <- dmpk_df[dmpk_df$Status == "Expanded (Potential Myotonic Dystrophy Risk)", c("Sample","Allele","Consensus.Size")]
cat("\nSamples with Potentially Pathogenic DMPK Expansions (>50 repeats):\n")
if (nrow(pathogenic) > 0) print(pathogenic) else cat("None found.\n")
