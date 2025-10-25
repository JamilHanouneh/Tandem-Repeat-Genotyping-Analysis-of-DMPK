# --- Libraries ---
suppressPackageStartupMessages({
  library(VariantAnnotation)  # Bioconductor
  library(ggplot2)
  library(dplyr)
  library(stringr)
})

# --- Step 1: Parse VCFs for repeat data ---
vcf_dir <- "C:/Users/jamel/Desktop/ALS/data/TRGT-VCF_files"
repeat_rows <- list()

if (dir.exists(vcf_dir)) {
  vcf_files <- list.files(vcf_dir, pattern = "\\.vcf(\\.gz)?$", full.names = TRUE)
  if (length(vcf_files) == 0) {
    message("No .vcf files found in: ", vcf_dir)
  }

  for (vf in vcf_files) {
    sample <- sub("\\.vcf(\\.gz)?$", "", basename(vf))
    message("Parsing: ", basename(vf))

    tryCatch({
      vcf <- readVcf(vf)

      nrec <- nrow(vcf)
      if (is.null(nrec) || nrec == 0) next

      # Access INFO fields; they can be CharacterList / List-like per row
      info_trid   <- info(vcf)$TRID
      info_motifs <- info(vcf)$MOTIFS

      # Access FORMAT field "AL" (allele lengths). This is a matrix [variants x samples]
      al_mat <- geno(vcf)$AL

      for (i in seq_len(nrec)) {
        # TRID
        trid <- tryCatch({
          x <- info_trid[[i]]
          if (length(x) == 0 || is.na(x)) "Unknown" else as.character(x)[1]
        }, error = function(e) "Unknown")

        # Keep only DMPK
        if (!identical(trid, "DMPK")) next

        # Skip if AL not present
        if (is.null(al_mat)) next

        # Use the first sample column, matching your Python logic
        al_raw <- tryCatch(al_mat[i, 1], error = function(e) NA)

        if (length(al_raw) == 0 || is.na(al_raw)) next

        # Coerce to integer vector (handles "30,140" or numeric)
        if (is.numeric(al_raw)) {
          lengths <- as.integer(al_raw)
        } else {
          parts <- unlist(strsplit(as.character(al_raw), ",", fixed = TRUE))
          lengths <- suppressWarnings(as.integer(parts[parts != ""]))
          lengths <- lengths[!is.na(lengths)]
        }
        if (length(lengths) == 0) next

        # MOTIF: take first value if multiple
        motif <- tryCatch({
          m <- info_motifs[[i]]
          if (length(m) == 0 || all(is.na(m))) "Unknown" else as.character(m)[1]
        }, error = function(e) "Unknown")

        # Build rows (Allele numbering starts at 0 to mirror Python enumerate)
        for (idx in seq_along(lengths)) {
          repeat_rows[[length(repeat_rows) + 1]] <- data.frame(
            Sample = sample,
            Gene   = trid,
            Length = lengths[idx],
            Motif  = motif,
            Allele = paste0("Allele ", idx - 1L),
            stringsAsFactors = FALSE
          )
        }
      }
    }, error = function(e) {
      message("Error parsing ", basename(vf), ": ", conditionMessage(e))
    })
  }
} else {
  message("VCF directory not found - check path or ensure TRGT-VCF_files exists.")
}

vcf_df <- if (length(repeat_rows)) bind_rows(repeat_rows) else data.frame()

# --- Step 2: Visualize DMPK repeats ---
if (nrow(vcf_df) > 0) {
  vcf_df <- vcf_df %>%
    mutate(SampleAllele = paste0(Sample, "_", Allele))

  threshold <- 50L

  p <- ggplot(vcf_df, aes(x = SampleAllele, y = Length)) +
    geom_col(fill = "salmon") +
    geom_hline(yintercept = threshold, linetype = "dashed", color = "red") +
    annotate("text", x = Inf, y = threshold, label = "Pathogenic Threshold (>50 repeats)",
             hjust = 1.05, vjust = -0.5, size = 3.2) +
    labs(
      x = "Sample and Allele",
      y = "DMPK Repeat Length (bp)",
      title = "DMPK Repeat Expansions in Coriell Samples"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(face = "bold")
    ) +
    geom_text(aes(label = Length), vjust = -0.3)

  print(p)
} else {
  message(
    "No DMPK data found in VCFs - check files or download TRGT VCFs from ",
    "https://downloads.pacbcloud.com/public/dataset/PureTargetRE/Coriell/PBMM2-BAM-Input-For-IGV-And-TRGT/."
  )
}

# --- Optional: Simulate 'novel' discovery by checking motif interruptions ---
if (nrow(vcf_df) > 0) {
  is_interrupted <- function(m) {
    if (is.na(m) || m == "Unknown") return("No")
    chars <- unlist(strsplit(m, "", fixed = TRUE))
    ifelse(length(unique(chars)) > 1, "Yes (Potential Novel Variant)", "No")
  }

  vcf_df <- vcf_df %>%
    mutate(Interrupted = vapply(Motif, is_interrupted, character(1)))

  cat("\nInterruptions in DMPK motifs (mimicking novel patterns):\n")
  print(vcf_df %>% select(Sample, Allele, Motif, Interrupted))
}
