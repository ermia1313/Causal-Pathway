#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
pway <- args[1]

file_in  <- file.path("data", "results", paste0("pway_2024jun_", pway, ".txt"))
file_out <- file_in

DF <- read.table(
  file_in,
  sep = "\t", header = TRUE, as.is = TRUE, stringsAsFactors = FALSE,
  comment.char = "", quote = ""
)

if (nrow(DF) > 0) {
  w <- which(DF[, "sum"] < 0)
  DF[w, "padj"] <- p.adjust(DF[w, "pval"], method = "fdr")
  w <- which(DF[, "sum"] > 0)
  DF[w, "padj"] <- p.adjust(DF[w, "pval"], method = "fdr")
  O <- order(DF[, "pval"], -abs(DF[, "sum"]))
  DF <- DF[O, ]
}

write.table(
  DF, file = file_out,
  sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE
)
