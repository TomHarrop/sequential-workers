#!/usr/bin/env Rscript

log <- file(snakemake@log[[1]],
            open = "wt")
sink(log,
     type = "message")
sink(log,
     type = "output",
     append = TRUE)

library(data.table)


vcf_file <- snakemake@input[["vcf"]]

# read the vcf
vcf_names <- c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", 
               "FORMAT")

message(paste(date(), "Parsing header"))
cn <- fread(cmd = paste("zcat", vcf_file, '| grep "^#[^#]" | head -n 1'),
      sep = "\t", header = FALSE)

message(paste(date(), "Reading VCF"))
vcf <- fread(cmd = paste("zcat", vcf_file, '| grep -v "^#"'),
      sep = "\t",
      header = FALSE,
      col.names = t(cn)[, 1])


# extract the AO and RO fields
message(paste(date(), "Converting to long data.table"))
indiv_cols <- names(vcf)[!names(vcf) %in% vcf_names]
gt_format <- vcf[1, strsplit(FORMAT, ":", fixed = TRUE)][, V1]
vcf_long <- melt(vcf,
     id.vars = c("#CHROM", "POS", "ID"),
     measure.vars = indiv_cols,
     variable.name = "indiv",
     value.name = "genotype_field")
message(paste(date(), "Parsing the GT field"))
vcf_split <- vcf_long[, tstrsplit(genotype_field,
                     split = ":",
                     fixed = TRUE),
                     by = c("indiv", "#CHROM", "POS", "ID")]
setnames(vcf_split, paste0("V", 1:7), gt_format)

# calculate MAF
message(paste(date(), "Calculating MAF"))
vcf_split[, AO := as.numeric(AO)]
vcf_split[, RO := as.numeric(RO)]
vcf_split[, 
    maf := min(AO, RO) / (AO + RO),
    by = c("indiv", "#CHROM", "POS", "ID")]

# convert for export
message(paste(date(), "Converting back to wide data.table"))
maf_dt <- dcast(vcf_split, `#CHROM` + POS ~ indiv,
                value.var = "maf")
maf_dt[, loc_id := paste(`#CHROM`, POS, sep = "_")]
maf_mat <- t(as.matrix(maf_dt[, c("loc_id", indiv_cols), with = FALSE],
          rownames = "loc_id"))
maf_dt[, loc_id := NULL]

message(paste(date(), "Writing MAF for", dim(maf_dt)[[1]], "loci"))
saveRDS(maf_mat, snakemake@output[["maf_mat"]])
fwrite(maf_dt, snakemake@output[["maf_dt"]])

sessionInfo()
