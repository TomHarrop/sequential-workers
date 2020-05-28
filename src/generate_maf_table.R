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

cn <- fread(cmd = paste("zcat", vcf_file, '| grep "^#[^#]" | head -n 1'),
      sep = "\t", header = FALSE)

vcf <- fread(cmd = paste("zcat", vcf_file, '| grep -v "^#"'),
      sep = "\t",
      header = FALSE,
      col.names = t(cn)[, 1])


# extract the AO and RO fields
indiv_cols <- names(vcf)[!names(vcf) %in% vcf_names]
gt_format <- vcf[1, strsplit(FORMAT, ":", fixed = TRUE)][, V1]
vcf_long <- melt(vcf,
     id.vars = c("#CHROM", "POS", "ID"),
     measure.vars = indiv_cols,
     variable.name = "indiv",
     value.name = "genotype_field")
vcf_split <- vcf_long[, tstrsplit(genotype_field,
                     split = ":",
                     fixed = TRUE),
                     by = c("indiv", "#CHROM", "POS", "ID")]
setnames(vcf_split, paste0("V", 1:7), gt_format)

# calculate MAF
vcf_split[, 
    maf := min(as.numeric(AO), as.numeric(RO)) / (as.numeric(AO) + as.numeric(RO)),
    by = c("indiv", "#CHROM", "POS", "ID")]

# convert for export
maf_dt <- dcast(vcf_split, `#CHROM` + POS ~ indiv,
                value.var = "maf")
maf_dt[, loc_id := paste(`#CHROM`, POS, sep = "_")]
maf_mat <- t(as.matrix(maf_dt[, c("loc_id", indiv_cols), with = FALSE],
          rownames = "loc_id"))
maf_dt[, loc_id := NULL]

saveRDS(maf_mat, snakemake@output[["maf_mat"]])
fwrite(maf_dt, snakemake@output[["maf_dt"]])

sessionInfo()
