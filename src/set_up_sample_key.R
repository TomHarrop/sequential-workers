#!/usr/bin/env Rscript

library(data.table)
library(seqinr)

RevComp <- function(x){
    toupper(c2s(rev(comp(s2c(x)))))
}


sample_sheet_raw <- fread("data/SampleSheet_MiSeq_E7600.csv", skip = 21)
final_pool <- fread("data/ws_libprep_layout-i79_12-20200122.csv")

# mung the barcodes in sample_sheet
munged_barcodes <- rbind(
    sample_sheet_raw[, .(index_id = I7_Index_ID,
                         index = index)],
    unique(sample_sheet_raw[, .(index_id = I5_Index_ID,
                                index = RevComp(index2)),
                            by = index2][, .(index_id, index)]))

# set up barcode sequences
samples_only <- final_pool[, .(
    sample = id,
    i7_primer,
    i5_primer)]

# fix the weird i7s
samples_only[grepl("^i7[[:digit:]]{3}", i7_primer),
             i7_primer := sub("i70", "i7", i7_primer)]

# merge
samples_i7 <- merge(samples_only,
                    munged_barcodes[, .(index_id, i7_index = index)],
                    by.x = "i7_primer", by.y = "index_id")
samples_both <- merge(samples_i7,
                      munged_barcodes[, .(index_id, i5_index = index)],
                      by.x = "i5_primer", by.y = "index_id")
samples_both[, barcode_sequence := paste(i7_index, i5_index, sep = "+")]

# set up file locations
samples_both[, r1_path := paste0("data/reads/",
                                 sample,
                                 "/",
                                 sample,
                                 "_R1.fq.gz")]
samples_both[, r2_path := paste0("data/reads/",
                                 sample,
                                 "/",
                                 sample,
                                 "_R2.fq.gz")]


# set up metadata
sample_info <- samples_both[, .(
    sample, barcode = barcode_sequence, r1_path, r2_path
)]

fwrite(sample_info, "data/sample_key.csv")
