setwd("C:\\old_E\\GeoTcgaData_work\\甲基化数据hg19转hg38")
ann <- minfi::getAnnotation(
    IlluminaHumanMethylation450kanno.ilmn12.hg19::IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann <- as.data.frame(ann)
ann <- ann[, c("chr", "pos")]
ann$end <-  ann$pos
ann$start <- ann$end-1
ann <- ann[, c("chr", "start", "end")]
aa <- paste(ann[, 1], ann[, 2], sep = ":")
bb <- paste(aa, ann[, 3], sep = "-")
write.table(bb, "hg19_1.bed", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(ann, "hg19_2.bed", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)