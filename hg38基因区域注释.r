setwd("E:\\GeoTcgaData_work\\甲基化数据hg19转hg38")
## 获取原版注释情况，备用
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann <- minfi::getAnnotation(
    IlluminaHumanMethylation450kanno.ilmn12.hg19::IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann <- as.data.frame(ann)
ann_anto <- ann[, c("chr", "pos", "UCSC_RefGene_Name", "UCSC_RefGene_Group")]
## 首先获取hg19和hg38都有的坐标
ncbi <- fread("report_hg19_2.txt")
class(ncbi) <- "data.frame"
ncbi <- ncbi[, c(1, 4, 8, 5, 13)]

## 将坐标对应的cpg编号标上去（可选）
ann_anto2 <- ann_anto
rownames(ann_anto2) <- paste(ann_anto2[, 1], ann_anto2[, 2], sep = "_")
ann_anto2$cpg <- rownames(ann_anto)
ncbi$cpg <- ann_anto2[ncbi[, 1], "cpg"]
ncbi <- na.omit(ncbi)

## 根据坐标，分别从hg19原版、ChIPseeker_19, ChIPseeker_38来注释
pos_19 <- ncbi[, c(2, 3)]
pos_19$end <- pos_19[, 2]
colnames(pos_19) <- c("chr", "start", "end")
pos_19$start <- pos_19$end - 1


pos_38 <- ncbi[, c(4, 5)]
pos_38$end <- pos_38[, 2]
colnames(pos_38) <- c("chr", "start", "end")
pos_38$start <- pos_38$end - 1

## 使用hg19的坐标，并用ChIPseeker注释到hg19基因组
require(ChIPseeker)
library(org.Hs.eg.db)
library(GenomicRanges)
pos <- pos_19
# pos_anno=as.data.frame(peakAnno)
require(ChIPseeker)
library(org.Hs.eg.db)
library(GenomicRanges)
peak <- GRanges(seqnames=Rle(pos[,1]),
                ranges=IRanges(pos[,2], pos[,3]), strand=rep(c("*"), nrow(pos)))

library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb_19=TxDb.Hsapiens.UCSC.hg19.knownGene
peakAnno_19 <- annotatePeak(peak, tssRegion=c(-2000, 200), level = "gene",
                         TxDb=txdb_19, annoDb="org.Hs.eg.db")
						 
peakAnno_19 <- as.data.frame(peakAnno_19)
peakAnno_19$annotation <- gsub("\\(.*\\)", "", peakAnno_19$annotation)
aa_19 <- table(peakAnno_19$annotation)



## 使用hg38的坐标，并用ChIPseeker注释到hg38基因组
require(ChIPseeker)
library(org.Hs.eg.db)
library(GenomicRanges)
pos <- pos_38
# pos_anno=as.data.frame(peakAnno)
require(ChIPseeker)
library(org.Hs.eg.db)
library(GenomicRanges)
peak <- GRanges(seqnames=Rle(pos[,1]),
                ranges=IRanges(pos[,2], pos[,3]), strand=rep(c("*"), nrow(pos)))

library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb_38=TxDb.Hsapiens.UCSC.hg38.knownGene
peakAnno_38 <- annotatePeak(peak, tssRegion=c(-2000, 200), level = "gene",
                         TxDb=txdb_38, annoDb="org.Hs.eg.db")
						 
peakAnno_38 <- as.data.frame(peakAnno_38)
peakAnno_38$annotation <- gsub("\\(.*\\)", "", peakAnno_38$annotation)
aa_38 <- table(peakAnno_38$annotation)

## 给peakAnno_19和peakAnno_38添加cpg位点
ncbi$hg19 <- ncbi[, 1]
ncbi$hg38 <- paste(ncbi[, "mapped_id"], ncbi[, "mapped_start"], sep = "_")

peakAnno_19$pos <- paste(peakAnno_19$seqnames,  peakAnno_19$end, sep ="_" )
peakAnno_38$pos <- paste(peakAnno_38$seqnames,  peakAnno_38$end, sep ="_" )
peakAnno_19$cpg <- ncbi[match(peakAnno_19$pos, ncbi$hg19), "cpg"]
peakAnno_38$cpg <- ncbi[match(peakAnno_38$pos, ncbi$hg38), "cpg"]

cpgAnno_19 <- peakAnno_19[, c("cpg", "geneId", "annotation", "SYMBOL")]
cpgAnno_38 <- peakAnno_38[, c("cpg", "geneId", "annotation", "SYMBOL")]
save(cpgAnno_19, cpgAnno_38, file = "cpg_annotation.Rdata")
write.table(cpgAnno_19, file = "cpgAnno_19.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(cpgAnno_38, file = "cpgAnno_38.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
