## Transforming raw RNA-seq data

## Installing necessary packages
install.packages()

## loading necessary packages
library(tidyverse)
library()

## load data
counts <- read.csv("data/2024_04_10/kctie2_rnaseq_counts.csv")


## load required packages
library(limma);library(edgeR); library(EMMREML);library(RColorBrewer);library(ggplot2)
library(gridExtra); library(grid); library(doParallel)

# load in data
load(url("https://github.com/nsmackler/status_genome_2016/blob/master/macaque_status_cell.RData?raw=true"))
# There are 5 R objects
# read_counts:  raw read counts for all 440 samples after removal of genes with a median RPKM < 2
dim(read_counts)
# read_count_per_cell: list of raw read counts for each cell type as its own dataframe in a list. Genes with a median RPKM < 2 within a cell type have been removed from the data frame for that cell type:
names(read_count_per_cell)
lapply(read_count_per_cell,dim)
# info: metadata for all 440 samples
dim(info)

## make the cell specific metadata and order it ("info_cell")
info_cell=lapply(unique(info$cell),function(x){return(info[as.character(rownames(info)[info$cell==x]),])})
names(info_cell)=unique(info$cell)
info_cell=lapply(info_cell,function(x){tmp=x;rownames(tmp)=tmp$sample; tmp=tmp[order(rownames(tmp)),];return(tmp)})

## Voom normalize and remove group ('batch') effects using limma: 
residuals=lapply(names(info_cell),function(x){
  design <- model.matrix(~0+info_cell[[x]]$group)
  dge <- DGEList(counts=read_count_per_cell[[x]])
  dge <- calcNormFactors(dge)
  v <- voom(dge,design,plot=FALSE)
  fit <-lmFit(v,design)
  return(residuals.MArrayLM(object=fit, v))
});names(residuals)=names(info_cell)

save(residuals, info_cell, file = "Week_6.RData")