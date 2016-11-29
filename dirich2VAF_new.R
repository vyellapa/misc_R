library(ggplot2)
library(sciClone)
setwd("~/Desktop/scratch/")
for(i in grep("xvfb",list.files(pattern = ".*_bestConsensusAssignments.bed", recursive = T), value=T)) {
  name=gsub("_bestConsensusAssignments.bed", "",strsplit(i,"/")[[1]][4])
  T1=name
  T2=gsub("-T1-","-T2-",T1)
  
  t1c=paste(name,".caveman.tsv.gz",sep=".*")
  t2c=paste(T2,".caveman.tsv.gz",sep=".*")
  
  t1c=list.files(pattern=t1c)
  t2c=list.files(pattern=t2c)
  
  caveman_t1c=read.table(gzfile(t1c), sep="\t", header=T, stringsAsFactors = FALSE)
  colnames(caveman_t1c)=paste("T1_CAVE",colnames(caveman_t1c),sep="_")
  caveman_t2c=read.table(gzfile(t2c), sep="\t", header=T, stringsAsFactors = FALSE)
  colnames(caveman_t2c)=paste("T2_CAVE",colnames(caveman_t2c),sep="_")
  caveman_t1c$merger=paste(caveman_t1c$T1_CAVE_CHR, caveman_t1c$T1_CAVE_START,sep=":")
  caveman_t2c$merger=paste(caveman_t2c$T2_CAVE_CHR, caveman_t2c$T2_CAVE_START,sep=":")
  
  cluster=read.table(sep="\t", header=T,i)
  colnames(cluster)=paste("clust",colnames(cluster),sep="_")
  cluster$merger=paste(cluster$clust_chr, cluster$clust_end,sep=":")
  
  merged=merge(cluster, caveman_t1c,by.x="merger", by.y="merger", all.x = T )
  merged=merge(merged, caveman_t2c,by.x="merger", by.y="merger", all.x = T )
  
  
  merged[is.na(merged$T1_CAVE_TARGET_VAF),]$T1_CAVE_TARGET_VAF=0.00
  merged[is.na(merged$T2_CAVE_TARGET_VAF),]$T2_CAVE_TARGET_VAF=0.00
  print(t1c)
  print(t2c)
}