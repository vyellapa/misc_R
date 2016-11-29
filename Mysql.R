#Library loading
library(RMySQL)


#Function that returns a table with exon coordinates (w.o. UTR)
#Example use : ensemblExons = GetExonsFromName(geneName)
GetExonsFromName <- function(name,genome="hg19")
{
  con = dbConnect(MySQL(),user="genome",dbname=genome,host="genome-mysql.cse.ucsc.edu", A="")
  gene = dbGetQuery(con,paste('select * from refGene where name2="',name,'"',sep=''))
  if(dim(gene)[1]>0)
  {
    exon.start = list()
    exon.end = list()
    refseq.list = list()
    chromosome.list = list()
    for(i in 1:dim(gene)[1])
    {
      chromosome = gene[i,3]
      cds.start = as.numeric(gene[i,7])
      cds.end = as.numeric(gene[i,8])
      gene.name = gene[i,13]
      gene.refseq = gene[i,2]
      
      
      exon.start.buffer = unlist(as.numeric(strsplit(gene[i,10],",",fixed=T)[[1]]))
      exon.end.buffer = unlist(as.numeric(strsplit(gene[i,11],",",fixed=T)[[1]]))
      
      exon.start.buffer = lapply(exon.start.buffer,function(start){return(max(cds.start,start))})
      exon.end.buffer = lapply(exon.end.buffer,function(start){return(max(cds.start,start))})
      exon.start.buffer = lapply(exon.start.buffer,function(end){return(min(cds.end,end))})
      exon.end.buffer = lapply(exon.end.buffer,function(end){return(min(cds.end,end))})
      exon.start = c(exon.start,exon.start.buffer)
      exon.end = c(exon.end,exon.end.buffer)
      chromosome.list = c(chromosome.list,rep(chromosome,length(exon.end.buffer)))
      refseq.list = c(refseq.list,rep(gene.refseq,length(exon.end.buffer)))
    }
    
    dbDisconnect(con)
    
    # Added to exclude UTR regions
    wUTR=cbind(refseq.list,chromosome.list,exon.start,exon.end,name)
    woUTR=list()
    
    for(i in 1:dim(wUTR)[1]) {
      if(wUTR[i,]$exon.start != wUTR[i,]$exon.end) {
        woUTR=rbind(woUTR,c(wUTR[i,]$refseq.list,wUTR[i,]$chromosome.list,wUTR[i,]$exon.start,wUTR[i,]$exon.end,wUTR[i,]$name))
      }
    }
    colnames(woUTR)=c("refseq","chromosome","exon.start","exon.end","name")
    exon.matrix=as.data.frame(woUTR)
    
    return(woUTR)	
  }
  else
  {
    dbDisconnect(con)
    cat("missing the following gene :",name,"\n")
    return(c("","","","",name))
  }
}

gnames=read.table("~/Desktop/MM_genes_final_1col.txt",header=F,sep="\t")
for(gene in gnames$V1) {exons=GetExonsFromName(gene); write.table(exons,file="USCS_woUTR_mmGenes.tsv", append=TRUE, quote=F,sep="\t", eol="\n", row.names=F, col.names=TRUE)}


g=GRanges(seqnames= Rle( as.character(aa$chr)), ranges = IRanges(start = aa$start, end=aa$end))
g=reduce(g)
lala=data.frame(seqnames=seqnames(g),starts=start(g)-1,ends=end(g))


