#Library loading
library(RMySQL)

#Function that returns a table with exon coordinates (w.o. UTR)
#Example use : ensemblExons = GetExonsFromName(geneName)

GetGeneCoordsFromName <- function(name,genome="hg19")
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
      tx.start = as.numeric(gene[i,5])
      tx.end = as.numeric(gene[i,6])
      gene.name = gene[i,13]
      gene.refseq = gene[i,2]
      
      
      
    }
    
    dbDisconnect(con)
    return(cat(chromosome,":", tx.start,"-", tx.end, gene.name, gene.refseq, sep="\t"))	
  }
  else
  {
    dbDisconnect(con)
    cat("missing the following gene :",name,"\n")
    return(c("","","","",name))
  }
}


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
        exon.uniq=unique(subset(exon.matrix,select=-c(refseq)))
        exons.sum=abs(sum(as.numeric(exon.uniq$exon.start)-as.numeric(exon.uniq$exon.end)))
        #exons.sum=exons.sum+nrow(exon.uniq)
        #woUTR.exon.sum=cbind(exon.matrix,exons.sum)
        
        return(woUTR)	
    }
    else
    {
        dbDisconnect(con)
        cat("missing the following gene :",name,"\n")
        return(c("","","","",name))
    }
}

#Function to "prune" the snps from Gunes's script
#Example use : pruned = pruneBaits(baitTable,minDistance=500,maxBaits=50,minBaits=10)
#Return a bait table with less than maxBaits baits, more than minBaits baits and which is going to remove baits with a distance inferior to minDistance in priority

pruneBaits = function(baitTable,minDistance=500,maxBaits=50,minBaits=10){
  toRemove = list()
  genes = lapply(1:dim(baitTable)[1],function(iter){
    gene = strsplit(baitTable[iter,4],'_',fixed=T)[[1]][1]
    gene = substr(gene,2,50)
    return(gene)
  })
  genes = unlist(genes)
  baitSplits = lapply(unique(genes),function(gene){
    buffer = baitTable[which(genes==gene),]
    return(buffer[sort(as.numeric(buffer[,2]),index.return=T)$ix,])
  })
  toRemove = lapply(baitSplits,function(geneBaits){
      distances = lapply(2:dim(geneBaits)[1],function(iter){
    return(as.numeric(geneBaits[iter,2])-as.numeric(geneBaits[iter-1,3]))
  })
      
      if(length(which(distances>minDistance))<maxBaits&&length(which(distances>minDistance))>minBaits){return(which(distances<minDistance)+1)}
      else{
        if(length(which(distances>minDistance))<minBaits){
          if(length(distances)>minBaits)
          {          return(sort(unlist(distances),decreasing=T,index.return=T)$ix[(minBaits+1):length(distances)])
          }
          else{return(list())}
        }
          else
            {
        return(sort(unlist(distances),decreasing=T,index.return=T)$ix[maxBaits:length(distances)])
            }
      }
      
  })
  results = lapply(1:length(toRemove),function(iter){
    if(length(toRemove[[iter]])==0){return(baitSplits[[iter]])}
    else{return(baitSplits[[iter]][-toRemove[[iter]],])}
  })
  results = do.call(rbind,results)
  return(results)
}


gnames=read.table("~/Desktop/gene_names",header=F,sep="\t")
for(gene in gnames$V1) {exons=GetExonsFromName(gene); print(exons)}
