## http://www.bioconductor.org/packages/release/bioc/vignettes/biomaRt/inst/doc/biomaRt.pdf

library("biomaRt")
listMarts(host="www.ensembl.org")
ensembl=useMart(host="www.ensembl.org", "ENSEMBL_MART_ENSEMBL" )
listDatasets(ensembl)
ensembl = useMart(host="www.ensembl.org","ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl")
listFilters(ensembl)
affyids=c("202763_at","209310_s_at","207500_at")
geneNames=c("TRAF3","TP53","NRAS","KRAS","FAM46C")
getBM(attributes=c('affy_hg_u133_plus_2', 'hgnc_symbol','ensembl_gene_id' ,'chromosome_name','start_position','end_position', 'band'), filters = 'affy_hg_u133_plus_2', values = affyids, mart = ensembl)
getBM(attributes=c('affy_hg_u133_plus_2', 'hgnc_symbol','ensembl_gene_id' ,'chromosome_name','start_position','end_position', 'band'),
      filters = 'hgnc_symbol', values = geneNames, mart = ensembl)


chr=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y","MT")
all=getBM(attributes=c('chromosome_name','start_position','end_position', 'hgnc_symbol','ensembl_gene_id','strand','wikigene_name','external_gene_name'),
          mart = ensembl)
all=all[all$chromosome_name %in% chr, ]


#Older versions
ensembl_74 = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="dec2013.archive.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
ensembl_75 = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="feb2014.archive.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
all74=getBM(attributes=c('chromosome_name','start_position','end_position','hgnc_symbol','ensembl_gene_id','strand'),
          mart = ensembl_74)

all75=getBM(attributes=c('chromosome_name','start_position','end_position','hgnc_symbol','ensembl_gene_id','strand'),
            filters = 'hgnc_symbol', values = geneNames, mart = ensembl_75)
all74=all74[all74$chromosome_name %in% chr, ]