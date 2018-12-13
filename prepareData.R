suppressMessages({
    library(readr)
    library(dplyr)
    library(tidyr)
    library(VariantAnnotation)
    library(rtracklayer)
})

args <- commandArgs(trailingOnly = TRUE)

# Input files: 

fileSNP = "Data/snpid_imp.tsv" # SNPs (snpid,chr,pos) 
fileRES = "Data/chr02.cis.tsv" # MatrixEQTL results
fileGTF = "Data/gencode.v12.annotationLiftover_02.gtf.gz"

gtf <- readGFF(fileGTF)
g <- gtf[gtf$type=='gene', c("gene_id","gene_name")]

snps = read_tsv(fileSNP)

### MatrixEQTL dataframe, columns used: snpid gene_name pValue
res = read_tsv(fileRES)
names(res) =c("snpid","gene_name", "beta","tStat","pValue","FDR")

# if gene_id use gene_name
if(sum(grepl('ENSG', res$gene_name)) > length(res$gene_name)/2){
    names(res)[2]='gene_id'
    res = res %>% left_join(g, by=c('gene_id'))
    res$gene_id <- NULL 
    res = res[!is.na(res$gene_name),]
}
# Only SNPs in results
snps = snps %>% filter(snpid %in% res$snpid)
