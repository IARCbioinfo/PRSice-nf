library(data.table)
library(parallel)
library(dplyr)

snp_file=as.character(commandArgs(TRUE)[1])  
info=as.character(commandArgs(TRUE)[2])  
chr=as.character(commandArgs(TRUE)[3])  

snp_list=fread(snp_file,header = F,stringsAsFactors = F)
colnames(snp_list)=c("SNP","CHR","BP","A1","A2","STAT","SE","PVAL")
info_file=fread(info,header = F,stringsAsFactors = F)
colnames(info_file)=c("alternate_ids","rsid", "chromosome","position","number_of_alleles","first_allele","alternative_alleles" )

custom_f=function(chr,pos,a1,a2){
  chr_pos=paste(c(chr,pos),collapse = ":")
  alleles=sort(c(a1,a2))
  return(paste(chr_pos,alleles[1],alleles[2],sep=":"))
}

snp_list[, new_ID := mcmapply(custom_f,CHR,BP,A1,A2,mc.cores=3L)]
info_file[, new_ID := mcmapply(custom_f,chromosome,position,first_allele,alternative_alleles,mc.cores=3L)]

snp_list_common=inner_join(snp_list,info_file[,.SD,.SDcols=c("new_ID","alternate_ids")],by="new_ID")
PRSice_input=snp_list_common[,c(10,2:8)]
colnames(PRSice_input)=c("SNP","CHR","BP","A1","A2","STAT","SE","PVAL")
write.table(PRSice_input,paste0("PRSice_input_chr_",chr,".txt"),row.names = F,col.names = T,quote = F,sep=" ")

# update fliped snps