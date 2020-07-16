library(data.table)
library(parallel)
library(dplyr)

snp_file=as.character(commandArgs(TRUE)[1])  
info=as.character(commandArgs(TRUE)[2])  
chr=as.character(commandArgs(TRUE)[3])  


# Functions ---------------------------------------------------------------

fliped_alleles=c("A","C","T","G")
names(fliped_alleles)=c("T","G","A","C")
strand_flipping=function(allele){
  allele_split=unlist(strsplit(allele,""))
  return(paste(fliped_alleles[allele_split],collapse = ""))
}
get_newID=function(chr,pos,ref_old,alt_old){
  alleles_1=sort(c(ref_old,alt_old))
  alleles_2=sort(c(strand_flipping(ref_old),strand_flipping(alt_old)))
  alleles=unlist(strsplit(sort(c(paste(alleles_1,collapse = ":"),paste(alleles_2,collapse = ":")))[1],":")) 
  return(paste(chr,pos,alleles[1],alleles[2],sep=":"))
}
get_swap=function(ref_old,alt_old,alt_new){
  swap=F
  if(alt_new==alt_old){
    swap=F
  }else if(alt_new!=alt_old & ref_old==alt_new){
    swap=T
  }else if(alt_new!=alt_old & ref_old!=alt_new){
    if(ref_old %in% c(alt_new,strand_flipping(alt_new)) ){
      swap=T
    }else{
      swap=F
    }
  }
  return(swap)
}
get_strand=function(ref_old,alt_old,alt_new){
  strand=F
  if(alt_new==alt_old){
    strand=F
  }else if(alt_new!=alt_old & ref_old==alt_new){
    strand=F
  }else if(alt_new!=alt_old & ref_old!=alt_new){
    strand=T
  }
  return(strand)
}


# Match IDs ---------------------------------------------------------------

snp_list=fread(snp_file,header = F,stringsAsFactors = F)
colnames(snp_list)=c("SNP","CHR","BP","A1","A2","STAT","SE","PVAL","QUAL")
info_file=fread(info,header = F,stringsAsFactors = F)
colnames(info_file)=c("alternate_ids","rsid", "chromosome","position","number_of_alleles","first_allele","alternative_alleles" )

snp_list[, new_ID := mcmapply(get_newID,CHR,BP,A1,A2,mc.cores=3L)]
info_file[, new_ID := mcmapply(get_newID,chromosome,position,first_allele,alternative_alleles,mc.cores=3L)]

snp_list_common=inner_join(snp_list,info_file[,.SD,.SDcols=c("new_ID","alternate_ids")],by="new_ID")
PRSice_input=snp_list_common[,c(11,2:9)]
colnames(PRSice_input)=c("SNP","CHR","BP","A1","A2","STAT","SE","PVAL","QUAL")
write.table(PRSice_input,paste0("PRSice_input_chr_",chr,".txt"),row.names = F,col.names = T,quote = F,sep=" ")

# update fliped snps