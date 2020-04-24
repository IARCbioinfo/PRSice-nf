rm(list=ls())

# Parameters -------------------------------------------------------------------
pheno_file=as.character(commandArgs(TRUE)[1])  
output_name=as.character(commandArgs(TRUE)[2])  

prs_scores_combined=read.table(pheno_file,header = F)
scores_files=list.files(path = ".",pattern = "all.score")
snps_files=list.files(path = ".",pattern = ".snp")
colnames(prs_scores_combined)=c("FID","IID","Status")

for(i in 1:length(scores_files)){
  file_name=scores_files[i]
  prs_scores=read.table(file_name,header = T,check.names = F)
  snp_file=read.table(snps_files[i],header = T,check.names = F)
  if(i==1){
    prs_scores=prs_scores[,-c(1:2),drop=F]
    prs_scores_combined=cbind(prs_scores_combined,prs_scores)
    snp_file_combined=snp_file
  }else{
    common_col=colnames(prs_scores)[colnames(prs_scores) %in% colnames(prs_scores_combined) & !(colnames(prs_scores) %in% c("FID","IID") )]
    unique_col=colnames(prs_scores)[!colnames(prs_scores) %in% colnames(prs_scores_combined) & !(colnames(prs_scores) %in% c("FID","IID") )]
    if(length(unique_col)!=0){
      prs_scores_combined=cbind(prs_scores_combined,prs_scores[,unique_col,drop=F])
    }
    if(length(common_col)!=0){
      prs_scores_combined[,common_col]=prs_scores_combined[,common_col]+prs_scores[,common_col]
    }
    snp_file_combined=rbind(snp_file_combined,snp_file)
  }
}

write.table(prs_scores_combined,file=paste0("prs_scores_",output_name,".txt"),quote = F,row.names = F,col.names = T)
write.table(snp_file_combined,file=paste0("snp_file_",output_name,".txt"),quote = F,row.names = F,col.names = T)