library(data.table)

snp_file=as.character(commandArgs(TRUE)[1])  

snp_list=fread(snp_file,header = T,stringsAsFactors = F)

#awk '{ if (($4=="T" && $5=="A")||($4=="A" && $5=="T")||($4=="C" && $5=="G")||($4=="G" && $5=="C"))"
ambiguous_snps=snp_list[which( (snp_list$A1=="A" & snp_list$A2=="T") | (snp_list$A1=="T" & snp_list$A2=="A") | (snp_list$A1=="G" & snp_list$A2=="C") | (snp_list$A1=="C" & snp_list$A2=="G")), ]
head(ambiguous_snps)
low_qual_snps=snp_list[which(snp_list$QUAL<0.3),]
head(low_qual_snps)
snp_list_filtered=snp_list[which(!snp_list$SNP %in% c(low_qual_snps$SNP,ambiguous_snps$SNP)),]
nrow(snp_list_filtered)
if(nrow(snp_list_filtered)==0){
  head(snp_list_filtered)
  write.table(snp_list_filtered,file=snp_file,row.names = F,col.names = T,sep="\t",quote = F)
}
