file_name<-list.files(pattern = '\\.peptide\\.csv$')
file<-read.csv(file_name,1)
file<-file[order(file["detectability"],decreasing = T),]
if (nrow(file)>80000) {
  index1<-1:80000
}else{
  index1<-1:nrow(file)
}

peptide_2charge<-file[index1,]

peptide_3charge<-file[index1,]

peptide_irt<-file[index1,]

write.csv(peptide_irt,paste0("irt/",format(Sys.time(), "%Y%m%d"),".peptide.csv"),row.names =F)
write.csv(peptide_2charge,paste0("charge2/",format(Sys.time(), "%Y%m%d"),"_charge2.peptide.csv"),row.names =F)
write.csv(peptide_3charge,paste0("charge3/",format(Sys.time(), "%Y%m%d"),"_charge3.peptide.csv"),row.names =F)
