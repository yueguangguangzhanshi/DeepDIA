file_name<-list.files(pattern = '\\.Peptides\\.csv$')
file<-read.csv(file_name,1)

if (nrow(file)>80000) {
  index1<-sample(1:nrow(file),size = 80000)
}else{
  index1<-1:nrow(file)
}

peptide_2charge<-file[index1,]
peptide_2charge<-peptide_2charge[unlist(lapply(peptide_2charge[,"sequence"],FUN = function(x){nchar(as.vector(x))[[1]][1]}))<=50,]
peptide_2charge = peptide_2charge[which(
  !grepl('[^ACDEFGHIKLMNPQRSTVWY]', peptide_2charge$sequence)
), ]

peptide_3charge<-file[index1,]
peptide_3charge<-peptide_3charge[unlist(lapply(peptide_3charge[,"sequence"],FUN = function(x){nchar(as.vector(x))[[1]][1]}))<=50,]
peptide_3charge = peptide_3charge[which(
  !grepl('[^ACDEFGHIKLMNPQRSTVWY]', peptide_3charge$sequence)
), ]

peptide_irt<-file[index1,]
peptide_irt<-peptide_irt[unlist(lapply(peptide_irt[,"sequence"],FUN = function(x){nchar(as.vector(x))[[1]][1]}))<=50,]
peptide_irt = peptide_irt[which(
  !grepl('[^ACDEFGHIKLMNPQRSTVWY]', peptide_irt$sequence)
), ]

write.csv(peptide_irt,paste0("irt/",format(Sys.time(), "%Y%m%d"),".peptide.csv"),row.names =F)
write.csv(peptide_2charge,paste0("charge2/",format(Sys.time(), "%Y%m%d"),"_charge2.peptide.csv"),row.names =F)
write.csv(peptide_3charge,paste0("charge3/",format(Sys.time(), "%Y%m%d"),"_charge3.peptide.csv"),row.names =F)
