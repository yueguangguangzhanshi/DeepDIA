library(cleaver)
library(dplyr)
library(data.table)
library(tibble)
fasta.file<-list.files(pattern = '*.fasta$')
proteins = read.fasta(fasta.file)
fasta<-list()
for (i in 1:length(proteins)) {
  fasta[[i]]<-sapply(proteins[[i]]$sequence,FUN = function(x)cleave(x,"trypsin",missedCleavages=0))
  names(fasta[[i]])<-proteins[[i]]$name
}

fasta2<-lapply(fasta, function(x){
  data<-as.data.frame(x,stringsAsFactors=F,check.names=F)
  data[,2]<-colnames(data)[1]
  return(data)
})
df <- rbindlist(fasta2,use.names=FALSE)%>%.[,c(2,1)]
names(df)<-c("Protein_Name","Sequence")
df<-filter(df,nchar(Sequence)>7&nchar(Sequence)<=50)%>%add_column(miss=0)

#write.table(df,sub('\\.fasta$', '_digested\\.txt', fasta.file),row.names = FALSE,sep = "\t")

digested=df
peptides = data.frame(
  protein = digested$Protein_Name,
  sequence = digested$Sequence,
  miss = digested$miss,
  stringsAsFactors = FALSE
)

proteins = data.frame(
  accession = sapply(proteins, function(x) x$name),
  sequence = sapply(proteins, function(x) x$sequence),
  description = sapply(proteins, function(x) sub('^>[^ ]* ', '', x$description)),
  stringsAsFactors = FALSE
)

cleavages = find.cleavageWindow(peptides, proteins)

peptides = cbind(peptides, cleavages)

peptides = peptides[which(
  !grepl('[^ACDEFGHIKLMNPQRSTVWY]', peptides$sequence) &
    !grepl('[^_ACDEFGHIKLMNPQRSTVWY]', peptides$nTerminal) &
    !grepl('[^_ACDEFGHIKLMNPQRSTVWY]', peptides$cTerminal)
), ]

write.csv(
  peptides,
  sub('(.*)\\.fasta$', '\\1.peptide.csv', fasta.file),
  row.names = FALSE
)