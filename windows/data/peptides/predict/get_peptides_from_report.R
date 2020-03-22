library(readr)
library(dplyr)
library(stringr)
source("../../../code/deepdetect/R/init.R")
proteins = read.fasta(list.files(pattern = '\\.fasta$'))
file_name<-list.files(pattern = 'PeptideReport\\.csv$')
file<-read.csv(file_name,1)

peptide<-file[,c("PG.ProteinAccessions","PEP.StrippedSequence","PEP.NrOfMissedCleavages")]
peptide[,"PG.ProteinAccessions"]<-str_split(file$PG.FASTAHeader," ",simplify = T)[,1]%>%str_sub(2)
peptide<-distinct(peptide,PEP.StrippedSequence,.keep_all = T)
names(peptide)<-c("Protein_Name","Sequence","miss")
digested<-peptide

peptides = data.frame(
  protein = as.character(digested$Protein_Name),
  sequence = as.character(digested$Sequence),
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

peptides = peptides[nchar(peptides$sequence) <= 50, ]

write.csv(
  peptides,
  sub('\\.PeptideReport\\.csv$', '\\1.peptide.csv', file_name),
  row.names = FALSE
)
