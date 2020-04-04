metr_pkgs<-c('rJava','PKI','digest','RSQLite','jsonlite','dplyr','DT','V8','readr','tidyverse','mixOmics','gplots','pheatmap', 'plotly', 'ggplotify', 'dragulaR', 'factoextra', 'FactoMineR', 'ggsci', 'seqinr', 'ggrepel', 'ggplot2','cleaver','data.table','rjson')

for(i in 1:length(metr_pkgs)){
  
  library(metr_pkgs[i], character.only = TRUE,quietly=TRUE)
  
}

args=commandArgs(T)

DeepFastaPredict<-function(FastaFile,testdir,models){
  source("../../../../code/deepdetect/R/init.R")
  proteins = read.fasta(FastaFile)
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
  df<-filter(df,nchar(Sequence)>=7&nchar(Sequence)<=50)%>%add_column(miss=0)
  
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
    paste0("../../../",models,"/detectability/validated.peptide.csv"),
    row.names = FALSE
  )
  setwd(paste0("../../../",models,"/detectability/"))
  system("~/anaconda3/envs/tensorflow/bin/python ../../../code/deepdetect/py/predict.py",wait = T)
  source("../../../code/deepdetect/R/filter_peptides_by_detectability.R")
  
  file_name<-list.files(pattern = '(.)*detectability(.)*\\.peptide\\.csv$')
  file<-read.csv(file_name,1)
  file<-file[order(file["detectability"],decreasing = T),]
  if (nrow(file)>80000) {
    #index1<-sample(1:nrow(file),size = 80000)
    # index1<-1:80000
    index1<-1:nrow(file)
  }else{
    index1<-1:nrow(file)
  }
  peptide_2charge<-file[index1,]
  peptide_3charge<-file[index1,]
  peptide_irt<-file[index1,]
  write.csv(peptide_irt,paste0("../irt/Report.peptide.csv"),row.names =F)
  write.csv(peptide_2charge,paste0("../charge2/Report_charge2.peptide.csv"),row.names =F)
  write.csv(peptide_3charge,paste0("../charge3/Report_charge3.peptide.csv"),row.names =F)
  setwd("../charge2/")
  system("~/anaconda3/envs/tensorflow/bin/python ../../../code/deepms2/py/predict.py")
  system(paste0("cp *.ions.json ../../library/",jobid,"/",testdir))
  setwd("../charge3/")
  system("~/anaconda3/envs/tensorflow/bin/python ../../../code/deepms2/py/predict.py")
  system(paste0("cp *.ions.json ../../library/",jobid,"/",testdir))
  setwd("../irt/")
  system("~/anaconda3/envs/tensorflow/bin/python ../../../code/deeprt/py/predict.py")
  system(paste0("cp *.csv ../../library/",jobid,"/",testdir))
  
  setwd(paste0("../../library/",jobid,"/",testdir))
  source("../../../../code/init.R")
  create.assays = function(peptides, ions, irt, precursorCharge) {
    assay = lapply(1:length(ions), function(i) {
      sequence = ions[[i]]$peptide
      precursorMz = (peptide.mass.1(sequence, carbamidomethyl = TRUE) + n.terminus.mass.monoisotopic['H'] * precursorCharge) / precursorCharge
      
      fragmentIntensity = c(
        ions[[i]]$ions$b1,
        ions[[i]]$ions$bn1,
        ions[[i]]$ions$bo1,
        ions[[i]]$ions$b2,
        ions[[i]]$ions$bn2,
        ions[[i]]$ions$bo2,
        ions[[i]]$ions$y1,
        ions[[i]]$ions$yn1,
        ions[[i]]$ions$yo1,
        ions[[i]]$ions$y2,
        ions[[i]]$ions$yn2,
        ions[[i]]$ions$yo2
      )
      fragmentType = c(
        rep('b', length(ions[[i]]$ions$b1)),
        rep('b', length(ions[[i]]$ions$bn1)),
        rep('b', length(ions[[i]]$ions$bo1)),
        rep('b', length(ions[[i]]$ions$b2)),
        rep('b', length(ions[[i]]$ions$bn2)),
        rep('b', length(ions[[i]]$ions$bo2)),
        rep('y', length(ions[[i]]$ions$y1)),
        rep('y', length(ions[[i]]$ions$yn1)),
        rep('y', length(ions[[i]]$ions$yo1)),
        rep('y', length(ions[[i]]$ions$y2)),
        rep('y', length(ions[[i]]$ions$yn2)),
        rep('y', length(ions[[i]]$ions$yo2))
      )
      fragmentNumber = c(
        1:length(ions[[i]]$ions$b1),
        1:length(ions[[i]]$ions$bn1),
        1:length(ions[[i]]$ions$bo1),
        1:length(ions[[i]]$ions$b2),
        1:length(ions[[i]]$ions$bn2),
        1:length(ions[[i]]$ions$bo2),
        1:length(ions[[i]]$ions$y1),
        1:length(ions[[i]]$ions$yn1),
        1:length(ions[[i]]$ions$yo1),
        1:length(ions[[i]]$ions$y2),
        1:length(ions[[i]]$ions$yn2),
        1:length(ions[[i]]$ions$yo2)
      )
      fragmentCharge = c(
        rep('1', length(ions[[i]]$ions$b1)),
        rep('1', length(ions[[i]]$ions$bn1)),
        rep('1', length(ions[[i]]$ions$bo1)),
        rep('2', length(ions[[i]]$ions$b2)),
        rep('2', length(ions[[i]]$ions$bn2)),
        rep('2', length(ions[[i]]$ions$bo2)),
        rep('1', length(ions[[i]]$ions$y1)),
        rep('1', length(ions[[i]]$ions$yn1)),
        rep('1', length(ions[[i]]$ions$yo1)),
        rep('2', length(ions[[i]]$ions$y2)),
        rep('2', length(ions[[i]]$ions$yn2)),
        rep('2', length(ions[[i]]$ions$yo2))
      )
      fragmentLossType = c(
        rep('noloss', length(ions[[i]]$ions$b1)),
        rep('NH3', length(ions[[i]]$ions$bn1)),
        rep('H2O', length(ions[[i]]$ions$bo1)),
        rep('noloss', length(ions[[i]]$ions$b2)),
        rep('NH3', length(ions[[i]]$ions$bn2)),
        rep('H2O', length(ions[[i]]$ions$bo2)),
        rep('noloss', length(ions[[i]]$ions$y1)),
        rep('NH3', length(ions[[i]]$ions$yn1)),
        rep('H2O', length(ions[[i]]$ions$yo1)),
        rep('noloss', length(ions[[i]]$ions$y2)),
        rep('NH3', length(ions[[i]]$ions$yn2)),
        rep('H2O', length(ions[[i]]$ions$yo2))
      )
      fragmentMz = fragment.ions.mz.1(
        sequence, 
        type = c('b', 'y'), 
        charge = c(1, 2),
        loss = c('NH3', 'H2O'),
        carbamidomethyl = TRUE
      )[paste0(
        fragmentType,
        c('noloss' = '', 'NH3' = '*', 'H2O' = 'o')[fragmentLossType],
        fragmentNumber,
        '^', fragmentCharge
      )]
      fragment = data.frame(
        fragmentType, fragmentNumber, fragmentCharge, fragmentLossType,
        fragmentMz, fragmentIntensity,
        stringsAsFactors = FALSE
      )
      fragment = fragment[fragment$fragmentIntensity > 0, ]
      if (nrow(fragment) == 0)
        return(NULL)
      
      
      iRT = irt$irt[match(sequence, irt$sequence)]
      
      # proteinId = report$PG.ProteinAccessions[match(sequence, peptides$PEP.StrippedSequence)]
      proteinId = peptides$protein[match(sequence, peptides$sequence)]
      
      data.frame(
        PrecursorMz = precursorMz,
        FragmentMz = fragment$fragmentMz,
        iRT = iRT,
        #iRTSourceSpecific,
        RelativeFragmentIntensity = fragment$fragmentIntensity,
        StrippedSequence = sequence,
        PrecursorCharge = precursorCharge,
        FragmentType = fragment$fragmentType,
        FragmentNumber = fragment$fragmentNumber,
        FragmentCharge = fragment$fragmentCharge,
        ProteinId = proteinId,
        #ProteinName,
        #ProteinDescription,
        ModifiedSequence = paste0('_', gsub('C', 'C[Carbamidomethyl (C)]', sequence), '_'),
        #LabelModifiedSequence,
        FragmentLossType = fragment$fragmentLossType,
        #ExcludeFromAssay,
        #IsProteotypic,
        stringsAsFactors = FALSE
      )
    })
  }
  
  file = list.files(pattern = '\\.peptide\\.csv$')[1]
  peptides = read_csv(file)
  irt = read_csv(sub('\\.peptide\\.csv$', '.prediction.irt.csv', file))
  
  charges = 2:3
  
  assays = lapply(charges, function(charge) {
    ions = rjson::fromJSON(file = sub('\\.peptide\\.csv$', paste0('_charge', charge, '.prediction.ions.json'), file))
    assay = create.assays(peptides, ions, irt, precursorCharge = charge)
  })
  #return(data)
  write.csv(
    do.call(rbind, do.call(c, assays)),
    sub('\\.peptide\\.csv$', '.prediction.library.csv', file),
    row.names = FALSE
  )
}

UpdateStatus <- function(status,id) {
  projects_list <- dbConnect(RSQLite::SQLite(), "~/shiny-server/project/db/projects_list.db", flags = SQLITE_RWC)	
  dbSendStatement(projects_list, paste("UPDATE Project SET status = '",status,"' WHERE id='",id,"'",sep=""))
  dbDisconnect(projects_list)
}


FastaFile<-list.files(pattern = "\\.fasta$")
testdir<-args[1]
jobid<-args[2]
models<-args[3]
system(paste0("rm -rf ../../../",models,"/charge2/*.csv"))
system(paste0("rm -rf ../../../",models,"/charge2/*.json"))
system(paste0("rm -rf ../../../",models,"/charge3/*.csv"))
system(paste0("rm -rf ../../../",models,"/charge3/*.json"))
system(paste0("rm -rf ../../../",models,"/detectability/*.csv"))
system(paste0("rm -rf ../../../",models,"/irt/*.csv"))
tryCatch({
  DeepFastaPredict(FastaFile,testdir,models)
  projects_list <- dbConnect(RSQLite::SQLite(), "~/shiny-server/project/db/projects_list.db")		
  project_status_table <- as.data.frame(dbGetQuery(projects_list, paste("SELECT * FROM Project WHERE user='","deepDIA","'",sep="")))
  dbDisconnect(projects_list)
  i=project_status_table[which(project_status_table[,"id"]==jobid),"status"]
  UpdateStatus(i-0.1,jobid)
},error=function(e){
  cat(conditionMessage(e),"\n","建库失败\n")
  UpdateStatus(3,jobid)
  quit()
})