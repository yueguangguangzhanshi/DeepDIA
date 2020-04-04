metr_pkgs<-c('rJava','PKI','digest','RSQLite','jsonlite','dplyr','DT','V8','readr','tidyverse','mixOmics','gplots','pheatmap', 'plotly', 'ggplotify', 'dragulaR', 'factoextra', 'FactoMineR', 'ggsci', 'seqinr', 'ggrepel', 'ggplot2','cleaver','data.table','rjson')

for(i in 1:length(metr_pkgs)){
  
  library(metr_pkgs[i], character.only = TRUE,quietly=TRUE)
  
}

# DeepReportPredict<-function(ReportFile){
#   file<-read.csv(ReportFile,1)
#   peptide<-file[,c("PG.ProteinAccessions","PEP.StrippedSequence","PP.Charge")]
#   index1<-sample(1:nrow(peptide),size = 80000)
#   peptide_2charge<-peptide[index1,]
#   peptide_2charge<-peptide_2charge[peptide_2charge[,"PP.Charge"]==2,c("PG.ProteinAccessions","PEP.StrippedSequence")]
#   names(peptide_2charge)<-c("protein","sequence")
#   peptide_2charge<-peptide_2charge[unlist(lapply(peptide_2charge[,"sequence"],FUN = function(x){nchar(as.vector(x))[[1]][1]}))<=50,]
#   peptide_2charge = peptide_2charge[which(
#     !grepl('[^ACDEFGHIKLMNPQRSTVWY]', peptide_2charge$sequence)
#   ), ]
#   
#   # index2<-sample(1:nrow(peptide),size = 80000)
#   # peptide_3charge<-peptide[index2,]
#   peptide_3charge<-peptide[index1,]
#   peptide_3charge<-peptide_3charge[peptide_3charge[,"PP.Charge"]==3,c("PG.ProteinAccessions","PEP.StrippedSequence")]
#   names(peptide_3charge)<-c("protein","sequence")
#   peptide_3charge<-peptide_3charge[unlist(lapply(peptide_3charge[,"sequence"],FUN = function(x){nchar(as.vector(x))[[1]][1]}))<=50,]
#   peptide_3charge = peptide_3charge[which(
#     !grepl('[^ACDEFGHIKLMNPQRSTVWY]', peptide_3charge$sequence)
#   ), ]
#   
#   # index3<-sample(1:nrow(peptide),size = 80000)
#   # peptide_irt<-peptide[index3,]
#   peptide_irt<-peptide[index1,]
#   peptide_irt<-peptide_irt[,c("PG.ProteinAccessions","PEP.StrippedSequence")]
#   names(peptide_irt)<-c("protein","sequence")
#   peptide_irt<-peptide_irt[unlist(lapply(peptide_irt[,"sequence"],FUN = function(x){nchar(as.vector(x))[[1]][1]}))<=50,]
#   peptide_irt = peptide_irt[which(
#     !grepl('[^ACDEFGHIKLMNPQRSTVWY]', peptide_irt$sequence)
#   ), ]
#   
#   write.csv(peptide_irt,"../data/deepdia/data/models/irt/Report.peptide.csv",row.names =F)
#   write.csv(peptide_2charge,"../data/deepdia/data/models/charge2/Report_charge2.peptide.csv",row.names =F)
#   write.csv(peptide_3charge,"../data/deepdia/data/models/charge3/Report_charge3.peptide.csv",row.names =F)
#   
#   setwd("../data/deepdia/data/models/charge2/")
#   system("~/anaconda3/envs/tensorflow/bin/python ../../../code/deepms2/py/predict.py")
#   system("cp *.ions.json ../../library/")
#   setwd("../charge3/")
#   system("~/anaconda3/envs/tensorflow/bin/python ../../../code/deepms2/py/predict.py")
#   system("cp *.ions.json ../../library/")
#   setwd("../irt/")
#   system("~/anaconda3/envs/tensorflow/bin/python ../../../code/deeprt/py/predict.py")
#   system("cp *.csv ../../library/")
#   
#   setwd("../../library/")
#   source("../../code/init.R")
#   create.assays = function(peptides, ions, irt, precursorCharge) {
#     assay = lapply(1:length(ions), function(i) {
#       sequence = ions[[i]]$peptide
#       precursorMz = (peptide.mass.1(sequence, carbamidomethyl = TRUE) + n.terminus.mass.monoisotopic['H'] * precursorCharge) / precursorCharge
#       
#       fragmentIntensity = c(
#         ions[[i]]$ions$b1,
#         ions[[i]]$ions$bn1,
#         ions[[i]]$ions$bo1,
#         ions[[i]]$ions$b2,
#         ions[[i]]$ions$bn2,
#         ions[[i]]$ions$bo2,
#         ions[[i]]$ions$y1,
#         ions[[i]]$ions$yn1,
#         ions[[i]]$ions$yo1,
#         ions[[i]]$ions$y2,
#         ions[[i]]$ions$yn2,
#         ions[[i]]$ions$yo2
#       )
#       fragmentType = c(
#         rep('b', length(ions[[i]]$ions$b1)),
#         rep('b', length(ions[[i]]$ions$bn1)),
#         rep('b', length(ions[[i]]$ions$bo1)),
#         rep('b', length(ions[[i]]$ions$b2)),
#         rep('b', length(ions[[i]]$ions$bn2)),
#         rep('b', length(ions[[i]]$ions$bo2)),
#         rep('y', length(ions[[i]]$ions$y1)),
#         rep('y', length(ions[[i]]$ions$yn1)),
#         rep('y', length(ions[[i]]$ions$yo1)),
#         rep('y', length(ions[[i]]$ions$y2)),
#         rep('y', length(ions[[i]]$ions$yn2)),
#         rep('y', length(ions[[i]]$ions$yo2))
#       )
#       fragmentNumber = c(
#         1:length(ions[[i]]$ions$b1),
#         1:length(ions[[i]]$ions$bn1),
#         1:length(ions[[i]]$ions$bo1),
#         1:length(ions[[i]]$ions$b2),
#         1:length(ions[[i]]$ions$bn2),
#         1:length(ions[[i]]$ions$bo2),
#         1:length(ions[[i]]$ions$y1),
#         1:length(ions[[i]]$ions$yn1),
#         1:length(ions[[i]]$ions$yo1),
#         1:length(ions[[i]]$ions$y2),
#         1:length(ions[[i]]$ions$yn2),
#         1:length(ions[[i]]$ions$yo2)
#       )
#       fragmentCharge = c(
#         rep('1', length(ions[[i]]$ions$b1)),
#         rep('1', length(ions[[i]]$ions$bn1)),
#         rep('1', length(ions[[i]]$ions$bo1)),
#         rep('2', length(ions[[i]]$ions$b2)),
#         rep('2', length(ions[[i]]$ions$bn2)),
#         rep('2', length(ions[[i]]$ions$bo2)),
#         rep('1', length(ions[[i]]$ions$y1)),
#         rep('1', length(ions[[i]]$ions$yn1)),
#         rep('1', length(ions[[i]]$ions$yo1)),
#         rep('2', length(ions[[i]]$ions$y2)),
#         rep('2', length(ions[[i]]$ions$yn2)),
#         rep('2', length(ions[[i]]$ions$yo2))
#       )
#       fragmentLossType = c(
#         rep('noloss', length(ions[[i]]$ions$b1)),
#         rep('NH3', length(ions[[i]]$ions$bn1)),
#         rep('H2O', length(ions[[i]]$ions$bo1)),
#         rep('noloss', length(ions[[i]]$ions$b2)),
#         rep('NH3', length(ions[[i]]$ions$bn2)),
#         rep('H2O', length(ions[[i]]$ions$bo2)),
#         rep('noloss', length(ions[[i]]$ions$y1)),
#         rep('NH3', length(ions[[i]]$ions$yn1)),
#         rep('H2O', length(ions[[i]]$ions$yo1)),
#         rep('noloss', length(ions[[i]]$ions$y2)),
#         rep('NH3', length(ions[[i]]$ions$yn2)),
#         rep('H2O', length(ions[[i]]$ions$yo2))
#       )
#       fragmentMz = fragment.ions.mz.1(
#         sequence, 
#         type = c('b', 'y'), 
#         charge = c(1, 2),
#         loss = c('NH3', 'H2O'),
#         carbamidomethyl = TRUE
#       )[paste0(
#         fragmentType,
#         c('noloss' = '', 'NH3' = '*', 'H2O' = 'o')[fragmentLossType],
#         fragmentNumber,
#         '^', fragmentCharge
#       )]
#       fragment = data.frame(
#         fragmentType, fragmentNumber, fragmentCharge, fragmentLossType,
#         fragmentMz, fragmentIntensity,
#         stringsAsFactors = FALSE
#       )
#       fragment = fragment[fragment$fragmentIntensity > 0, ]
#       if (nrow(fragment) == 0)
#         return(NULL)
#       
#       
#       iRT = irt$irt[match(sequence, irt$sequence)]
#       
#       # proteinId = report$PG.ProteinAccessions[match(sequence, peptides$PEP.StrippedSequence)]
#       proteinId = peptides$protein[match(sequence, peptides$sequence)]
#       
#       data.frame(
#         PrecursorMz = precursorMz,
#         FragmentMz = fragment$fragmentMz,
#         iRT = iRT,
#         #iRTSourceSpecific,
#         RelativeFragmentIntensity = fragment$fragmentIntensity,
#         StrippedSequence = sequence,
#         PrecursorCharge = precursorCharge,
#         FragmentType = fragment$fragmentType,
#         FragmentNumber = fragment$fragmentNumber,
#         FragmentCharge = fragment$fragmentCharge,
#         ProteinId = proteinId,
#         #ProteinName,
#         #ProteinDescription,
#         ModifiedSequence = paste0('_', gsub('C', 'C[Carbamidomethyl (C)]', sequence), '_'),
#         #LabelModifiedSequence,
#         FragmentLossType = fragment$fragmentLossType,
#         #ExcludeFromAssay,
#         #IsProteotypic,
#         stringsAsFactors = FALSE
#       )
#     })
#   }
#   
#   file = list.files(pattern = '\\.peptide\\.csv$')[1]
#   peptides = read_csv(file)
#   irt = read_csv(sub('\\.peptide\\.csv$', '.prediction.irt.csv', file))
#   
#   charges = 2:3
#   
#   assays = lapply(charges, function(charge) {
#     ions = rjson::fromJSON(file = sub('\\.peptide\\.csv$', paste0('_charge', charge, '.prediction.ions.json'), file))
#     assay = create.assays(peptides, ions, irt, precursorCharge = charge)
#   })
#   data<-do.call(rbind, do.call(c, assays))
#   dir.create(jobid)
#   #return(data)
#   setwd(jobid)
#   write.csv(
#     do.call(rbind, do.call(c, assays)),
#     sub('\\.peptide\\.csv$', '.prediction.library.csv', file),
#     row.names = FALSE
#   )
# }

UpdateStatus <- function(status,id) {
  projects_list <- dbConnect(RSQLite::SQLite(), "~/shiny-server/project/db/projects_list.db", flags = SQLITE_RWC)	
  dbSendStatement(projects_list, paste("UPDATE Project SET status = '",status,"' WHERE id='",id,"'",sep=""))
  dbDisconnect(projects_list)
}

setwd("~/shiny-server/project/deepDIA")
projects_list <- dbConnect(RSQLite::SQLite(), "../db/projects_list.db")		
project_status_table <- as.data.frame(dbGetQuery(projects_list, paste("SELECT * FROM Project WHERE user='","deepDIA","'",sep="")))
dbDisconnect(projects_list)

#在运行项目检查
processing_projects<-project_status_table[which(project_status_table$status<=2&project_status_table$status>1),]
if (nrow(processing_projects)>0) {
  cat("有正在运行的fasta预测任务","\n")
  quit()
}

#unfinished_train_projects <- project_status_table[which((project_status_table$status==0)&(project_status_table$name=="train")),]
unfinished_fastapredict_projects <- project_status_table[which((project_status_table$status==0)&(grepl("^fasta_predict_",project_status_table$name))),]

#modelid<<-unfinished_train_projects[which(unfinished_train_projects$time==min(unfinished_train_projects$time)),"id"]
tryCatch({
  jobid<<-unfinished_fastapredict_projects[which(unfinished_fastapredict_projects$time==min(unfinished_fastapredict_projects$time)),"id"]
  modelid<<-unfinished_fastapredict_projects[which(unfinished_fastapredict_projects$time==min(unfinished_fastapredict_projects$time)),"name"]%>%str_split("\\_")%>%.[[1]]%>%.[3]
},
  error=function(e){
    cat(conditionMessage(e),"\n","没有fasta预测任务\n")
    quit()
    }
  #finally={print("有fasta预测任务")}
)
if (length(modelid)>0&length(jobid)>0) {
  t1=proc.time()
  UpdateStatus(2,jobid)
  setwd(paste0("../data/deepdia/",modelid,"/library/",jobid))
  test_index<-list.dirs()%>%sub("[\\./test]+","",.)%>%as.numeric%>%na.omit%>%max
  for (i in 1:test_index) {
    setwd(paste0("test",i))
    system(paste("Rscript ../../../../code/fasta_predict.R",paste0("test",i),jobid,paste0("models",i)),wait = F)
    setwd("../")
  }
  i=2
  while (i>(2-0.1*test_index)) {
    projects_list <- dbConnect(RSQLite::SQLite(), "~/shiny-server/project/db/projects_list.db")		
    project_status_table <- as.data.frame(dbGetQuery(projects_list, paste("SELECT * FROM Project WHERE user='","deepDIA","'",sep="")))
    dbDisconnect(projects_list)
    i=project_status_table[which(project_status_table[,"id"]==jobid),"status"]
    if (i==3) {
      stop
    }
    Sys.sleep(300)
  }
  library_partly<-list()
  library_partly<-lapply(1:test_index,function(x){
    library_partly[[x]]<-read_csv(paste0("test",x,"/Report.prediction.library.csv"))
  })
  library<-do.call(rbind,library_partly)
  write_csv(library,"Report.prediction.library.csv")
  system("zip -r Report.prediction.library.zip *.prediction.library.csv")
  file.create("success.txt")
  setwd("~/shiny-server/project/deepDIA")
  if (file.exists(paste0("../data/deepdia/",modelid,"/library/",jobid,"/success.txt"))) {
    UpdateStatus(1,jobid)
  }else{
    UpdateStatus(3,jobid)
    stop
  }
  t2=proc.time()
  t=t2-t1
  print(paste0('执行时间：',t[3][[1]]/60/60,' hours'))
}else{
  stop
}
