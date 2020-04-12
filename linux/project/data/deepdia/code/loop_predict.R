metr_pkgs<-c('rJava','PKI','digest','RSQLite','jsonlite','dplyr','DT','V8','readr','tidyverse','mixOmics','gplots','pheatmap', 'plotly', 'ggplotify', 'dragulaR', 'factoextra', 'FactoMineR', 'ggsci', 'seqinr', 'ggrepel', 'ggplot2','cleaver','data.table','rjson')

for(i in 1:length(metr_pkgs)){
  
  library(metr_pkgs[i], character.only = TRUE,quietly=TRUE)
  
}

UpdateStatus <- function(status,id) {
  projects_list <- dbConnect(RSQLite::SQLite(), "~/shiny-server/project/db/projects_list.db", flags = SQLITE_RWC)	
  dbSendStatement(projects_list, paste("UPDATE Project SET status = '",status,"' WHERE id='",id,"'",sep=""))
  dbDisconnect(projects_list)
}

setwd("~/shiny-server/project/deepDIA")
projects_list <- dbConnect(RSQLite::SQLite(), "../db/projects_list.db")		
project_status_table <- as.data.frame(dbGetQuery(projects_list, paste("SELECT * FROM Project WHERE name GLOB'","fasta_predict*","'",sep="")))
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
  modelid<<-unfinished_fastapredict_projects[which(unfinished_fastapredict_projects$time==min(unfinished_fastapredict_projects$time)),"name"]%>%str_split("\\_",n=3)%>%.[[1]]%>%.[3]
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
    system(paste("/home/deepdia/anaconda3/lib/R/bin/Rscript ../../../../code/fasta_predict.R",paste0("test",i),jobid,paste0("models",i)),wait = F)
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
