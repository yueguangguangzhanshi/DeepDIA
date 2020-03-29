metr_pkgs<-c('shinyjs','shiny','rJava','PKI','digest','shinydashboard','shinycssloaders','RSQLite','jsonlite','dplyr','DT','V8','shinyBS','readr','tidyverse','mixOmics','gplots','pheatmap', 'plotly', 'ggplotify', 'dragulaR', 'factoextra', 'FactoMineR', 'ggsci', 'seqinr', 'ggrepel', 'ggplot2','mailR','shinyAce','cleaver','data.table','rjson')

for(i in 1:length(metr_pkgs)){
  
  library(metr_pkgs[i], character.only = TRUE)
  
}

DeepDetectTrain<-function(ProteinFile,PeptideFile,FastaFile){
  source("../code/deepdetect/R/init.R")
  peptideReport = read_csv(PeptideFile)
  
  peptideQuantity = local({
    runs = unique(peptideReport$R.FileName)
    peptideQuantity = Reduce(
      function(x, y) list(
        run = '',
        data = merge(
          x$data, y$data, 
          by = 'sequence', all = TRUE
        )
      ), 
      lapply(runs, function(run) {
        peptideReport = peptideReport[peptideReport$R.FileName == run, ]
        list(
          run = run,
          data = data.frame(
            sequence = peptideReport$PEP.StrippedSequence,
            intensity = peptideReport$`PEP.Label-Free Quant`,
            stringsAsFactors = FALSE
          )
        )
      })
    )$data
    colnames(peptideQuantity)[-1] = runs
    peptideQuantity
  })
  
  peptides = local({
    sequence = peptideQuantity$sequence
    intensity = apply(peptideQuantity[, -1], 1, function(x) median(x, na.rm = TRUE))
    
    rowIndexes = match(sequence, peptideReport$PEP.StrippedSequence)
    missCleavages = peptideReport$PEP.NrOfMissedCleavages[rowIndexes]
    qValue = peptideReport$PEP.QValue[rowIndexes]
    protein = sapply(strsplit(peptideReport$PG.ProteinAccessions[rowIndexes], ';'), function(x) x[1])
    start = sapply(strsplit(peptideReport$PEP.StartingPositions[rowIndexes], ';'), function(x) as.integer(gsub('^\\(([0-9]+)\\).*', '\\1', x[1])))
    end = sapply(strsplit(peptideReport$PEP.EndingPositions[rowIndexes], ';'), function(x) as.integer(gsub('^\\(([0-9]+)\\).*', '\\1', x[1])))
    
    data.frame(
      protein = protein,
      sequence = sequence,
      start = start,
      end = end,
      missCleavages = missCleavages,
      intensity = intensity,
      stringsAsFactors = FALSE
    )
  })
  
  peptides = local({
    peptides = peptides[order(peptides$protein), ]
    protein.indexes = which(!duplicated(peptides$protein))
    peptides$relativeIntensity = do.call(c, lapply(1:length(protein.indexes), function(i) {
      if (i < length(protein.indexes)) {
        rowIndexes = protein.indexes[i]:(protein.indexes[i + 1] - 1)
      }
      else {
        rowIndexes = protein.indexes[i]:nrow(peptides)
      }
      
      intensity = peptides$intensity[rowIndexes]
      relativeIntensity = intensity / max(intensity)
    }))
    peptides
  })
  
  peptides$detectability = pmax((5 + log10(peptides$relativeIntensity)) / 5 * 0.5, 0) + 0.5
  
  proteinReport = read_csv(ProteinFile)
  
  proteinAccession = local({
    coverageThreshold = 25
    coverage = sapply(strsplit(proteinReport$PG.Coverage, ';'), function(x) as.numeric(sub('%', '', x[1])))
    unique(sapply(strsplit(proteinReport$PG.ProteinAccessions[which(
      proteinReport$PG.UniquePeptides > 1 &
        coverage >= coverageThreshold
    )], ';'), function(x) x[1]))
  })
  
  writeLines(
    proteinAccession, 
    'train/detect/train_excludeSingleHit_coverage25.proteinAccession.txt'
  )
  
  peptides = peptides[peptides$protein %in% proteinAccession, ]
  
  peptides = local({
    proteinHit = table(peptides$protein)
    proteins = names(proteinHit)[proteinHit > 1]
    peptides[peptides$protein %in% proteins, ]
  })
  
  write.csv(peptides,
            'train/detect/train_excludeSingleHit_coverage25.detectability.csv',
            row.names = FALSE
  )
  setwd("train/detect/")
  fastaFile2<-'train_excludeSingleHit_coverage25.proteinAccession.txt'
  system(paste0("pwsh -File 'Filter-Fasta.ps1' ",FastaFile," ",fastaFile2))
  system("mv *.filtered.fasta train_excludeSingleHit_coverage25.fasta")
  proteins = read.fasta("train_excludeSingleHit_coverage25.fasta")
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
  digestedPeptides=df
  
  # detectability.file = list.files(pattern = '\\.detectability\\.csv$')[1]
  # peptides = read_csv(detectability.file)
  peptides.negative = local({
    digestedPeptides = digestedPeptides[!grepl('[^ACDEFGHIKLMNPQRSTVWY]', digestedPeptides$Sequence), ]
    digestedPeptides = digestedPeptides[nchar(digestedPeptides$Sequence) <= 50, ]
    
    digestedPeptides = digestedPeptides[!duplicated(digestedPeptides$Sequence), ]
    
    indexes = which(!(digestedPeptides$Sequence %in% peptides$sequence))
    
    data.frame(
      protein = digestedPeptides$Protein_Name[indexes],
      sequence = digestedPeptides$Sequence[indexes],
      detectability = 0,
      stringsAsFactors = FALSE
    )
  })
  
  peptides.negative$protein = sub('^[A-Za-z0-9]+\\|([A-Z0-9]+)\\|.*', '\\1', peptides.negative$protein)
  
  write.csv(
    peptides.negative, 
    'train_excludeSingleHit_coverage25_negative.detectability.csv',
    row.names = FALSE
  )
  source("../../../code/deepdetect/R/get_cleavage_window.R")
  system("~/anaconda3/envs/tensorflow/bin/python ../../../code/deepdetect/py/train_hard_negative.py")
  if (file.exists("training_0/models/last_epoch.hdf5")) {
    train_score = rjson::fromJSON(file = "./training.json")
    x<-do.call(rbind,lapply(train_score$rounds,function(x){
      unlist(x$evaluate)
    }))%>%as.data.frame
    if (nrow(x)==0) {
      y<-0
    }else{
      x['score']=(x[,1]*(1/3))-(x[,2]*(1/3))+(x[,3]*(1/3))
      y<-which(x[,4]==max(x[,4]))-1
    }
    num<-str_split(dir(paste0("training_",y,"/models/")),"\\.",simplify = T)[,1]%>%str_extract_all("\\d+",simplify = T)%>%max()
    DetectModel<-paste0("epoch_",num,".hdf5")
    file.copy(paste0("training_",y,"/models/",DetectModel),"../../models/detectability/models/")
    file.create("success.txt")
  }
}

DeepMsmsTrain<-function(FragmentFile){
  setwd("~/shiny-server/project/deepDIA")
  setwd(paste0("../data/deepdia/",modelid))
  fragmentReport = read_csv(FragmentFile)
  psm.indexes = which(sapply(1:nrow(fragmentReport), function(i) {
    if (i == 1) {
      return(TRUE)
    }
    
    if (fragmentReport$PSM.MS2ScanNumber[i] != fragmentReport$PSM.MS2ScanNumber[i - 1] ||
        fragmentReport$PEP.StrippedSequence[i] != fragmentReport$PEP.StrippedSequence[i - 1] ||
        fragmentReport$PP.Charge[i] != fragmentReport$PP.Charge[i - 1] ||
        fragmentReport$PG.ProteinAccessions[i] != fragmentReport$PG.ProteinAccessions[i - 1] ||
        fragmentReport$R.FileName[i] != fragmentReport$R.FileName[i - 1]) {
      TRUE
    }
    else {
      FALSE
    }
  }))
  
  lapply(2:3, function(charge) {
    extracted.ions = lapply(1:length(psm.indexes), function(i) {
      if (fragmentReport$PP.Charge[psm.indexes[i]] != charge) {
        return(NULL)
      }
      
      if (i < length(psm.indexes)) {
        rowIndexes = psm.indexes[i]:psm.indexes[i + 1]
      }
      else {
        rowIndexes = psm.indexes[i]:nrow(fragmentReport)
      }
      
      sequence = fragmentReport$PEP.StrippedSequence[rowIndexes[1]]
      charge = fragmentReport$PP.Charge[rowIndexes[1]]
      
      ions = list(
        b1 = rep(0, nchar(sequence) - 1),
        bn1 = rep(0, nchar(sequence) - 1),
        bo1 = rep(0, nchar(sequence) - 1),
        b2 = rep(0, nchar(sequence) - 1),
        bn2 = rep(0, nchar(sequence) - 1),
        bo2 = rep(0, nchar(sequence) - 1),
        y1 = rep(0, nchar(sequence) - 1),
        yn1 = rep(0, nchar(sequence) - 1),
        yo1 = rep(0, nchar(sequence) - 1),
        y2 = rep(0, nchar(sequence) - 1),
        yn2 = rep(0, nchar(sequence) - 1),
        yo2 = rep(0, nchar(sequence) - 1)
      )
      
      lapply(rowIndexes, function(row) {
        type = paste0(
          fragmentReport$FI.FrgType[row],
          ifelse(
            fragmentReport$FI.LossType[row] == 'noloss', '',
            ifelse(
              fragmentReport$FI.LossType[row] == 'NH3', 'n',
              ifelse(
                fragmentReport$FI.LossType[row] == 'H2O', 'o', 
                ''
              )
            )
          ),
          fragmentReport$FI.Charge[row]
        )
        num = as.integer(fragmentReport$FI.FrgNum[row])
        intensity = as.numeric(fragmentReport$FI.Intensity[row])
        if (type %in% names(ions) && num > 0 && num <= nchar(sequence) - 1 && !is.na(intensity)) {
          ions[[type]][num] <<- intensity
        }
      })
      
      qvalue = fragmentReport$PSM.Qvalue[rowIndexes[1]]
      
      list(
        peptide = sequence,
        charge = charge,
        ions = ions,
        qvalue = qvalue
      )
    })
    
    extracted.ions = extracted.ions[!sapply(extracted.ions, is.null)]
    save(extracted.ions, file = paste0('train/msms/train_charge', charge, '.ions.RData'))
  })
  
  setwd("train/msms/")
  source("../../../code/deepms2/R/remove_redundant_ions.R")
  system("mv train_charge2.ions.json charge2")
  system("mv train_charge3.ions.json charge3")
  setwd("charge2")
  system("~/anaconda3/envs/tensorflow/bin/python ../../../../code/deepms2/py/train.py")
  if (file.exists("models/last_epoch.hdf5")) {
    num<-str_split(dir("models/"),"\\.",simplify = T)[,1]%>%str_extract_all("\\d+",simplify = T)%>%max()
    Charge2Model<-paste0("epoch_",num,".hdf5")
    file.copy(paste0("models/",Charge2Model),"../../../models/charge2/models/")
    file.create("success.txt")
  }
  
  setwd("../charge3")
  system("~/anaconda3/envs/tensorflow/bin/python ../../../../code/deepms2/py/train.py")
  if (file.exists("models/last_epoch.hdf5")) {
    num<-str_split(dir("models/"),"\\.",simplify = T)[,1]%>%str_extract_all("\\d+",simplify = T)%>%max()
    Charge3Model<-paste0("epoch_",num,".hdf5")
    file.copy(paste0("models/",Charge3Model),"../../../models/charge3/models/")
    file.create("success.txt")
  }
  
  setwd("../irt")
  psmReport = fragmentReport
  irt = local({
    rank = order(psmReport$PSM.Qvalue, decreasing = FALSE)
    dup = duplicated(psmReport$PEP.StrippedSequence[rank])
    uni = rank[!dup]
    indexes = uni[order(uni)]
    
    irt = data.frame(
      sequence = psmReport$PEP.StrippedSequence[indexes],
      irt = psmReport$PP.iRTEmpirical[indexes],
      rt = psmReport$PP.EmpiricalRT[indexes],
      stringsAsFactors = FALSE
    )
  })
  write.csv(
    irt, 
    'train.irt.csv',
    row.names = FALSE
  )
  system("~/anaconda3/envs/tensorflow/bin/python ../../../../code/deeprt/py/train.py")
  if (file.exists("models/last_epoch.hdf5")) {
    num<-str_split(dir("models/"),"\\.",simplify = T)[,1]%>%str_extract_all("\\d+",simplify = T)%>%max()
    iRTModel<-paste0("epoch_",num,".hdf5")
    file.copy(paste0("models/",iRTModel),"../../../models/irt/models/")
    file.create("success.txt")
  }
}

UpdateStatus <- function(status,id) {
  projects_list <- dbConnect(RSQLite::SQLite(), "../db/projects_list.db", flags = SQLITE_RWC)	
  dbSendStatement(projects_list, paste("UPDATE Project SET status = '",status,"' WHERE id='",id,"'",sep=""))
  dbDisconnect(projects_list)
}

setwd("~/shiny-server/project/deepDIA")
projects_list <- dbConnect(RSQLite::SQLite(), "../db/projects_list.db")		
project_status_table <- as.data.frame(dbGetQuery(projects_list, paste("SELECT * FROM Project WHERE user='","deepDIA","'",sep="")))
dbDisconnect(projects_list)

#在运行项目检查
processing_projects<-project_status_table[which(project_status_table$status==2),]
if (nrow(processing_projects)>0) {
  cat("有正在运行的train任务","\n")
  quit()
}

tryCatch({
  unfinished_train_projects <- project_status_table[which((project_status_table$status==0)&(project_status_table$name=="train")),]
  #unfinished_fastapredict_projects <- project_status_table[which((project_status_table$status==0)&(project_status_table$name=="fasta_predict")),]
  
  modelid<<-unfinished_train_projects[which(unfinished_train_projects$time==min(unfinished_train_projects$time)),"id"]
  #jobid<<-unfinished_fastapredict_projects[which(unfinished_fastapredict_projects$time==min(unfinished_fastapredict_projects$time)),"id"]
},
error=function(e){
  cat(conditionMessage(e),"\n","没有train任务\n")
  quit()
}
)

if (length(modelid)>0) {
  t1=proc.time()
  UpdateStatus(2,modelid)
  setwd(paste0("../data/deepdia/",modelid))
  ProteinFile<-list.files(pattern = "\\.ProteinReport\\.csv$")
  PeptideFile<-list.files(pattern = "\\.PeptideReport\\.csv$")
  FastaFile<-list.files("train/detect/",pattern = "\\.fasta$")
  FragmentFile<-list.files(pattern = "\\.FragmentReport\\.csv$")
  DeepDetectTrain(ProteinFile,PeptideFile,FastaFile)
  DeepMsmsTrain(FragmentFile)
  setwd("~/shiny-server/project/deepDIA")
  if (file.exists(paste0("../data/deepdia/",modelid,"/train/msms/irt/success.txt"))) {
    UpdateStatus(1,modelid)
  }else{
    UpdateStatus(3,modelid)
  }
  t2=proc.time()
  t=t2-t1
  print(paste0('执行时间：',t[3][[1]]/60/60,' hours'))
}else{
  stop
}
