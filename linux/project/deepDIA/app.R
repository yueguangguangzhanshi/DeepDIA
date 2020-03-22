metr_pkgs<-c('shinyjs','shiny','rJava','PKI','digest','shinydashboard','shinycssloaders','RSQLite','jsonlite','dplyr','DT','V8','shinyBS','readr','tidyverse','mixOmics','gplots','pheatmap', 'plotly', 'ggplotify', 'dragulaR', 'factoextra', 'FactoMineR', 'ggsci', 'seqinr', 'ggrepel', 'ggplot2','mailR','shinyAce','cleaver','data.table','rjson')

for(i in 1:length(metr_pkgs)){
  
  library(metr_pkgs[i], character.only = TRUE)
  
}

options(shiny.maxRequestSize=26*1000*1024^2)

ui<-dashboardPage(
  dashboardHeader(title = "DeepDIA Analysis", 
                  dropdownMenuOutput("messageMenu"),
                  tags$li(class = "dropdown", style = "font-weight:bold;", tags$a(textOutput("usertext"), href="../dashboard")),
                  tags$li(class = "dropdown", style = "font-weight:bold;", tags$a(icon("sign-out"), "Logout", onclick ="Shiny.onInputChange('logout_switch', new Date().toLocaleString());", href="#"))
  ),
  dashboardSidebar(
    sidebarMenuOutput("menu"),
    uiOutput("test")
    # tags$div(id = "id-progress-grey-out", class = "progress-grey-out",
    #          ""
    # )
  ),
  dashboardBody(
    shinyjs::useShinyjs(),
    tags$style("#shiny-notification-panel {position: absolute; top: 40% !important;left: 40% !important; zoom: 150%; margin-left: -3%;}"),
    tags$style(".progress-message {zoom: 85%; font-family: Arial,'Helvetica Neue',Helvetica,sans-serif; font-size: 12px !important; font-weight:bold !important;}"),
    tags$style(".shiny-notification {z-index:99; background-color: #e8e8e8; color: #333; border: 0.5px solid #ccc; border-radius: 5px; opacity: 0.95; padding: 10px 8px 10px 10px; margin: 2px; text-align: center;}"),
    tags$style(".progress-grey-out {z-index:97; background-color:rgba(195,195,195,0.85); position: absolute; top: 0; left: 0;}"),
    tags$script(
      '$(document).on("hidden.bs.modal", function (e) { if (e.target.id === "Peptidesmodal") { x = new Date().toLocaleString();
			Shiny.onInputChange("last_modal_close",x); } });
			$(document).ready(function() {$(\'#login_panel\').css(\'margin-top\',$(window).height()*0.2);});'
    ),
    extendShinyjs(text = "shinyjs.resetClick = function() { Shiny.onInputChange('select_button', 'null'); }"),
    column(9,
           fluidRow(
             uiOutput("dashboard_panel")
             )
           ),
    # column(3,
    #        fluidRow(
    #          uiOutput("server_announcement")
    #        ),
    #        fluidRow(
    #          uiOutput("server_status")
    #        )
    # ),
    div(style="margin:2em 0 2em 0;clear:both;","")
  )
)

server <- function(input, output, session) {
  
  myFun <- function(n = 5000) {
    a <- do.call(paste0, replicate(5, sample(LETTERS, n, TRUE), FALSE))
    paste0(round(as.numeric(Sys.time()),0), a, sprintf("%04d", sample(9999, n, TRUE)), sample(LETTERS, n, TRUE))
  }
  
  # DeepFastaPredict<-function(FastaFile,TrainModelId){
  #   setwd(paste0("../data/deepdia/",TrainModelId,"/library/",jobid))
  #   source("../../../code/deepdetect/R/init.R")
  #   proteins = read.fasta(FastaFile)
  #   fasta<-list()
  #   for (i in 1:length(proteins)) {
  #     fasta[[i]]<-sapply(proteins[[i]]$sequence,FUN = function(x)cleave(x,"trypsin",missedCleavages=0))
  #     names(fasta[[i]])<-proteins[[i]]$name
  #   }
  #   
  #   fasta2<-lapply(fasta, function(x){
  #     data<-as.data.frame(x,stringsAsFactors=F,check.names=F)
  #     data[,2]<-colnames(data)[1]
  #     return(data)
  #   })
  #   df <- rbindlist(fasta2,use.names=FALSE)%>%.[,c(2,1)]
  #   names(df)<-c("Protein_Name","Sequence")
  #   df<-filter(df,nchar(Sequence)>=7&nchar(Sequence)<=50)%>%add_column(miss=0)
  #   
  #   digested=df
  #   peptides = data.frame(
  #     protein = digested$Protein_Name,
  #     sequence = digested$Sequence,
  #     miss = digested$miss,
  #     stringsAsFactors = FALSE
  #   )
  #   
  #   proteins = data.frame(
  #     accession = sapply(proteins, function(x) x$name),
  #     sequence = sapply(proteins, function(x) x$sequence),
  #     description = sapply(proteins, function(x) sub('^>[^ ]* ', '', x$description)),
  #     stringsAsFactors = FALSE
  #   )
  #   
  #   cleavages = find.cleavageWindow(peptides, proteins)
  #   
  #   peptides = cbind(peptides, cleavages)
  #   
  #   peptides = peptides[which(
  #     !grepl('[^ACDEFGHIKLMNPQRSTVWY]', peptides$sequence) &
  #       !grepl('[^_ACDEFGHIKLMNPQRSTVWY]', peptides$nTerminal) &
  #       !grepl('[^_ACDEFGHIKLMNPQRSTVWY]', peptides$cTerminal)
  #   ), ]
  #   
  #   write.csv(
  #     peptides,
  #     "../../models/detectability/validated.peptide.csv",
  #     row.names = FALSE
  #   )
  #   setwd("../../models/detectability/")
  #   system("~/anaconda3/envs/tensorflow/bin/python ../../../code/deepdetect/py/predict.py")
  #   source("../../../code/deepdetect/R/filter_peptides_by_detectability.R")
  #   
  #   file_name<-list.files(pattern = '(.)*detectability(.)*\\.peptide\\.csv$')
  #   file<-read.csv(file_name,1)
  #   if (nrow(file)>50000) {
  #     index1<-sample(1:nrow(file),size = 50000)
  #   }else{
  #     index1<-1:nrow(file)
  #   }
  #   peptide_2charge<-file[index1,]
  #   peptide_2charge<-peptide_2charge[unlist(lapply(peptide_2charge[,"sequence"],FUN = function(x){nchar(as.vector(x))[[1]][1]}))<=50,]
  #   peptide_2charge = peptide_2charge[which(
  #     !grepl('[^ACDEFGHIKLMNPQRSTVWY]', peptide_2charge$sequence)
  #   ), ]
  #   
  #   # index2<-sample(1:nrow(file),size = 80000)
  #   # peptide_3charge<-file[index2,]
  #   peptide_3charge<-file[index1,]
  #   peptide_3charge<-peptide_3charge[unlist(lapply(peptide_3charge[,"sequence"],FUN = function(x){nchar(as.vector(x))[[1]][1]}))<=50,]
  #   peptide_3charge = peptide_3charge[which(
  #     !grepl('[^ACDEFGHIKLMNPQRSTVWY]', peptide_3charge$sequence)
  #   ), ]
  #   
  #   # index3<-sample(1:nrow(file),size = 80000)
  #   # peptide_irt<-file[index3,]
  #   peptide_irt<-file[index1,]
  #   peptide_irt<-peptide_irt[unlist(lapply(peptide_irt[,"sequence"],FUN = function(x){nchar(as.vector(x))[[1]][1]}))<=50,]
  #   peptide_irt = peptide_irt[which(
  #     !grepl('[^ACDEFGHIKLMNPQRSTVWY]', peptide_irt$sequence)
  #   ), ]
  #   
  #   write.csv(peptide_irt,paste0("../irt/Report.peptide.csv"),row.names =F)
  #   write.csv(peptide_2charge,paste0("../charge2/Report_charge2.peptide.csv"),row.names =F)
  #   write.csv(peptide_3charge,paste0("../charge3/Report_charge3.peptide.csv"),row.names =F)
  #   setwd("../charge2/")
  #   system("~/anaconda3/envs/tensorflow/bin/python ../../../code/deepms2/py/predict.py")
  #   system(paste0("cp *.ions.json ../../library/",jobid))
  #   setwd("../charge3/")
  #   system("~/anaconda3/envs/tensorflow/bin/python ../../../code/deepms2/py/predict.py")
  #   system(paste0("cp *.ions.json ../../library/",jobid))
  #   setwd("../irt/")
  #   system("~/anaconda3/envs/tensorflow/bin/python ../../../code/deeprt/py/predict.py")
  #   system(paste0("cp *.csv ../../library/",jobid))
  #   
  #   setwd(paste0("../../library/",jobid))
  #   source("../../../code/init.R")
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
  #   #return(data)
  #   write.csv(
  #     do.call(rbind, do.call(c, assays)),
  #     sub('\\.peptide\\.csv$', '.prediction.library.csv', file),
  #     row.names = FALSE
  #   )
  #   file.create("success.txt")
  # }
  # 
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
  # 
  # DeepDetectTrain<-function(ProteinFile,PeptideFile,FastaFile){
  #   setwd(paste0("../data/deepdia/",modelid))
  #   
  #   source("../code/deepdetect/R/init.R")
  #   peptideReport = read_csv(PeptideFile)
  #   peptideReport.HPRP = peptideReport[grepl('HPRP', peptideReport$R.FileName), ]
  #   peptideReport = peptideReport[!grepl('HPRP', peptideReport$R.FileName), ]
  #   
  #   peptideQuantity = local({
  #     runs = unique(peptideReport$R.FileName)
  #     peptideQuantity = Reduce(
  #       function(x, y) list(
  #         run = '',
  #         data = merge(
  #           x$data, y$data, 
  #           by = 'sequence', all = TRUE
  #         )
  #       ), 
  #       lapply(runs, function(run) {
  #         peptideReport = peptideReport[peptideReport$R.FileName == run, ]
  #         list(
  #           run = run,
  #           data = data.frame(
  #             sequence = peptideReport$PEP.StrippedSequence,
  #             intensity = peptideReport$`PEP.Label-Free Quant`,
  #             stringsAsFactors = FALSE
  #           )
  #         )
  #       })
  #     )$data
  #     colnames(peptideQuantity)[-1] = runs
  #     peptideQuantity
  #   })
  #   
  #   peptides = local({
  #     sequence = peptideQuantity$sequence
  #     intensity = apply(peptideQuantity[, -1], 1, function(x) median(x, na.rm = TRUE))
  #     
  #     rowIndexes = match(sequence, peptideReport$PEP.StrippedSequence)
  #     missCleavages = peptideReport$PEP.NrOfMissedCleavages[rowIndexes]
  #     qValue = peptideReport$PEP.QValue[rowIndexes]
  #     protein = sapply(strsplit(peptideReport$PG.ProteinAccessions[rowIndexes], ';'), function(x) x[1])
  #     start = sapply(strsplit(peptideReport$PEP.StartingPositions[rowIndexes], ';'), function(x) as.integer(gsub('^\\(([0-9]+)\\).*', '\\1', x[1])))
  #     end = sapply(strsplit(peptideReport$PEP.EndingPositions[rowIndexes], ';'), function(x) as.integer(gsub('^\\(([0-9]+)\\).*', '\\1', x[1])))
  #     
  #     data.frame(
  #       protein = protein,
  #       sequence = sequence,
  #       start = start,
  #       end = end,
  #       missCleavages = missCleavages,
  #       intensity = intensity,
  #       stringsAsFactors = FALSE
  #     )
  #   })
  #   
  #   peptides = local({
  #     peptides = peptides[order(peptides$protein), ]
  #     protein.indexes = which(!duplicated(peptides$protein))
  #     peptides$relativeIntensity = do.call(c, lapply(1:length(protein.indexes), function(i) {
  #       if (i < length(protein.indexes)) {
  #         rowIndexes = protein.indexes[i]:(protein.indexes[i + 1] - 1)
  #       }
  #       else {
  #         rowIndexes = protein.indexes[i]:nrow(peptides)
  #       }
  #       
  #       intensity = peptides$intensity[rowIndexes]
  #       relativeIntensity = intensity / max(intensity)
  #     }))
  #     peptides
  #   })
  #   peptides$detectability = pmax((5 + log10(peptides$relativeIntensity)) / 5 * 0.5, 0) + 0.5
  #   
  #   if (nrow(peptideReport.HPRP)>0) {
  #     peptides.HPRP = local({
  #       peptideReport = peptideReport.HPRP
  #       
  #       sequence = setdiff(peptideReport.HPRP$PEP.StrippedSequence, peptides$sequence)
  #       
  #       rowIndexes = match(sequence, peptideReport$PEP.StrippedSequence)
  #       missCleavages = peptideReport$PEP.NrOfMissedCleavages[rowIndexes]
  #       qValue = peptideReport$PEP.QValue[rowIndexes]
  #       protein = sapply(strsplit(peptideReport$PG.ProteinAccessions[rowIndexes], ';'), function(x) x[1])
  #       start = sapply(strsplit(peptideReport$PEP.StartingPositions[rowIndexes], ';'), function(x) as.integer(gsub('^\\(([0-9]+)\\).*', '\\1', x[1])))
  #       end = sapply(strsplit(peptideReport$PEP.EndingPositions[rowIndexes], ';'), function(x) as.integer(gsub('^\\(([0-9]+)\\).*', '\\1', x[1])))
  #       
  #       value = min(peptides$relativeIntensity) / 2
  #       detectability = pmax((5 + log10(value)) / 5 * 0.5, 0) + 0.5
  #       
  #       data.frame(
  #         protein = protein,
  #         sequence = sequence,
  #         start = start,
  #         end = end,
  #         missCleavages = missCleavages,
  #         intensity = NA,
  #         relativeIntensity = NA,
  #         detectability = detectability,
  #         stringsAsFactors = FALSE
  #       )
  #     })
  #   }
  #   proteinReport = read_csv(ProteinFile)
  #   
  #   proteinAccession = local({
  #     coverageThreshold = 25
  #     coverage = sapply(strsplit(proteinReport$PG.Coverage, ';'), function(x) as.numeric(sub('%', '', x[1])))
  #     unique(sapply(strsplit(proteinReport$PG.ProteinAccessions[which(
  #       proteinReport$PG.UniquePeptides > 1 &
  #         coverage >= coverageThreshold
  #     )], ';'), function(x) x[1]))
  #   })
  #   
  #   writeLines(
  #     proteinAccession, 
  #     'train/detect/train_excludeSingleHit_coverage25.proteinAccession.txt'
  #   )
  #   
  #   peptides = peptides[peptides$protein %in% proteinAccession, ]
  #   
  #   if (nrow(peptideReport.HPRP)>0) {
  #     peptides.HPRP = peptides.HPRP[peptides.HPRP$protein %in% proteinAccession, ]
  #     peptides = rbind(peptides, peptides.HPRP)
  #   }
  #   
  #   peptides = local({
  #     proteinHit = table(peptides$protein)
  #     proteins = names(proteinHit)[proteinHit > 1]
  #     peptides[peptides$protein %in% proteins, ]
  #   })
  #   
  #   write.csv(peptides,
  #             'train/detect/train_excludeSingleHit_coverage25.detectability.csv',
  #             row.names = FALSE
  #   )
  #   setwd("train/detect/")
  #   fastaFile2<-'train_excludeSingleHit_coverage25.proteinAccession.txt'
  #   system(paste0("pwsh -File 'Filter-Fasta.ps1' ",FastaFile," ",fastaFile2))
  #   system("mv *.filtered.fasta train_excludeSingleHit_coverage25.fasta")
  #   proteins = read.fasta("train_excludeSingleHit_coverage25.fasta")
  #   fasta<-list()
  #   for (i in 1:length(proteins)) {
  #     fasta[[i]]<-sapply(proteins[[i]]$sequence,FUN = function(x)cleave(x,"trypsin",missedCleavages=0))
  #     names(fasta[[i]])<-proteins[[i]]$name
  #   }
  #   
  #   fasta2<-lapply(fasta, function(x){
  #     data<-as.data.frame(x,stringsAsFactors=F,check.names=F)
  #     data[,2]<-colnames(data)[1]
  #     return(data)
  #   })
  #   df <- rbindlist(fasta2,use.names=FALSE)%>%.[,c(2,1)]
  #   names(df)<-c("Protein_Name","Sequence")
  #   df<-filter(df,nchar(Sequence)>=7&nchar(Sequence)<=50)%>%add_column(miss=0)
  #   digestedPeptides=df
  #   
  #   # detectability.file = list.files(pattern = '\\.detectability\\.csv$')[1]
  #   # peptides = read_csv(detectability.file)
  #   peptides.negative = local({
  #     digestedPeptides = digestedPeptides[!grepl('[^ACDEFGHIKLMNPQRSTVWY]', digestedPeptides$Sequence), ]
  #     digestedPeptides = digestedPeptides[nchar(digestedPeptides$Sequence) <= 50, ]
  #     
  #     digestedPeptides = digestedPeptides[!duplicated(digestedPeptides$Sequence), ]
  #     
  #     indexes = which(!(digestedPeptides$Sequence %in% peptides$sequence))
  #     
  #     data.frame(
  #       protein = digestedPeptides$Protein_Name[indexes],
  #       sequence = digestedPeptides$Sequence[indexes],
  #       detectability = 0,
  #       stringsAsFactors = FALSE
  #     )
  #   })
  #   
  #   peptides.negative$protein = sub('^[A-Za-z0-9]+\\|([A-Z0-9]+)\\|.*', '\\1', peptides.negative$protein)
  #   
  #   write.csv(
  #     peptides.negative, 
  #     'train_excludeSingleHit_coverage25_negative.detectability.csv',
  #     row.names = FALSE
  #   )
  #   source("../../../code/deepdetect/R/get_cleavage_window.R")
  #   system("~/anaconda3/envs/tensorflow/bin/python ../../../code/deepdetect/py/train_hard_negative.py")
  #   if (file.exists("training_0/models/last_epoch.hdf5")) {
  #     num<-str_split(dir("training_0/models/"),"\\.",simplify = T)[,1]%>%str_extract_all("\\d+",simplify = T)%>%max()
  #     DetectModel<-paste0("epoch_",num,".hdf5")
  #     file.copy(paste0("training_0/models/",DetectModel),"../../models/detectability/models/")
  #     file.create("success.txt")
  #   }
  #   
  # }
  # 
  # DeepMsmsTrain<-function(FragmentFile){
  #   setwd(paste0("../data/deepdia/",modelid))
  #   fragmentReport = read_csv(FragmentFile)
  #   psm.indexes = which(sapply(1:nrow(fragmentReport), function(i) {
  #     if (i == 1) {
  #       return(TRUE)
  #     }
  #     
  #     if (fragmentReport$PSM.MS2ScanNumber[i] != fragmentReport$PSM.MS2ScanNumber[i - 1] ||
  #         fragmentReport$PEP.StrippedSequence[i] != fragmentReport$PEP.StrippedSequence[i - 1] ||
  #         fragmentReport$PP.Charge[i] != fragmentReport$PP.Charge[i - 1] ||
  #         fragmentReport$PG.ProteinAccessions[i] != fragmentReport$PG.ProteinAccessions[i - 1] ||
  #         fragmentReport$R.FileName[i] != fragmentReport$R.FileName[i - 1]) {
  #       TRUE
  #     }
  #     else {
  #       FALSE
  #     }
  #   }))
  #   
  #   lapply(2:3, function(charge) {
  #     extracted.ions = lapply(1:length(psm.indexes), function(i) {
  #       if (fragmentReport$PP.Charge[psm.indexes[i]] != charge) {
  #         return(NULL)
  #       }
  #       
  #       if (i < length(psm.indexes)) {
  #         rowIndexes = psm.indexes[i]:psm.indexes[i + 1]
  #       }
  #       else {
  #         rowIndexes = psm.indexes[i]:nrow(fragmentReport)
  #       }
  #       
  #       sequence = fragmentReport$PEP.StrippedSequence[rowIndexes[1]]
  #       charge = fragmentReport$PP.Charge[rowIndexes[1]]
  #       
  #       ions = list(
  #         b1 = rep(0, nchar(sequence) - 1),
  #         bn1 = rep(0, nchar(sequence) - 1),
  #         bo1 = rep(0, nchar(sequence) - 1),
  #         b2 = rep(0, nchar(sequence) - 1),
  #         bn2 = rep(0, nchar(sequence) - 1),
  #         bo2 = rep(0, nchar(sequence) - 1),
  #         y1 = rep(0, nchar(sequence) - 1),
  #         yn1 = rep(0, nchar(sequence) - 1),
  #         yo1 = rep(0, nchar(sequence) - 1),
  #         y2 = rep(0, nchar(sequence) - 1),
  #         yn2 = rep(0, nchar(sequence) - 1),
  #         yo2 = rep(0, nchar(sequence) - 1)
  #       )
  #       
  #       lapply(rowIndexes, function(row) {
  #         type = paste0(
  #           fragmentReport$FI.FrgType[row],
  #           ifelse(
  #             fragmentReport$FI.LossType[row] == 'noloss', '',
  #             ifelse(
  #               fragmentReport$FI.LossType[row] == 'NH3', 'n',
  #               ifelse(
  #                 fragmentReport$FI.LossType[row] == 'H2O', 'o', 
  #                 ''
  #               )
  #             )
  #           ),
  #           fragmentReport$FI.Charge[row]
  #         )
  #         num = as.integer(fragmentReport$FI.FrgNum[row])
  #         intensity = as.numeric(fragmentReport$FI.Intensity[row])
  #         if (type %in% names(ions) && num > 0 && num <= nchar(sequence) - 1 && !is.na(intensity)) {
  #           ions[[type]][num] <<- intensity
  #         }
  #       })
  #       
  #       qvalue = fragmentReport$PSM.Qvalue[rowIndexes[1]]
  #       
  #       list(
  #         peptide = sequence,
  #         charge = charge,
  #         ions = ions,
  #         qvalue = qvalue
  #       )
  #     })
  #     
  #     extracted.ions = extracted.ions[!sapply(extracted.ions, is.null)]
  #     save(extracted.ions, file = paste0('train/msms/train_charge', charge, '.ions.RData'))
  #   })
  #   
  #   setwd("train/msms/")
  #   source("../../../code/deepms2/R/remove_redundant_ions.R")
  #   system("mv train_charge2.ions.json charge2")
  #   system("mv train_charge3.ions.json charge3")
  #   setwd("charge2")
  #   system("~/anaconda3/envs/tensorflow/bin/python ../../../../code/deepms2/py/train.py")
  #   if (file.exists("models/last_epoch.hdf5")) {
  #     num<-str_split(dir("models/"),"\\.",simplify = T)[,1]%>%str_extract_all("\\d+",simplify = T)%>%max()
  #     Charge2Model<-paste0("epoch_",num,".hdf5")
  #     file.copy(paste0("models/",Charge2Model),"../../../models/charge2/models/")
  #     file.create("success.txt")
  #   }
  #   
  #   setwd("../charge3")
  #   system("~/anaconda3/envs/tensorflow/bin/python ../../../../code/deepms2/py/train.py")
  #   if (file.exists("models/last_epoch.hdf5")) {
  #     num<-str_split(dir("models/"),"\\.",simplify = T)[,1]%>%str_extract_all("\\d+",simplify = T)%>%max()
  #     Charge3Model<-paste0("epoch_",num,".hdf5")
  #     file.copy(paste0("models/",Charge3Model),"../../../models/charge3/models/")
  #     file.create("success.txt")
  #   }
  #   
  #   setwd("../irt")
  #   psmReport = fragmentReport
  #   irt = local({
  #     rank = order(psmReport$PSM.Qvalue, decreasing = FALSE)
  #     dup = duplicated(psmReport$PEP.StrippedSequence[rank])
  #     uni = rank[!dup]
  #     indexes = uni[order(uni)]
  #     
  #     irt = data.frame(
  #       sequence = psmReport$PEP.StrippedSequence[indexes],
  #       irt = psmReport$PP.iRTEmpirical[indexes],
  #       rt = psmReport$PP.EmpiricalRT[indexes],
  #       stringsAsFactors = FALSE
  #     )
  #   })
  #   write.csv(
  #     irt, 
  #     'train.irt.csv',
  #     row.names = FALSE
  #   )
  #   system("~/anaconda3/envs/tensorflow/bin/python ../../../../code/deeprt/py/train.py")
  #   if (file.exists("models/last_epoch.hdf5")) {
  #     num<-str_split(dir("models/"),"\\.",simplify = T)[,1]%>%str_extract_all("\\d+",simplify = T)%>%max()
  #     iRTModel<-paste0("epoch_",num,".hdf5")
  #     file.copy(paste0("models/",iRTModel),"../../../models/irt/models/")
  #     file.create("success.txt")
  #   }
  # }
  
  # DeepGetModel=function(){
  #   DetectModel<-""
  #   Charge2Model<-""
  #   Charge3Model<-""
  #   iRTModel<-""
  #   
  #   if (length(dir("../data/deepdia/data/train/detect/training_0/models/"))>0) {
  #     num<-str_split(dir("../data/deepdia/data/train/detect/training_0/models/"),"\\.",simplify = T)[,1]%>%str_extract_all("\\d+",simplify = T)%>%max()
  #     DetectModel<-paste0("epoch_",num,".hdf5")
  #     #file.copy(paste0("../data/deepdia/data/train/detect/training_0/models/",DetectModel),"../data/deepdia/data/train/model/detectability/")
  #   }
  #   if (length(dir("../data/deepdia/data/train/msms/charge2/models/"))>0) {
  #     num<-str_split(dir("../data/deepdia/data/train/msms/charge2/models/"),"\\.",simplify = T)[,1]%>%str_extract_all("\\d+",simplify = T)%>%max()
  #     Charge2Model<-paste0("epoch_",num,".hdf5")
  #   }
  #   if (length(dir("../data/deepdia/data/train/msms/charge3/models/"))>0) {
  #     num<-str_split(dir("../data/deepdia/data/train/msms/charge3/models/"),"\\.",simplify = T)[,1]%>%str_extract_all("\\d+",simplify = T)%>%max()
  #     Charge3Model<-paste0("epoch_",num,".hdf5")
  #   }
  #   if (length(dir("../data/deepdia/data/train/msms/irt/models/"))>0) {
  #     num<-str_split(dir("../data/deepdia/data/train/msms/irt/models/"),"\\.",simplify = T)[,1]%>%str_extract_all("\\d+",simplify = T)%>%max()
  #     iRTModel<-paste0("epoch_",num,".hdf5")
  #   }
  #   return(c(DetectModel,Charge2Model,Charge3Model,iRTModel))
  # }

  # SendLibrary<-function(data){
  #   data<-read.csv(da ta)
  #   write.csv(data,"test.csv")
  #   return(data)
  # }
  
  output$dashboard_panel <-renderUI({
    mainPanel(
      tabsetPanel(type = "tabs",
                  tabPanel("Model Train",
                           fileInput("TrainProtein", "Choose Protein csv File",
                                     multiple = FALSE,
                                     accept = c(".csv")),
                           fileInput("TrainPeptide", "Choose Peptide csv File",
                                     multiple = FALSE,
                                     accept = c(".csv")),
                           selectInput("TrainFasta", "Choose fasta:",
                                       c("human" = "human.fasta",
                                         "mouse" = "mouse.fasta")),
                           # fileInput("TrainFasta", "Choose Fasta File",
                           #           multiple = FALSE,
                           #           accept = c(".fasta")),
                           fileInput("TrainFragment", "Choose Fragment csv File",
                                     multiple = FALSE,
                                     accept = c(".csv")),
                           actionButton("TrainConfirm", "确认")
                           # downloadButton("DownloadModel", "Download Models")
                  ),
                  tabPanel("Fasta Predict",
                           fileInput("PredictFasta", "Choose Fasta File",
                                     multiple = FALSE,
                                     accept = c(".fasta")),
                           textInput("TrainId1", "输入Train Models Task ID"),
                           actionButton("PredictFastaConfirm", "确认")
                           # downloadButton("DownloadFastaLibrary", "Download Library")
                  ),
                  # tabPanel("Report Predict",
                  #          fileInput("PredictFasta2", "Choose Report fasta File",
                  #                    multiple = FALSE,
                  #                    accept = c(".fasta")),
                  #          fileInput("PredictPeptide", "Choose Peptide csv File",
                  #                    multiple = FALSE,
                  #                    accept = c(".csv")),
                  #          actionButton("PredictReportConfirm", "确认")
                  #          # downloadButton("DownloadReportLibrary", "Download Library")
                  # ),
                  tabPanel("Library Download",
                           textInput("TrainId2", "输入Train Models Task ID"),
                           textInput("LibraryId", "输入Library Task ID"),
                           actionButton("LibraryConfirm", "确认"),
                           fluidRow(
                             uiOutput("LibraryDownload")
                           )
                           # downloadButton("DownloadReportLibrary", "Download Library")
                  )
      )
    )
  })
  
  observeEvent(input$PredictFastaConfirm,{
    setwd("~/shiny-server/project/deepDIA/")
    if (input$TrainId1 %in% dir("../data/deepdia/")) {
      jobid<<-myFun(1)
      showModal(modalDialog(
        title = "Important message",
        paste0("This is your Predict LibraryTask ID: ",jobid),
        easyClose = FALSE,
        footer = NULL
      ))
      dir.create(paste0("../data/deepdia/",input$TrainId1,"/library/",jobid))
      projects_list <- dbConnect(RSQLite::SQLite(), "../db/projects_list.db", flags = SQLITE_RWC)
      dbSendQuery(projects_list, paste0("INSERT INTO Project VALUES ('", jobid, "','", paste0("fasta_predict_",trimws(input$TrainId1)), "', '", "deepDIA", "', '0', '", as.character(Sys.time()),"')"))
      dbDisconnect(projects_list)
      file.copy(input$PredictFasta$datapath,paste0("../data/deepdia/",input$TrainId1,"/library/",jobid))
      #DeepFastaPredict(input$PredictFasta$datapath,input$TrainId1)
      #DeepFastaLibrary<-DeepFastaPredict(input$PredictFasta$datapath)
    }else{
      showModal(modalDialog(
        title = "Important message",
        "你输入的ID未找到models",
        easyClose = TRUE
      ))
    }
  })
  
  observeEvent(input$PredictReportConfirm,{
    jobid<<-myFun(1)
    showModal(modalDialog(
      title = "Important message",
      paste0("This is your Library Task ID: ",jobid),
      easyClose = FALSE
    ))
    DeepReportPredict(input$PredictPeptide$datapath)
    #DeepReportLibrary<-DeepReportPredict(input$PredictPeptide$datapath)
  })
  
  observeEvent(input$TrainConfirm,{
    setwd("~/shiny-server/project/deepDIA/")
    modelid<<-myFun(1)
    showModal(modalDialog(
      title = "Important message",
      paste0("This is your Train Models Task ID: ",modelid),
      easyClose = FALSE
    ))
    dir.create(paste0("../data/deepdia/",modelid))
    projects_list <- dbConnect(RSQLite::SQLite(), "../db/projects_list.db", flags = SQLITE_RWC)
    dbSendQuery(projects_list, paste0("INSERT INTO Project VALUES ('", modelid, "','", "train", "', '", "deepDIA", "', '0', '", as.character(Sys.time()),"')"))
    dbDisconnect(projects_list)
    file.copy(input$TrainProtein$datapath,paste0("../data/deepdia/",modelid))
    system(paste0("mv ../data/deepdia/",modelid,"/0.csv ../data/deepdia/",modelid,"/train.ProteinReport.csv"))
    file.copy(input$TrainPeptide$datapath,paste0("../data/deepdia/",modelid))
    system(paste0("mv ../data/deepdia/",modelid,"/0.csv ../data/deepdia/",modelid,"/train.PeptideReport.csv"))
    file.copy(input$TrainFragment$datapath,paste0("../data/deepdia/",modelid))
    system(paste0("mv ../data/deepdia/",modelid,"/0.csv ../data/deepdia/",modelid,"/train.FragmentReport.csv"))
    system(paste0("cp -rf ../data/deepdia/data/* ","../data/deepdia/",modelid))
    system(paste0("cp -rf ../data/deepdia/database/",input$TrainFasta," ../data/deepdia/",modelid,"/train/detect"))
    #DeepDetectTrain(input$TrainProtein$datapath,input$TrainPeptide$datapath,input$TrainFasta)
    #DeepMsmsTrain(input$TrainFragment$datapath)
    # DeepModel=DeepGetModel()
  })
  
  observeEvent(input$LibraryConfirm,{
    setwd("~/shiny-server/project/deepDIA/")
    if (file.exists(paste0("../data/deepdia/",input$TrainId2,"/library/",input$LibraryId,"/success.txt"))) {
      # if (input$LibraryId=="") {
      #   showModal(modalDialog(
      #     title = "Important message",
      #     "请输入Library Task ID",
      #     easyClose = TRUE
      #   ))
      # }else
        output$LibraryDownload<-renderUI({
          downloadButton("DownloadLibrary", "Download Library")
        })
        output$DownloadLibrary <- downloadHandler(
          filename = function() {
            paste("library_fasta", Sys.Date(), ".csv", sep="")
          },
          content = function(file) {
            #file=read_csv(paste0("../data/deepdia/",input$TrainId2,"/library/",input$LibraryId,"/Report.prediction.library.csv"))
            file.copy(paste0("../data/deepdia/",input$TrainId2,"/library/",input$LibraryId,"/Report.prediction.library.csv"), file)
          }
        )
    }else{
      if ((input$TrainId2 %in% dir("../data/deepdia/"))&(input$LibraryId %in% dir(paste0("../data/deepdia/",input$TrainId2,"/library")))) {
        showModal(modalDialog(
          title = "Important message",
          "你的Library还没好，请等待",
          easyClose = TRUE
        ))
      }else{
        showModal(modalDialog(
          title = "Important message",
          "请输入正确ID",
          easyClose = TRUE
        ))
      }
    }
  })

  output$test <- renderUI({
    
    paste0(input$TrainId2, ',', input$LibraryId, ',', getwd()) %>% print
    
  })
  
  # output$DownloadModel <- downloadHandler(
  #   filename = function() {
  #     paste("Model", Sys.Date(), ".zip", sep="")
  #   },
  #   content = function(file) {
  #     tmpdir <- tempdir()
  #     setwd(tmpdir)
  #     dir.create("models")
  #     setwd(paste0(tmpdir,"/models"))
  #     dir.create("detectability")
  #     dir.create("charge2")
  #     dir.create("charge3")
  #     dir.create("iRT")
  #     file.copy(paste0("/home/saitoasuka/shiny-server/project/data/deepdia/data/train/detect/training_0/models/",DeepModel[1]),"./detectability")
  #     file.copy(paste0("/home/saitoasuka/shiny-server/project/data/deepdia/data/train/msms/charge2/models/",DeepModel[2]),"./charge2")
  #     file.copy(paste0("/home/saitoasuka/shiny-server/project/data/deepdia/data/train/msms/charge3/models/",DeepModel[3]),"./charge3")
  #     file.copy(paste0("/home/saitoasuka/shiny-server/project/data/deepdia/data/train/msms/irt/models/",DeepModel[4]),"./iRT")
  #     zip(zipfile = file, files = paste0(tmpdir,"/models"))
  #     #zip(zipfile = file, files ="./")
  #   }
  # )
}

shinyApp(ui, server)