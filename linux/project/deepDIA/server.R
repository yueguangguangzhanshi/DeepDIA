metr_pkgs<-c('shinyjs','shiny','rJava','PKI','digest','shinydashboard','shinycssloaders','RSQLite','jsonlite','dplyr','DT','V8','shinyBS','readr','tidyverse','mixOmics','gplots','pheatmap', 'plotly', 'ggplotify', 'dragulaR', 'factoextra', 'FactoMineR', 'ggsci', 'seqinr', 'ggrepel', 'ggplot2','mailR','shinyAce','cleaver','data.table','rjson')

for(i in 1:length(metr_pkgs)){
  
  library(metr_pkgs[i], character.only = TRUE)
  
}

options(shiny.maxRequestSize=70*1024^2)

shinyServer(function(input, output, session){
  
  DeepFastaPredict<-function(FastaFile){
    source("../data/deepdia/code/deepdetect/R/init.R")
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
    df<-filter(df,nchar(Sequence)>7&nchar(Sequence)<=50)%>%add_column(miss=0)
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
      "../data/deepdia/data/models/detectability/validated.peptide.csv",
      row.names = FALSE
    )
    setwd("../data/deepdia/data/models/detectability/")
    system("~/anaconda3/envs/tensorflow/bin/python ../../../code/deepdetect/py/predict.py")
    source("../../../code/deepdetect/R/filter_peptides_by_detectability.R")
    
    file_name<-list.files(pattern = '(.)*detectability(.)*\\.peptide\\.csv$')
    file<-read.csv(file_name,1)
    
    index1<-sample(1:nrow(file),size = 80000)
    peptide_2charge<-file[index1,]
    peptide_2charge<-peptide_2charge[unlist(lapply(peptide_2charge[,"sequence"],FUN = function(x){nchar(as.vector(x))[[1]][1]}))<=50,]
    peptide_2charge = peptide_2charge[which(
      !grepl('[^ACDEFGHIKLMNPQRSTVWY]', peptide_2charge$sequence)
    ), ]
    
    index2<-sample(1:nrow(file),size = 80000)
    peptide_3charge<-file[index2,]
    peptide_3charge<-peptide_3charge[unlist(lapply(peptide_3charge[,"sequence"],FUN = function(x){nchar(as.vector(x))[[1]][1]}))<=50,]
    peptide_3charge = peptide_3charge[which(
      !grepl('[^ACDEFGHIKLMNPQRSTVWY]', peptide_3charge$sequence)
    ), ]
    
    index3<-sample(1:nrow(file),size = 80000)
    peptide_irt<-file[index3,]
    peptide_irt<-peptide_irt[unlist(lapply(peptide_irt[,"sequence"],FUN = function(x){nchar(as.vector(x))[[1]][1]}))<=50,]
    peptide_irt = peptide_irt[which(
      !grepl('[^ACDEFGHIKLMNPQRSTVWY]', peptide_irt$sequence)
    ), ]
    
    write.csv(peptide_irt,paste0("../irt/Report.peptide.csv"),row.names =F)
    write.csv(peptide_2charge,paste0("../charge2/Report_charge2.peptide.csv"),row.names =F)
    write.csv(peptide_3charge,paste0("../charge3/Report_charge3.peptide.csv"),row.names =F)
    setwd("../charge2/")
    system("~/anaconda3/envs/tensorflow/bin/python ../../../code/deepms2/py/predict.py")
    system("cp *.ions.json ../../library/")
    setwd("../charge3/")
    system("~/anaconda3/envs/tensorflow/bin/python ../../../code/deepms2/py/predict.py")
    system("cp *.ions.json ../../library/")
    setwd("../irt/")
    system("~/anaconda3/envs/tensorflow/bin/python ../../../code/deeprt/py/predict.py")
    system("cp *.csv ../../library/")
    
    setwd("../../library/")
    source("../../code/init.R")
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
    data<-do.call(rbind, do.call(c, assays))
    # write.csv(
    #   do.call(rbind, do.call(c, assays)),
    #   sub('\\.peptide\\.csv$', '.prediction.library.csv', file), 
    #   row.names = FALSE
    # )
    return(data)
  }
  
  DeepReportPredict<-function(ReportFile){
    file<-read.csv(ReportFile,1)
    peptide<-file[,c("PG.ProteinAccessions","PEP.StrippedSequence","PP.Charge")]
    index1<-sample(1:nrow(peptide),size = 80000)
    peptide_2charge<-peptide[index1,]
    peptide_2charge<-peptide_2charge[peptide_2charge[,"PP.Charge"]==2,c("PG.ProteinAccessions","PEP.StrippedSequence")]
    names(peptide_2charge)<-c("protein","sequence")
    peptide_2charge<-peptide_2charge[unlist(lapply(peptide_2charge[,"sequence"],FUN = function(x){nchar(as.vector(x))[[1]][1]}))<=50,]
    peptide_2charge = peptide_2charge[which(
      !grepl('[^ACDEFGHIKLMNPQRSTVWY]', peptide_2charge$sequence)
    ), ]
    
    index2<-sample(1:nrow(peptide),size = 80000)
    peptide_3charge<-peptide[index2,]
    peptide_3charge<-peptide_3charge[peptide_3charge[,"PP.Charge"]==3,c("PG.ProteinAccessions","PEP.StrippedSequence")]
    names(peptide_3charge)<-c("protein","sequence")
    peptide_3charge<-peptide_3charge[unlist(lapply(peptide_3charge[,"sequence"],FUN = function(x){nchar(as.vector(x))[[1]][1]}))<=50,]
    peptide_3charge = peptide_3charge[which(
      !grepl('[^ACDEFGHIKLMNPQRSTVWY]', peptide_3charge$sequence)
    ), ]
    
    index3<-sample(1:nrow(peptide),size = 80000)
    peptide_irt<-peptide[index3,]
    peptide_irt<-peptide_irt[,c("PG.ProteinAccessions","PEP.StrippedSequence")]
    names(peptide_irt)<-c("protein","sequence")
    peptide_irt<-peptide_irt[unlist(lapply(peptide_irt[,"sequence"],FUN = function(x){nchar(as.vector(x))[[1]][1]}))<=50,]
    peptide_irt = peptide_irt[which(
      !grepl('[^ACDEFGHIKLMNPQRSTVWY]', peptide_irt$sequence)
    ), ]
    
    write.csv(peptide_irt,"../data/deepdia/data/models/irt/Report.peptide.csv",row.names =F)
    write.csv(peptide_2charge,"../data/deepdia/data/models/charge2/Report_charge2.peptide.csv",row.names =F)
    write.csv(peptide_3charge,"../data/deepdia/data/models/charge3/Report_charge3.peptide.csv",row.names =F)
    
    setwd("../data/deepdia/data/models/charge2/")
    system("~/anaconda3/envs/tensorflow/bin/python ../../../code/deepms2/py/predict.py")
    system("cp *.ions.json ../../library/")
    setwd("../charge3/")
    system("~/anaconda3/envs/tensorflow/bin/python ../../../code/deepms2/py/predict.py")
    system("cp *.ions.json ../../library/")
    setwd("../irt/")
    system("~/anaconda3/envs/tensorflow/bin/python ../../../code/deeprt/py/predict.py")
    system("cp *.csv ../../library/")
    
    setwd("../../library/")
    source("../../code/init.R")
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
    data<-do.call(rbind, do.call(c, assays))
    # write.csv(
    #   do.call(rbind, do.call(c, assays)),
    #   sub('\\.peptide\\.csv$', '.prediction.library.csv', file), 
    #   row.names = FALSE
    # )
    return(data)
  }
  
  DeepDetectTrain<-function(ProteinFile,PeptideFile,FastaFile){
    source("../data/deepdia/code/deepdetect/R/init.R")
    peptideReport = read_csv(PeptideFile)
    peptideReport.HPRP = peptideReport[grepl('HPRP', peptideReport$R.FileName), ]
    peptideReport = peptideReport[!grepl('HPRP', peptideReport$R.FileName), ]
    
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
    
    if (nrow(peptideReport.HPRP)>0) {
      peptides.HPRP = local({
        peptideReport = peptideReport.HPRP
        
        sequence = setdiff(peptideReport.HPRP$PEP.StrippedSequence, peptides$sequence)
        
        rowIndexes = match(sequence, peptideReport$PEP.StrippedSequence)
        missCleavages = peptideReport$PEP.NrOfMissedCleavages[rowIndexes]
        qValue = peptideReport$PEP.QValue[rowIndexes]
        protein = sapply(strsplit(peptideReport$PG.ProteinAccessions[rowIndexes], ';'), function(x) x[1])
        start = sapply(strsplit(peptideReport$PEP.StartingPositions[rowIndexes], ';'), function(x) as.integer(gsub('^\\(([0-9]+)\\).*', '\\1', x[1])))
        end = sapply(strsplit(peptideReport$PEP.EndingPositions[rowIndexes], ';'), function(x) as.integer(gsub('^\\(([0-9]+)\\).*', '\\1', x[1])))
        
        value = min(peptides$relativeIntensity) / 2
        detectability = pmax((5 + log10(value)) / 5 * 0.5, 0) + 0.5
        
        data.frame(
          protein = protein,
          sequence = sequence,
          start = start,
          end = end,
          missCleavages = missCleavages,
          intensity = NA,
          relativeIntensity = NA,
          detectability = detectability,
          stringsAsFactors = FALSE
        )
      })
    }
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
      '../data/deepdia/data/train/detect/train_excludeSingleHit_coverage25.proteinAccession.txt'
    )
    
    peptides = peptides[peptides$protein %in% proteinAccession, ]
    
    if (nrow(peptideReport.HPRP)>0) {
      peptides.HPRP = peptides.HPRP[peptides.HPRP$protein %in% proteinAccession, ]
      peptides = rbind(peptides, peptides.HPRP)
    }
    
    peptides = local({
      proteinHit = table(peptides$protein)
      proteins = names(proteinHit)[proteinHit > 1]
      peptides[peptides$protein %in% proteins, ]
    })
    
    write.csv(peptides,
              '../data/deepdia/data/train/detect/train_excludeSingleHit_coverage25.detectability.csv',
              row.names = FALSE
    )
    setwd("../data/deepdia/data/train/detect/")
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
    df<-filter(df,nchar(Sequence)>7&nchar(Sequence)<=50)%>%add_column(miss=0)
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
  }
  
  DeepMsmsTrain<-function(FragmentFile){
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
      save(extracted.ions, file = paste0('../data/deepdia/data/train/msms/train_charge', charge, '.ions.RData'))
    })
    
    setwd("../data/deepdia/data/train/msms/")
    source("../../../code/deepms2/R/remove_redundant_ions.R")
    system("mv train_charge2.ions.json charge2")
    system("mv train_charge3.ions.json charge3")
    setwd("charge2")
    system("~/anaconda3/envs/tensorflow/bin/python ../../../../code/deepms2/py/train.py")
    setwd("../charge3")
    system("~/anaconda3/envs/tensorflow/bin/python ../../../../code/deepms2/py/train.py")
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
  }
  
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
                  tabPanel("Fasta Predict",
                           fileInput("PredictFasta", "Choose Fasta File",
                                     multiple = FALSE,
                                     accept = c(".fasta")),
                           actionButton("PredictFastaConfirm", "确认"),
                           downloadButton("DownloadFastaLibrary", "Download Library")
                  ),
                  tabPanel("Report Predict",
                           # fileInput("PredictProtein", "Choose Protein csv File",
                           #           multiple = FALSE,
                           #           accept = c(".csv")),
                           fileInput("PredictPeptide", "Choose Peptide csv File",
                                     multiple = FALSE,
                                     accept = c(".csv")),
                           actionButton("PredictReportConfirm", "确认"),
                           downloadButton("DownloadReportLibrary", "Download Library")
                  ),
                  tabPanel("Model Train",
                           fileInput("TrainProtein", "Choose Protein csv File",
                                     multiple = FALSE,
                                     accept = c(".csv")),
                           fileInput("TrainPeptide", "Choose Peptide csv File",
                                     multiple = FALSE,
                                     accept = c(".csv")),
                           fileInput("TrainFasta", "Choose Fasta File",
                                     multiple = FALSE,
                                     accept = c(".fasta")),
                           fileInput("TrainFragment", "Choose Fragment csv File",
                                     multiple = FALSE,
                                     accept = c(".csv")),
                           actionButton("TrainConfirm", "确认")
                           # downloadButton("DownloadModel", "Download Models")
                  )
      )
    )
  })
  
  observeEvent(input$PredictFastaConfirm,{
    DeepFastaLibrary<-DeepFastaPredict(input$PredictFasta$datapath)
    # DeepFastaLibrary<-SendLibrary(input$PredictFasta$datapath)
  })
  
  observeEvent(input$PredictReportConfirm,{
    DeepReportLibrary<-DeepReportPredict(input$PredictPeptide$datapath)
  })
  
  observeEvent(input$TrainConfirm,{
    DeepDetectTrain(input$TrainProtein$datapath,input$TrainPeptide$datapath,input$TrainFasta$datapath)
    DeepMsmsTrain(input$TrainFragment$datapath)
    # DeepModel=DeepGetModel()
  })
  
  output$DownloadFastaLibrary <- downloadHandler(
    filename = function() {
      paste("library_fasta", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      file=DeepFastaLibrary
      #file.copy('./test.csv', file)
	  
    }
  )
  
  output$DownloadReportLibrary <- downloadHandler(
    filename = function() {
      paste("library_report", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      file=DeepReportLibrary
    }
  )
  
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
})


