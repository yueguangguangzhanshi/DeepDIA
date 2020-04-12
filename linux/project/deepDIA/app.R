metr_pkgs<-c('shinythemes','shinyjs','shiny','rJava','PKI','digest','shinydashboard','shinycssloaders','RSQLite','jsonlite','dplyr','DT','V8','shinyBS','readr','tidyverse','mixOmics','gplots','pheatmap', 'plotly', 'ggplotify', 'dragulaR', 'factoextra', 'FactoMineR', 'ggsci', 'seqinr', 'ggrepel', 'ggplot2','mailR','shinyAce','cleaver','data.table','rjson')

for(i in 1:length(metr_pkgs)){
  
  library(metr_pkgs[i], character.only = TRUE)
  
}

options(shiny.maxRequestSize=30*1000*1024^2)

ui<-fluidPage(
  theme = shinythemes::shinytheme("simplex"),
  titlePanel(title = "DeepDIA Analysis"),
  sidebarLayout(
    sidebarPanel(
      width =0,
      sidebarMenuOutput("menu"),
      uiOutput("test")
      # tags$div(id = "id-progress-grey-out", class = "progress-grey-out",
      #          ""
      # )
    ),
    mainPanel(
      tags$head(
        tags$link(rel = "stylesheet", type = "text/css", href = "custom.css")
      ),
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
      navlistPanel(
        tabPanel("网站说明",
                 fluidRow(style = 'margin:2em 2em 2em 0;',uiOutput("explain_en")),
                 fluidRow(style = 'margin:2em 2em 2em 0;',uiOutput("explain_ch"))
                 #img(src="../data/DIA-ApexRT_8w.png",width="1280",height="1280")
                 ),
        tabPanel("建模预测分析",
                 fluidRow(style = 'margin-bottom: 80px;',
                          uiOutput("dashboard_panel")
                          )
                 ),
        widths = c(3, 9)
      ),
      div(style="margin:1em 0 1em 0;clear:both;","")
    )
  )
)

server <- function(input, output, session) {
  
  myFun <- function(n = 5000) {
    a <- do.call(paste0, replicate(5, sample(LETTERS, n, TRUE), FALSE))
    paste0(round(as.numeric(Sys.time()),0), a, sprintf("%04d", sample(9999, n, TRUE)), sample(LETTERS, n, TRUE))
  }
  
  output$explain_en<-renderText({
    "DeepDIA offers high quality MS2 train/predicted spectra for any organism fasta as well as iRT train/prediction. DeepDIA is part of the Omicsolution (www.omicsolution.com) shiny project and its models was trained on the project's/lab's High resolution mass spectrometry raw data or DDA Spectronaut search results(protein/peptide/fragment report). When using deepDIA is helpful for your research, please cite : Yang, Y., Liu, X., Shen, C. et al. In silico spectral libraries by deep learning facilitate data-independent acquisition proteomics. Nat Commun 11, 146 (2020). https://doi.org/10.1038/s41467-019-13866-z"
  })
  output$explain_ch<-renderText({
    "复旦大学化学系乔亮课题组与生物医学研究院杨芃原课题组、上海易算生物科技有限公司沈诚频团队、澳大利亚国立大学林宇课题组合作，在《自然-通讯》上在线发表了题为“In silico spectral libraries by deep learning facilitate data-independent acquisition proteomics”的研究论文。该研究利用深度学习技术开发了从肽段或蛋白质序列构建谱图库的工具deepDIA，实现了不依赖于DDA的DIA数据直接分析。使用deepDIA构建专用于实验室或项目的特定仪器的模型，并且预测生成高质量的二级谱图结果和保留时间结果，效果接近DDA构建的谱图库。除此以外，我们还设计了预测肽段在质谱中的可检测性的模型，预测蛋白的理论酶切肽段的可检测性，筛选可检测性分数达到一定阈值的肽段来构建谱图库，最后可以直接将预测的数据库导入Spectronaut软件进行DIA分析
  .\n如果deepDIA对您的研究有帮助，请引用此内容：Yang, Y., Liu, X., Shen, C. et al. In silico spectral libraries by deep learning facilitate data-independent acquisition proteomics. Nat Commun 11, 146 (2020).  https://doi.org/10.1038/s41467-019-13866-z"
  })
  output$dashboard_panel <-renderUI({
    mainPanel(
      tabsetPanel(type = "tabs",id="tabs",
                  tabPanel("Model Train",
                           br(),
                           helpText("注意:请将spectronaut搜库结果以如下名称命名：",
                                    "*.ProteinReport.csv , *.PeptideReport.csv , *.FragmentReport.csv",
                                    "例如HeLa.ProteinReport.csv , HeLa.PeptideReport.csv , HeLa.FragmentReport.csv",
                                    "然后进行zip压缩，上传至下面栏框。",style="color:black;"),
                           # p("请将spectronaut搜库结果以如下名称命名："),
                           # p("*.ProteinReport.csv , *.PeptideReport.csv , *.FragmentReport.csv"),
                           # p("例如HeLa.ProteinReport.csv , HeLa.PeptideReport.csv , HeLa.FragmentReport.csv"),
                           # p("然后进行zip压缩，上传至下面栏框。"),
                           fileInput("TrainProtein", "Upload your Train File",
                                     multiple = FALSE,
                                     accept = c(".zip")),
                           # fileInput("TrainPeptide", "Choose Peptide csv File",
                           #           multiple = FALSE,
                           #           accept = c(".csv")),
                           selectInput("TrainFasta", "Choose fasta:",
                                       c("human" = "human.fasta",
                                         "mouse" = "mouse.fasta")),
                           # fileInput("TrainFasta", "Choose Fasta File",
                           #           multiple = FALSE,
                           #           accept = c(".fasta")),
                           # fileInput("TrainFragment", "Choose Fragment csv File",
                           #           multiple = FALSE,
                           #           accept = c(".csv")),
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
  
  
  observeEvent(input$TrainConfirm,{
    if (is.null(input$TrainProtein)) {
      showModal(modalDialog(
        title = "Important message",
        paste0("请上传文件"),
        easyClose = T
      ))
    }else{
      shinyjs::removeCssClass(id = "id-progress-grey-out", class = "hide")
      shinyjs::addCssClass(id = "id-progress-grey-out", class = "show")
      runjs('$(\'.progress-grey-out\').css(\'height\',$(document).height()); $(\'.progress-grey-out\').css(\'width\',$(document).width()*1.2);')
      withProgress(message = '提交分析...', style = "notification", value = 0.01, {
        # #将文件复制到临时文件解压验证
        # tmp<-myFun(1)
        # dir.create(paste(tempdir(),tmp,sep = "/"))
        # setwd(paste(tempdir(),tmp,sep = "/"))
        # file.copy(from = input$TrainProtein$datapath,to = "./")
        # setProgress(1/3, message = '开始解压文件...')
        # system("unzip *.zip")
        setProgress(1/3, message = '开始检查文件...')
        x<-system(paste0("zip -sf ",input$TrainProtein$datapath),intern = T)
        
        if ((str_detect(trimws(x),pattern = ".*\\.ProteinReport\\.csv$")%>%sum)==1&(str_detect(trimws(x),pattern = ".*\\.PeptideReport\\.csv$")%>%sum)==1&(str_detect(trimws(x),pattern = ".*\\.FragmentReport\\.csv$")%>%sum)==1) {
          #证明成功即将压缩包复制到项目目录解压
          setProgress(2/3, message = '文件名格式检查正确，进行解压...')
          setwd("~/shiny-server/project/deepDIA/")
          modelid<<-myFun(1)
          dir.create(paste0("~/shiny-server/project/data/deepdia/",modelid))
          #在SQLite中加入项目信息状态为0
          projects_list <- dbConnect(RSQLite::SQLite(), "../db/projects_list.db", flags = SQLITE_RWC)
          dbSendQuery(projects_list, paste0("INSERT INTO Project VALUES ('", modelid, "','", "train", "', '", "deepDIA", "', '0', '", as.character(Sys.time()),"')"))
          dbDisconnect(projects_list)
          #进行文件解压和必要文件复制
          file.copy(input$TrainProtein$datapath,paste0("../data/deepdia/",modelid))
          system(paste0("unzip ../data/deepdia/",modelid,"/*.zip -d ../data/deepdia/",modelid))
          system(paste0("mv ../data/deepdia/",modelid,"/*.ProteinReport.csv ../data/deepdia/",modelid,"/train.ProteinReport.csv"))
          # file.copy(input$TrainPeptide$datapath,paste0("../data/deepdia/",modelid))
          system(paste0("mv ../data/deepdia/",modelid,"/*.PeptideReport.csv ../data/deepdia/",modelid,"/train.PeptideReport.csv"))
          # file.copy(input$TrainFragment$datapath,paste0("../data/deepdia/",modelid))
          system(paste0("mv ../data/deepdia/",modelid,"/*.FragmentReport.csv ../data/deepdia/",modelid,"/train.FragmentReport.csv"))
          system(paste0("cp -rf ../data/deepdia/data/* ","../data/deepdia/",modelid))
          system(paste0("cp -rf ../data/deepdia/database/",input$TrainFasta," ../data/deepdia/",modelid,"/train/detect"))
          #DeepDetectTrain(input$TrainProtein$datapath,input$TrainPeptide$datapath,input$TrainFasta)
          #DeepMsmsTrain(input$TrainFragment$datapath)
          setProgress(1, message = '解压完成，开始进行模型训练')
          showModal(modalDialog(
            title = "Important message",
            paste0("This is your Train Models Task ID: ",modelid),
            easyClose = FALSE,
            footer = NULL
          ))
        }else{
          showModal(modalDialog(
            title = "Important message",
            "上传失败，请检查文件或文件名完整性",
            easyClose = TRUE
          ))
        }
      })
      shinyjs::removeCssClass(id = "id-progress-grey-out", class = "show")
      shinyjs::addCssClass(id = "id-progress-grey-out", class = "hide")
      updateTabItems(session, "tabs", "Model Train")
    }
  })
  
  
  observeEvent(input$PredictFastaConfirm,{
    setwd("~/shiny-server/project/deepDIA/")
    if (is.null(input$PredictFasta)) {
      showModal(modalDialog(
        title = "Important message",
        "请上传fasta文件",
        easyClose = TRUE
      ))
    }else{
      if (file.exists(paste0("../data/deepdia/",input$TrainId1,"/train/msms/irt/success.txt"))) {
        shinyjs::removeCssClass(id = "id-progress-grey-out", class = "hide")
        shinyjs::addCssClass(id = "id-progress-grey-out", class = "show")
        runjs('$(\'.progress-grey-out\').css(\'height\',$(document).height()); $(\'.progress-grey-out\').css(\'width\',$(document).width()*1.2);')
        withProgress(message = '提交分析...', style = "notification", value = 0.01, {
          setProgress(1/4, message = '开始检查文件...')
          jobid<<-myFun(1)
          dir.create(paste0("../data/deepdia/",input$TrainId1,"/library/",jobid))
          projects_list <- dbConnect(RSQLite::SQLite(), "../db/projects_list.db", flags = SQLITE_RWC)
          dbSendQuery(projects_list, paste0("INSERT INTO Project VALUES ('", jobid, "','", paste0("fasta_predict_",trimws(input$TrainId1)), "', '", "deepDIA", "', '0', '", as.character(Sys.time()),"')"))
          dbDisconnect(projects_list)
          file.copy(input$PredictFasta$datapath,paste0("../data/deepdia/",input$TrainId1,"/library/",jobid))
          setwd(paste0("../data/deepdia/",input$TrainId1,"/library/",jobid))
          setProgress(2/3, message = '创建分布式任务...')
          seq <- seqinr::read.fasta("0.fasta",seqtype="AA")
          annotation <- getAnnot(seq)
          sequence <- getSequence(seq)
          index<-list()
          splitlength<-as.integer(length(annotation)/8)
          for (i in 1:7) {
            index[[i]]<-(splitlength*(i-1)+1):(i*splitlength)
          }
          index[[8]]<-(7*splitlength):length(annotation)
          lapply(1:8, function(x){
            dir.create(paste0("test",x))
            write.fasta(sequence[index[[x]]], annotation[index[[x]]], paste0('test',x,'.fasta'))
            system(paste("mv",paste0('test',x,'.fasta'),paste0("test",x)))
            if (dir.exists(paste0("../../models",x))) {
              print("file exist")
            }else{
              dir.create(paste0("../../models",x))
              system(paste("cp -rf","../../models/*",paste0("../../models",x)))
            }
          })
          setProgress(1, message = '任务创建成功，请保留您的建库ID！')
          setwd("~/shiny-server/project/deepDIA/")
        })
        shinyjs::removeCssClass(id = "id-progress-grey-out", class = "show")
        shinyjs::addCssClass(id = "id-progress-grey-out", class = "hide")
        updateTabItems(session, "tabs", "Fasta Predict")
        showModal(modalDialog(
          title = "Important message",
          paste0("This is your Predict LibraryTask ID: ",jobid),
          easyClose = FALSE,
          footer = NULL
        ))
      }else{
        if (input$TrainId1 %in% dir("../data/deepdia/")) {
          showModal(modalDialog(
            title = "Important message",
            "您的建模任务在运行中，请等待",
            easyClose = TRUE
          ))
        }else{
          showModal(modalDialog(
            title = "Important message",
            "请输入正确的建模任务ID",
            easyClose = TRUE
          ))
        }
      }
    }
  }
  )
  
  # observeEvent(input$PredictReportConfirm,{
  #   jobid<<-myFun(1)
  #   showModal(modalDialog(
  #     title = "Important message",
  #     paste0("This is your Library Task ID: ",jobid),
  #     easyClose = FALSE
  #   ))
  #   DeepReportPredict(input$PredictPeptide$datapath)
  #   #DeepReportLibrary<-DeepReportPredict(input$PredictPeptide$datapath)
  # })
  

  
  observeEvent(input$LibraryConfirm,{
    setwd("~/shiny-server/project/deepDIA/")
    if (file.exists(paste0("../data/deepdia/",input$TrainId2,"/library/",input$LibraryId,"/success.txt"))) {
        output$LibraryDownload<-renderUI({
          downloadButton("DownloadLibrary", "Download Library")
        })
        output$DownloadLibrary <- downloadHandler(
          filename = function() {
            paste("library_fasta_", format(Sys.Date(),"%Y%m%d"), ".zip", sep="")
          },
          content = function(file) {
            #file=read_csv(paste0("../data/deepdia/",input$TrainId2,"/library/",input$LibraryId,"/Report.prediction.library.csv"))
            file.copy(paste0("../data/deepdia/",input$TrainId2,"/library/",input$LibraryId,"/Report.prediction.library.zip"), file)
          }
        )
    }else{
      if ((input$TrainId2 %in% dir("../data/deepdia/"))&(input$LibraryId %in% dir(paste0("../data/deepdia/",input$TrainId2,"/library")))) {
        showModal(modalDialog(
          title = "Important message",
          "您的建库任务在运行中，请等待",
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

  #output$test <- renderUI({
    
    #paste0(input$TrainId2, ',', input$LibraryId, ',', getwd()) %>% print
    
  #})
}

shinyApp(ui, server)
