metr_pkgs<-c('shinyjs','shiny','rJava','PKI','digest','shinydashboard','shinycssloaders','RSQLite','jsonlite','dplyr','DT','V8','shinyBS','readr','tidyverse','mixOmics','gplots','pheatmap', 'plotly', 'ggplotify', 'dragulaR', 'factoextra', 'FactoMineR', 'ggsci', 'seqinr', 'ggrepel', 'ggplot2','mailR','shinyAce','cleaver','data.table')
for(i in 1:length(metr_pkgs)){
  library(metr_pkgs[i],character.only = TRUE)
}

shinyUI(
  ui<-dashboardPage(
    dashboardHeader(title = "Project System", 
                    dropdownMenuOutput("messageMenu"),
                    tags$li(class = "dropdown", style = "font-weight:bold;", tags$a(textOutput("usertext"), href="../dashboard")),
                    tags$li(class = "dropdown", style = "font-weight:bold;", tags$a(icon("sign-out"), "Logout", onclick ="Shiny.onInputChange('logout_switch', new Date().toLocaleString());", href="#"))
    ),
    dashboardSidebar(
      sidebarMenuOutput("menu")
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
             uiOutput("dashboard_panel")
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
)
