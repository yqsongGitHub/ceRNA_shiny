#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
options(shiny.maxRequestSize=70*1024^2)
options(stringsAsFactors = F)

library(shiny)
library(shinythemes)
library(reshape2)
library(edgeR)
library(DESeq2)
library(limma)
library(pheatmap)
library(ggplot2)
library(DT)
library(yulab.utils)
library(rvcheck)
library(KEGG.db)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(igraph)
# library(openxlsx)
source("./script/all_function.R")

# Define UI for application 
ui <- fluidPage(
  theme = shinytheme('readable'),
  #shinythemes::themeSelector()
  #theme = shinytheme('sandstone'),
  navbarPage("ceRNA", 
             tabPanel("Step by Step",
                      fluidRow(
                        column(2,
                               wellPanel(
                                 h4(strong("Input your expression file")),
                                 selectInput("file1",label= "Choose an example or your own data", 
                                             choices = c("Example_array" = "Example1",
                                                         "Example_RNAseq" = "Example1_RNAseq",
                                                         "Your own data" = "load_my_own1")),
                                 conditionalPanel("input.file1 == 'Example1'",
                                                  downloadButton('downloadEx1', 'Download example')),
                                 conditionalPanel("input.file1 == 'Example1_RNAseq'",
                                                  downloadButton('downloadEx1_RNAseq', 'Download example')),
                                 conditionalPanel("input.file1 == 'load_my_own1'",
                                                  fileInput('loadfile1', 'Choose xlsx File', 
                                                            accept=c('.xlsx','text/csv', 'text/comma-separated-values,text/plain', '.csv', '.txt')))
                               ),
                               wellPanel(
                                 h4(strong("Input your sample information file")),
                                 selectInput("file2",label= "Choose an example or your own data", 
                                             choices = c("Example_array"="Example2",
                                                         "Example_RNAseq"="Example2_RNAseq",
                                                         "Your own data" = "load_my_own2")),
                                 conditionalPanel("input.file2 == 'Example2'",
                                                  downloadButton('downloadEx2', 'Download example')),
                                 conditionalPanel("input.file2 == 'Example2_RNAseq'",
                                                  downloadButton('downloadEx2_RNAseq', 'Download example')),
                                 conditionalPanel("input.file2 == 'load_my_own2'",
                                                  fileInput('loadfile2', 'Choose xlsx File', 
                                                            accept=c('.xlsx','text/csv', 'text/comma-separated-values,text/plain', '.csv', '.txt')))
                               ),
                               wellPanel(
                                 h4(strong("Input your gene ID information file")),
                                 selectInput("file3",label= "Choose an example or your own data", 
                                             choices = c("Example_array"="Example3", 
                                                         "Example_RNAseq"="Example3_RNAseq",
                                                         "Your own data" = "load_my_own3")),
                                 conditionalPanel("input.file3 == 'Example3'",
                                                  downloadButton('downloadEx3', 'Download example')),
                                 conditionalPanel("input.file3 == 'Example3_RNAseq'",
                                                  downloadButton('downloadEx3_RNAseq', 'Download example')),
                                 conditionalPanel("input.file3 == 'load_my_own3'",
                                                  fileInput('loadfile3', 'Choose xlsx File', 
                                                            accept=c('.xlsx','text/csv', 'text/comma-separated-values,text/plain', '.csv', '.txt')))
                               ),
   
                               conditionalPanel("input.cPanels1 == 2",
                                                wellPanel(
                                                  h4(strong("Heatmap")),
                                                  selectInput("select_heatmap_filetype", 
                                                              label= "The file type:",
                                                              choices=c("array", "RNAseq"),
                                                              selected = "array",
                                                              multiple = FALSE),
                                                  # sliderInput("select_heatmap_expression_threshold", 
                                                  #             label = "The cutoff of gene expression:",
                                                  #             min = 0, max = 10, value = 2, step = 1),
                                                  sliderInput("select_heatmap_pvalue", 
                                                              label = "The cutoff of pvalue:",
                                                              min = 0, max = 1, value = 0.05, step = 0.01),
                                                  sliderInput("select_heatmap_fc", 
                                                              label = "The cutoff of fold change:",
                                                              min = 0, max = 5, value = 1, step = 0.1),
                                                  selectInput("select_heatmap_Scale", 
                                                              label= "scale",
                                                              choices=c("row", "column", "none"),
                                                              selected = "row",
                                                              multiple = FALSE),
                                                  checkboxInput("select_heatmap_Clusterrows", 
                                                                label= "Clustered by rows",
                                                                TRUE),
                                                  checkboxInput("select_heatmap_Clustercols", 
                                                                label= "Clustered by columns",
                                                                TRUE),
                                                  checkboxInput("select_Annotation_legend", 
                                                                label= "Show annatation legend",
                                                                TRUE),
                                                  checkboxInput("select_Showrownames", 
                                                                label= "Show row names",
                                                                FALSE),
                                                  checkboxInput("select_heatmap_inputemiss", 
                                                                label= "Inpute missing value",
                                                                FALSE)
                                                )),
                               
                               conditionalPanel("input.cPanels1 == 2",
                                                h4(strong("Download")),
                                                wellPanel(
                                                  textInput("fnameHeatmap", "filename", value = "Heatmap"),
                                                  downloadButton('DownloadHeatmap', 'Download Heatmap'),
                                                  br(),
                                                  br(),
                                                  downloadButton('DownloadDEGheatmaptable', 'Download DEG table')
                                                )),
                               
                               
                               conditionalPanel("input.cPanels1 == 3",
                                                wellPanel(
                                                  h4(strong("Volcano plot")),
                                                  selectInput("select_volcano_filetype", 
                                                              label= "The file type:",
                                                              choices=c("array", "RNAseq"),
                                                              selected = "array",
                                                              multiple = FALSE),
                                                  # sliderInput("select_volcano_expression_threshold", 
                                                  #             label = "The cutoff of gene expression:",
                                                  #             min = 0, max = 10, value = 2, step = 1),
                                                  sliderInput("select_volcano_pvalue", 
                                                              label = "The cutoff of pvalue:",
                                                              min = 0, max = 1, value = 0.05, step = 0.01),
                                                  sliderInput("select_volcano_fc", 
                                                              label = "The cutoff of fold change:",
                                                              min = 0, max = 5, value = 1.5, step = 0.1),
                                                  checkboxInput("select_volcano_inputemiss", 
                                                                label= "Inpute missing value",
                                                                FALSE),
                                                  textInput("select_volcano_Title", "Title name", value = "PD vs NC")
                                                )),
                               
                               conditionalPanel("input.cPanels1 == 3",
                                                h4(strong("Download")),
                                                wellPanel(
                                                  textInput("fnamevolcano", "filename", value = "Volcano"),
                                                  downloadButton('Downloadvolcano', 'Download volcano plot'),
                                                  br(),
                                                  br(),
                                                  downloadButton('DownloadDEGvolcanotable', 'Download DEG table')
                                                )),

                               conditionalPanel("input.cPanels1 == 4",
                                                wellPanel(
                                                  h4(strong("Enrichment")),
                                                  selectInput("select_enrichment_filetype", 
                                                              label= "The file type:",
                                                              choices=c("array", "RNAseq"),
                                                              selected = "array",
                                                              multiple = FALSE),
                                                  selectInput("select_enrichment_genesymbol", 
                                                              label= "The up or down regulated differential gene:",
                                                              choices=c("up", "down"),
                                                              selected = "up",
                                                              multiple = FALSE),
                                                  selectInput("select_enrichment_methods", 
                                                              label= "The methods of enrichment:",
                                                              choices=c("KEGG","GO","Reactome"),
                                                              selected = "KEGG",
                                                              multiple = FALSE),
                                                  sliderInput("select_enrichment_pvaluecutoff", 
                                                              label = "The pvalue of enrichment:",
                                                              min = 0, max = 1, value = 1, step = 0.01),
                                                  sliderInput("select_enrichment_qvaluecutoff", 
                                                              label = "The qvalue of enrichment:",
                                                              min = 0, max = 1, value = 1, step = 0.01),
                                                  sliderInput("select_enrichment_showCategory", 
                                                              label = "The number of category:",
                                                              min = 1, max = 10, value = 5, step = 1)
                                                )),
                               
                               conditionalPanel("input.cPanels1 == 4",
                                                h4(strong("Download")),
                                                wellPanel(
                                                  textInput("fnameenrichment", "filename", value = "Enrichment"),
                                                  downloadButton('Downloadenrichment', 'Download enrichment plot'),
                                                  br(),
                                                  br(),
                                                  downloadButton('Downloadenrichmenttable', 'Download enrichment table')
                                                )),
                               
                               conditionalPanel("input.cPanels1 == 5",
                                                wellPanel(
                                                  h4(strong("Classification")),
                                                  selectInput("select_classification_filetype", 
                                                              label= "The file type:",
                                                              choices=c("array", "RNAseq"),
                                                              selected = "array",
                                                              multiple = FALSE),
                                                  selectInput("select_classification_type", 
                                                              label= "Show the type of  differential genes:",
                                                              choices=c("All", "mRNA","miRNA","lncRNA"),
                                                              selected = "mRNA",
                                                              multiple = FALSE),
                                                  sliderInput("select_classification_pvalue", 
                                                              label = "The pvalue of differential genes:",
                                                              min = 0, max = 1, value = 0.05, step = 0.01),
                                                  sliderInput("select_classification_fc", 
                                                              label = "The fold change of differential genes:",
                                                              min = 0, max = 5, value = 1.5, step = 0.1)
                                                )),
                               
                               conditionalPanel("input.cPanels1 == 5",
                                                h4(strong("Download")),
                                                wellPanel(
                                                  textInput("fnameclassification", "filename", value = "classification"),
                                                  downloadButton('Downloadclassificationtable', 'Download classification table')
                                                )),
                               
                               # conditionalPanel("input.cPanels1 == 6",
                               #                  wellPanel(
                               #                    h4(strong("DEG Correlation")),
                               #                    selectInput("select_correlation_filetype", 
                               #                                label= "The file type:",
                               #                                choices=c("array", "RNAseq"),
                               #                                selected = "array",
                               #                                multiple = FALSE),
                               #                    sliderInput("select_calculate_cor_share_miRNAs", 
                               #                                label = "The number of share miRNAs:",
                               #                                min = 0, max = 10, value = 3, step = 1),
                               #                    sliderInput("select_calculate_cor_padj", 
                               #                                label = "The p adjusted value of positive correlation:",
                               #                                min = 0, max = 1, value = 0.05, step = 0.01),
                               #                    sliderInput("select_calculate_cor_sc", 
                               #                                label = "The sensitivity partial pearson correlation:",
                               #                                min = 0, max = 1, value = 0.3, step = 0.01),
                               #                    sliderInput("select_calculate_cor_corvalue", 
                               #                                label = "The value of correlation:",
                               #                                min = -1, max = 1, value =0, step = 0.1),
                               #                    # enrichment
                               #                    selectInput("select_cor_enrichment_methods", 
                               #                                label= "The methods of enrichment:",
                               #                                choices=c("KEGG"),
                               #                                selected = "KEGG",
                               #                                multiple = FALSE),
                               #                    sliderInput("select_cor_enrichment_pvaluecutoff", 
                               #                                label = "The pvalue of enrichment:",
                               #                                min = 0, max = 1, value = 1, step = 0.01),
                               #                    sliderInput("select_cor_enrichment_qvaluecutoff", 
                               #                                label = "The qvalue of enrichment:",
                               #                                min = 0, max = 1, value = 1, step = 0.01)
                               #                  )),
                               # 
                               # conditionalPanel("input.cPanels1 == 6",
                               #                  h4(strong("Download")),
                               #                  wellPanel(
                               #                    textInput("fnamecorrelation", "filename", value = "correlation"),
                               #                    downloadButton('Downloadcorrelationtable', 'Download correlation table'),
                               #                    br(),
                               #                    br(),
                               #                    downloadButton('Downloadcorenrichment', 'Download enrichment plot'),
                               #                    br(),
                               #                    br(),
                               #                    downloadButton('Downloadecorenrichmenttable', 'Download enrichment table')
                               #                  )),
                               conditionalPanel("input.cPanels1 == 6",
                                                wellPanel(
                                                  h4(strong("DEG Correlation")),
                                                  selectInput("select_correlation_filetype", 
                                                              label= "The file type:",
                                                              choices=c("array", "RNAseq"),
                                                              selected = "array",
                                                              multiple = FALSE),
                                                  selectInput("select_correlation_algorithm", 
                                                              label= "The algorithm:",
                                                              choices=c("PCC", "SPCC","PC"),
                                                              selected = "PCC",
                                                              multiple = FALSE),
                                                  conditionalPanel("input.select_correlation_algorithm == 'PCC'",
                                                                   sliderInput("select_calculate_cor_pvalue",
                                                                               label = "The pvalue of correlation:",
                                                                               min = 0, max = 1, value = 0.05, step = 0.01),
                                                                   sliderInput("select_calculate_correlation",
                                                                               label = "The value of correlation:",
                                                                               min = 0, max = 1, value = 0.4, step = 0.01)),
                                                  conditionalPanel("input.select_correlation_algorithm == 'SPCC'",
                                                                   sliderInput("select_calculate_cor_SPCC",
                                                                               label = "The value of correlation:",
                                                                               min = 0, max = 1, value = 0.1, step = 0.01),
                                                                   sliderInput("select_calculate_cor_share_miRNAs_SPCC", 
                                                                               label = "The number of share miRNAs:",
                                                                               min = 0, max = 10, value = 2, step = 1),
                                                                   sliderInput("select_calculate_cor_padj_SPCC", 
                                                                               label = "The p adjusted value of positive correlation:",
                                                                               min = 0, max = 1, value = 0.05, step = 0.01),
                                                                   sliderInput("select_calculate_cor_sc_SPCC", 
                                                                               label = "The sensitivity partial pearson correlation:",
                                                                               min = 0, max = 1, value = 0.3, step = 0.01)),
                                                  conditionalPanel("input.select_correlation_algorithm == 'PC'",
                                                                   sliderInput("select_calculate_cor_PC",
                                                                               label = "The value of correlation:",
                                                                               min = 0, max = 1, value = 0.1, step = 0.01),
                                                                   sliderInput("select_calculate_cor_share_miRNAs_PC", 
                                                                               label = "The number of share miRNAs:",
                                                                               min = 0, max = 10, value = 2, step = 1),
                                                                   sliderInput("select_calculate_cor_sc_PC", 
                                                                               label = "The p adjusted value of positive correlation:",
                                                                               min = 0, max = 1, value = 0.3, step = 0.01))

                                                )),
                               
                               conditionalPanel("input.cPanels1 == 6",
                                                h4(strong("Download")),
                                                wellPanel(
                                                  textInput("fnamecorrelation", "filename", value = "correlation"),
                                                  downloadButton('Downloadcorrelationtable', 'Download correlation table'),
                                                  br(),
                                                  br(),
                                                  downloadButton('Downloadcorenrichment', 'Download enrichment plot'),
                                                  br(),
                                                  br(),
                                                  downloadButton('Downloadecorenrichmenttable', 'Download enrichment table')
                                                )),
                               
                               conditionalPanel("input.cPanels1 == 7",
                                                wellPanel(
                                                  h4(strong("Prediction")),
                                                  selectInput("select_prediction_filetype", 
                                                              label= "The file type:",
                                                              choices=c("array", "RNAseq"),
                                                              selected = "array",
                                                              multiple = FALSE),
                                                  selectInput("select_prediction_type", 
                                                              label= "Show the type of prediction:",
                                                              choices=c("miRNA-lnc","miRNA-mRNA","miRNA-lnc-mRNA"),
                                                              selected = "miRNA-mRNA",
                                                              multiple = FALSE),
                                                  selectInput("select_prediction_database", 
                                                              label= "The type of database:",
                                                              choices=c("starBase v2.0"),
                                                              selected = "starBase v2.0",
                                                              multiple = FALSE),
                                                  # enrichment
                                                  selectInput("select_prediction_enrichment_methods", 
                                                              label= "The methods of enrichment:",
                                                              choices=c("KEGG","GO"),
                                                              selected = "KEGG",
                                                              multiple = FALSE),
                                                  sliderInput("select_prediction_enrichment_pvaluecutoff", 
                                                              label = "The pvalue of enrichment:",
                                                              min = 0, max = 1, value = 1, step = 0.01),
                                                  sliderInput("select_prediction_enrichment_qvaluecutoff", 
                                                              label = "The qvalue of enrichment:",
                                                              min = 0, max = 1, value = 1, step = 0.01)
                                                )),
                               
                               conditionalPanel("input.cPanels1 == 7",
                                                h4(strong("Download")),
                                                wellPanel(
                                                  textInput("fnameprediction", "filename", value = "prediction"),
                                                  downloadButton('Downloadpredictiontable', 'Download prediction table'),
                                                  br(),
                                                  br(),
                                                  downloadButton('Downloadpreenrichment', 'Download enrichment plot'),
                                                  br(),
                                                  br(),
                                                  downloadButton('Downloadepreenrichmenttable', 'Download enrichment table')
                                                )),
                               
                               conditionalPanel("input.cPanels1 == 8",
                                                wellPanel(
                                                  h4(strong("Network")),
                                                  selectInput("select_ce_filetype", 
                                                              label= "The file type:",
                                                              choices=c("array", "RNAseq"),
                                                              selected = "array",
                                                              multiple = FALSE),
                                                  selectInput("select_layout_ce", 
                                                              label= "The algorithm of layout:",
                                                              choices=c("fruchterman.reingold",
                                                                        "circle",
                                                                        "reingold.tilford",
                                                                        "star",
                                                                        "with_drl",
                                                                        "random"),
                                                              selected = "fruchterman.reingold",
                                                              multiple = FALSE),
                                                  selectInput("select_vertex.shape_ce", 
                                                              label= "The vertex shape:",
                                                              choices=c("circle","rectangle","none"),
                                                              selected = "circle",
                                                              multiple = FALSE),
                                                  # sliderInput("select_network_DEG_p", 
                                                  #             label = "The pvalue of differential genes:",
                                                  #             min = 0, max = 1, value = 0.05, step = 0.01),
                                                  # sliderInput("select_network_DEG_fc", 
                                                  #             label = "The fold change of differential genes:",
                                                  #             min = 0, max = 5, value = 1, step = 0.1),
                                                  # sliderInput("select_network_cor_pvalue", 
                                                  #             label = "The pvalue of correlation:",
                                                  #             min = 0, max = 1, value = 0.05, step = 0.01),
                                                  # sliderInput("select_network_cor_corvalue", 
                                                  #             label = "The value of correlation:",
                                                  #             min = -1, max = 1, value = -0.4, step = 0.01),
                                                  sliderInput("select_vertex.size_ce", 
                                                              label = "The value of vertex size:",
                                                              min = 1, max = 10, value = 7, step = 1),
                                                  sliderInput("select_vertex.label.cex_ce", 
                                                              label = "The value of vertex label size:",
                                                              min = 0.1, max = 1.5, value = 0.9, step = 0.1),
                                                  sliderInput("select_vertex.label.dist_ce", 
                                                              label = "The value of vertex label dist:",
                                                              min = 0, max = 2, value = 0, step = 0.1),
                                                  sliderInput("select_edge.width_ce", 
                                                              label = "The value of edge width:",
                                                              min = 0, max = 2, value = 0.2, step = 0.1),
                                                  # enrichment
                                                  selectInput("select_network_enrichment_methods", 
                                                              label= "The methods of enrichment:",
                                                              choices=c("KEGG","GO"),
                                                              selected = "KEGG",
                                                              multiple = FALSE),
                                                  sliderInput("select_network_enrichment_pvaluecutoff", 
                                                              label = "The pvalue of enrichment:",
                                                              min = 0, max = 1, value = 1, step = 0.01),
                                                  sliderInput("select_network_enrichment_qvaluecutoff", 
                                                              label = "The qvalue of enrichment:",
                                                              min = 0, max = 1, value = 1, step = 0.01)
                                                )),
                               
                               conditionalPanel("input.cPanels1 == 8",
                                                h4(strong("Download")),
                                                wellPanel(
                                                  textInput("fnamenetwork", "filename", value = "Network"),
                                                  downloadButton('Downloadnetwork', 'Download network plot'),
                                                  br(),
                                                  br(),
                                                  downloadButton('Downloadceenrichment', 'Download enrichment'),
                                                  br(),
                                                  br(),
                                                  downloadButton('Downloadceenrichmentable', 'Download enrichment plot'),
                                                ))
                               
                        ),
                        
                        column(10,
                               tabsetPanel(
                                 tabPanel("Manual", htmlOutput("ReadMe1"), value =1),
                                 tabPanel("Heatmap plot", plotOutput("heatmap1", height= 800, width =800), DT::dataTableOutput("DEGheatmaptable",width = 800),value = 2),
                                 tabPanel("Volcano plot", plotOutput("volcano", height= 800, width = 1000), DT::dataTableOutput("DEGvolcanotable",width = 800),value = 3),
                                 tabPanel("Enrichment", plotOutput("enrichment", height= 600, width = 1000),DT::dataTableOutput("enrichmenttable",width = 800),value = 4),
                                 tabPanel("RNA classification", DT::dataTableOutput("classificationtable",width = 800),value = 5),
                                 tabPanel("ceRNA network1", DT::dataTableOutput("correlationtable",width = 500),plotOutput("correlationEnrichment",height= 500,width = 1000),
                                          #DT::dataTableOutput("corenrichmenttable",height= 400,width =600),
                                          value = 6),                                 
                                 tabPanel("ceRNA network2", DT::dataTableOutput("predictiontable",width = 1000),plotOutput("predictionEnrichment", height= 500, width = 1000),value = 7),
                                 tabPanel("ceRNA network3",  textOutput("networktext"),plotOutput("network", height=800, width = 800),
                                          plotOutput("networkEnrichment", height= 500, width =1000),value = 8),
                                 id = "cPanels1"
                               )                
                               
                        ),
                        
                        column(12,
                               tags$head(tags$style(type="text/css", "
                                                    #loadmessage {
                                                    position: fixed;
                                                    bottom: 0px;
                                                    right: 0px;
                                                    width: 100%;
                                                    padding: 5px 0px 5px 0px;
                                                    text-align: center; 
                                                    font-weight: bold;
                                                    font-size: 100%;
                                                    color: #000000;
                                                    background-color: #b8b8b8;
                                                    z-index: 105;
                                                    }
                                                    ")),
                               conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                                                tags$div("Loading...",id="loadmessage"))
                        )
                      )
             ),
             # another Panel
  )
)

# Define server
server <- function(input, output, session) {
  data_input1 <- reactive({
    if(input$file1 == 'Example1'){
      d1 <- read.csv("./example/GSE7621.csv",row.names = 1)
    }else if(input$file1 == 'Example1_RNAseq'){
      d1 <- read.csv("./example/GSE136666.csv",row.names = 1)
    }else if(input$file1 == 'load_my_own1'){
      inFile <- input$loadfile1
      if (is.null(inFile))
        return(NULL)
      else if(grepl(".xlsx", inFile[1])) { 
        d1 <-  read.xlsx(as.character(inFile$datapath), colNames = TRUE, rowNames = F)
      }
      else if(grepl(".csv", inFile[1])) { 
        d1 <-  read.csv(as.character(inFile$datapath), header = TRUE, sep = ",", stringsAsFactors = F, as.is = T, fill = T) 
      }
      else if(grepl(".txt", inFile[1])) { 
        d1 <-  read.table(as.character(inFile$datapath), header = TRUE, sep = "\t", stringsAsFactors = F, as.is = T, fill = T)
      }
    }
    else 
      return(NULL)
    Dataset1 <- data.frame(d1)
    return(as.data.frame(Dataset1))
  })
  
  output$downloadEx1 <- downloadHandler( 
    filename <- function() {
      paste0('GSE7621','.csv')
    },
    content <- function(file) {
      ds1 <- data_input1()
      write.csv(ds1, file)
    }
  )
  
  output$downloadEx1_RNAseq <- downloadHandler( 
    filename <- function() {
      paste0('GSE136666','.csv')
    },
    content <- function(file) {
      ds1 <- data_input1()
      write.csv(ds1, file)
    }
  )
  
  data_input2 <- reactive({
    if(input$file2 == 'Example2'){
      d2 <- read.csv("./example/samplelist.csv")
    }else if(input$file2 == 'Example2_RNAseq'){
      d2 <- read.csv("./example/samplelist_RNAseq.csv")
    }else if(input$file2 == 'load_my_own2'){
      inFile <- input$loadfile2
      if (is.null(inFile))
        return(NULL)
      else if(grepl(".xlsx", inFile[1])) { 
        d2 <-  read.xlsx(as.character(inFile$datapath), colNames = TRUE, rowNames = F) 
      }
      else if(grepl(".csv", inFile[1])) { 
        d2 <-  read.csv(as.character(inFile$datapath), header = TRUE, sep = ",", stringsAsFactors = F, as.is = T, fill = T) 
      }
      else if(grepl(".txt", inFile[1])) { 
        d2 <-  read.table(as.character(inFile$datapath), header = TRUE, sep = "\t", stringsAsFactors = F, as.is = T, fill = T)
      }
    }
    else 
      return(NULL)
    Dataset2 <- data.frame(d2)
    return(as.data.frame(Dataset2))
  })
  
  output$downloadEx2 <- downloadHandler( 
    filename <- function() {
      paste0("samplelist",".csv")
    },
    content <- function(file) {
      ds2 <- data_input2()
      write.csv(ds2, file, row.names = F)
    }
  )
  
  output$downloadEx2_RNAseq <- downloadHandler( 
    filename <- function() {
      paste0("samplelist",".csv")
    },
    content <- function(file) {
      ds2 <- data_input2()
      write.csv(ds2, file, row.names = F)
    }
  )
  
  data_input3 <- reactive({
    if(input$file3 == 'Example3'){
      d3 <- read.csv("./data/GPL570.csv")
    }else if(input$file3 == 'Example3_RNAseq'){
      d3 <- read.csv("./data/Novaseq6000.csv")
    }else if(input$file3 == 'load_my_own3'){
      inFile <- input$loadfile3
      if (is.null(inFile))
        return(NULL)
      else if(grepl(".xlsx", inFile[1])) { 
        d3 <-  read.xlsx(as.character(inFile$datapath), colNames = TRUE, rowNames = F)
      }
      else if(grepl(".csv", inFile[1])) { 
        d3 <-  read.csv(as.character(inFile$datapath), header = TRUE, sep = ",", stringsAsFactors = F, as.is = T, fill = T) 
      }
      else if(grepl(".txt", inFile[1])) { 
        d3 = read.table(as.character(inFile$datapath), header = TRUE, sep = "\t", stringsAsFactors = F, as.is = T, fill = T) 
      }
    }
    else 
      return(NULL)
    Dataset3 <- data.frame(d3)
    return(as.data.frame(Dataset3))
  })
  
  output$downloadEx3 <- downloadHandler( 
    filename <- function() {
      paste0("GPL570",".csv")
    },
    content <- function(file) {
      ds3 <- data_input3()
      write.csv(ds3, file, row.names = F)
    }
  )
  
  output$downloadEx3_RNAseq <- downloadHandler( 
    filename <- function() {
      paste0("Novaseq6000",".csv")
    },
    content <- function(file) {
      ds3 <- data_input3()
      write.csv(ds3, file, row.names = F)
    }
  )
  
  RNAmap_input  <- reactive({read.csv("./data/RNAmap.csv")})
  #miRNA_trans_input <- reactive({read.csv("./data/miRNA_database.csv",row.names = 1,header = T)})
  #mi_lnc_input <- reactive({read.csv("./data/database-miRNA-LncRNA.csv",header = T)})
  #mi_m_input <- reactive({read.csv("./data/database_mRNA_miRNAN.csv",header = T)})
  
  observe({
    dsnames1 <- colnames(data_input1())
  })
  
  #pheatmap
  heatmapforuse0 <- reactive(
    if(input$select_heatmap_filetype=="array"){
      if(!input$select_heatmap_inputemiss){
        heatmap.data <- read.csv("./data/preprocessing/DEG/DEG_array.csv",row.names = 1)
      }else{
        heatmap.data <- read.csv("./data/preprocessing/DEG/DEG_array_inputmissvalue.csv",row.names = 1)
      }
    }else{
      if(!input$select_heatmap_inputemiss){
        heatmap.data <- read.csv("./data/preprocessing/DEG/DEG_RNA.csv",row.names = 1)
      }else{
        heatmap.data <- read.csv("./data/preprocessing/DEG/DEG_RNA_inputmissvalue.csv",row.names = 1)
      }
    }
  )
  
  pheat.data1 <- reactive(
    if(input$select_heatmap_filetype=="array"){
      data1 <- read.csv("./example/GSE7621.csv",row.names = 1)
    }else{
      data1 <- read.csv("./example/GSE136666.csv",row.names = 1)
    }
  )
  pheat.data2 <- reactive(
    if(input$select_heatmap_filetype=="array"){
      data2 <- read.csv("./example/samplelist.csv")
    }else{
      data2 <- read.csv("./example/samplelist_RNAseq.csv")
    }
  )
  
  heatmapforuse <- function(){
    data1 <- pheat.data1()
    #data1 <- data1[rowSums(data1>1) >input$select_heatmap_expression_threshold,]
    data1 <- data1[rowSums(data1>1) >2,]
    data2 <-  pheat.data2()
    DEGG <- heatmapforuse0()
    DEGG <- DEGG[DEGG[["P.Value"]]<input$select_heatmap_pvalue,]
    DEGG <- DEGG[DEGG[["logFC"]]<log(input$select_heatmap_fc),]
    DEGpheatmap(data1,
                data2,
                DEGG,
                Scale=input$select_heatmap_Scale,
                Clusterrows=input$select_heatmap_Clusterrows,
                Clustercols=input$select_heatmap_Clustercols,
                Annotation_legend=input$select_Annotation_legend,
                Showrownames=input$select_Showrownames
    )
  }

  output$heatmap1 <- renderPlot({
    print(heatmapforuse())
  })

  output$DownloadHeatmap <- downloadHandler(
    filename <- function() {
      pdf_file <<- as.character(input$fnameHeatmap)
      paste(pdf_file,".pdf", sep='')
    },
    content <- function(file) {
      pdf(file, height= 10, width=12,onefile = FALSE)
      print(heatmapforuse())
      dev.off()
      #file.copy(paste(pdf_file,'.pdf', sep='') ,file, overwrite=TRUE)
    }, contentType = 'image/pdf')

  output$DEGheatmaptable <- renderDataTable({
    data1 <- pheat.data1()
    #data1 <- data1[rowSums(data1>1) >input$select_heatmap_expression_threshold,]
    data1 <- data1[rowSums(data1>1) >2,]
    DEGG <- heatmapforuse0()
    DEGG <- DEGG[DEGG[["P.Value"]]<input$select_heatmap_pvalue,]
    DEGG <- DEGG[DEGG[["logFC"]]<log(input$select_heatmap_fc),]
    DEGG <- DEGG[rownames(DEGG) %in% rownames(data1),]
    DEGG
  })

  output$DownloadDEGheatmaptable <- downloadHandler(
    filename <- function() {
      pdf_file <<- as.character(input$fnameHeatmap)
      paste0(pdf_file,'.csv')
    },
    content <- function(file) {
      data1 <- pheat.data1()
      #data1 <- data1[rowSums(data1>1) >input$select_heatmap_expression_threshold,]
      data1 <- data1[rowSums(data1>1) >2,]
      DEGG <- heatmapforuse0()
      DEGG <- DEGG[DEGG[["P.Value"]]<input$select_heatmap_pvalue,]
      DEGG <- DEGG[DEGG[["logFC"]]<log(input$select_heatmap_fc),]
      DEGG <- DEGG[rownames(DEGG) %in% rownames(data1),]
      DEGG
      write.csv(DEGG,file)
    })
  
  ##volcano
  volcanoforuse0 <- reactive(
    if(input$select_volcano_filetype=="array"){
      if(!input$select_volcano_inputemiss){
        vol.data <- read.csv("./data/preprocessing/DEG/DEG_array.csv",row.names = 1)
      }else{
        vol.data <- read.csv("./data/preprocessing/DEG/DEG_array_inputmissvalue.csv",row.names = 1)
      }
    }else{
      if(!input$select_volcano_inputemiss){
        vol.data <- read.csv("./data/preprocessing/DEG/DEG_RNA.csv",row.names = 1)
      }else{
        vol.data <- read.csv("./data/preprocessing/DEG/DEG_RNA_inputmissvalue.csv",row.names = 1)
      }
    }
  )
  vol.data1 <- reactive(
    if(input$select_volcano_filetype=="array"){
      data1 <- read.csv("./example/GSE7621.csv",row.names = 1)
    }else{
      data1 <- read.csv("./example/GSE136666.csv",row.names = 1)
    }
  )
  
  volcanoforuse <- function(){
    data1 <- vol.data1()
    #data1 <- data1[rowSums(data1>1) >input$select_volcano_expression_threshold,]
    data1 <- data1[rowSums(data1>1) >2,]
    DEGG <- volcanoforuse0()
    DEGG <- DEGG[rownames(DEGG) %in% rownames(data1),]
    DEGvolcano(DEGG,
               pvalue=input$select_volcano_pvalue,
               fc=input$select_volcano_fc,
               Title=input$select_volcano_Title)
  }
  
  output$volcano <- renderPlot({
    volcanoforuse()
  })
  
  output$Downloadvolcano <- downloadHandler(
    filename <- function() {
      pdf_file <<- as.character(input$fnamevolcano)
      paste(pdf_file,".pdf", sep='')
    },
    content <- function(file) {
      pdf(file, height= 10, width=12,onefile = FALSE)
      print(volcanoforuse())
      dev.off()
      #file.copy(paste(pdf_file,'.pdf', sep='') ,file, overwrite=TRUE)
    }, contentType = 'image/pdf')
  
  output$DEGvolcanotable <- renderDataTable({
    data1 <- vol.data1()
    #data1 <- data1[rowSums(data1>1) >input$select_volcano_expression_threshold,]
    data1 <- data1[rowSums(data1>1) >2,]
    DEGG <- volcanoforuse0()
    DEGG <- DEGG[rownames(DEGG) %in% rownames(data1),]
    DEGG
  })
  
  output$DownloadDEGvolcanotable <- downloadHandler(
    filename <- function() {
      pdf_file <<- as.character(input$fnamevolcano)
      paste0(pdf_file,'.csv')
    },
    content <- function(file) {
      data1 <- vol.data1()
      #data1 <- data1[rowSums(data1>1) >input$select_volcano_expression_threshold,]
      data1 <- data1[rowSums(data1>1) >2,]
      DEGG <- volcanoforuse0()
      DEGG <- DEGG[rownames(DEGG) %in% rownames(data1),]
      write.csv(DEGG,file)
    })

  # enrichment
  enrichmentforuse0  <- reactive(
    if(input$select_enrichment_filetype == 'array'){
      if(input$select_enrichment_genesymbol=="up"){
        if(input$select_enrichment_methods=="KEGG"){
          enr <- readRDS("./data/preprocessing/enrichment/DEG_array_up_kegg.rds")
        }else if(input$select_enrichment_methods=="GO"){
          enr<- readRDS("./data/preprocessing/enrichment/DEG_array_up_GO.rds")
        }else if(input$select_enrichment_methods=="Reactome"){
          enr<- readRDS("./data/preprocessing/enrichment/DEG_array_up_Reactome.rds")
        }
      }else{
        if(input$select_enrichment_methods=="KEGG"){
          enr <- readRDS("./data/preprocessing/enrichment/DEG_array_down_kegg.rds")
        }else if(input$select_enrichment_methods=="GO"){
          enr<- readRDS("./data/preprocessing/enrichment/DEG_array_down_GO.rds")
        }else if(input$select_enrichment_methods=="Reactome"){
          enr<- readRDS("./data/preprocessing/enrichment/DEG_array_down_Reactome.rds")
        }
      }
    }else if(input$select_enrichment_filetype == 'RNAseq'){
      if(input$select_enrichment_genesymbol=="up"){
        if(input$select_enrichment_methods=="KEGG"){
          enr <- readRDS("./data/preprocessing/enrichment/DEG_RNA_up_kegg.rds")
        }else if(input$select_enrichment_methods=="GO"){
          enr<- readRDS("./data/preprocessing/enrichment/DEG_RNA_up_GO.rds")
        }else if(input$select_enrichment_methods=="Reactome"){
          enr<- readRDS("./data/preprocessing/enrichment/DEG_RNA_up_Reactome.rds")
        }
      }else{
        if(input$select_enrichment_methods=="KEGG"){
          enr <- readRDS("./data/preprocessing/enrichment/DEG_RNA_down_kegg.rds")
        }else if(input$select_enrichment_methods=="GO"){
          enr<- readRDS("./data/preprocessing/enrichment/DEG_RNA_down_GO.rds")
        }else if(input$select_enrichment_methods=="Reactome"){
          enr<- readRDS("./data/preprocessing/enrichment/DEG_RNA_down_Reactome.rds")
        }
      }
    }
  )
  
  enrichmentforuse <- function(){
    enr.dat <- enrichmentforuse0()
    enr.dat1 <- data.frame(enr.dat)
    enr.dat1 <- enr.dat1[enr.dat1[["pvalue"]]< input$select_enrichment_pvaluecutoff,]
    enr.dat1 <- enr.dat1[enr.dat1[["qvalue"]]< input$select_enrichment_qvaluecutoff,]
  }

  output$enrichment <- renderPlot({
    if(input$select_enrichment_methods=="KEGG"){
      p <- dotplot(enrichmentforuse0(),color="pvalue",orderBy = "x",showCategory=input$select_enrichment_showCategory)
    }else if(input$select_enrichment_methods=="Reactome"){
      p <- dotplot(enrichmentforuse0(),color="pvalue",orderBy = "x",showCategory=input$select_enrichment_showCategory)
    }else if(input$select_enrichment_methods=="GO"){
      p <- dotplot(all,color="pvalue",orderBy = "x",split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free")
    }
    print(p)
  })
  
  output$Downloadenrichment <- downloadHandler(
    filename <- function() {
      pdf_file <<- as.character(input$fnameenrichment)
      paste(pdf_file,".pdf", sep='')
    },
    content <- function(file) {
      pdf(file, height= 10, width=12,onefile = FALSE)
      if(input$select_enrichment_methods=="KEGG"){
        p <- dotplot(enrichmentforuse0(),color="pvalue",orderBy = "x",showCategory=input$select_enrichment_showCategory)
      }else if(input$select_enrichment_methods=="Reactome"){
        p <- dotplot(enrichmentforuse0(),color="pvalue",orderBy = "x",showCategory=input$select_enrichment_showCategory)
      }else if(input$select_enrichment_methods=="GO"){
        p <- dotplot(all,color="pvalue",orderBy = "x", split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free")
      }
      print(p)
      dev.off()
      #file.copy(paste(pdf_file,'.pdf', sep='') ,file, overwrite=TRUE)
    }, contentType = 'image/pdf')
  
  output$enrichmenttable <- renderDataTable({
    data.frame(enrichmentforuse())
  })
  
  output$Downloadenrichmenttable <- downloadHandler(
    filename <- function() {
      pdf_file <<- as.character(input$fnameenrichment)
      paste0(pdf_file,'.csv')
    },
    content <- function(file) {
      enrichment_table <- data.frame(enrichmentforuse())
      write.csv(enrichment_table,file,row.names = F)
    })
  
  # Classification
  classficationforuse0 <- reactive(
    if(input$select_classification_filetype== 'array'){
      dat.class <- read.csv("./data/preprocessing/DEG/DEG_array.csv",row.names = 1)
    }else if(input$select_classification_filetype== 'RNAseq'){
      dat.class <- read.csv("./data/preprocessing/DEG/DEG_RNA.csv",row.names = 1)
    }
  )
  classficationforuse <- function(){
    dat.class1 <- classficationforuse0()
    dat.class1 <- dat.class1[dat.class1[["P.Value"]]<input$select_classification_pvalue,]
    dat.class1 <- dat.class1[dat.class1[["logFC"]]>log(input$select_classification_fc),]
    RNAmap  <- RNAmap_input()
    data_all <- merge(dat.class1,RNAmap,by.x="GeneSymbol",by.y="gene_name")
    if (input$select_classification_type=="All"){
      return(data_all)
    }else if(input$select_classification_type=="mRNA"){
      data_mRNA <- data_all[data_all$type=="mRNA",]
      return(data_mRNA)
    }else if(input$select_classification_type=="miRNA"){
      data_miRNA <- data_all[data_all$type=="miRNA",]
      return(data_miRNA)
    }else if(input$select_classification_type=="lncRNA"){
      data_lnc <-data_all[data_all$type=="lncRNA",]
      return(data_lnc)
    }
  }
  
  output$classificationtable <- renderDataTable({
    out <- classficationforuse()
    datatable(out,
              options = list(
                "pageLength" =20))
  })
  
  output$Downloadclassificationtable <- downloadHandler(
    filename <- function() {
      pdf_file <<- as.character(input$fnameclassification)
      paste0(pdf_file,'.csv')
    },
    content <- function(file) {
      fnameclassification_table <- classficationforuse()
      write.csv(fnameclassification_table,file,row.names = F)
    })
  
  # correlation
  data_input_cor  <- reactive(
    if(input$select_correlation_filetype == 'array'){
      if(input$select_correlation_algorithm=="PCC"){
        data_input_cor <- read.csv("./data/preprocessing/correlation/DEG.array.pcc.csv")
      }else if(input$select_correlation_algorithm=="SPCC"){
        data_input_cor <- read.csv("./data/preprocessing/correlation/DEG.pcc.data.array.sponge.sppc.csv")
      }else if(input$select_correlation_algorithm=="PC"){
        data_input_cor <- read.csv("./data/preprocessing/correlation/DEG.pcc.data.array.sponge.pc.csv")
      }
    }else if(input$select_correlation_filetype == 'RNAseq') {
      if(input$select_correlation_algorithm=="PCC"){
        data_input_cor <- read.csv("./data/preprocessing/correlation/DEG.pcc.csv")
      }else if(input$select_correlation_algorithm=="SPCC"){
        data_input_cor <- read.csv("./data/preprocessing/correlation/DEG.pcc.data.sponge.sppc.csv")
      }else if(input$select_correlation_algorithm=="PC"){
        data_input_cor <- read.csv("./data/preprocessing/correlation/DEG.pcc.data.sponge.pc.csv")
      }
    }
  )
  
  output$correlationtable <- renderDataTable({
    outcor <- data_input_cor()
    if(input$select_correlation_algorithm=="PCC"){
      outcor <- outcor[outcor$pvalue <  input$select_calculate_cor_pvalue,]
      outcor <- outcor[abs(outcor$cor) >  input$select_calculate_correlation,]
    }else if(input$select_correlation_algorithm=="SPCC"){
      outcor <- outcor[outcor$correlation > input$select_calculate_cor_SPCC,]
      outcor <- outcor[outcor[,3] >= input$select_calculate_cor_share_miRNAs_SPCC,]
      outcor <- outcor[outcor[,6] >= input$select_calculate_cor_padj_SPCC,]
      outcor <- outcor[outcor[,7] >= input$select_calculate_cor_sc_SPCC,]
    }else if(input$select_correlation_algorithm=="PC"){
      outcor <- outcor[outcor$correlation >= input$select_calculate_cor_PC,]
      outcor <- outcor[outcor[,3] >= input$select_calculate_cor_share_miRNAs_PC,]
      outcor <- outcor[outcor[,6] >= input$select_calculate_cor_sc_PC,]
    }
    return(outcor)
  })
  
  data_input_cor_enrichment_foruse <- function(){
    if(input$select_correlation_filetype == 'array'){
      if(input$select_correlation_algorithm=="PCC"){
        data_input_cor1 <- readRDS("./data/preprocessing/correlation/enrichment/DEG.array.enrich.pcc.rds")
      }else if(input$select_correlation_algorithm=="SPCC"){
        data_input_cor1 <- readRDS("./data/preprocessing/correlation/enrichment/DEG.array.enrich.spcc.rds")
      }else if(input$select_correlation_algorithm=="PC"){
        data_input_cor1 <- readRDS("./data/preprocessing/correlation/enrichment/DEG.array.enrich.pc.rds")
      }
    }else if(input$select_correlation_filetype == 'RNAseq') {
      if(input$select_correlation_algorithm=="PCC"){
        data_input_cor1 <- readRDS("./data/preprocessing/correlation/enrichment/DEG.enrich.pcc.rds")
      }else if(input$select_correlation_algorithm=="SPCC"){
        data_input_cor1 <- readRDS("./data/preprocessing/correlation/enrichment/DEG.enrich.spcc.rds")
      }else if(input$select_correlation_algorithm=="PC"){
        data_input_cor1 <- readRDS("./data/preprocessing/correlation/enrichment/DEG.enrich.pc.rds")
      }
    }
    return(data_input_cor1)
  }
  
  output$correlationEnrichment <- renderPlot({
    dotplot(data_input_cor_enrichment_foruse() ,color="pvalue",orderBy = "x",showCategory=5)
  })

  output$Downloadcorrelationtable <- downloadHandler(
    filename <- function() {
      pdf_file <<- as.character(input$fnamecorrelation)
      paste0(pdf_file,'.csv')
    },
    content <- function(file) {
      fnamecorrelation_table <- data_input_cor()
      write.csv(fnamecorrelation_table,file,row.names = F)
    })

  output$Downloadcorenrichment <- downloadHandler(
    filename <- function() {
      pdf_file <<- as.character(input$fnamecorrelation)
      paste(pdf_file,"_enrichment.pdf", sep='')
    },
    content <- function(file) {
      pdf(file, height= 10, width=12,onefile = FALSE)
      dotplot(data_input_cor_enrichment_foruse() ,color="pvalue",orderBy = "x",showCategory=5)
      dev.off()
      #file.copy(paste(pdf_file,'.pdf', sep='') ,file, overwrite=TRUE)
    }, contentType = 'image/pdf')

  output$Downloadecorenrichmenttable <- downloadHandler(
    filename <- function() {
      pdf_file <<- as.character(input$fnamecorrelation)
      paste0(pdf_file,'_enrichment.csv')
    },
    content <- function(file) {
      enrichment_table <- data.frame(data_input_cor_enrichment_foruse())
      write.csv(enrichment_table,file,row.names = F)
    })
  
  # prediction
  predictionforuse0 <- reactive(
    if(input$select_prediction_filetype=="array"){
      mi_mRNA_lnc <- readRDS("./data/preprocessing/prediction/mi_mRNA_lnc_array.rds")
    }else if(input$select_prediction_filetype=="RNAseq"){
      mi_mRNA_lnc <- readRDS("./data/preprocessing/prediction/mi_mRNA_lnc_RNA.rds")
    }
  )
  
  data_input_pre_enrichment_foruse0 <- reactive(
    if(input$select_prediction_filetype=="array"){
      if(input$select_prediction_enrichment_methods=="KEGG"){
        dar.pre.enr <- readRDS("./data/preprocessing/prediction/mi_mRNA_lnc_array_kegg.rds")
      }else if(input$select_prediction_enrichment_methods=="GO"){
        dar.pre.enr <- readRDS("./data/preprocessing/prediction/mi_mRNA_lnc_array_GO.rds")
      }
    }else if(input$select_prediction_filetype=="RNAseq"){
      if(input$select_prediction_enrichment_methods=="KEGG"){
        dar.pre.enr <- readRDS("./data/preprocessing/prediction/mi_mRNA_lnc_RNA_kegg.rds")
      }else if(input$select_prediction_enrichment_methods=="GO"){
        dar.pre.enr <- readRDS("./data/preprocessing/prediction/mi_mRNA_lnc_RNA_GO.rds")
      }
    }
  )
  
  data_input_pre_enrichment_foruse <- function(){
    dar.pre.enr <- data_input_pre_enrichment_foruse0()
    if(input$select_prediction_enrichment_methods=="KEGG"){
      p <- dotplot(dar.pre.enr ,color="pvalue",orderBy = "x",showCategory=5)
    }else if(input$select_prediction_enrichment_methods=="GO"){
      p <- dotplot(dar.pre.enr, color="pvalue",orderBy = "x",split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free")
    }
    print(p)
  }

  output$predictiontable <- renderDataTable({
    outcor <- predictionforuse0()
    if (input$select_prediction_type=="miRNA-lnc"){
      return(outcor[["miRNA-lnc"]])
    }else if(input$select_prediction_type=="miRNA-mRNA"){
      return(outcor[["miRNA-mRNA"]])
    }else if(input$select_prediction_type=="miRNA-lnc-mRNA"){
      return(outcor[["miRNA-lnc-mRNA"]])
    }
  })
  
  output$predictionEnrichment <- renderPlot({
    data_input_pre_enrichment_foruse()
  })
  
  output$Downloadpredictiontable <- downloadHandler(
    filename <- function() {
      pdf_file <<- as.character(input$fnameprediction)
      paste0(pdf_file,'.csv')
    },
    content <- function(file) {
      fnameprediction_table1 <-  predictionforuse0()
      if (input$select_prediction_type=="miRNA-lnc"){
        fnameprediction_table <- fnameprediction_table1[["miRNA-lnc"]]
      }else if(input$select_prediction_type=="miRNA-mRNA"){
        fnameprediction_table <- fnameprediction_table1[["miRNA-mRNA"]]
      }else if(input$select_prediction_type=="miRNA-lnc-mRNA"){
        fnameprediction_table <- fnameprediction_table1[["miRNA-lnc-mRNA"]]
      }
      write.csv(fnameprediction_table,file,row.names = F)
    })
  
  output$Downloadpreenrichment<- downloadHandler(
    filename <- function() {
      pdf_file <<- as.character(input$fnameprediction)
      paste(pdf_file,"_enrichment.pdf", sep='')
    },
    content <- function(file) {
      pdf(file, height= 10, width=12,onefile = FALSE)
      data_input_pre_enrichment_foruse()
      dev.off()
      #file.copy(paste(pdf_file,'.pdf', sep='') ,file, overwrite=TRUE)
    }, contentType = 'image/pdf')
  
  output$Downloadepreenrichmenttable <- downloadHandler(
    filename <- function() {
      pdf_file <<- as.character(input$fnameprediction)
      paste0(pdf_file,'_enrichment.csv')
    },
    content <- function(file) {
      enrichment_table <- data.frame(data_input_pre_enrichment_foruse0())
      enrichment_table <- enrichment_table[enrichment_table[["pvalue"]]<input$select_prediction_enrichment_pvaluecutoff,]
      enrichment_table <- enrichment_table[enrichment_table[["qvalue"]]<input$select_prediction_enrichment_qvaluecutoff,]
      write.csv(enrichment_table,file,row.names = F)
    })

  # network
  networkforuse0 <- reactive(
    if(input$select_ce_filetype=="array"){
      mi_m_lnc <- read.csv("./data/preprocessing/ceRNA/ceRNA.csv")
    }else if(input$select_ce_filetype=="RNAseq"){
      mi_m_lnc <- read.csv("./data/preprocessing/ceRNA/ceRNA_RNA.csv")
    }
  )
  
  networkenrichmentforuse <- function(){
    if(input$select_ce_filetype=="array"){
      if(input$select_network_enrichment_methods=="KEGG"){
        pp <- readRDS("./data/preprocessing/ceRNA/ceRNA_kegg.rds")
      }else if(input$select_network_enrichment_methods=="GO"){
        pp <- readRDS("./data/preprocessing/ceRNA/ceRNA_GO.rds")
      }
    }
  }
  
  networkplot <- function(){
    mi_m_lnc <- networkforuse0()
    if(nrow(mi_m_lnc)==0){
      print("There is no result, please change your CUTOFF")
    }else{
      if (input$select_layout_ce=="fruchterman.reingold"){
        draw_ceRNA(mi_m_lnc, 
                   layout_ce=layout.fruchterman.reingold,
                   vertex.size_ce=input$select_vertex.size_ce,
                   vertex.shape_ce=input$select_vertex.shape_ce,
                   vertex.label.cex_ce=input$select_vertex.label.cex_ce,
                   vertex.label.dist_ce=input$select_vertex.label.dist_ce,
                   edge.width_ce=input$select_edge.width_ce,
                   edge.color_ce="gray",
                   vertex.label.color_ce="black")
      }else if(input$select_layout_ce=="circle"){
        draw_ceRNA(mi_m_lnc, 
                   layout_ce=layout.circle,
                   vertex.size_ce=input$select_vertex.size_ce,
                   vertex.shape_ce=input$select_vertex.shape_ce,
                   vertex.label.cex_ce=input$select_vertex.label.cex_ce,
                   vertex.label.dist_ce=input$select_vertex.label.dist_ce,
                   edge.width_ce=input$select_edge.width_ce,
                   edge.color_ce="gray",
                   vertex.label.color_ce="black")
      }else if(input$select_layout_ce=="reingold.tilford"){
        draw_ceRNA(mi_m_lnc, 
                   layout_ce=layout.reingold.tilford,
                   vertex.size_ce=input$select_vertex.size_ce,
                   vertex.shape_ce=input$select_vertex.shape_ce,
                   vertex.label.cex_ce=input$select_vertex.label.cex_ce,
                   vertex.label.dist_ce=input$select_vertex.label.dist_ce,
                   edge.width_ce=input$select_edge.width_ce,
                   edge.color_ce="gray",
                   vertex.label.color_ce="black")
      }else if(input$select_layout_ce=="star"){
        draw_ceRNA(mi_m_lnc, 
                   layout_ce=layout_as_star,
                   vertex.size_ce=input$select_vertex.size_ce,
                   vertex.shape_ce=input$select_vertex.shape_ce,
                   vertex.label.cex_ce=input$select_vertex.label.cex_ce,
                   vertex.label.dist_ce=input$select_vertex.label.dist_ce,
                   edge.width_ce=input$select_edge.width_ce,
                   edge.color_ce="gray",
                   vertex.label.color_ce="black")
      }else if(input$select_layout_ce=="with_drl"){
        draw_ceRNA(mi_m_lnc, 
                   layout_ce=layout_with_drl,
                   vertex.size_ce=input$select_vertex.size_ce,
                   vertex.shape_ce=input$select_vertex.shape_ce,
                   vertex.label.cex_ce=input$select_vertex.label.cex_ce,
                   vertex.label.dist_ce=input$select_vertex.label.dist_ce,
                   edge.width_ce=input$select_edge.width_ce,
                   edge.color_ce="gray",
                   vertex.label.color_ce="black")
      }else if(input$select_layout_ce=="random"){
        draw_ceRNA(mi_m_lnc, 
                   layout_ce=layout.random,
                   vertex.size_ce=input$select_vertex.size_ce,
                   vertex.shape_ce=input$select_vertex.shape_ce,
                   vertex.label.cex_ce=input$select_vertex.label.cex_ce,
                   vertex.label.dist_ce=input$select_vertex.label.dist_ce,
                   edge.width_ce=input$select_edge.width_ce,
                   edge.color_ce="gray",
                   vertex.label.color_ce="black")
      }
    }
  }

  output$networktext <- renderText({
    mi_m_lnc <- networkforuse0()
    if(nrow(mi_m_lnc)==0){
     info <- "There is no result, please change your CUTOFF"
     info
    }
  })
  
  output$network <- renderPlot({
    networkplot()
  })
  output$networkEnrichment <- renderPlot({
    if(input$select_ce_filetype=="array"){
      pp <- networkenrichmentforuse()
      if(input$select_network_enrichment_methods=="KEGG"){
        pp1 <- dotplot(pp,color="pvalue",orderBy = "x",showCategory=5)
      }else if(input$select_network_enrichment_methods=="GO"){
        pp1 <-dotplot(all,color="pvalue",orderBy = "x", split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free")
      }
      print(pp1)
    }else if(input$select_ce_filetype=="RNAseq"){
      print("There is no result, please change your CUTOFF")
    }
  })
  
  output$Downloadnetwork <- downloadHandler(
    filename <- function() {
      pdf_file <<- as.character(input$fnamenetwork)
      paste(pdf_file,".pdf", sep='')
    },
    content <- function(file) {
      pdf(file, height= 10, width=12,onefile = FALSE)
      networkplot()
      dev.off()
    }, contentType = 'image/pdf')
  
  output$Downloadceenrichment <- downloadHandler(
    filename <- function() {
      pdf_file <<- as.character(input$fnamenetwork)
      paste(pdf_file,"_enrichment.pdf", sep='')
    },
    content <- function(file) {
      pdf(file, height= 10, width=12,onefile = FALSE)
      pp <- networkenrichmentforuse()
      if(input$select_network_enrichment_methods=="KEGG"){
        pp1 <- dotplot(pp,color="pvalue",orderBy = "x",showCategory=5)
      }else if(input$select_network_enrichment_methods=="GO"){
        pp1 <-dotplot(all, color="pvalue",orderBy = "x",split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free")
      }
      print(pp1)
      dev.off()
      #file.copy(paste(pdf_file,'.pdf', sep='') ,file, overwrite=TRUE)
    }, contentType = 'image/pdf')
  
  output$Downloadceenrichmentable <- downloadHandler(
    filename <- function() {
      pdf_file <<- as.character(input$fnamenetwork)
      paste0(pdf_file,'_enrichment.csv')
    },
    content <- function(file) {
      enrichment_table <- data.frame(networkenrichmentforuse())
      enrichment_table <- enrichment_table[enrichment_table[["pvalue"]]<select_network_enrichment_pvaluecutoff,]
      enrichment_table <- enrichment_table[enrichment_table[["qvalue"]]<select_network_enrichment_qvaluecutoff,]
      write.csv(enrichment_table,file,row.names = F)
    })
  
  ## ReadMe
  output$ReadMe1 <- renderUI({
    str00 <- paste("&emsp;")
    str0 <- paste("Example")
    str1 <- paste("&emsp; 1. Input data: GSE7621 or GSE136666.csv")
    str2 <- paste("&emsp; 2. Easily click the buttons from table panel to get corresponding results")
    str3 <- paste("&emsp; 3. NOTICE: Please check that your file types are the same as the ones you clicked on, e.g. all are Array types")

    str21 <- paste("Format of your files")
    str211 <- paste("Array Data")
    str22 <- paste("&emsp; 1. Upload your file: expression matrix, group list and annotation platform information")
    str23 <- paste("&emsp; 2. Row names of the expression matrix are probe names")
    str24 <- paste("&emsp; 3. Column names are sample lists")
    str212 <- paste("High-throughput RNA sequencing data")
    str222 <- paste("&emsp; 1. Upload your file: expression matrix and group list")
    str232 <- paste("&emsp; 2. Row names of the expression matrix are Ensembl numbers")
    str242 <- paste("&emsp; 3. Column names are sample lists")
    
    str31 <- paste("Module")
    str32 <- paste("&emsp; 1. Heatmap plot &emsp;")
    str33 <- paste("&emsp; 2. Volcanol plot &emsp;")
    str34 <- paste("&emsp; 3. Enrichment &emsp;")
    str35 <- paste("&emsp; 4. RNA classification &emsp;")
    str36 <- paste("&emsp; 5. ceRNA network1: ceRNA network based on PCC, SPCC and PC &emsp;")
    str37 <- paste("&emsp; 6. ceRNA network2: ceRNA network based on databases &emsp;")
    str38 <- paste("&emsp; 7. ceRNA network3: ceRNA network based on PCC and databases &emsp;")

    HTML(paste(str00,h4(strong(str0)), str1, str2, str3,
               str00,h4(strong(str21)),strong(str211),str22,str23,str24,str00,strong(str212),str222,str232,str242,
               str00,h4(strong(str31)),str32,str33,str34,str35,str36,str37,str38,
               sep = '<br/>'))
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
