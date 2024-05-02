#library(ISLR)
library(ggplot2)
# library(plotly)
library(stringr)
library(stringi)
library(shinythemes)
library(tensorflow)
library(keras)
#library(corrplot)
library(tidyr)
#library(psych)
library(dplyr)
library(DT)

library("protr")
library("bio3d")
source("generateFeatures.R")

library(shiny)
library(shinyjs)

library(seqinr)
#library(Biostrings)
library(phylotools)
#library(shinyFeedback)

css <- "
.busy { 
  position: fixed;
  z-index: 1000;
  top: 50%;
  left: 50%;
  margin-top: -100px;
  margin-left: -50px;
  display: none;
  background-color: rgba(230,230,230,.8);
  text-align: center;
  padding-top: 20px;
  padding-left: 30px;
  padding-bottom: 40px;
  padding-right: 30px;
  border-radius: 5px;
}"

amp.imp.features <- c("ARGP820103.hydrophobicity.nf.g1","ARGP820103.hydrophobicity.nf.g2","ARGP820103.hydrophobicity.nf.g3","ARGP820103.vanderwa.nf.g1","ARGP820103.vanderwa.nf.g2","ARGP820103.vanderwa.nf.g3","ARGP820103.charge.nf.g1","ARGP820103.charge.nf.g2","ARGP820103.charge.nf.g3","ARGP820103.polarity.nf.g1","ARGP820103.polarity.nf.g2","ARGP820103.polarity.nf.g3","ARGP820103.polarizability.nf.g1","ARGP820103.polarizability.nf.g2","ARGP820103.polarizability.nf.g3","ARGP820103.second.struct.nf.g1","ARGP820103.second.struct.nf.g2","ARGP820103.second.struct.nf.g3","ARGP820103.solv.acc.nf.g1","ARGP820103.solv.acc.nf.g2","ARGP820103.solv.acc.nf.g3","CHAM830105.hydrophobicity.nf.g1","CHAM830105.hydrophobicity.nf.g2","CHAM830105.hydrophobicity.nf.g3","CHAM830105.vanderwa.nf.g1","CHAM830105.vanderwa.nf.g2","CHAM830105.vanderwa.nf.g3","CHAM830105.charge.nf.g1","CHAM830105.charge.nf.g2","CHAM830105.charge.nf.g3","CHAM830105.polarity.nf.g1","CHAM830105.polarity.nf.g2","CHAM830105.polarity.nf.g3","CHAM830105.polarizability.nf.g1","CHAM830105.polarizability.nf.g2","CHAM830105.polarizability.nf.g3","CHAM830105.second.struct.nf.g1","CHAM830105.second.struct.nf.g2","CHAM830105.second.struct.nf.g3","CHAM830105.solv.acc.nf.g1","CHAM830105.solv.acc.nf.g2","CHAM830105.solv.acc.nf.g3","SUEM840102.hydrophobicity.nf.g1","SUEM840102.hydrophobicity.nf.g2","SUEM840102.hydrophobicity.nf.g3","SUEM840102.vanderwa.nf.g1","SUEM840102.vanderwa.nf.g2","SUEM840102.vanderwa.nf.g3","SUEM840102.charge.nf.g1","SUEM840102.charge.nf.g2","SUEM840102.charge.nf.g3","SUEM840102.polarity.nf.g1","SUEM840102.polarity.nf.g2","SUEM840102.polarity.nf.g3","SUEM840102.polarizability.nf.g1","SUEM840102.polarizability.nf.g2","SUEM840102.polarizability.nf.g3","SUEM840102.second.struct.nf.g1","SUEM840102.second.struct.nf.g2","SUEM840102.second.struct.nf.g3","SUEM840102.solv.acc.nf.g1","SUEM840102.solv.acc.nf.g2","SUEM840102.solv.acc.nf.g3","HUTJ700101.hydrophobicity.nf.g1","HUTJ700101.hydrophobicity.nf.g2","HUTJ700101.hydrophobicity.nf.g3","HUTJ700101.vanderwa.nf.g1","HUTJ700101.vanderwa.nf.g2","HUTJ700101.vanderwa.nf.g3","HUTJ700101.charge.nf.g1","HUTJ700101.charge.nf.g2","HUTJ700101.charge.nf.g3","HUTJ700101.polarity.nf.g1","HUTJ700101.polarity.nf.g2","HUTJ700101.polarity.nf.g3","HUTJ700101.polarizability.nf.g1","HUTJ700101.polarizability.nf.g2","HUTJ700101.polarizability.nf.g3","HUTJ700101.second.struct.nf.g1","HUTJ700101.second.struct.nf.g2","HUTJ700101.second.struct.nf.g3","HUTJ700101.solv.acc.nf.g1","HUTJ700101.solv.acc.nf.g2","HUTJ700101.solv.acc.nf.g3","GEIM800104.hydrophobicity.nf.g1","GEIM800104.hydrophobicity.nf.g2","GEIM800104.hydrophobicity.nf.g3","GEIM800104.vanderwa.nf.g1","GEIM800104.vanderwa.nf.g2","GEIM800104.vanderwa.nf.g3","GEIM800104.charge.nf.g1","GEIM800104.charge.nf.g2","GEIM800104.charge.nf.g3","GEIM800104.polarity.nf.g1","GEIM800104.polarity.nf.g2","GEIM800104.polarity.nf.g3","GEIM800104.polarizability.nf.g1","GEIM800104.polarizability.nf.g2","GEIM800104.polarizability.nf.g3","GEIM800104.second.struct.nf.g1","GEIM800104.second.struct.nf.g2","GEIM800104.second.struct.nf.g3","GEIM800104.solv.acc.nf.g1","GEIM800104.solv.acc.nf.g2","GEIM800104.solv.acc.nf.g3")
abp.imp.features <- c("SUEM840102.hydrophobicity.nf.g1","SUEM840102.hydrophobicity.nf.g2","SUEM840102.hydrophobicity.nf.g3","SUEM840102.vanderwa.nf.g1","SUEM840102.vanderwa.nf.g2","SUEM840102.vanderwa.nf.g3","SUEM840102.charge.nf.g1","SUEM840102.charge.nf.g2","SUEM840102.charge.nf.g3","SUEM840102.polarity.nf.g1","SUEM840102.polarity.nf.g2","SUEM840102.polarity.nf.g3","SUEM840102.polarizability.nf.g1","SUEM840102.polarizability.nf.g2","SUEM840102.polarizability.nf.g3","SUEM840102.second.struct.nf.g1","SUEM840102.second.struct.nf.g2","SUEM840102.second.struct.nf.g3","SUEM840102.solv.acc.nf.g1","SUEM840102.solv.acc.nf.g2","SUEM840102.solv.acc.nf.g3","QIAN880128.hydrophobicity.nf.g1","QIAN880128.hydrophobicity.nf.g2","QIAN880128.hydrophobicity.nf.g3","QIAN880128.vanderwa.nf.g1","QIAN880128.vanderwa.nf.g2","QIAN880128.vanderwa.nf.g3","QIAN880128.charge.nf.g1","QIAN880128.charge.nf.g2","QIAN880128.charge.nf.g3","QIAN880128.polarity.nf.g1","QIAN880128.polarity.nf.g2","QIAN880128.polarity.nf.g3","QIAN880128.polarizability.nf.g1","QIAN880128.polarizability.nf.g2","QIAN880128.polarizability.nf.g3","QIAN880128.second.struct.nf.g1","QIAN880128.second.struct.nf.g2","QIAN880128.second.struct.nf.g3","QIAN880128.solv.acc.nf.g1","QIAN880128.solv.acc.nf.g2","QIAN880128.solv.acc.nf.g3","HUTJ700101.hydrophobicity.nf.g1","HUTJ700101.hydrophobicity.nf.g2","HUTJ700101.hydrophobicity.nf.g3","HUTJ700101.vanderwa.nf.g1","HUTJ700101.vanderwa.nf.g2","HUTJ700101.vanderwa.nf.g3","HUTJ700101.charge.nf.g1","HUTJ700101.charge.nf.g2","HUTJ700101.charge.nf.g3","HUTJ700101.polarity.nf.g1","HUTJ700101.polarity.nf.g2","HUTJ700101.polarity.nf.g3","HUTJ700101.polarizability.nf.g1","HUTJ700101.polarizability.nf.g2","HUTJ700101.polarizability.nf.g3","HUTJ700101.second.struct.nf.g1","HUTJ700101.second.struct.nf.g2","HUTJ700101.second.struct.nf.g3","HUTJ700101.solv.acc.nf.g1","HUTJ700101.solv.acc.nf.g2","HUTJ700101.solv.acc.nf.g3","GEIM800104.hydrophobicity.nf.g1","GEIM800104.hydrophobicity.nf.g2","GEIM800104.hydrophobicity.nf.g3","GEIM800104.vanderwa.nf.g1","GEIM800104.vanderwa.nf.g2","GEIM800104.vanderwa.nf.g3","GEIM800104.charge.nf.g1","GEIM800104.charge.nf.g2","GEIM800104.charge.nf.g3","GEIM800104.polarity.nf.g1","GEIM800104.polarity.nf.g2","GEIM800104.polarity.nf.g3","GEIM800104.polarizability.nf.g1","GEIM800104.polarizability.nf.g2","GEIM800104.polarizability.nf.g3","GEIM800104.second.struct.nf.g1","GEIM800104.second.struct.nf.g2","GEIM800104.second.struct.nf.g3","GEIM800104.solv.acc.nf.g1","GEIM800104.solv.acc.nf.g2","GEIM800104.solv.acc.nf.g3")    
avp.imp.features <- c("ARGP820102.hydrophobicity.nf.g1","ARGP820102.hydrophobicity.nf.g2","ARGP820102.hydrophobicity.nf.g3","ARGP820102.vanderwa.nf.g1","ARGP820102.vanderwa.nf.g2","ARGP820102.vanderwa.nf.g3","ARGP820102.charge.nf.g1","ARGP820102.charge.nf.g2","ARGP820102.charge.nf.g3","ARGP820102.polarity.nf.g1","ARGP820102.polarity.nf.g2","ARGP820102.polarity.nf.g3","ARGP820102.polarizability.nf.g1","ARGP820102.polarizability.nf.g2","ARGP820102.polarizability.nf.g3","ARGP820102.second.struct.nf.g1","ARGP820102.second.struct.nf.g2","ARGP820102.second.struct.nf.g3","ARGP820102.solv.acc.nf.g1","ARGP820102.solv.acc.nf.g2","ARGP820102.solv.acc.nf.g3","CHAM830105.hydrophobicity.nf.g1","CHAM830105.hydrophobicity.nf.g2","CHAM830105.hydrophobicity.nf.g3","CHAM830105.vanderwa.nf.g1","CHAM830105.vanderwa.nf.g2","CHAM830105.vanderwa.nf.g3","CHAM830105.charge.nf.g1","CHAM830105.charge.nf.g2","CHAM830105.charge.nf.g3","CHAM830105.polarity.nf.g1","CHAM830105.polarity.nf.g2","CHAM830105.polarity.nf.g3","CHAM830105.polarizability.nf.g1","CHAM830105.polarizability.nf.g2","CHAM830105.polarizability.nf.g3","CHAM830105.second.struct.nf.g1","CHAM830105.second.struct.nf.g2","CHAM830105.second.struct.nf.g3","CHAM830105.solv.acc.nf.g1","CHAM830105.solv.acc.nf.g2","CHAM830105.solv.acc.nf.g3","HUTJ700101.hydrophobicity.nf.g1","HUTJ700101.hydrophobicity.nf.g2","HUTJ700101.hydrophobicity.nf.g3","HUTJ700101.vanderwa.nf.g1","HUTJ700101.vanderwa.nf.g2","HUTJ700101.vanderwa.nf.g3","HUTJ700101.charge.nf.g1","HUTJ700101.charge.nf.g2","HUTJ700101.charge.nf.g3","HUTJ700101.polarity.nf.g1","HUTJ700101.polarity.nf.g2","HUTJ700101.polarity.nf.g3","HUTJ700101.polarizability.nf.g1","HUTJ700101.polarizability.nf.g2","HUTJ700101.polarizability.nf.g3","HUTJ700101.second.struct.nf.g1","HUTJ700101.second.struct.nf.g2","HUTJ700101.second.struct.nf.g3","HUTJ700101.solv.acc.nf.g1","HUTJ700101.solv.acc.nf.g2","HUTJ700101.solv.acc.nf.g3")
afp.imp.features <- c("COSI940101.hydrophobicity.nf.g1","COSI940101.hydrophobicity.nf.g2","COSI940101.hydrophobicity.nf.g3","COSI940101.vanderwa.nf.g1","COSI940101.vanderwa.nf.g2","COSI940101.vanderwa.nf.g3","COSI940101.charge.nf.g1","COSI940101.charge.nf.g2","COSI940101.charge.nf.g3","COSI940101.polarity.nf.g1","COSI940101.polarity.nf.g2","COSI940101.polarity.nf.g3","COSI940101.polarizability.nf.g1","COSI940101.polarizability.nf.g2","COSI940101.polarizability.nf.g3","COSI940101.second.struct.nf.g1","COSI940101.second.struct.nf.g2","COSI940101.second.struct.nf.g3","COSI940101.solv.acc.nf.g1","COSI940101.solv.acc.nf.g2","COSI940101.solv.acc.nf.g3","QIAN880128.hydrophobicity.nf.g1","QIAN880128.hydrophobicity.nf.g2","QIAN880128.hydrophobicity.nf.g3","QIAN880128.vanderwa.nf.g1","QIAN880128.vanderwa.nf.g2","QIAN880128.vanderwa.nf.g3","QIAN880128.charge.nf.g1","QIAN880128.charge.nf.g2","QIAN880128.charge.nf.g3","QIAN880128.polarity.nf.g1","QIAN880128.polarity.nf.g2","QIAN880128.polarity.nf.g3","QIAN880128.polarizability.nf.g1","QIAN880128.polarizability.nf.g2","QIAN880128.polarizability.nf.g3","QIAN880128.second.struct.nf.g1","QIAN880128.second.struct.nf.g2","QIAN880128.second.struct.nf.g3","QIAN880128.solv.acc.nf.g1","QIAN880128.solv.acc.nf.g2","QIAN880128.solv.acc.nf.g3","PALJ810109.hydrophobicity.nf.g1","PALJ810109.hydrophobicity.nf.g2","PALJ810109.hydrophobicity.nf.g3","PALJ810109.vanderwa.nf.g1","PALJ810109.vanderwa.nf.g2","PALJ810109.vanderwa.nf.g3","PALJ810109.charge.nf.g1","PALJ810109.charge.nf.g2","PALJ810109.charge.nf.g3","PALJ810109.polarity.nf.g1","PALJ810109.polarity.nf.g2","PALJ810109.polarity.nf.g3","PALJ810109.polarizability.nf.g1","PALJ810109.polarizability.nf.g2","PALJ810109.polarizability.nf.g3","PALJ810109.second.struct.nf.g1","PALJ810109.second.struct.nf.g2","PALJ810109.second.struct.nf.g3","PALJ810109.solv.acc.nf.g1","PALJ810109.solv.acc.nf.g2","PALJ810109.solv.acc.nf.g3","HUTJ700101.hydrophobicity.nf.g1","HUTJ700101.hydrophobicity.nf.g2","HUTJ700101.hydrophobicity.nf.g3","HUTJ700101.vanderwa.nf.g1","HUTJ700101.vanderwa.nf.g2","HUTJ700101.vanderwa.nf.g3","HUTJ700101.charge.nf.g1","HUTJ700101.charge.nf.g2","HUTJ700101.charge.nf.g3","HUTJ700101.polarity.nf.g1","HUTJ700101.polarity.nf.g2","HUTJ700101.polarity.nf.g3","HUTJ700101.polarizability.nf.g1","HUTJ700101.polarizability.nf.g2","HUTJ700101.polarizability.nf.g3","HUTJ700101.second.struct.nf.g1","HUTJ700101.second.struct.nf.g2","HUTJ700101.second.struct.nf.g3","HUTJ700101.solv.acc.nf.g1","HUTJ700101.solv.acc.nf.g2","HUTJ700101.solv.acc.nf.g3")
amp.model <- load_model_hdf5('amp-best-apq-index.h5')
abp.model <- load_model_hdf5('abp-best-apq-index.h5')
avp.model <- load_model_hdf5('avp-best-apq-index.h5')
afp.model <- load_model_hdf5('afp-best-apq-index.h5')

impIndexGroup <- c("HUTJ700101","SUEM840102","GEIM800104","CHAM830105","ARGP820103", "QIAN880128","ARGP820102","COSI940101","PALJ810109")
ampImpIndexGroup <- c("HUTJ700101","SUEM840102","GEIM800104","CHAM830105","ARGP820103")
abpImpIndexGroup <- c("HUTJ700101","GEIM800104","QIAN880128","SUEM840102")
avpImpIndexGroup <- c("HUTJ700101","ARGP820102","CHAM830105")
afpImpIndexGroup <- c("HUTJ700101","QIAN880128", "COSI940101", "PALJ810109")

ui <- tagList(
    useShinyjs(),
    tags$head(tags$style(css)),
    tags$div(class = "busy", 
             tags$img(src = "https://loading.io/spinners/comets/lg.comet-spinner.gif")),
    shinythemes::themeSelector(),
    navbarPage(
    
        "AMPs Recognition App",
        tabPanel("Type", 
                 mainPanel(
                     textAreaInput("ampSeq", "Input sequences in FASTA format (length between 5 and 50) (Example)", "", width = "1300px", height = "170px"),
                     p("Example:"),
                     p(">SEQ1",tags$br(), "IRGYKGGYCKGAFKQTCKCY",tags$br(),">SEQ2",tags$br(),"YGAMMMKKKDDD"),
                     actionButton("predictBtn", "Predict", class = "btn-primary"),
                     actionButton("clearBtn", "Clear", class = "btn-danger"),
                     br(),
                     br(),
                     span(textOutput("response_msg"), style="color:red"),
                     # downloadButton("download1", "Download results in csv format"),
                     tabsetPanel(
                         tabPanel("Classification",
                                  DT::dataTableOutput("rsTable")
                         ),
                         tabPanel("Score",
                                  DT::dataTableOutput("scoreTable")
                         )
                     ),
                     shinycssloaders::withSpinner(
                         plotOutput("plot")
                     )
                         
                 )
        ),
        tabPanel("Upload", 
                 mainPanel(
                   fileInput("file1", "Choose Fasta File",
                             multiple = FALSE,
                             accept = c("text/fasta",
                                        ".fasta"),
                             ),
                   actionButton("predictBtn2", "Predict", class = "btn-primary"),
                   actionButton("clearBtn2", "Clear", class = "btn-danger"),
                   br(),
                   br(),
                   span(textOutput("response_msg2"), style="color:red"),
                   tabsetPanel(
                     tabPanel("Classification",
                              DT::dataTableOutput("rsTable2")
                     ),
                     tabPanel("Score",
                              DT::dataTableOutput("scoreTable2")
                     )
                   ),
                   shinycssloaders::withSpinner(
                     plotOutput("plot2")
                   )
                 )
        )#,
        # tabPanel("Generate AAIW", 
        #          mainPanel(
        #            fileInput("file2", "Choose Fasta File",
        #                      multiple = FALSE,
        #                      accept = c("text/fasta",
        #                                 ".fasta"),
        #            ),
        #            actionButton("generateBtn", "Generate", class = "btn-primary"),
        #            actionButton("clearGenBtn", "Clear", class = "btn-danger"),
        #            br(),
        #            br(),
        #            span(textOutput("response_msg3"), style="color:red"),
        #            downloadButton("downloadData", "Download"),
        #            tabsetPanel(
        #              tabPanel("Classification",
        #                       DT::dataTableOutput("rsTable3")
        #              )#,
        #              # tabPanel("Score",
        #              #          DT::dataTableOutput("scoreTable3")
        #              # )
        #            ),
        #            shinycssloaders::withSpinner(
        #              plotOutput("plot3")
        #            )
        #          )
        # )
        
    )
)

server <- function(input, output) {
    is.odd <- function(x) x %% 2 != 0
    df.amp.rs1 <- data.frame()
    df.amp.rs2 <- data.frame()
    observeEvent(input$predictBtn, {
        if(input$ampSeq != ""){
          word_splits <- strsplit(input$ampSeq, "\n")
          word_splits <- word_splits[[1]]
          seq_id <- c()
          seq_name <- c()
          valid_format <- 1
          for (i in 1: length(word_splits)){
            if (is.odd(i)){
              id.cl <- str_trim(word_splits[i]) 
              if(substr(id.cl,1,1)!='>'){
                valid_format <- 0
                break
              }
              seq_id <- append(seq_id, substr(id.cl, 2, str_length(id.cl)))
            }
            else
              seq_name <- append(seq_name,str_trim(word_splits[i]))
          }
          
          df.amp <- data.frame(seq_id=seq_id,seq_name=seq_name)
          df.amp$seq_length <- str_length(df.amp$seq_name)
          
          if(valid_format==1){
            df.amp.rs1 <- predictDf(df.amp)
            if(!is.null(df.amp.rs1)){
              output$rsTable = DT::renderDataTable({
                df.amp.rs1 %>% select(seq_id, seq_name, amp, abp, avp, afp)
                
              })
              output$scoreTable = DT::renderDataTable({
                df.amp.rs1 %>% select(seq_id, seq_name, ampscore,abpscore,avpscore,afpscore)
              })
              
              output$response_msg <- renderText("")
            }else{
              output$response_msg <- renderText("Error: Amino acid invalid")
              output$rsTable = DT::renderDataTable(NULL)
              output$scoreTable = DT::renderDataTable(NULL)
              message("There was an error message.")
            }
          }else{
            output$response_msg <- renderText("Please input fasta format")
            output$rsTable = DT::renderDataTable(NULL)
            output$scoreTable = DT::renderDataTable(NULL)
          }
        }else{
          output$response_msg <- renderText("Please input sequence to predict")
        }
        
    })
    
    output$plot <- renderPlot({
        input$predictBtn
        Sys.sleep(1.3)
    })
    
    output$plot2 <- renderPlot({
      input$predictBtn2
      Sys.sleep(1.3)
    })
    
    output$plot3 <- renderPlot({
      input$generateBtn
      Sys.sleep(1.3)
    })
    
    observeEvent(input$clearBtn, {
        reset("ampSeq")
        output$response_msg <- renderText("")
        output$rsTable = DT::renderDataTable(NULL)
        output$scoreTable = DT::renderDataTable(NULL)
    })
    
    observeEvent(input$clearBtn2, {
      reset("file1")
      output$response_msg <- renderText("")
      output$rsTable2 = DT::renderDataTable(NULL)
      output$scoreTable2 = DT::renderDataTable(NULL)
      
    })
    
    predictDf <- function(df.amp){
      tryCatch(                      
        expr = {                      
            amp.indtest.apq <- generate_features(df.amp)
            origin.features <- colnames(amp.indtest.apq)
            amp.indtest.index <- generate_impFeatures(df.amp, impIndexGroup)
            amp.indtest <- cbind(amp.indtest.apq, amp.indtest.index)
            
            df.amp.encode <- amp.indtest %>% select(origin.features, amp.imp.features)
            df.amp.mat <- df.amp.encode %>% as.matrix
            amp.score <- round(predict(amp.model, df.amp.mat),2)
            amp.predict <- factor(round(amp.score))
            df.amp$ampscore <- amp.score
            df.amp$amp <- amp.predict
            
            df.abp.encode <- amp.indtest %>% select(origin.features, abp.imp.features)
            df.abp.mat <- df.abp.encode %>% as.matrix
            abp.score <- round(predict(abp.model, df.abp.mat),2)
            abp.predict <- factor(round(abp.score))
            df.amp$abpscore <- abp.score
            df.amp$abp <- abp.predict
            
            df.avp.encode <- amp.indtest %>% select(origin.features, avp.imp.features)
            df.avp.mat <- df.avp.encode %>% as.matrix
            avp.score <- round(predict(avp.model, df.avp.mat),2)
            avp.predict <- factor(round(avp.score))
            df.amp$avpscore <- avp.score
            df.amp$avp <- avp.predict
            
            df.afp.encode <- amp.indtest %>% select(origin.features, afp.imp.features)
            df.afp.mat <- df.afp.encode %>% as.matrix
            afp.score <- round(predict(afp.model, df.afp.mat),2)
            afp.predict <- factor(round(afp.score))
            df.amp$afpscore <- afp.score
            df.amp$afp <- afp.predict
        },
        
        error = function(e){
          df.amp <- NULL
        },
        warning = function(w){
          
          # output$response_msg <- renderText("Amino acid invalid")
          # message("There was a warning message.")
          df.amp <- NULL
        }
      )
      df.amp
    }
    generateAAIWFeature <- function(df.amp){
      df.rs <- data.frame()
      tryCatch(                      
        expr = {                      
          amp.indtest.index <- generate_impFeatures(df.amp, impIndexGroup)
          df.rs <- amp.indtest.index
          # df.amp.encode <- amp.indtest %>% select(origin.features, amp.imp.features)
          # df.amp.mat <- df.amp.encode %>% as.matrix
          # amp.score <- round(predict(amp.model, df.amp.mat),2)
          # amp.predict <- factor(round(amp.score))
          # df.amp$ampscore <- amp.score
          # df.amp$amp <- amp.predict
          # 
          # df.abp.encode <- amp.indtest %>% select(origin.features, abp.imp.features)
          # df.abp.mat <- df.abp.encode %>% as.matrix
          # abp.score <- round(predict(abp.model, df.abp.mat),2)
          # abp.predict <- factor(round(abp.score))
          # df.amp$abpscore <- abp.score
          # df.amp$abp <- abp.predict
          # 
          # df.avp.encode <- amp.indtest %>% select(origin.features, avp.imp.features)
          # df.avp.mat <- df.avp.encode %>% as.matrix
          # avp.score <- round(predict(avp.model, df.avp.mat),2)
          # avp.predict <- factor(round(avp.score))
          # df.amp$avpscore <- avp.score
          # df.amp$avp <- avp.predict
          # 
          # df.afp.encode <- amp.indtest %>% select(origin.features, afp.imp.features)
          # df.afp.mat <- df.afp.encode %>% as.matrix
          # afp.score <- round(predict(afp.model, df.afp.mat),2)
          # afp.predict <- factor(round(afp.score))
          # df.amp$afpscore <- afp.score
          # df.amp$afp <- afp.predict
        },
        
        error = function(e){
          df.rs <- NULL
        },
        warning = function(w){
          
          # output$response_msg <- renderText("Amino acid invalid")
          # message("There was a warning message.")
          df.rs <- NULL
        }
      )
      df.rs
    }
    ############## for file upload ############
    observeEvent(input$predictBtn2, {
      if(!is.null(input$file1)){
        df2 <- read.fasta(input$file1$datapath) #colname: seq.name, seq.text
        names(df2)[1]<- 'seq_id'
        names(df2)[2]<- 'seq_name'
        df2$seq_length <- str_length(df2$seq_name)
        
        df.amp.rs <- predictDf(df2)
        if(!is.null(df.amp.rs)){
          output$rsTable2 = DT::renderDataTable({
            df.amp.rs %>% select(seq_id, seq_name, amp, abp, avp, afp)
          })
          output$scoreTable2 = DT::renderDataTable({
            df.amp.rs %>% select(seq_id, seq_name, ampscore,abpscore,avpscore,afpscore)
          })
          output$response_msg2 <- renderText("")
        }else{
          output$response_msg2 <- renderText("Error: Amino acid invalid")
          output$rsTable2 = DT::renderDataTable(NULL)
          output$scoreTable2 = DT::renderDataTable(NULL)
          message("There was an error message.")
        }
      }else{
        output$response_msg2 <- renderText("Please insert fasta file")
      }
      # datatt <- mtcars
      output$download <- downloadHandler(
          filename = function () {
            paste("PredictedAMPs.csv", sep = "")
            # paste("data-", Sys.Date(), ".csv", sep="")
          },
          
          content = function(file) {
            write.csv(datatt, file)
          }
        )
    })
  ################ Generate AAIW features ######################
    observeEvent(input$generateBtn, {
      if(!is.null(input$file2)){
        df2 <- read.fasta(input$file2$datapath) #colname: seq.name, seq.text
        names(df2)[1]<- 'seq_id'
        names(df2)[2]<- 'seq_name'
        df2$seq_length <- str_length(df2$seq_name)
        
        df.amp.rs <- generateAAIWFeature(df2)
        if(!is.null(df.amp.rs)){
          output$rsTable3 = DT::renderDataTable({
            df.amp.rs 
          })
          # output$scoreTable2 = DT::renderDataTable({
          #   df.amp.rs %>% select(seq_id, seq_name, ampscore,abpscore,avpscore,afpscore)
          # })
          output$response_msg3 <- renderText("")
        }else{
          output$response_msg3 <- renderText("Error: Amino acid invalid")
          output$rsTable3 = DT::renderDataTable(NULL)
          # output$scoreTable2 = DT::renderDataTable(NULL)
          message("There was an error message.")
        }
      }else{
        output$response_msg3 <- renderText("Please insert fasta file")
      }
      # datatt <- mtcars
      output$download <- downloadHandler(
        filename = function () {
          paste("PredictedAMPs.csv", sep = "")
          # paste("data-", Sys.Date(), ".csv", sep="")
        },
        
        content = function(file) {
          write.csv(datatt, file)
        }
      )
      
      # Our dataset
      data <- df.amp.rs#mtcars
      
      output$downloadData <- downloadHandler(
        filename = function() {
          paste("data-", Sys.Date(), ".csv", sep="")
        },
        content = function(file) {
          write.csv(data, file)
        }
      )
    })
}

# Run the application 
shinyApp(ui = ui, server = server)