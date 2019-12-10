library(shiny)
# Define UI for application that draws a histogram
library(shiny)
library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(plotly)
library(shinythemes)
#============================================================================
#============================================================================
## TITLE PANEL ##
#============================================================================
#============================================================================
ui <- fluidPage(shinytheme("united"),
    titlePanel(# div(#column(width = 1, h1("scMelon"),h3("Seurat Explorer")), 
        h4("BlockChain Dashboard",align="center"),
        h5("Moffitt Lab - Stony Brook University",align="center")
    
    ),
    fluidRow(
        # tags$img(src = "scmelon.jpg"),
        column(3, 
               wellPanel(
                #   fileInput(inputId="file_name", label="Select Saved Seurat Object File"),
                   # verbatimTextOutput("console1"),
                 #  radioButtons(inputId = "cellIdent", label="Set Ident", choices = c("CellType","Phase","Patient","Condition","seurat_clusters"), selected="CellType")
                   textInput(inputId = "Market", label="Enter Market Name", value="Binance"),
                   textInput(inputId = "Crypto", label="Enter Crypto Name", value="ENJ")
                   #      tags$hr(),
               #    sliderInput(inputId = "dotSize", label = "Pointsize", value=1, min=0, max=3),
                #   sliderInput(inputId = "labelSize", label = "Label size", value=4, min=0.5, max=10)
               )),
        column(8,offset=0,
               # verbatimTextOutput("console1"),
               tabsetPanel(
                   tabPanel("Crypto PCA", plotlyOutput("UMAP")),
                   tabPanel("BiPlot", plotlyOutput("biplot"))
                 #  tabPanel("tSNE", plotOutput("featureTSNE")),
                #   tabPanel("PCA", plotOutput("pcaPlot")),
                 #  tabPanel("Violin", plotOutput("vlnPlot")
                            ,type="pills")
               )),
        column(3,offset=0,
               verbatimTextOutput("console1"),
               wellPanel(
                   textInput(inputId = "Crypto Name", label="Enter Crypto Name", value="ENJ"),
                   actionButton("goButton", "Analyze"),
                   radioButtons(inputId = "Color", label="Color By", choices = c("Red","Blue"), selected="Red")
               ))
    )

#============================================================================
#============================================================================

#wellPanel(
#textInput(inputId = "geneName", label="Enter Gene Name", value="EPCAM"),
# actionButton("goButton", "Query Gene")
#)


#============================================================================
#============================================================================
#============================================================================
## SERVER ##
#============================================================================
#============================================================================

server <- function(input, output) {
    market_today <- get_crypto_listings()
    # increase the max upload file size
    options(shiny.maxRequestSize=1000*1024^2)
    df1 <- na.omit(market_today[,c('symbol','market_cap_usd', 'percent_change_24h', 'price_btc', 'percent_change_7d', 'total_supply','X24h_volume_usd', 'available_supply')])
    
    #as numeric
    df1$market_cap_usd <- as.numeric(df1$market_cap_usd)
    df1$percent_change_24h <- as.numeric(df1$percent_change_24h)
    df1$price_btc <- as.numeric(df1$price_btc)
    df1$percent_change_7d <- as.numeric(df1$percent_change_7d)
    df1$total_supply <- as.numeric(df1$total_supply)
    df1$X24h_volume_usd <- as.numeric(df1$X24h_volume_usd)
    df1$available_supply <- as.numeric(df1$available_supply)
    
    #df1$formatted_market_cap <-  paste0(df1$id,'\n','$',format(df1$market_cap_usd,big.mark = ',',scientific = F, trim = T))
    
    cryptolist <- df1[1:50, 2:8]
    
    res.pca <- prcomp(cryptolist, scale = TRUE)
   # fviz_eig(res.pca)
    # set demo dataset file here
    #demo = "powersintegrated.RData"
    demo ="pengN3C3B3_Powers_integrated.RData"
    # demo = "~/Documents/GitHub/DCA/Baron_DCA.RData"
    #============================================================================
    SeuratObject <- reactive({
        if(is.null(input$file_name)){
            # Load a default demo dataset
            load(file=demo)
        } else {
            # Load user-defined dataset
            load(input$file_name$datapath)
            #  DefaultAssay(SeuratObject())<-"RNA"
        }
    })
    #============================================================================
    output$contents <- renderText({
        if (!is.null(input$file_name)){
            "Currently displaying 8.5k Pancreas dataset as demo"
        } else {
            load(file=inFile$name)
        }
    })
    #============================================================================  
    output$console1 <- renderPrint({
        if (is.null(input$file_name)){
            get(load(demo))
        } else {
            # should be a more elegant way of doing this, but this works...
            get(load(input$file_name$name))
        }
    })
    #============================================================================ 
    # Display TSNE Plot as overview toggled options for display legends/labels and setting size. Default displays demo dataset
    output$UMAP <- renderPlotly({
        fviz_pca_ind(res.pca,
                     col.ind = "cos2", # Color by the quality of representation
                     gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                #     repel = TRUE,     # Avoid text overlapping
                    label="ind",
                      axes = c(1,2)
        )
        # DefaultAssay(SeuratObject())<-"RNA"
        #    df<-get(SeuratObject())
        # Idents(df)<-"CellType"
        #    df <- subset(get(SeuratObject()), downsample = 1000)
      #  DimPlot(subset(get(SeuratObject()), downsample = 500), group.by=input$cellIdent, 
      #          do.label = input$labelBoolean, 
     #           pt.size = input$dotSize, 
     ##           label.size = input$labelSize, 
     #           no.legend=input$legendBoolean)
        #DefaultAssay(SeuratObject())<-"RNA"
        
        #   HoverLocator(plot = df, information = FetchData(object = get(SeuratObject()), vars = c("ident", "CellType", "Patient","Condition","Dataset")))
    })
    output$biplot <- renderPlotly({
        fviz_pca_biplot(res.pca, repel = TRUE,
                        col.var = "#2E9FDF", # Variables color
                        col.ind = "#696969"  ,# Individuals color
                        axes = c(1,2)
        )
        # DefaultAssay(SeuratObject())<-"RNA"
        #    df<-get(SeuratObject())
        # Idents(df)<-"CellType"
        #    df <- subset(get(SeuratObject()), downsample = 1000)
        #  DimPlot(subset(get(SeuratObject()), downsample = 500), group.by=input$cellIdent, 
        #          do.label = input$labelBoolean, 
        #           pt.size = input$dotSize, 
        ##           label.size = input$labelSize, 
        #           no.legend=input$legendBoolean)
        #DefaultAssay(SeuratObject())<-"RNA"
        
        #   HoverLocator(plot = df, information = FetchData(object = get(SeuratObject()), vars = c("ident", "CellType", "Patient","Condition","Dataset")))
    })    
    #Biplot
 

    #============================================================================
    vPlot <- eventReactive(input$goButton, {
        if (is.null(input$file_name)){ 
            load(file=demo)
        } else {
            load(input$file_name$datapath)
        }
        # DefaultAssay(SeuratObject())<-"RNA"
        # VlnPlot(get(SeuratObject()), features = input$geneName, pt.size = input$dotSize, group.by=input$cellIdent,sort=T,assay="RNA")+ggtitle(input$geneName)
        VlnPlot(subset(get(SeuratObject()), downsample = 500), group.by=input$cellIdent, 
                do.label = input$labelBoolean, 
                features = input$geneName,
                pt.size = input$dotSize, 
                label.size = input$labelSize, 
                no.legend=input$legendBoolean,sort=T,assay="RNA")+ggtitle(input$geneName)
    })
    #============================================================================

    #============================================================================
    #============================================================================
    pPlot <- eventReactive(input$goButton, {
        if (is.null(input$file_name)){ 
            load(file=demo)
        } else {
            load(input$file_name$datapath)
        }
        #DefaultAssay(SeuratObject())<-"RNA"
        DimPlot(get(SeuratObject()), features=input$geneName,pt.size = input$dotSize,reduction="pca")+ggtitle(input$geneName)
    })
    #============================================================================
    uPlot2 <- eventReactive(input$goButton, {
        if (is.null(input$file_name)){ 
            load(file=demo)
        } else {
            load(input$file_name$datapath)
        }
        #DefaultAssay(SeuratObject())<-"RNA"
        df<-get(SeuratObject())
        Idents(df)<-"CellType"
        df <- subset(df, downsample = 3000)
        FeaturePlot(df, features=input$geneName,pt.size = input$dotSize,reduction="pca")+ggtitle(input$geneName)
    })
    
    output$featureTSNE <- renderPlot(fPlot())
    output$vlnPlot <- renderPlot(vPlot())
    output$umapPlot <- renderPlot(uPlot())
    output$pcaPlot <- renderPlot(pPlot())
    
}

shinyApp(ui = ui, server = server)
