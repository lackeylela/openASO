library(shiny)
library(dplyr)
library(ggplot2)
library(DT)
library(plotly)
tif (interactive()) {
  
  ui <- fluidPage(
    headerPanel("openASO"),
    sidebarLayout(
      sidebarPanel(
        fileInput("file1", "Choose Bed File",
                  accept = c('text/csv',
                             'text/comma-separated-values',
                             'text/tab-separated-values',
                             'text/plain','.csv','.tsv')
        ),
        tags$hr(),
        checkboxInput("header", "Header", TRUE)
        
      ),
      mainPanel(
        tabsetPanel(
          tabPanel("ASO Mean Plot", plotOutput("graph2")),
          tabPanel("Interactive Genomic Plot", plotlyOutput("graph1")),
          tabPanel("RNA Structure Score Plot", plotlyOutput("graph3"))
        ),
        DTOutput("table")
        
      )
    )
    
  )
  
  server <- function(input, output) {
    
    observe({
      inFile <- input$file1
      if (is.null(inFile))
        return(NULL)
      raw <- read.csv(inFile$datapath, header = TRUE)
      colnames(raw) <- c("ASO_Sequence","Chromosome", "RNAstructScore1","RNAstructScore2","RNAstructScore3","RNAstructScore4",
                         "RNAstructScore5","RNAstructScore6","RNAstructScore7","RNAstructScore8","RNAstructScore9",
                         "RNAstructScore10","RNAstructScore11","RNAstructScore12","RNAstructScore13",
                         "RNAstructScore14","RNAstructScore15","RNAstructScore16","RNAstructScore17",
                         "RNAstructScore18","Gene", "OriginalGene", "Effectivity", "HGNC", "GeneName", 
                         "Sequence", "ASO_Sequence_Reverse", "TranscritID", "Transcript_Start_Position", 
                         "Transcript_End_Position")
      
      output$table <- renderDT(raw %>% dplyr::select(Gene, Effectivity, ASO_Sequence, Transcript_Start_Position, Transcript_End_Position), filter="top", options=list(pageLength=10))
      
      output$graph1 <- renderPlotly({
        inFile <- input$file1
        if (is.null(inFile))
          return(NULL)
        
        if (!is.null(input$table_rows_all))
          (raw <- raw[input$table_rows_all, ])
        
        graphs <- ggplot(raw, aes(text = paste(
          "Gene: ", Gene, "\n",
          "Sequence: ", ASO_Sequence, "\n",
          "Score: ", Effectivity, "\n",
          sep = ""), color=Effectivity)) +
          geom_point( aes(x=Effectivity, y=Transcript_Start_Position), size=3, shape=15 ) +
          coord_flip()+
          theme(legend.position = "left") +
          xlab("Score") +
          ylab("Genomic Coordinates")+
          scale_color_gradient(low="blue", high="red")
        
        ggplotly(graphs, tooltip = "text")
      })
      
      output$graph2 <- renderPlot({
        inFile <- input$file1
        if (is.null(inFile))
          return(NULL)
        if (!is.null(input$table_rows_all))
          (raw <- raw[input$table_rows_all, ])
        Means <- raw %>% group_by(Gene) %>% add_tally() %>% summarize(Means = mean(Effectivity, na.rm =TRUE), stdev = sd(Effectivity, na.rm = TRUE), count = unique(n)) %>% mutate(upper = Means + stdev, lower = Means - stdev)
        Means$Gene = with(Means, reorder(Gene, Means, median))
        Means %>% filter(count > 10) %>% ggplot() +
          geom_errorbar(aes(x=Gene, ymin=lower, ymax=upper), width=0.1, size=1, color="black") +
          geom_point(aes(Gene, Means, colour = factor(Means)), size = 6) + theme_bw() +
          theme(axis.text.x = element_text(angle = 90, hjust = 1)) + theme(legend.position = "none") +
          theme(axis.text.x = element_text(size = 10)) + theme(axis.text.x = element_text(vjust = 0.5)) +ylim(0,1.2) + theme(plot.title = element_blank(), axis.title.x = element_blank()) + ylab("ASO Effectivity Score")
      })
      
      output$graph3 <- renderPlotly({
        inFile <- input$file1
        if (is.null(inFile))
          return(NULL)
        
        if (!is.null(input$table_rows_all))
          (raw <- raw[input$table_rows_all, ])
        
        bed <- raw %>% 
          dplyr::select(ASO_Sequence, GeneName, "1" = RNAstructScore1, "2" = RNAstructScore2, "3" = RNAstructScore3, 
                        "4" = RNAstructScore4,"5" = RNAstructScore5, "6" = RNAstructScore6, 
                        "7" = RNAstructScore7, "8" = RNAstructScore8, "9" = RNAstructScore9, 
                        "10" = RNAstructScore10, "11" = RNAstructScore11, "12" = RNAstructScore12,
                        "13" =  RNAstructScore13, "14" = RNAstructScore14, "15" = RNAstructScore15, 
                        "16" = RNAstructScore16, "17" = RNAstructScore17, "18" = RNAstructScore18) %>%
          group_by(GeneName)
        
        new.bed <- gather(bed, "BasePos", "RNAStructureScore", c(3:20)) %>% mutate(BasePos = as.numeric(BasePos))
        
        graphs <- ggplot(data=new.bed, aes(x=BasePos, y=RNAStructureScore, group=ASO_Sequence)) +
          geom_line(aes(color=ASO_Sequence))+
          geom_point(aes(color=ASO_Sequence))
        
        ggplotly(graphs, tooltip = "text")
      })
    })
  
    
  }
  
  shinyApp(ui, server)
}
