# app.R
#
# data files are generated from the tandem_get_data and intensity_plot function
# r shiny plot function for visualizing the location of the PTM modification

library(shiny)

# Define UI for sample app
myUi <- fluidPage(
  
  titlePanel("Post-translational Modification Plot"),
  colourInput("col", "Select colour", "purple"),
  sidebarLayout( position = "right",
                 sidebarPanel(
                   sliderInput(inputId = "Seq_index",
                               label = "#th sequence:",
                               min = 1,
                               max = 6,
                               value = 2)
                 ),
                 
                 
                 mainPanel(
                   textOutput(outputId = "uid")
                   ,
                   plotOutput(outputId = "bar")
                 )
  )
)

# Define server logic to subject samples to large numbers of small changes,
# subject to a left-hand bound
# read data from csv files to get uid and peptide id and modification of the sequence
myServer <- function(input, output) {
  
  output$bar <- renderPlot({
    #read data from cvs files
    data_peptide<-  read.table(paste("../extdata/data_tandem.csv"), header=TRUE, sep=",")
    data_PTM <- read.table(paste("../extdata/data_PTM_",1,".csv",sep=""), header=TRUE, sep=",")
    #plot the graph
    title<-paste("Protein sequence:",data_peptide$Uid[[input$Seq_index]] ," with peptide id: ",
                 data$Pep_ids[[input$Seq_index]],"vs PTM")
    p<-ggplot(data=data_PTM, aes(x=seq, y=indicators)) + 
      geom_bar(stat="identity", fill=input$col)+
      theme_minimal() +
      labs(title=title, y="Post-translational Modification", x="Protein sequence")
    print(p)
    
  })
  output$uid <- renderText({
    data_peptide<-  read.table(paste("../extdata/data_tandem.csv"), header=TRUE, sep=",")
    sprintf("Uid: %f", data_peptide$Uid[[input$Seq_index]])
  })
  
}

shinyApp(ui = myUi, server = myServer)

# [END]
