#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
if (!"devtools" %in% installed.packages()){
    install.packages("devtools", repos = "https://stat.ethz.ch/CRAN/")
}
devtools::install_github("lldelisle/usefulLDfunctions")
library(usefulLDfunctions)
safelyLoadAPackageInCRANorBioconductor("GenomicRanges")
safelyLoadAPackageInCRANorBioconductor("plotrix")
safelyLoadAPackageInCRANorBioconductor("ggplot2")
safelyLoadAPackageInCRANorBioconductor("reshape")
source("helpers.R")
ui <- fluidPage(
    # Application title
    titlePanel("Statistics on your nanopore run"),
    h3("Inputs"),
    # flowLayout(
    # fluidRow(
    # column(width = 2,
    fileInput("filesMS_det", "Choose all mappingStats_det.txt",
              accept = c(
                  "text/csv",
                  "text/comma-separated-values,text/plain"),
              multiple = T
    ),
    textOutput("nbOfInputFiles"),
    textInput("expeName", "Line to add to all plot titles"),
    h3("Size distribution without categories"),
    # ),
    fluidRow(
        column(width = 2,
               uiOutput("sliderSizeRange")
        ),
        column(width = 2,
               sliderInput("bins",
                           "Number of bins:",
                           min = 1,
                           max = 50,
                           value = 30)
        ),
        column(width = 2,
               uiOutput("mySampleSelector")
        )
    ),
    fluidRow(
        column(width = 1,
               uiOutput("maxValueForPlotPerRead"),
               downloadButton("downloadPlotPerRead", "Export Plot")
        ),
        column(width = 11,
               plotOutput("plotPerRead")
        )
    ),
    fluidRow(
        column(width = 1,
               uiOutput("maxValueForPlotPerBase"),
               downloadButton("downloadPlotPerBase", "Export Plot")
        ),
        column(width = 11,
               plotOutput("plotPerBase")
        )
    ),
    h3("Size distribution with categories"),
    uiOutput("myIntermediateCategoriesSelector"),
    fluidRow(
        column(width = 1,
               #        uiOutput("maxValueForPlotPerRead"),
               downloadButton("downloadIndividualMappingStats", "Export Plot")
        ),
        column(width = 11,
               plotOutput("plotIndividualMappingStats")
        )
    ),
    fluidRow(
        column(width = 1,
               #        uiOutput("maxValueForPlotPerRead"),
               downloadButton("downloadIndividualMappingStatsProportion",
                              "Export Plot")
        ),
        column(width = 11,
               plotOutput("plotIndividualMappingStatsProportion")
        )
    ),
    fluidRow(
        column(width = 1,
               #        uiOutput("maxValueForPlotPerRead"),
               downloadButton("downloadReadvsBase", "Export Plot")
        ),
        column(width = 11,
               plotOutput("plotReadvsBase")
        )
    )
    # )
)



# Define server logic required to draw a histogram
server <- function(input, output) {
    createFullSizeResData <- reactive({
        if (is.null(input$filesMS_det)){
            return(NULL)
        }
        generateFullSizeResData(input$filesMS_det)
    })
    output$nbOfInputFiles <- renderText({
        if (is.null(input$filesMS_det)){
            return("No file uploaded.")
        }
        paste("You uploaded ", nrow(input$filesMS_det), "files.")
    })
    output$maxFragmentSize <- renderText({
        fullSizeResData <- createFullSizeResData()
        if (is.null(fullSizeResData)){
            return("0")
        }
        return(max(fullSizeResData$length))
    })
    output$sliderSizeRange <- renderUI({
        fullSizeResData <- createFullSizeResData()
        if (is.null(fullSizeResData)){
            myMax <- 100
        } else {
            myMax <- max(fullSizeResData$length) + 1
        }
        sliderInput("sizeRange",
                    "Range of size to display",
                    min = 0,
                    max = myMax,
                    value = c(0, myMax))
    })
    output$mySampleSelector <- renderUI({
        fullSizeResData <- createFullSizeResData()
        selectInput("sampleToDisplay",
                    "Which sample do you want to display:",
                    choices = c("all", unique(fullSizeResData$sample)),
                    selected = "all",
                    multiple = F)
    })
    output$maxValueForPlotPerRead <- renderUI({
        fullSizeResData <- createFullSizeResData()
        if (is.null(fullSizeResData)){
            myMax <- 100
        } else {
            w <- with(fullSizeResData,
                      weighted.hist(length, number, breaks = bins(), plot = F))
            myMax <- max(w$counts)
        }
        numericInput("ymaxPerRead",
                     "Maximum value for y axis",
                     min = 0,
                     max = myMax,
                     value = myMax)
    })
    bins <- reactive({
        seq(input$sizeRange[1], input$sizeRange[2], length.out = input$bins + 1)
    })
    firstLineTitle <- reactive({
        if (input$expeName != ""){
            return(paste0(input$expeName, "\n"))
        } else {
            return("")
        }
    })
    plotFragmentSizePerRead <- function(){
        fullSizeResData <- createFullSizeResData()
        if (!is.null(fullSizeResData)){
            if (input$sampleToDisplay != "all"){
                fullSizeResData <- subset(fullSizeResData,
                                          sample %in% input$sampleToDisplay)
            }
            bins <- bins()
            with(fullSizeResData,
                 weighted.hist(length, number,
                               main = paste0(firstLineTitle(),
                                             input$sampleToDisplay,
                                             "\nby read number"),
                               breaks = bins, ylim = c(0, input$ymaxPerRead)))
            with(fullSizeResData,
                 weighted.hist(length, nMapped, breaks = bins,
                               col = 2, add = T))
            legend("topright", legend = c("Unmapped", "Mapped"),
                   fill = c("white", 2))
        }
    }
    output$plotPerRead <- renderPlot({
        # generate bins based on input$bins from ui.R
        plotFragmentSizePerRead()
    })
    output$downloadPlotPerRead <- downloadHandler(
        filename = function() {
            "nbOfReadPerSize.pdf"
        },
        content = function(file) {
            pdf(file, title = paste(input$expeName, "nbOfReadPerSize"))
            plotFragmentSizePerRead()
            dev.off()
        })
    output$maxValueForPlotPerBase <- renderUI({
        fullSizeResData <- createFullSizeResData()
        if (is.null(fullSizeResData)){
            myMax <- 100
        } else {
            w <- with(fullSizeResData,
                      weighted.hist(length, number * length,
                                    breaks = bins(), plot = F))
            myMax <- max(w$counts)
        }
        numericInput("ymaxPerBase",
                     "Maximum value for y axis",
                     min = 0,
                     max = myMax,
                     value = myMax)
    })
    plotFragmentSizePerBase <- function(){
        fullSizeResData <- createFullSizeResData()
        if (!is.null(fullSizeResData)){
            if (input$sampleToDisplay != "all"){
                fullSizeResData <- subset(fullSizeResData,
                                          sample %in% input$sampleToDisplay)
            }
            bins <- bins()
            with(fullSizeResData,
                 weighted.hist(length, number * length,
                               main = paste0(firstLineTitle(),
                                             input$sampleToDisplay,
                                             "\nby base number"),
                               breaks = bins, ylim = c(0, input$ymaxPerBase)))
            with(fullSizeResData,
                 weighted.hist(length, Coverage,
                               breaks = bins, col = 2, add = T))
            legend("topright", legend = c("Unmapped", "Mapped"),
                   fill = c("white", 2))
        }
    }
    output$plotPerBase <- renderPlot({
        # generate bins based on input$bins from ui.R
        plotFragmentSizePerBase()
    })
    output$downloadPlotPerBase <- downloadHandler(
        filename = function(){
            "nbOfBasesPerBase.pdf"
        },
        content = function(file) {
            pdf(file, title = paste(input$expeName, "nbOfReadPerBase"))
            plotFragmentSizePerBase()
            dev.off()
        })
    output$myIntermediateCategoriesSelector <- renderUI({
        fullSizeResData <- createFullSizeResData()
        allChoices <- c(0.5, 1, 5, 10, 20, 30, 40, 50,
                        100, 200, 500) * 1e3
        if (! is.null(fullSizeResData)){
            possibleChoices <- allChoices[allChoices
                                          < max(fullSizeResData$length)]
        } else {
            possibleChoices <- allChoices
        }
        selectInput("intermediateCat",
                    "Which are the intermediate categories you want to use:",
                    choices = possibleChoices,
                    selected = intersect(possibleChoices,
                                         c(1, 5, 10, 50, 100) * 1e3),
                    multiple = T)
    })
    intervalls <- reactive({
        fullSizeResData <- createFullSizeResData()
        if (! is.null(fullSizeResData)){
            c(0, as.numeric(input$intermediateCat),
              max(fullSizeResData$length) + 1)
        } else {
            c(0, as.numeric(input$intermediateCat))
        }
    })
    individualMatSizePerBase <- reactive({
        fullSizeResData <- createFullSizeResData()
        if (! is.null(fullSizeResData)){
            getIndividualMatSizePerBase(getDfIntervalsTotal(
                getMergedSizeRes(intervalls(), fullSizeResData)))
        } else {
            return(NULL)
        }
    })
    plotIndividualMappingStats <- function(){
        matSizePerBase <- individualMatSizePerBase()
        if (! is.null(matSizePerBase)){
            colorUsed <- rainbow(100)[seq(1, 100,
                                          length.out = nrow(matSizePerBase) + 1)
                                      [1:nrow(matSizePerBase)]]
            par(mar = c(10, 4, 4, 2) + 0.1)
            tryCatch(barplot(matSizePerBase,
                             main = paste0(firstLineTitle(),
                                           "Number of bases ",
                                           "in each size category"),
                             legend = T,
                             ylim = c(0,
                                      max(apply(matSizePerBase, 2, sum)) *
                                          (1 + nrow(matSizePerBase) / 10)),
                             col = colorUsed, las = 2),
                     warning = function(w){
                         barplot(matSizePerBase / 1e6,
                                 main = paste0(firstLineTitle(),
                                               "Number of million of bases",
                                               " in each size category"),
                                 legend = T,
                                 ylim = c(0,
                                          max(apply(matSizePerBase / 1e6,
                                                    2, sum)) *
                                              (1 + nrow(matSizePerBase) / 10)),
                                 col = colorUsed, las = 2)
                     })
            par(mar = c(5, 4, 4, 2) + 0.1)
        }
    }
    output$plotIndividualMappingStats <- renderPlot({
        plotIndividualMappingStats()
    })
    output$downloadIndividualMappingStats <- downloadHandler(
        filename = function(){
            "individualMappingStats.pdf"
            },
        content = function(file) {
            pdf(file, title = paste(input$expeName, "individualMappingStats"))
            plotIndividualMappingStats()
            dev.off()
        })
    plotIndividualMappingStatsProportion <- function(){
        matSizePerBase <- individualMatSizePerBase()
        if (! is.null(matSizePerBase)){
            colorUsed <- rainbow(100)[seq(1,
                                          100,
                                          length.out = nrow(matSizePerBase) + 1)
                                      [1:nrow(matSizePerBase)]]
            par(mar = c(10, 4, 4, 2) + 0.1)
            barplot(prop.table(matSizePerBase, margin = 2),
                    main = paste0(firstLineTitle(),
                                  "Proportion of bases",
                                  " in each size category"),
                    legend = T,
                    ylim = c(0, 1 * (1 + nrow(matSizePerBase) / 10)),
                    col = colorUsed, las = 2)
            par(mar = c(5, 4, 4, 2) + 0.1)
        }
    }
    output$plotIndividualMappingStatsProportion <- renderPlot({
        plotIndividualMappingStatsProportion()
    })
    output$downloadIndividualMappingStatsProportion <- downloadHandler(
        filename = function(){
            "individualMappingStatsProportion.pdf"
            },
        content = function(file) {
            pdf(file, title = paste(input$expeName,
                                    "individualMappingStatsProportion"))
            plotIndividualMappingStatsProportion()
            dev.off()
        })
    output$plotReadvsBase <- renderPlot({
        fullSizeResData <- createFullSizeResData()
        if (! is.null(fullSizeResData)){
            plotReadVsBase(intervalls(), fullSizeResData, firstLineTitle())
        }
    })
    output$downloadReadvsBase <- downloadHandler(
        filename = function(){
            "readvsBase.pdf"
            },
        content = function(file) {
            pdf(file, title = paste(input$expeName, "readvsBase"))
            plotReadVsBase(intervalls(), createFullSizeResData(), 
                           firstLineTitle())
            dev.off()
        })
}

# Run the application
shinyApp(ui = ui, server = server)
