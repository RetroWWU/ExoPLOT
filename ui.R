library(stringr)
library(rmarkdown)
library(plotly)

# list for exonization tab
list.te = readRDS("database/exonization/list.rds.gz")
# exonization: list of ref
ref.list = split(list.te$ref.list[,1], 
                 factor(list.te$ref.list$V2))[list.te$ref.list$V2]


fluidPage(
  titlePanel("Exonization Plots"),
  
  fluidRow(column(width = 12, 
                  br(),
                  div("Real-time screening for reads of alternative exonizations based on a multi-organ and -developmental stages from 313 human transcriptome RNAseq data. A collaboration with ",
                      a("the Henrik Kaessmann lab.", href = "https://www.zmbh.uni-heidelberg.de/kaessmann/")),
                  br()
  )),
  
  sidebarLayout(
    sidebarPanel(
      tabsetPanel(
        tabPanel(title = "Exonization", 
                 br(),
                 # return string ["ref.0", "ref.1", ...]
                 selectInput(inputId = "refSelect",
                             label = "Select a TE exon list",
                             choices = ref.list
                 ),
                 uiOutput("ui.teSelect"),
                 uiOutput("ui.locSelect"),
                 actionButton(inputId = "rExonization",
                              label = "Start!"),
                 helpText("Exonization shown cases were screened prior on the server, the whole process will finish in seconds.")
        ),
        tabPanel(title = "Customized",
                 br(),
                 selectInput(inputId = "genomeSelect",
                             label = "Select search genome",
                             choices = list("", # empty string as first
                                            "GRCh38 - hg38" = "hg38",
                                            "GRCh37 - hg19" = "hg19")
                 ),
                 textInput(inputId = "coordInput",
                           label = "Fill in a pair of coordinates"),
                 textInput(inputId = "locInput",
                           label = "Fill in a name for the locus (without spaces)",
                           value = "loc"),
                 actionButton(inputId = "rCustomized",
                              label = "Start!"),
                 helpText("The customized approach screens 313 RNAseq datasets for your input regions. This process takes about 1 min."),
                 helpText("For the best results, key regions of alternative splicing should be carefully selected (see below)"),
                 a("How to select the fitting set of coordinates?", href = "explainAS.html", target = "_blank")
        )
      )
    ),
    mainPanel(
      h3("Main plot"),
      a("How to use the tool?", href = "how2use.html", target = "_blank"),
      # more options for cssloader see `https://daattali.com/shiny/shinycssloaders-demo/`
      shinycssloaders::withSpinner(plotlyOutput(outputId = "mainPlot", width = "100%", height = "400px"),
                                   type = 5, size = 1.5),
      #shinycssloaders::withSpinner(plotlyOutput(outputId = "threeDPlot", width = "100%"),
      #                             type = 5, size = 1.5),
      br(),
      conditionalPanel(condition = "(input.rExonization != 0 || input.rCustomized != 0)",
                       plotOutput(outputId = "heatmapPlot", width = "98%" ,height = "175px"),
                       br(),
                       # conditionalPanel maybe
                       downloadButton("downloadTSV", "Save cpm table as tsv"),
                       downloadButton("downloadPlot", "Save main plot as pdf"),
                       uiOutput("link2ucscExon"), uiOutput("link2ucscGene"), 
                       downloadButton("customeTrack.ucsc", "Save UCSC custome track"),
                       uiOutput("link2vastdb")
                       )
    )
  ),
  fluidRow(
    br(), br(),
    column(width = 4,
    ),
    column(width = 8,
           h5("Further Reading:"),
           helpText("ExoPLOT: a web-based tool to visualize expression of exonized transposons from multi-organ and developmental stage RNA-seq data. Zhang F., Moreira M.C., Kaessmann H., Schmitz J."),
           helpText("Cardoso-Moreira, M., Halbert, J., Valloton, D., Velten, B., Chen, C., Shao, Y., … & Mazin, P. V. (2019). Gene expression across mammalian organ development. *Nature*, *571*(7766), 505-509.")
    )
  )
)
