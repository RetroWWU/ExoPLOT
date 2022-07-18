function(input, output) {
  library(shinycustomloader)
  library(stringr)
  library(plotly)
  source("cpmPlot.R")
  
  # list for exonization tab
  list.te = readRDS("database/exonization/list.rds.gz")
  
  # chain file for liftover from hg38 back to hg19
  f.chain = "database/liftOver/chain.hg38ToHg19.rds.gz"

  # assign folder to keep records from user
  ## each job will create a dir with `run-time-randomNum`
  ## given directory should be relative to /srv/shiny-server/runningShiny/
  dir.user = "/srv/shiny-server/exz-plot-d/.data"
  
  # ui for ref list
   output$ui.teSelect = renderUI({
     tes = names(list.te[[input$refSelect]])
     names(tes) = tools::toTitleCase(tes)
     selectInput(inputId = "teSelect",
                 label = "Select a TE class",
                 choices = c(list(""), as.list(tes))
     )
   })

  # ui for select loc reactive to selected TE class
   output$ui.locSelect = renderUI({
     loci = list.te[[input$refSelect]][[input$teSelect]]
     names(loci) = loci
     loci[names(loci)] = str_extract(loci, pattern = "^\\S*")
     selectInput(inputId = "locSelect",
                 label = "Select a loc from list",
                 choices = c(list(""), as.list(loci))
     )
   })

  # running triggered by input$rExonization
  observeEvent(eventExpr = input$rExonization, handlerExpr = {
    
    if (input$locSelect == "" | is.null(input$locSelect)){
      list.reactive$list.parameters = list("error" = "locus is not selected.")
      return()
    }
    
    # create temp directory
    list.reactive$dir.month = paste0(dir.user, "/", format(Sys.time(), "%Y-%m"))
    dir.create(list.reactive$dir.month, showWarnings = FALSE)
    list.reactive$id.run = paste(format(Sys.time(), "%Y%m%d%H%M%S"),
                                 as.character(floor(runif(1)*100000000)),
                                 sep = "-")
    list.reactive$dir.temp = paste0(list.reactive$dir.month, "/run-",
                                    list.reactive$id.run)
    dir.create(list.reactive$dir.temp)

    # prepare count table: copy from pre-generated tables
    list.reactive$ct.file.ori = paste0("database/exonization/", 
                                       input$refSelect, "/", 
                                       input$teSelect, "/",
                                       input$locSelect, "_edgeRtable.txt")
    list.reactive$ct.file = paste0(list.reactive$dir.temp, "/",
                                   input$locSelect, "_edgeRtable.txt")
    file.copy(list.reactive$ct.file.ori, list.reactive$ct.file)

    # other parameters
    list.reactive$loc = input$locSelect
    list.reactive$file.query = readLines(paste0("database/exonization/",
                                                input$refSelect, "/list.exons/",
                                                input$teSelect, "_exon.txt"))
    list.reactive$line.query = grep(input$locSelect, list.reactive$file.query)
    list.reactive$string.query = list.reactive$file.query[list.reactive$line.query]
    list.reactive$query = gsub(pattern = "\t", replacement = " ",
                               x = list.reactive$string.query)
    
    # cpm running 
    message(paste(list.reactive$ct.file, list.reactive$query))
    list.reactive$list.parameters = cpmPlot(list.reactive$ct.file, 
                                            "database/lib.info.rds",
                                            "database/edgeR-result.rds", 
                                            "database/gr.rds",
                                            query = list.reactive$query)
    
    # UCSC custome track file
    output$customeTrack.ucsc = downloadHandler(
      filename = function(){paste0("customeTrack", list.reactive$loc, ".txt")},
      content = function(file) {
        data = paste0(
          "track name=teExon description='Exonization' color=0,0,255,\n",
          stringr::str_remove(list.reactive$query, "^\\S+\\s"), "\n"
        )
        write(data, file)
      }
    )
    
    ## store cpm result
    saveRDS(list.reactive$list.parameters, 
            paste0(list.reactive$dir.temp, "/temp.rds.gz"))
  })

  # running triggered by input$rCustomized
  observeEvent(eventExpr = input$rCustomized, handlerExpr = {
    
    if (input$coordInput == "" | is.null(input$coordInput)){
      list.reactive$list.parameters = list("error" = "Coordinates are not provided.")
      return()
    }
    
    # create temp directory
    list.reactive$dir.month = paste0(dir.user, "/", format(Sys.time(), "%Y-%m"))
    dir.create(list.reactive$dir.month, showWarnings = FALSE)
    list.reactive$id.run = paste(format(Sys.time(), "%Y%m%d%H%M%S"),
                                 as.character(floor(runif(1)*100000000)),
                                 sep = "-")
    list.reactive$dir.temp = paste0(list.reactive$dir.month, "/run-",
                                    list.reactive$id.run)
    dir.create(list.reactive$dir.temp)

    # prepare count table: run perl
    ## prepare input for perl
    list.reactive$file.input = paste0(list.reactive$dir.temp, "/input.txt")
    list.reactive$loc = str_remove_all(str_trim(input$locInput), 
                                       pattern = " ")
    ## formating coordinates input$coordInput
    list.reactive$coordInput =
      tolower(
        str_replace_all(
          str_remove_all(str_trim(input$coordInput), 
                         pattern = ","), 
          pattern = "[:-]", replacement = " ")
      )
    ## transform hg38 to hg19
    if (input$genomeSelect == "hg38") {
       library(rtracklayer)
       library(GenomicRanges)
       chain = readRDS(f.chain)
       ### turn the input$coordInput into a `GRange` object.size
       list.reactive$coordInput = str_split_fixed(list.reactive$coordInput, 
                                                  pattern = " ", n = 3)
       list.reactive$gr.hg38 = GRanges(
         seqnames = Rle(as.character(list.reactive$coordInput[,1])),
         ranges = IRanges(start = as.integer(list.reactive$coordInput[,2]),
                          end = as.integer(list.reactive$coordInput[,3]))
       )
       list.reactive$gr.hg19 = liftOver(list.reactive$gr.hg38, chain)[[1]]
       list.reactive$coordInput = paste(as.character(seqnames(list.reactive$gr.hg19)),
                                        start(list.reactive$gr.hg19),
                                        end(list.reactive$gr.hg19))
    }
    list.reactive$query = paste0(list.reactive$loc, " ", list.reactive$coordInput)
    list.reactive$fh = file(list.reactive$file.input)
    write.table(list.reactive$query, list.reactive$fh,
                col.names = FALSE, row.names = FALSE, quote = FALSE)
    unlink(list.reactive$fh)
    ## run perl
    list.reactive$ct.file <- paste0(paste0(list.reactive$dir.temp, "/",
                                           list.reactive$loc),
                                    "_edgeRtable.txt")
    system(paste("perl", "count.pl",
                 "-i", list.reactive$file.input,
                 "-o", list.reactive$ct.file,
                 "-l", "database/hg19-junc-coord"))
    
    # cpm running 
    message(paste(list.reactive$ct.file, list.reactive$query))
    list.reactive$list.parameters = cpmPlot(list.reactive$ct.file, 
                                            "database/lib.info.rds",
                                            "database/edgeR-result.rds", 
                                            "database/gr.rds",
                                            query = list.reactive$query)
    
    # UCSC custome track file
    output$customeTrack.ucsc = downloadHandler(
      filename = function(){paste0("customeTrack", list.reactive$loc, ".txt")},
      content = function(file) {
        data = paste0(
          "track name=teExon description='Exonization' color=0,0,255,\n",
          stringr::str_remove(list.reactive$query, "^\\S+\\s"), "\n"
        )
        write(data, file)
      }
    )
    
    ## store cpm result
    saveRDS(list.reactive$list.parameters, 
            paste0(list.reactive$dir.temp, "/temp.rds.gz"))
  })

  # reactive control
  list.reactive = reactiveValues()
  
  # main plot
  output$mainPlot = renderPlotly({
    if (input$rExonization == 0 & input$rCustomized == 0) {return()}
    # check traces to hide in style(): plotly_json(p)
    # traces in plotly_json(p) is JS indexing, 21:34 -> 22:35 in R
    # in R style: 1:7 -> vline, 22:35 -> error bar
    ggplotly(graphPlot(list.reactive$list.parameters)) %>%
      style(hoverinfo = "none", traces = c(1:7,22:35)) %>%
      config(scrollZoom = TRUE) %>%
      layout(yaxis = list(fixedrange = TRUE),
             showlegend = FALSE,
             margin = list(t=150))
  })
  
  # 3d plot
#  output$threeDPlot = renderPlotly({
#    if (input$rExonization == 0 & input$rCustomized == 0) {return()}
#    list.reactive$p.3d
#  })

  # heatmap
  output$heatmapPlot = renderPlot({
	  if (input$rExonization == 0 & input$rCustomized == 0){return()}
	  heatmapPlot(list.reactive$list.parameters)
  })
  
  
  # download main plot
  output$downloadPlot = downloadHandler(
    filename = function(){paste0("plot-", list.reactive$loc, ".pdf")},
    content = function(file) {
      ggsave(file,width=13, height=5, 
             plot = graphPlot(list.reactive$list.parameters))
    }
  )
  
  # download table 
  output$downloadTSV = downloadHandler(
    filename = function(){paste0("summary-", list.reactive$loc,".tsv")},
    content = function(file){
      if (names(list.reactive$list.parameters) == "names") {
        write.table(list.reactive$list.parameters[["error"]], file, 
                    row.names = FALSE, quote = FALSE)
      } else {
        write.table(list.reactive$list.parameters[["df.dl"]], file, 
                    sep = "\t", row.names = FALSE, quote = FALSE)
      }
    }
  )
  
  # external link
  output$link2ucscExon = renderUI(a(
    href = paste0("https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=", 
                  str_replace(str_replace(str_remove(list.reactive$query, 
                                                     "^\\S* "), 
                                          " ", "%3A"), " ", "%2D")),
    "Link to exon in UCSC genome broswer (hg19)",
    target = "_blank"
  ))
  
  output$link2ucscGene = renderUI(a(
    href = paste0("https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=", list.reactive$list.parameters$geneSy.approved),
    "Link to search gene in UCSC genome broswer (hg19)",
    target = "_blank"
  ))
  
  output$link2vastdb = renderUI(a(
    href = paste0("https://vastdb.crg.eu/wiki/Genes?query=",
                  list.reactive$list.parameters$geneSy.approved, 
                  "&genome=hg19"),
    "Link to VastDB search (hg19)",
    target = "_blank"
  ))
}
