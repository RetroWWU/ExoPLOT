
library(stringr)
library(GenomicRanges)
library(HGNChelper)
library(ggplot2)
library(plotly)
library(viridis)

std.error <- function(x) {
  y <- sd(x)/sqrt(sum(!is.na(x)))
  return(y)
}
v.max <- function(x) {
  y <- max(na.omit(x))
  return(y)
}

cpmPlot <- function(f_ct, f_lib, f_edgeR, f_gene, query) {
  
  list.parameters = list()
  
  loc <- str_replace(basename(f_ct), "_.+", "")
  
  # read count file
  df <- read.table(f_ct, row.names = 1)
  df <- t(df)
  colnames(df)[1] <- "LibraryID"
  # read library information
  lib.info <- readRDS(f_lib)
  df <- cbind(df, lib.info[match(df[, 1], lib.info[, 1]), 2:5])
  
  # calculate PSI for cpm
  cpm <- df[, 2:3] 
  if (max(cpm) < 10){
    list.parameters[["error"]] = "Maxium read counts is less than 10."
    return(list.parameters) 
  }
  colnames(cpm) <- c("cpm_hit", "cpm_nonhit")
  cpm.all = cpm$cpm_hit + cpm$cpm_nonhit
  cpm = cpm / cpm.all
  cpm[is.nan(cpm[,1]),1] = 0 # remove NaN
  cpm[is.nan(cpm[,2]),2] = 0
  ## this `cpm` matrix stores raw PSI matrix temporarily
  
  # read normalized data
  ct.lib.all <- readRDS(f_edgeR)
  gr.gene <- readRDS(f_gene)
  # look for gene expression in normalized data
  # stored into ct.lib as matrix
  coord_query = str_split(query, pattern = " ")
  index.chr = coord_query[[1]][2]
  index.chr = as.character(substr(index.chr, 4, nchar(index.chr)))
  start.query = as.integer(coord_query[[1]][3])
  end.query = as.integer(coord_query[[1]][4])
  gr.query = GRanges(Rle(index.chr),
                     IRanges(start = start.query, end = end.query))
  id.gene = findOverlaps(gr.query, gr.gene)
  id.gene = subjectHits(id.gene)
  if (length(id.gene) == 0) {
    list.parameters[["error"]] = "No overlapping Ensembl gene found (in current version)."
    return(list.parameters) 
  }
  id.embl.gene = gr.gene[id.gene,]$gene_id
  # store the current HGNC names as geneSy.approved
  # the alternative HGNC names as geneSy.other
  for (i in 1:length(id.embl.gene)){
    if (id.embl.gene[i] %in% rownames(ct.lib.all$counts)){
      ct.lib = ct.lib.all$counts[id.embl.gene[i],]
	  geneSy.approved = gr.gene[id.gene[i],]$gene_name
	  gene.biotype = gr.gene[id.gene[i],]$gene_biotype
	  geneSy.other = paste(hgnc.table[hgnc.table$Approved.Symbol==geneSy.approved,"Symbol"], collapse = " ")
      break
    }
    if (i == length(id.embl.gene)) {
      # means didnt find any exon overlapping in the whole list
      list.parameters[["error"]] = "No overlapping Ensembl gene found (in current version)."
      return(list.parameters) 
    }
  }
  
  # normlize PSI matrix based on size factor
  cpm = cpm * ct.lib.all$sample[as.character(df$LibraryID), "norm.factors"]
  # integrate TMM values
  cpm = cpm * ct.lib[as.character(df$LibraryID)] * 10^6/df[, 4] 
  
  df <- cbind(df, cpm)
  
  tmp <- order(df[, "order"])  # major sort by stage
  df <- df[tmp, ]
  tissue <- str_extract(df[, "Grouptag"], "\\w+")
  stage <- str_replace(df[, "Grouptag"], "\\w+\\.", "")
  df <- cbind(df, tissue, stage)
  tmp <- order(df[, "tissue"])  # minor sort by tissue
  df <- df[tmp, ]
  
  # paired colors
  color.list <- c("#3399CC", "#A6CEE3", "#33CCFF", "#99E6FF", "#CC0000", "#FB9A99",
                  "#CC9900", "#ECD89F", "#339900", "#B2DF8A", "#CC3399", "#EAADD6", "#FF6600",
                  "#FDBF6F")
  
  # combine same stage lib for main plot
  group <- unique(df[, "Grouptag"])
  df.plot <- data.frame(matrix(ncol = 10, nrow = 0))
  colnames(df.plot) <- c("tissue", "stage", "ct.mean.hit", "ct.conv.std.err.hit", "ct.mean.nonhit",
                     "ct.conv.std.err.nonhit", "cpm.mean.hit", "cpm.conv.std.err.hit", "cpm.mean.nonhit",
                     "cpm.conv.std.err.nonhit")
  for (i in 1:length(group)) {
    tmp <- df[df[, "Grouptag"] == group[i], ]
    ct.hit <- tmp[, 2]
    ct.nonhit <- tmp[, 3]
    cpm.hit <- tmp[, "cpm_hit"]
    cpm.nonhit <- tmp[, "cpm_nonhit"]
    tmp <- data.frame(tissue = as.character(unique(tmp[["tissue"]])), stage = as.character(unique(tmp[["stage"]])),
                      ct.mean.hit = mean(ct.hit), ct.conv.std.err.hit = std.error(ct.hit), ct.mean.nonhit = mean(ct.nonhit),
                      ct.conv.std.err.nonhit = std.error(ct.nonhit), cpm.mean.hit = mean(cpm.hit),
                      cpm.conv.std.err.hit = std.error(cpm.hit), cpm.mean.nonhit = mean(cpm.nonhit),
                      cpm.conv.std.err.nonhit = std.error(cpm.nonhit))
    df.plot <- rbind(df.plot, tmp)
  }
  df.plot[is.na(df.plot)] <- 0
    
  ## count table for download
  df.dl = df.plot[,c("stage","tissue", "cpm.mean.hit","cpm.mean.nonhit")]
  colnames(df.dl) <- c("stage","tissue", "cpm.withTE","cpm.withoutTE")
  
  ### other parameters: tissue list
  level.tissue <- unique(df.plot[, "tissue"])
  level.tissue[c(1, 2)] <- level.tissue[c(2, 1)]  # change order of forbrain and hindbrain
  
  ### other parameters: x,y lim
  cpm.min <- 0
  cpm.max <- max(c(v.max(df.plot[, "cpm.mean.hit"]) + v.max(df.plot[, "cpm.conv.std.err.hit"]),
                   v.max(df.plot[, "cpm.mean.nonhit"]) + v.max(df.plot[, "cpm.conv.std.err.nonhit"])))
#  cpm.max = cpm.max * 1.1
  xlim = c(0, nrow(df.plot) + length(level.tissue))  # set lim for locus graph
  ylim = c(cpm.min, cpm.max)
  
  list.parameters[["geneSy.approved"]] = geneSy.approved
  list.parameters[["gene.biotype"]] = gene.biotype
  list.parameters[["geneSy.other"]] = geneSy.other
  list.parameters[["df"]] = df
  list.parameters[["df.plot"]] = df.plot
  list.parameters[["df.dl"]] = df.dl
#  list.parameters[["p.3d"]] = p.3d
  list.parameters[["loc"]] = loc
  list.parameters[["xlim"]] = xlim
  list.parameters[["ylim"]] = ylim
  list.parameters[["level.tissue"]] = level.tissue
  list.parameters[["color.list"]] = color.list
  
  return(list.parameters)
  
}

graphPlot <- function(list.parameters){
  
  # check error
  if (names(list.parameters) == "error") {
    p = ggplot() +
      geom_text(aes(x=0, y=0, label = list.parameters[["error"]]),
                colour = "red", size = 6) +
      theme_void()
    return(p) 
  }
  
  df = list.parameters[["df.plot"]]
  loc = list.parameters[["loc"]]
  geneNames = paste0("Submit name:", loc, 
                     ", Current gene symbol:", list.parameters[["geneSy.approved"]], 
                     ",", list.parameters[["gene.biotype"]]
  )
  geneNames.2 = gsub(pattern = '(.{1,70})(\\s|$)',
                     replacement = '\\1\n',
                     x = paste0("Available gene names: ",
                                list.parameters[["geneSy.other"]],
                                collapse = " ")
  )
  xlim = list.parameters[["xlim"]]
  ylim = list.parameters[["ylim"]]
  level.tissue = list.parameters[["level.tissue"]]
  color.list = list.parameters[["color.list"]]
  
  # ggplot 
  # reform df to df.ggplot
  df.mean = reshape2::melt(data = df[, !stringr::str_detect(colnames(df), "^ct")], 
                           id.var = c("tissue", "stage"),
                           measure.vars = c("cpm.mean.hit", "cpm.mean.nonhit"),
                           variable.name = "TE", value.name = "cpm.mean")
  df.mean$TE = ifelse(stringr::str_detect(df.mean$TE, "nonhit"), "without", "with")
  df.std.error = reshape2::melt(data = df[, !stringr::str_detect(colnames(df), "^ct")], 
                                id.var = c("tissue", "stage"),
                                measure.vars = c("cpm.conv.std.err.hit", "cpm.conv.std.err.nonhit"),
                                variable.name = "TE", value.name = "std.error")
  ## df.std.error$TE = ifelse(stringr::str_detect(df.std.error$TE, "nonhit"), "without", "with")
  df.ggplot = cbind(df.mean, df.std.error$std.error)
  colnames(df.ggplot)[ncol(df.ggplot)] = "std.error"
  df.ggplot = df.ggplot[order(df.ggplot$tissue),]
  case = paste0(df.ggplot$tissue,"-",df.ggplot$TE,"TE")
  df.ggplot$case = case
  df.color = data.frame(case = unique(case), color.list)
  df.ggplot = dplyr::left_join(df.ggplot,df.color, by = c("case"))
  colnames(df.ggplot)[ncol(df.ggplot)] = "color"
  # sort stage from df
  levels.stage = df.ggplot[df.ggplot$tissue=="Liver" & df.ggplot$TE=="with","stage"]
  levels.stage = append(levels.stage, "1y", after = 16)
  df.ggplot$stage = factor(df.ggplot$stage, levels = levels.stage)
  df.ggplot$case = factor(df.ggplot$case, levels = unique(df.ggplot$case))
  df.ggplot$index = match(df.ggplot$stage, levels.stage)
  
  p = ggplot(df.ggplot, aes(x = index, y = cpm.mean, group = TE,
                            colour = case, stage = stage)) +
    geom_vline(xintercept = 15, color="grey") +
    geom_point(size = 0.6) +
    geom_line(size = 0.2) + 
    geom_segment(aes(x=index, xend=index,
                     y=cpm.mean-std.error, yend=cpm.mean+std.error),
                 size = 0.2) +
    scale_x_continuous(breaks = 1:length(levels.stage), labels = levels.stage) +
    scale_color_manual(values = color.list, name = "") +
    xlab("stage") +
    facet_grid(~tissue, scales = 'free') +
    ggtitle(label = geneNames, subtitle = geneNames.2) +
    theme_classic() +
    theme(axis.text.x = element_text(size=4, angle=90, hjust=1),
          plot.title = element_text(vjust = 1),
          panel.spacing = unit(0.05, "cm"),
          panel.border = element_rect(color = "black", fill = NA, size = 1))
  
  return(p)
  
}

heatmapPlot <- function(list.parameters){
	
	if (names(list.parameters) == "error") {
		plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
		text(x = 0.5, y = 0.5, list.parameters[["error"]], cex = 1.6, col = "red")
		return() 
    }
	
	df = list.parameters[["df.dl"]]
	loc = list.parameters[["loc"]]
	xlim = list.parameters[["xlim"]]
	level.tissue = list.parameters[["level.tissue"]]
	level.stage = c("prenatal", "postnatal")
	
	df$stage = str_replace(df$stage, ".*wpc", level.stage[1])
	df$stage = str_replace(df$stage, "[^(prenatal)].*", level.stage[2])
	
	data = data.frame()
	
	for (tissue in level.tissue) {
	  for (stage in level.stage) {
	    new.row = data.frame(stage, tissue, colMeans(df[(df$stage==stage & df$tissue==tissue),c("cpm.withTE", "cpm.withoutTE")]))
	    colnames(new.row)[3] = "cpm"
	    new.row$TE = str_remove(rownames(new.row), "cpm.")
	    rownames(new.row) = NULL
	    data = rbind(data, new.row)
	  }
	}
  data$stage = factor(data$stage, levels = c("prenatal", "postnatal"))
	data$TE = factor(data$TE, levels = c("withTE", "withoutTE"))
  data = data[!is.nan(data$cpm),]
	
	p = ggplot(data, aes(stage, TE, fill = cpm, text = cpm)) + geom_tile() + 
	  scale_fill_viridis(discrete=FALSE) + facet_wrap(~tissue, nrow = 1) +
    theme_minimal() +
	  theme(legend.position = "bottom",
	        panel.grid = element_blank(),
	        axis.title = element_blank(),
	        axis.ticks= element_blank(),
	        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.9),
	        panel.background = element_blank())

	return(p)
}
