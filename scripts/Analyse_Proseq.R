#### Packages ####
library(stringr)
library(glue)

#### Config ####
setwd("/Users/Pomato/mrc/project/sita_gropro-seq")

#### Functions ####
analyzeProseq <- function(samples, bgCov, geneCov, counts, agg, condition, batch, case, control) {
  require(DESeq2)
  
  # Background signal filtering
  bgCov$length <- bgCov$end - bgCov$start

  normBg <- (bgCov[, c(4:(length(samples) + 3))] * 10) / bgCov$length
  colnames(normBg) <- samples
  cutoff <- colMeans(normBg)
  

  # Long gene end signal for size factor generation  
  geneCov2 <- geneCov[(geneCov[, c(7:(length(samples) + 6))] / (geneCov$end - geneCov$start) > cutoff), c(4, 7:(length(samples) + 6))]
  colnames(geneCov2) <- c("gene", samples)
  geneCov2 <- geneCov2[!is.na(geneCov2$gene), ]
  print(nrow(geneCov2))
  
  sizeFact <- colSums(geneCov2[, c(2:ncol(geneCov2))]) / mean(colSums(geneCov2[, c(2:ncol(geneCov2))]))
  
  print(sizeFact)

  countTable <- NULL
  
  #For Pol3 genes/ TEs: aggregate counts per type
  if (agg) {
    counts$V9 <- gsub(".*gene_name (.*); transcript_type.*","\\1", counts$V9)
    counts$V9 <- gsub(".*transcript_id (.*); family_id.*","\\1", counts$V9)
    counts$V9 <- gsub("(RNU6ATAC).*","\\1", counts$V9)
    counts$V9 <- gsub("(RNU6)-.*","\\1", counts$V9)
    countsAgg <- aggregate( .~V9, data=counts[, c(9:ncol(counts))], FUN=sum)
    countTable <- countsAgg[2:ncol(countsAgg)]
    rownames(countTable) <- countsAgg$V9
    #rownames(countTable) = gsub(".*gene_name (.*); transcript_type.*","\\1",rownames(countTable))
    #rownames(countTable) = gsub("transcript_id (.*); family_id.*","\\1",rownames(countTable))
  } else {
    countTable <- counts[, c(7:(length(samples) + 6))]
    rownames(countTable) <- counts$name
  }
  
  colData <- data.frame(row.names=colnames(countTable), condition=condition, batch=batch)
  ddsProseq <- DESeqDataSetFromMatrix(countData=countTable,
                                      colData=colData,
                                      design=~condition)
  
  ddsProseq$condition <- relevel(ddsProseq$condition, ref=control)
  
  sizeFactors(ddsProseq) <- sizeFact
  ddsProseq <- estimateDispersions(ddsProseq)
  ddsProseq <- nbinomWaldTest(ddsProseq)
  
  result <- as.data.frame(results(ddsProseq, 
                                  contrast=c("condition", case, control), 
                                  alpha=0.05))
  
  if (!agg) {
    vsd <- vst(ddsProseq)
    
    p1 <- plotPCA(vsd, intgroup=c("condition", "batch"))
    print(p1)
  }

  return(result)
}

plotMAProseq <- function (df, cutoff, label, control, treatment) {
  require(ggplot2)
  require(ggrepel)
  require(dplyr)
  require(glue)
  
  df$color <- ifelse(is.na(df$padj), "grey", ifelse(df$padj < 0.05, "red", "grey"))
  df$shape <- ifelse(df$log2FoldChange > cutoff , 24, ifelse(df$log2FoldChange < -cutoff, 25, 21))
  df$log2FoldChange[df$log2FoldChange > cutoff] <- cutoff
  df$log2FoldChange[df$log2FoldChange < -cutoff] <- -cutoff
  
  nUp <- nrow(filter(df, !is.na(padj) & padj < 0.05 & log2FoldChange > 0))
  nDown <- nrow(filter(df, !is.na(padj) & padj < 0.05 & log2FoldChange < 0))
  ratioUpDown <- nUp / nDown
  message("Ratio of up:down regulated genes: ", ratioUpDown)
  
  if (label) {
    df$label <- ifelse(!is.na(df$padj) & df$padj < 0.05, rownames(df), "")
    gg <- ggplot(data=df, aes(x=log10(baseMean), y=log2FoldChange, label=label)) + 
      theme_classic() + 
      scale_shape_identity() + 
      geom_point(alpha=1, size=2, shape=df$shape, fill=df$color) + 
      geom_text_repel(size=3.5) + 
      scale_colour_manual(values=c('grey', 'red')) + 
      ylim(-cutoff, cutoff) + 
      theme(legend.position="none")
    
    print(gg)
    
  } else {
    gg <- ggplot(data=df, aes(x=log10(baseMean), y=log2FoldChange)) + 
      labs(title=glue("{treatment} vs. {control}"),
           caption=glue("UP={nUp}, DOWN={nDown}")) +
      theme_classic() + 
      scale_shape_identity() + 
      geom_point(alpha=1, size=2, shape=df$shape, fill=df$color) + 
      scale_colour_manual(values=c('grey', 'red')) + 
      ylim(-cutoff, cutoff) + 
      theme(legend.position="none")
    
    print(gg)
  }
}

format_condition <- function (colnames) {
  replace <- replace <- c("_[0-9]*$", "_rep[0-9]*$", "^[A-Z]{3}[0-9]+_", "^[0-9]+_")
  
  for (r in replace) {
    colnames <- gsub(r, "", colnames)
  }
  
  return(colnames)
}

#### Code ####
# Load files
proj <- "GSE112379"
treatment_time <- 90

background <- read.table(glue("./processed/{proj}_{treatment_time}min_GeneDesertCounts.txt"),
                         sep="\t", header=T)
longend <- read.table(glue("./processed/{proj}_{treatment_time}min_LongGeneEndCounts.txt"),
                      sep="\t", header=T)
protcoding <- read.table(glue("./processed/{proj}_{treatment_time}min_GeneBodyCounts.txt"),
                         sep="\t", header=T)

# Format conditions and replicates
samples <- colnames(protcoding)[c(7:ncol(protcoding))]
conditions <- format_condition(colnames(protcoding)[c(7:ncol(protcoding))])
rep <- str_extract(conditions, "[^_][A-Za-z0-9]+$")

control <- "HEK_NHS"
case <- "HEK_HS"

# PRO-seq analysis
proteinCoding_PROseq <- analyzeProseq(samples=samples, 
                                   bgCov=background, geneCov=longend, counts=protcoding, 
                                   agg=FALSE,
                                   condition=conditions, batch=rep,
                                   case=case, control=control)

# MA plot Pro-seq HS
plotMAProseq(df=proteinCoding_PROseq, control=control, treatment=case,
             cutoff=4, label=FALSE)
