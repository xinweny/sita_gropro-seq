analyzeProseq <- function(samples, bgCov, geneCov, counts, agg, condition, batch, cond1, cond2) {
  require(DESeq2)
  
  #Background signal filtering
  bgCov$length <- bgCov$V3 - bgCov$V2
  normBg <- (bgCov[, c(4:(length(samples) + 3))] * 10) / bgCov$length
  colnames(normBg) <- samples
  cutoff <- colMeans(normBg)
  
  #Long gene ens signal for size factor generation  
  geneCov2 <- geneCov[(geneCov[, c(7:(length(samples) + 6))] / (geneCov$V3 - geneCov$V2) > cutoff), c(4, 7:(length(samples) + 6))]
  colnames(geneCov2) <- c("gene", samples)
  geneCov2 <- geneCov2[ !is.na(geneCov2$gene),]
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
    rownames(countTable) = countsAgg$V9
    #rownames(countTable) = gsub(".*gene_name (.*); transcript_type.*","\\1",rownames(countTable))
    #rownames(countTable) = gsub("transcript_id (.*); family_id.*","\\1",rownames(countTable))
  } else {
    countTable <- counts[, c(7:(length(samples) + 6))]
    rownames(countTable) <- counts$V4
  }
  
  colData <- data.frame(row.names=colnames(countTable), condition=condition, batch=batch)
  ddsProseq <- DESeqDataSetFromMatrix(countData=countTable,
                                      colData=colData,
                                      design=~batch + condition)
  
  sizeFactors(ddsProseq) = sizeFact
  ddsProseq <- estimateDispersions(ddsProseq)
  ddsProseq <- nbinomWaldTest(ddsProseq)
  
  result <- as.data.frame(results(ddsProseq, contrast=c("condition", cond2, cond1), 
                                 alpha=0.05))
  
  if (!agg) {
    vsd <- vst(ddsProseq)
    
    p1 = plotPCA(vsd, intgroup=c("condition", "batch"))
    print(p1)
  }

  return(result)
}

plotMAProseq <- function (df, cutoff, label) {
  require(ggplot2)
  require(ggrepel)
  
  df$color <- ifelse(is.na(df$padj), "grey", ifelse(df$padj < 5e-2, "red", "grey"))
  df$shape <- ifelse(df$log2FoldChange > cutoff , 24, ifelse(df$log2FoldChange < -cutoff, 25, 21))
  df$log2FoldChange[df$log2FoldChange > cutoff] <- cutoff
  df$log2FoldChange[df$log2FoldChange < -cutoff] <- -cutoff
  
  if (label) {
    df$label <- ifelse(!is.na(df$padj) & df$padj < 0.05, rownames(df), "")
    gg <- ggplot(data=df, aes(x=log10(df$baseMean), y=log2FoldChange, label=label)) + 
      theme_classic() + 
      scale_shape_identity() + 
      geom_point(alpha=1, size=2, shape=df$shape, fill=df$color) + 
      geom_text_repel(size=3.5) + 
      scale_colour_manual(values=c('grey', 'red')) + 
      ylim(-cutoff, cutoff) + 
      theme(legend.position="none")
    
    print(gg)
    
  } else {
    gg <- ggplot(data=df, aes(x=log10(df$baseMean), y=log2FoldChange)) + 
      theme_classic() + 
      scale_shape_identity() + 
      geom_point(alpha=1, size=2, shape=df$shape, fill=df$color) + 
      scale_colour_manual(values=c('grey', 'red')) + 
      ylim(-cutoff, cutoff) + 
      theme(legend.position="none")
    
    print(gg)
  }
}

#MAplot Proseq K562 HS
Proseq_K562_Background <- read.table("/data/sawarkar/hummel/externalDataAnalysis/2018_10_09_K562_Proseq_HS_GSE89230/GeneDeserts500kb_Hg38_Gencode_Rel26.counts.txt", sep="\t", header=F)
Proseq_K562_LongEnd <- read.table("/data/sawarkar/hummel/externalDataAnalysis/2018_10_09_K562_Proseq_HS_GSE89230/LongGenes_200k_GRCh38_gencode_rel26.counts.txt", sep="\t", header=F)
Proseq_K562_ProtCoding <- read.table("/data/sawarkar/hummel/externalDataAnalysis/2018_10_09_K562_Proseq_HS_GSE89230/Genes_Filtered_Peaks_Rpm.sorted.distance.length.geneBody.counts.txt", sep="\t", header=F)

proteinCodingK562 <- analyzeProseq(c("NHS1","NHS2","HS1","HS2"), Proseq_K562_Background, Proseq_K562_LongEnd, Proseq_K562_ProtCoding, FALSE, c("NHS","NHS","HS","HS"), c("Rep1","Rep2","Rep1","Rep2"), "NHS", "HS")
plotMAProseq(proteinCodingK562, 4, FALSE)