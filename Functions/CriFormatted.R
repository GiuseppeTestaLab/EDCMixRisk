#' dotplot.multintersect.mod
#' 
#' Modified version of dotplot.multintersect
#' 
#' @param th threshold of significance for enrichment to be plotted 

dotplot.multintersect.mod <- function (m, sizeRange = c(0, 20), colors = c(low="blue", mid="grey", 
                                                                           high = "yellow"), forceDivergent=FALSE, th=0.05) 

  # Added threshold parameter: it should plot the dot only for enrichment that are above the specified significance threshold
  # The filtering is done on 'size', that is calculated as -log10(m$prob)
  
  {
  ml <- list(fill = log2(m$enr), val = m$m, size = -log10(m$prob))
  ml <- lapply(ml, FUN = function(x) {
    if (nrow(x) == ncol(x) && all(row.names(x) == colnames(x)) && 
        ncol(x) > 2) {
      x[lower.tri(x, diag = TRUE)] <- NA
      x <- x[-nrow(x), -1]
    }
    x <- t(x)
    x[nrow(x):1, ]
  })
  
  for (f in names(ml)) ml[[f]] <- melt(ml[[f]], value.name = f)
  d <- cbind(ml[[1]], do.call(cbind, lapply(ml[-1], FUN = function(x) x[, 
                                                                        3, drop = FALSE])))
  thlog <- -log10(th)
  p <- ggplot(d, aes(Var1, Var2, label = val)) + 
    geom_point(aes(size=ifelse(size >= thlog, size, 0), colour = fill)) + 
    geom_text() + 
    scale_size_continuous(range = sizeRange) + 
    theme_bw() + # added white background to increase readability
    theme(axis.line = element_blank(), axis.title = element_blank(), 
          axis.text.x = element_text(angle = 45, hjust = 1)) + 
    labs(colour = "log2(Enrichment)", size = "-log10(p-value)")
  
  if (any(d$fill < 0, na.rm = TRUE) || forceDivergent) {
    p <- p + scale_color_gradient2(low = colors["low"], mid = colors["mid"], 
                                   high = colors["high"], na.value = "white")
  }
  
  else {
    p <- p + scale_color_gradient(low = colors["low"], high = colors["high"], 
                                  na.value = "white")
  }
  
  p
}



#' compareResultsFCNew
#' 
#' Compare the results, in term of fold-changes, from two differential expression analyses (DEA) producing a scatterplot.
#' 
#' @param resultA Object organized as a topTag$table for the first DEA.
#' @param resultB Object organized as a topTag$table for the second DEA.
#' @param FDRth Numeric threshold on FDR to select DEGs; default 0.05.
#' @param FCth Numeric threshold on FC (in linear scale) to select DEGs; default 2.
#' @param logCPMth Numeric threshold on logCPM (log scale) to select genes to be compared. It is applied only on the first dataset. Default is null, so the filter is not applied. 
#' @param FDRceil Numeric ceiling on FDR; default 1e-30. 
#' @param FCceil Numeric ceiling on the FC for data visualization; default 8. 
#' @param title A character string for plot title; if NULL (default), a general title is put. 
#' @param geneLabel Logical; whether to put gene labels in the plot.
#' @param colA String specifying color for DEGs from the first DEA.
#' @param colB String specifying color for DEGs from the second DEA.
#' @param colShar String specifying color for DEGs shared between the two analyses.   
#' @param corMethod String specifying the method for calculating the correlation coefficient; default is spearman.
#' @param topLab Number of genes for which the gene name is shown; using a negative value, the top-n most significant genes across those significant in both datasets are selected according to FDR in dataset A.
#' 
#' @return A list containing the combined results of the two analyses, the scatter plot and the correlation results.

compareResultsFCNew <- function(resultA, resultB, FDRth=0.05, FCth=2, logCPMth=NULL, FDRceil=1e-30, FCceil=8, title=NULL, geneLabel=TRUE, 
                             colA='#5EB2F280', colB='#725BFF80', colShar='#ffa05c', corMethod='spearman', topLab=-12) {
  
  ## 1. Title definition
  name1 = deparse(substitute(resultA))
  name2 = deparse(substitute(resultB))
  
  if (is.null(title)) {
    title = paste0(name1, '_versus_', name2)
  } 
  
  
  ### 2. Pre-processing of results objects   
  # Starting from topTag objects, I generate a data frame by inner_join dplyr function 
  ## A. Selection of results of statistical test and order by gene 
  message(paste('The number of genes in the first data set is',  dim(resultA)[1]))
  message(paste('The number of genes in the second data set is',  dim(resultB)[1]))
  AllRes <- dplyr::inner_join(resultA, resultB, by='genes', suffix=c('_A', '_B'))
  # Discard genes with LogCPM lower than threshold if setted
  # The filtering is done on the first dataset
  if(!is.null(logCPMth)){
    AllRes <- dplyr::filter(AllRes, logCPM_A > logCPMth & logCPM_B > logCPMth)
  }
  message(paste('The number of common genes that will be examined is',  dim(AllRes)[1]))
  
  ## B.Log10 transformation for FDR and ceiling for FDR and FC
  AllRes <- mutate(AllRes, logFDRA = ifelse(FDR_A < FDRceil, -log10(FDRceil), -log10(FDR_A))) %>%
    mutate(logFDRB = ifelse(FDR_B < FDRceil, -log10(FDRceil), -log10(FDR_B))) %>%
    mutate(log2FCA = ifelse(abs(logFC_A) < log2(FCceil), logFC_A, ifelse(logFC_A > log2(FCceil), log2(FCceil), -log2(FCceil)))) %>%
    mutate(log2FCB = ifelse(abs(logFC_B) < log2(FCceil), logFC_B, ifelse(logFC_B > log2(FCceil), log2(FCceil), -log2(FCceil))))
  
  ## C. Vector of status for differential expression results
  AllRes <- mutate(AllRes, Status = ifelse(logFDRA > -log10(FDRth) & logFDRB > -log10(FDRth) & abs(log2FCA) > log2(FCth) & abs(log2FCB) > log2(FCth), "SignificantA&B",  
                                           ifelse(logFDRA > -log10(FDRth) & abs(log2FCA) > log2(FCth), "SignificantA",
                                                  ifelse(logFDRB > -log10(FDRth) & abs(log2FCB) > log2(FCth), "SignificantB", "NotSignificant"))))
  
  Significant <- dplyr::filter(AllRes, Status=='SignificantA&B' | Status=='SignificantA' | Status=='SignificantB' )
  
  
  ### 3. Correlation test
  CorTest <- cor.test(AllRes$logFC_A, AllRes$logFC_B, method=corMethod)
  
  
  ### 4. Plot for FC 
  P1 <- ggplot(data=AllRes, aes(x=log2FCA, y=log2FCB)) + 
    geom_point(aes(color=Status), fill=colShar, size=ifelse(AllRes$Status=='NotSignificant', 0.5, 
                                                            ifelse(AllRes$Status=='SignificantA&B', 2.5, 1.5)), 
               shape=ifelse(AllRes$Status=='SignificantA&B', 21, 19)) +
    
    geom_text_repel(data=(dplyr::filter(AllRes, Status=='SignificantA&B') %>% dplyr::top_n(topLab, FDR_A)), aes(x=log2FCA, y=log2FCB, label=genes),
              size=2.5, alpha=ifelse(geneLabel==TRUE, 1, 0), max.overlaps=30) +
    annotate("text", x=-log2(FCceil)+0.2, y=log2(FCceil), 
             label= paste(corMethod, '\n Coef:', round(CorTest$estimate, 2), '\n PValue:', round(CorTest$p.value, 8))) +
    geom_vline(xintercept=c(-log2(FCth), log2(FCth)), col='darkred', lty='longdash') +
    geom_hline(yintercept=c(-log2(FCth), log2(FCth)), col='darkred', lty='longdash') + 
    xlab(paste(name1, "log2FC")) + ylab(paste(name2, "log2FC")) +
    scale_colour_manual(values=c('NotSignificant'="#BEBEBE60", 'SignificantA'=colA,
                                 'SignificantB'=colB, 'SignificantA&B'=colShar)) + 
    ggtitle(paste('Fold Change: comparison between', name1, 'and', name2)) +
    xlim(-log2(FCceil), log2(FCceil)) + ylim(-log2(FCceil), log2(FCceil)) +
    #scale_y_continuous(limits=c(-log2(FCceil)-1, log2(FCceil))+0.5) +
    #scale_x_continuous(limits=c(-log2(FCceil)-1, log2(FCceil))+0.5) +
    theme_bw() +
    theme(plot.title = element_text(face='bold', colour='darkred', size=12, hjust=0.5), 
          legend.text=element_text(size=12), legend.position='bottom') 
  
  
  ### 5. Store results in list
  Results <- list()
  Results$AllRes <- AllRes
  Results$Significant <- Significant
  Results$Scatter <- P1
  Results$Correlation$Test <- corMethod
  Results$Correlation$Coefficient <- CorTest$estimate
  Results$Correlation$PVal <- CorTest$p.value  
  return(Results)
} 




