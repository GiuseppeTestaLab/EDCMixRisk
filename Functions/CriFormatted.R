#' compareResultsFC
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
#' 
#' @return A list containing the combined results of the two analyses, the scatter plot and the correlation results.

compareResultsFC <- function(resultA, resultB, FDRth=0.05, FCth=2, logCPMth=NULL, FDRceil=1e-30, FCceil=8, title=NULL, geneLabel=TRUE, 
                             colA='#5EB2F280', colB='#725BFF80', colShar='#FFE369', corMethod='spearman') {
  
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
                                                            ifelse(AllRes$Status=='SignificantA&B', 3, 1.5)), 
               shape=ifelse(AllRes$Status=='SignificantA&B', 21, 19)) +
    geom_text(data=dplyr::filter(AllRes, Status=='SignificantA&B'), aes(x=log2FCA, y=log2FCB, label=genes),
              size=2.2, vjust='top', nudge_y=-0.075, alpha=ifelse(geneLabel==TRUE, 1, 0)) +
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
    title = paste('Fold Change: comparison between', name1, 'and', name2)
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
    ggtitle(title) +
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





#' rankGeneVector
#' 
#' Retrieves ranked gene lists from toptag$table object
#' 
#' @param topTagsTable An object structured as topTagTable produced by edgeR.
#' @param Order A character string specifying the column of the topTagTable according to which the gene vector is ordered. Default: PValue.
#' @param PValSel If different from NULL (default), a selection is performed according to this threshold.
#' 
#' @return An ordered vector of genes (as row.names of the topTags). 

rankGeneVector <- function(topTagsTable, Order='PValue', PValSel=NULL){
  # Check that ordering is properly dealt with
  if(!(Order %in% c('PValue', 'logFC', 'StatSign'))){
    stop('Column for ordering not recognized!')
  }
  # Perform a selection based on p-value if PValSel is set
  topTagsTable <- topTagsTable %>% dplyr::mutate(Gene=row.names(topTagsTable))
  if(is.null(PValSel) == FALSE){
    topTagsTable <- topTagsTable %>% dplyr::filter(PValue < PValSel)
  }
  Res <- topTagsTable %>%  dplyr::arrange(!! rlang::sym(Order))
  # to use a variable in arrange I have to adopt the syntax written above
  GeneVector <- dplyr::pull(Res, which(names(Res)==Order))
  names(GeneVector) <- dplyr::pull(Res, Gene)
  # If the order is according to fold-change, I reverse so that top-genes are up-regulated
  if(Order=='logFC' | Order=='StatSign'){
    GeneVector <- rev(GeneVector)
  }
  return(GeneVector)
}



#' topGOResults
#' 
#' Generates TopGO object, performs TopGO enrichment analysis and return the result table
#' 
#' @param Genes Named vector: names are genes in the universe (e.g. gene symbol); 1/0 value for DEGs/nonDEGs.
#' @param gene2GO 
#' @param ontology String for ontology domain: BP (default), MF or CC
#' @param description
#' @param nodeSize topGO node size; default 8.
#' @param algorith String for topGO algorithm; default weight01.
#' @param statistic String for topGo statistical test; default fisher.
#' @param PvalTh Numeric: threshold on the pvalue for the selection of enriched terms.
#' @param EnTh Numeric: threshold on the enrichment, calculated as significant/expected
#' @param minTerms Numeric: minimum number of terms to be reported. If higher than significant terms, also not significant are reported.
#' 
#' @return A list containing the complete results and the GO categories selected as significantly enriched.

topGOResults <- function(Genes, gene2GO, ontology='BP', description=NULL, nodeSize=8, algorithm='weight01', statistic='fisher', PvalTh=0.01, EnTh=1.5, minTerms=15) {
  
  # From RNASeqDEFunctionsGenecode, version 6 February 2019
  
  # 1. Setting description if null
  if (is.null(description)) {
    desc <- paste(ontology, deparse(substitute(Genes)))
  }
  
  # 2. Create list in which useful results are stored
  
  TopGORes <- list()
  
  # 3. TopGO object 
  TopGORes$GOdata <- new('topGOdata', ontology=ontology, allGenes=Genes, annotat=annFUN.gene2GO, gene2GO=gene2GO, description=desc, nodeSize=nodeSize)
  
  # 4. Statistical test
  TopGORes$Test <- topGO::runTest(TopGORes$GOdata, algorithm=algorithm, statistic=statistic)
  # message(Test@geneData)
  
  # 5. Generation of result table
  TopGORes$ResAll <- topGO::GenTable(TopGORes$GOdata, Statistics=TopGORes$Test, topNodes=length(TopGORes$Test@score)) 
  # Selection based on enrichment threshold
  
  ESel <- TopGORes$ResAll %>% dplyr::mutate(ER=round(Significant/Expected,2)) %>% dplyr::filter(ER > EnTh)
  if(dim(ESel)[1] < minTerms){
    TopGORes$ResSel <- ESel
    print('No enough enriched terms for barplot production')
    return(TopGORes)
  }
  n <- max(dim(dplyr::filter(ESel, as.numeric(Statistics) <= PvalTh))[1], minTerms)
  TopGORes$ResSel <- ESel[1:n,]
  # n for the selection of the results is set to the minimal number of terms if the selection
  # on pval gives a number of terms lower than the minimum number
  return(TopGORes)
}



#' topGOBarplot
#' 
#' Visualizes the results of TopGO analysis as a p-value based barplot.
#' 
#' @param TopGORes Object (list) produced by topGOResults function and containing TopGO enrichment analysis results.
#' @param terms Number of ontology terms to be represented in the barplot; default 15.
#' @param pvalTh Threshold on p-value to be plotted as a line; default 0.01.
#' @param title Plot title. 
#' @param palette Color palette; if NULL (default), yellow-green palette is used.
#'
#' @return A barplot.

topGOBarplot <- function(TopGORes, terms=15, pvalTh=0.01, title=NULL, palette=NULL, fill=1) {
  
  # From RNASeqDEFunctionsGenecode, version 6 February 2019

  # 1. Definition of the color palette
  if (is.null(palette)) {
    palette <- colorRampPalette(c('gold', 'forestgreen'))(terms)
  }
  
  # 2. Title definition
  if (is.null(title)) {
    title = deparse(substitute(TopGORes))
  }
  
  # 3. Dataframe reorder
  # first I check for non numeric (<1e-30) values and put a ceiling at -30
  TopGORes$Statistics <- ifelse(grepl('<', TopGORes$Statistics), 1e-30, TopGORes$Statistics)
  # then I order the results
  ResOrdered <- base::transform(TopGORes, GO.ID=reorder(GO.ID, -as.numeric(Statistics)))[1:terms,]
  
  # 4. x-axis limit definition
  MaxVal <- round(max(-log10(as.numeric(ResOrdered$Statistics)), na.rm=TRUE), 0) +1
  
  # 4. BarPlot empty
  if (is.null(fill)){
    TopGOBarplot <- ggplot(data=ResOrdered, aes(x=GO.ID, y=-log10(as.numeric(Statistics)), fill=GO.ID)) + 
      geom_bar(stat='identity', aes(alpha=0.75)) +
      #geom_text(aes(y=0), label=ResOrdered$Term, hjust=0) + 
      scale_y_continuous(breaks=seq(0,MaxVal,2), labels=abs(seq(0, MaxVal, 2)), limits=c(0,MaxVal), expand=c(0.025, 0.025)) +
      geom_hline(yintercept=-log10(pvalTh), col='darkred', lty='longdash') +
      coord_flip() + 
      scale_fill_manual(values=palette) +
      ylab('-log10 PValue') + xlab('') +
      ggtitle(paste('TopGO Enrichment results: ', title)) +
      theme_bw() +
      theme(legend.position='none', axis.title.x = element_text(face = 'bold', colour = 'grey30', size=12), 
            plot.title= element_text(face='bold', colour='darkred', size=12))
  }
  else{
  # 4. BarPlot with text
  TopGOBarplot <- ggplot(data=ResOrdered, aes(x=GO.ID, y=-log10(as.numeric(Statistics)), fill=GO.ID)) + 
    geom_bar(stat='identity', aes(alpha=0.75)) +
    geom_text(aes(y=0), label=ResOrdered$Term, hjust=0) + 
    scale_y_continuous(breaks=seq(0,MaxVal,2), labels=abs(seq(0, MaxVal, 2)), limits=c(0,MaxVal), expand=c(0.025, 0.025)) +
    geom_hline(yintercept=-log10(pvalTh), col='darkred', lty='longdash') +
    coord_flip() + 
    scale_fill_manual(values=palette) +
    ylab('-log10 PValue') + xlab('') +
    ggtitle(paste('TopGO Enrichment results: ', title)) +
    theme_bw() +
    theme(legend.position='none', axis.title.x = element_text(face = 'bold', colour = 'grey30', size=12), 
          plot.title= element_text(face='bold', colour='darkred', size=12))
  }
}



#' topGOBarplotAll
#' 
#' From a DEA analysis, produces 3 barplots for enrichment analysis on all DEGs or splitted in up-regulated and down-regulated.
#' For all and up the generic function (topGOBarplot) is applied; for down-regulated genes the function is re-organized to work properly to produce a symmetric plot.
#' 
#' @param TopGOResAll Object (list) produced by topGOResults on all DEGs.
#' @param TopGOResDown Object (list) produced by topGOResults on down-regulated genes.
#' @param TopGOResUp Object (list) produced by topGOResults on up-regulated genes.
#' @param terms Number of ontology terms to be represented in the barplot; default 15.
#' @param pvalTh Threshold on p-value to be plotted as a line; default 0.01.
#' @param title Plot title. 
#'
#' @return A panel with 3 barplots.
#' 

topGOBarplotAll <- function(TopGOResAll, TopGOResDown, TopGOResUp, terms=15, pvalTh=0.01, title=NULL) {

  # From RNASeqDEFunctionsGenecode, version 6 February 2019
  
  # 1. Definition of the color palette
  YG <- colorRampPalette(c('gold', 'forestgreen'))(terms)
  YR <- colorRampPalette(c('gold', 'red'))(terms)
  YB <- colorRampPalette(c('gold', 'blue'))(terms)
  
  # 2. Barplot for all genes
  BarAll <- topGOBarplot(TopGOResAll, terms=terms, pvalTh=pvalTh, title=title, palette=YG)
  
  # 3. Barplot for up genes
  BarUp <- topGOBarplot(TopGOResUp, terms=terms, pvalTh=pvalTh, title=title, palette=YR)
  
  # 4. Barplot for down genes
  TopGOResDown$Statistics <- ifelse(grepl('<', TopGOResDown$Statistics), 1e-30, TopGOResDown$Statistics)
  ResOrdered <- base::transform(TopGOResDown, GO.ID=reorder(GO.ID, -as.numeric(Statistics)))[1:terms,]
  MaxVal <- round(max(-log10(as.numeric(ResOrdered$Statistics)), na.rm=TRUE), 0) +1
  
  if (is.null(title)) {
    title = deparse(substitute(TopGOResDown))
  }
  
  
  BarDown <- ggplot(data=ResOrdered, aes(x=GO.ID, y=log10(as.numeric(Statistics)), fill=GO.ID)) + 
    geom_bar(stat='identity', aes(alpha=0.75)) +
    geom_text(aes(y=0), label=ResOrdered$Term, hjust=1) + 
    scale_y_continuous(breaks=seq(-MaxVal,0,2), labels=abs(seq(-MaxVal,0, 2)), limits=c(-MaxVal,0), expand=c(0.025, 0.025)) +
    geom_hline(yintercept=log10(pvalTh), col='darkred', lty='longdash') +
    coord_flip() + 
    scale_fill_manual(values=YB) +
    ylab('-log10 PValue') + xlab('') +
    ggtitle(paste('TopGO Enrichment results: ', title)) +
    theme_bw() +
    theme(legend.position='none', axis.title.x = element_text(face = 'bold', colour = 'grey30', size=12), 
          plot.title= element_text(face='bold', colour='darkred', size=12))
  
  # 5. Baplots
  grid.arrange(BarAll, BarDown, BarUp, ncol=3)
}



#' geneStripPairEDCMix
#' 
#' Starting from an SE and a GeneSet, it produces a strip-plot stratifying expression levels for exposure. 
#' 
#' @param SE Object SE; stratification occurs according to EXPO in SE metadata.
#' @param GeneSet Vector of gene symbols to be plotted.
#' @param Key String to specify if genes are in a different nomenclature.
#' @param title Plot title.
#' @param SampleColors Color palette; if null (default) .
#' @param relev If different from NULL (default), vector of leves to re-order exposure.
#' @param printExp Logical; whether expression data should be returned.
#'
#' @return A stripchart.
#'


geneStripPairEDCMix <- function(SE, GeneSet, Key=NULL, title=NULL, SampleColors=NULL, relev=NULL, printExp=TRUE, ...){
  
  # Starting from SE and a GeneSet (vector of gene symbols) it produced a stripplot with the expression levels of the genes
  # stratified for exposure, as specified in SE metadata (from colData(SE)$EXPO)
  
  # 1. Setting of plot title
  if (is.null(title)) {
    title = deparse(substitute(SE))
  }
  
  # 2. Check for missing genes
  # If missing genes are present, a message is printet and then they are discarded from the GeneSet
  Missing <- GeneSet[!(GeneSet %in% row.names(assays(SE)$corrected))]
  
  if(length(Missing)>0){
    print(paste('Expression values are not available for the following genes:', 
                paste(Missing, collapse=' ')))
    GeneSet <- GeneSet[(GeneSet %in% row.names(assays(SE)$corrected))]
  }
  
  # 3. Data selection 
  # If key null, Gene is excpected to be in the same annotation as the DGE$genes row.names
  # Alternative keys still to be implemented
  
  for (i in 1:length(GeneSet)){
    
    Gene <- GeneSet[i]
    
    if (is.null(Key) == TRUE){
      GeneSelection <- row.names(assays(SE)$corrected) %in% Gene
    } else{
      stop('Key different from gene symbol not implemented yet for SE')
    }
    
    GeneCorrectedExp <- assays(SE)$corrected[GeneSelection, ]
    DataTemp <- data.frame(Sample=names(GeneCorrectedExp), Exp=GeneCorrectedExp, Group=colData(SE)$EXPO, Gene=Gene)
    
    if(i==1){
      Data <- DataTemp
    }else{
      Data <- rbind(Data, DataTemp)
    }
    
    if (! is.null(relev)){
      Data$Group <- relevel(Data$Group, ref=relev)
    } 
  } # Closure of cycle on GeneSet
  
  # 4. Print  expression data 
  if(printExp==TRUE){
    print(Data)
  }
  
  # 5. Color Scale
  if (is.null(SampleColors)){
    SampleColors <- scales::hue_pal()(length(levels(Data$Group)))
  } 
  
  else if (SampleColors=="Default"){
    SampleColors <- c("CNT"="#0000FF", "DMSO"="#2900D5", "0.1X"="#5500AA", "1X"="#7E0080", "10X"="#AA0054", "100X"="#D3002B", "1000X"="#FF0000", "BPA0.04X"="#117733", "BPA1X"="#999933", "Triclosan3nM"="#117744", "Triclosan100nM"="#999944", "VITC"="lightgrey", "T3"="yellow",  "T3MixN"="orange")
  } 
  
  # 6. Stripchart
  SC1 <- ggplot(data=Data, aes(y=Exp, x=Gene, fill=Group)) +
    #geom_jitter(position=position_dodge(0.8))
    geom_jitter(position=position_dodge(0.6), size=5, pch=21) + 
    stat_summary(fun.y=mean, fun.ymin=mean, fun.ymax=mean,
                 geom="crossbar", width=0.5, col='gray35', position=position_dodge(0.6)) + 
    ggtitle(paste(title, '\n', 'Gene Levels')) + ylab('SVA corrected expression') +
    scale_fill_manual(values=SampleColors) +
    theme_bw() + xlab('') +
    theme(plot.title = element_text(face='bold', colour='darkred', size=18, hjust=0.5), 
          axis.title=element_text(size=14), axis.text=element_text(size=12.5, angle=45, hjust=1)) 
  
  # 7. Plot 
  SC1
}

