#' sortRows
#'
#' @param x A numeric matrix or data.frame.
#' @param z Whether to scale rows for the purpose of calculating order (default FALSE).
#' @param na.rm Wheter to remove missing values and invariant rows (default FALSE).
#'
#' @return A reordered matrix or data.frame.
sortRows <- function(x, z=F, toporder=NULL, na.rm=F){
  library(seriation)
  if(na.rm){
    x <- x[which( apply(x, 1, FUN = function(y){ !any(is.na(y)) }) |
                    !(apply(x, 1, na.rm=T, FUN=sd) > 0) ), ]
  }
  y <- x
  if(z) y <- t(scale(t(x)))
  ss <- seriate(dist(y), method = "MDS_angle")
  if(is.null(toporder)) return(x[get_order(ss), ])
  o1 <- get_order(ss)
  oa <- aggregate(o1,by=list(top=as.factor(toporder)),FUN=median)
  oa <- oa[order(oa[,2]),1]
  toporder <- factor(as.character(toporder), levels=as.character(oa))
  x[order(as.numeric(toporder),o1),]
}


#' dround
#'
#' Trim to a certain number of digits (equivalent to `format(...,digits=digits)`, except that the output is numeric)-
#'
#' @param x A vector of numeric values
#' @param digits The number of digits to keep
#' @param roundGreaterThan1 Whether to trim also numbers greater than 1 (default FALSE)
#'
#' @return A numeric vector of the same length as `x`
#' @export
#'
#' @examples
#' dround( c(0.00002345, 554356, 12.56) )
dround <- function(x, digits=3, roundGreaterThan1=FALSE){
  if(is.matrix(x) || is.data.frame(x)){
    for(i in 1:ncol(x)){
      if(is.numeric(x[,i])){
        tryCatch({
          x[,i] <- dround(x[,i], digits, roundGreaterThan1)
        }, error=function(e) warning(e))
      }
    }
    return(x)
  }  
  if(roundGreaterThan1){
    w <- 1:length(x)
  }else{
    w <- which(abs(x)<1)
  }
  if(length(w)==0) return(x)
  e <- ceiling(-log10(abs(x[w])))
  x[w] <- round(10^e*x[w],digits-1)/10^e
  x
}


#' Transcription-factor activity analysis
#'
#' A wrapper around a viper-based TFA.
#' Filter low-counts genes before running. Consider SVA/VST.
#'
#' @param se An object of class SummarizedExperiment, or an expression matrix
#' @param dea A DEA table as produced by edgeR or coerced through `as.DEA()`
#' @param design A formula (if `se` is a `SummarizedExperiment`) or model matrix
#' @param regulon A regulon object, or the path to a 3-column network format
#' @param pleiotropy Logical; passed to `viper` (default TRUE)
#' @param testCoef Coefficient(s) of `design` to test for differential transcript activity
#' @param assayName name of the assay to use. By default, will take the first available in:
#' "vst", "corrected", "imputed", "logcpm", "lognorm"
#' @param degs Vector of differentially-expressed genes (if omitted, taken from `dea`)
#' @param ... passed to `viper`
#'
#' @return A `SummarizedExperiment` with TF-level activity
#' @export
TFA <- function(se, dea, design, regulon, pleiotropy=TRUE, testCoef=NULL, assayName=NULL, degs=NULL, ...){
  suppressPackageStartupMessages({
    library(SummarizedExperiment)  
    library(viper)
    library(edgeR)
  })
  if(!is(se, "SummarizedExperiment") && is(design, "formula")) stop("`design` can only be a formula if `se` is a SummarizedExperiment")
  mm <- switch(class(design),
               formula=model.matrix(design, data=as.data.frame(colData(se))),
               data.frame=design,
               matrix=design,
               character=model.matrix(as.formula(paste0("~",design)), data=as.data.frame(colData(se))),
               stop("Unknown `design`") )
  if(is.null(testCoef)){
    if(is.character(design)){
      testCoef <- design
    }else{
      testCoef <- colnames(mm)[ncol(mm)]
    }
    message("Testing ", testCoef)
  }
  if(is(se, "SummarizedExperiment")){
    if(is.null(assayName)) assayName <- intersect(assayNames(se), c("vst", "corrected", "imputed", "logcpm", "lognorm"))[1]
    e <- assays(se)[[assayName]]
  }else{
    e <- se
  }
  if(is.character(regulon) & length(regulon)==1){
    regulon <- aracne2regulon(regulon, e, verbose=FALSE)
  }
  tfg <- strsplit(names(regulon),"/",fixed=TRUE)
  names(tfg) <- names(regulon)
  w <- which(sapply(tfg, FUN=function(x){ any(x %in% row.names(dea)) }))
  regulon <- regulon[w]
  tfg <- tfg[w]
  vi1 <- viper(e, regulon, pleiotropy=pleiotropy, verbose=FALSE, ...)
  vi2 <- SummarizedExperiment( list(logFC=vi1-rowMeans(vi1), viper=vi1) )
  if(is(se, "SummarizedExperiment")) colData(vi2) <- colData(se)
  res1 <- topTable(eBayes(lmFit(vi1, mm)), testCoef, Inf)
  res1 <- res1[,grep("^t$|^F$|^B$",colnames(res1),invert=TRUE)]
  colnames(res1)[(ncol(res1)-3)+1:3] <- c("meanActivity", "PValue", "FDR")
  maxfc <- function(x){ x[order(abs(x), decreasing=TRUE)[1]] }
  if(ncol(res1)>4) res1 <- cbind(max.logFC=apply(res1[,-1*((ncol(res1)-3)+1:3)],1,FUN=maxfc), res1)
  
  vi2 <- vi2[row.names(res1),]
  rd <- dround(res1)
  colnames(rd) <- paste0("activity.", colnames(rd))
  tfg <- sapply(tfg[row.names(rd)], FUN=function(x){ intersect(row.names(dea), x)[1] })
  dea <- dea[order(dea$FDR, dea$PValue),]
  if(!("logFC" %in% colnames(dea)) && any(grepl("logFC|log2FC",colnames(dea)))){
    dea$logFC <- apply(dea[,grep("logFC|log2FC", colnames(dea)),drop=FALSE],1,FUN=maxfc)
  }
  rd2 <- dround(dea[tfg,intersect(c("baseMean","logFC","PValue","FDR"),colnames(dea))])
  row.names(rd2) <- row.names(rd)
  colnames(rd2) <- paste("expression",colnames(rd2),sep=".")
  
  tfs <- lapply(regulon, FUN=function(x){ names(x$tfmode) })
  if(is.null(degs)) degs <- row.names(dea)[which(dea$FDR < 0.05)]
  #bg <- unique(c(intersect(unlist(tfs),row.names(dea)),degs))
  bg <- row.names(dea)
  rd <- cbind(rd,rd2)
  try({
    ea  <- ORA(degs, bg, tfs, 4)
    row.names(ea) <- ea[,1]
    ea <- ea[row.names(rd),c(-1,-5)]
    colnames(ea) <- c("nbTargets", "nbTargetsInDEGs", "targetEnrichment", "targetEnrFDR", "DEtargets")
    rd <- cbind(rd,ea)
  }, silent=FALSE)
  rowData(vi2) <- rd
  vi2
}

ORA <- function (x, bg, gsets, minSize = 5){
  gsets <- lapply(gsets, y = bg, FUN = intersect)
  gsets <- gsets[which(sapply(gsets, length) >= minSize)]
  res <- data.frame(Term = names(gsets), setSize = sapply(gsets, 
                                                          length))
  res$overlap <- sapply(gsets, y = x, FUN = function(x, y) {
    length(intersect(x, y))
  })
  expected <- sapply(gsets, length) * length(x)/length(bg)
  res$Enrichment <- round(log2(res$overlap/expected), 2)
  res$PValue <- dround(sapply(gsets, set2 = x, universe = bg, 
                              FUN = overlap.prob))
  res <- res[which(res$overlap > 0), ]
  res$FDR <- dround(p.adjust(res$PValue))
  res$genes <- sapply(gsets[row.names(res)], y = x, FUN = function(x, 
                                                                   y) {
    paste(intersect(x, y), collapse = ", ")
  })
  row.names(res) <- NULL
  res[order(res$PValue), ]
}


#' regulon.filter
#'
#' Filters a viper regulon by likelihood
#'
#' @param regulon A regulon object
#' @param min.likelihood Minimum likelihood (default 0.15)
#'
#' @return The filtered regulon
#' @export
regulon.filter <- function(regulon, min.likelihood=0.15){
  regulon <- lapply(regulon, FUN=function(x){
    w <- which(x$likelihood>=min.likelihood)
    if(length(w)==0) return(NULL)
    lapply(x,FUN=function(y) y[w])
  })
  regulon[!sapply(regulon, is.null)]
}

#' plotTFA
#'
#' Plots a volcano-like overview of the results of `TFA`
#'
#' @param x A `SummarizedExperiment` produced by the `TFA` function
#' @param color formula to use as color
#' @param size formula to use as size
#'
#' @return A plotly widget
#' @export
plotTFA <- function(x, color=~expression.logFC, size=~-log10(expression.FDR)){
  x <- as.data.frame(rowData(x))
  suppressPackageStartupMessages(library(plotly))
  tx <- paste0(row.names(x),"\nmean activity: ", x$activity.meanActivity, "\nActivity log2FC: ")
  if(!("activity.logFC" %in% colnames(x))){
    x$activity.logFC <- x$activity.max.logFC
  }
  tx <- paste0(tx, x$activity.logFC)
  if("expression.logFC" %in% colnames(x)){
    tx <- paste0(tx, "\nexpression log2FC: ", x$expression.logFC)
  }
  if("targetEnrichment" %in% colnames(x)){
    tx <- paste0(tx, "\ntarget enrichment: ", x$targetEnrichment," (", x$nbTargets," targets)")
    tx <- paste0(tx, "\ntarget ORA FDR: ", x$targetEnrFDR)
  }
  
  p <- plot_ly( data=x,
                x=~activity.max.logFC,
                y=~-log10(activity.FDR),
                type="scatter",
                mode="markers",
                text=tx,
                hoverinfo="text",
                color=color,
                size=size
  )
  layout(p, title="TF activity analysis")
}

#' categoricalComparison
#'
#' Perform differential expression analysis between the different concentration of EDC exposure, considering each condition as a catogorical variable
#'
#' @param counts The matrix with the gene counts.
#' @param design The design matrix with the info about each sample
#' @param correctFor The variable to `correct for` in the GLM (default `line`, i.e. ~line+concentration).
#'
#' @export
categoricalComparison <- function(counts, design, correctFor="line"){
  cc <- relevel(as.factor(as.numeric(as.character(design$concentration))),"0")
  if(!is.na(correctFor)){
    correctFor <- as.character(design[[correctFor]])
    mm <- model.matrix(~correctFor+cc)
  }else{
    mm <- model.matrix(~cc)
  }
  dds <- calcNormFactors(DGEList(counts))
  dds <- estimateDisp(dds,mm)
  fit <- glmFit(dds,mm)
  lrt <- glmLRT(fit, coef=colnames(mm)[grep("^cc",colnames(mm))])
  return(as.data.frame(topTags(lrt, nrow(counts))))
}

#' combinedRank
#'
#' Returns gene ranks based on both foldchange and p-value.
#'
#' @param logFC a vector of log-foldchanges
#' @param significance a vector of p-values corresponding to each foldchange
#' @param signed logical; whether to consider the direction of the foldchange for ranking (default TRUE). If TRUE, negative foldchanges will be lower and positive foldchanges will be higher.
#' @param df degrees of freedom of the smoothing (default 10)
#'
#' @export
combinedRank <- function(logFC, significance, signed=TRUE, df=10){
  if(signed){
    wD <- which(logFC<0)
    wU <- which(logFC>=0)
    rD <-combinedRank(logFC[wD],significance[wD],signed=FALSE,df=df)
    rU <-combinedRank(logFC[wU],significance[wU],signed=FALSE,df=df)
    r <- vector(mode="integer",length=length(logFC))
    r[wD] <- -(max(rD)-rD-1)
    r[wU] <- max(rU)-rU
    return(r)
  }
  x <- -log10(significance)
  y <- abs(logFC)
  s <- smooth.spline(y~0+x,df=df)
  y2 <- as.data.frame(predict(s,newdata=data.frame(x=x)))
  row.names(y2) <- y2[,1]
  (length(x)+1)-rank(y2[as.character(x),2])
}


#' getDEGsFromRes
#'
#' Extract the list of genes that are significantly differentially expressed from the res object (output of EdgeR)
#'
#' @param res Output of the differential expression analysis.
#' @param threshold The minimal FDR value to be considered as significant (default=0.05)
#' @param minLogFC The minimal value of LogFoldchange to be extracted for each significant gene (default=0).
#'
#' @export
getDEGsFromRes <- function(res, threshold = 0.05, minLogFC=0, logCPM=0){
  if("list" %in% class(res))	return(lapply(res,threshold=threshold,minLogFC=minLogFC,logCPM=logCPM,FUN=getDEGsFromRes))
  fccols <- grep("logFC",colnames(res))
  if(length(fccols)>1)
    row.names(res)[which(res$FDR < threshold & (res$logCPM > logCPM) & apply(res[,fccols],1,minLogFC=minLogFC,FUN=function(x,minLogFC){ any(abs(x) > minLogFC) }))]
  else
    row.names(res)[which((res$FDR < threshold) & (abs(res$logFC) > minLogFC) & (res$logCPM > logCPM))]
}

#' doDEA
#'
#' A wrapper to perform differential expression analysis.
#' 
#' @param counts The matrix with the gene counts.
#' @param design The design matrix with the info about each sample.
#' @param independentVar The independent variable, default `concentration`.
#' @param correctFor The variable to `correct for` in the GLM (default `line`).
#' @param robust logical; whether to force using robust GLM (defaults to TRUE if there are more than 6 samples, FALSE otherwise).
#'
#' @return A data.frame.
#'
#' @export
doDEA <- function(counts, design, independentVar="concentration", correctFor="line", robust=NULL){
  design$line <- as.character(design$line)
  design$type <- as.character(design$type)
  dd <- design
  dd$concentration <- dd[[independentVar]]
  if(length(unique(dd$concentration))==2){
    dd$concentration <- as.character(dd$concentration)
  }else{
    if(independentVar!="line")	concentration <- as.numeric(dd$concentration)	
  }
  if(!is.null(correctFor)){
    if(correctFor %in% colnames(design) & length(unique(design[[correctFor]]))>1){
      dd$correctFor <- design[[correctFor]]
    }else{
      if(!is.na(correctFor))	warning("correctFor is not a valid column name; will be ignored.")
      correctFor <- NULL
    }
  }
  if(is.null(robust))	robust <- !is.null(correctFor) & (nrow(dd) > 6)
  if(is.null(correctFor) && length(unique(dd$concentration))==2){
    dds <- DGEList(counts, group=as.character(dd$concentration))
  }else{
    dds <- DGEList(counts)
  }
  dds <- calcNormFactors(dds)
  if(is.null(correctFor) && length(unique(dd$concentration))==2){
    dds <- estimateDisp(dds)
    res <- exactTest(dds)
  }else{
    if(is.null(correctFor)){
      mm <- model.matrix(~concentration, data=dd)
    }else{
      mm <- model.matrix(~correctFor+concentration, data=dd)
    }
    if(robust){
      dds <- estimateGLMRobustDisp(dds,mm)
    }else{
      dds <- estimateDisp(dds,mm)
    }
    res <- glmLRT(glmFit(dds,mm))
  }
  as.data.frame(topTags(res,nrow(counts)))
}

#' getClusterXvalues
#'
#' Models the mean smoothed foldchange values of a cluster of genes at each concentration.
#' 
#' @param e The matrix with the gene counts.
#' @param design The design matrix with the info about each sample
#' @param geneClusters A vector of cluster assignment for all genes (values are cluster numbers, names are gene names).
#' @param cluster Integer; the cluster for which to extract the smoothed mean foldchanges.
#' @param correctFor The variable to `correct for` in the GLM (default `line`).
#' @param showPoints Logical; whether to show ploints in the plot (default TRUE).
#' @param splitLines Logical; whether to plot cell lines separately (default FALSE).
#' @param spar Smoothing parameter. Default values were calibrated on sample data and depend on the xtype value: for 'ordinal' spar=0.33, otherwise spar=0.7.
#' @param interpolate Integer >= 0. Number of points to interpolate between actual points (default 2).
#' @param yuse Variable to use for the y-axis and clustering. Either 'FC' (foldchange, default), 'absFC' (absolute foldchange), or 'z' (z-score).
#' @param xtype Treatment of the x-axis (concentration) variable. Either 'ordinal' (concentration ranks, default), 'sqrt' (square-root concentration), or 'log10' (log10 or 1+concentration).
#' @param main Plot title.
#' @param col Line color(s).
#' @param shadingCol Shading color(s).
#'
#' @return Produces a plot and returns the modeled values.
#'
#' @export
getClusterXvalues <- function(e,design,geneClusters=NULL,cluster=NULL,correctFor="line",showPoints=TRUE,splitLines=FALSE,spar=NULL,interpolate=2,yuse="FC", xtype="ordinal", main=NULL,col="#E7943D", shadingCol="lightgrey"){
  yuse <- match.arg(yuse, c("absFC","FC","z"))
  xtype <- match.arg(xtype, c("ordinal","sqrt","log10"))
  if(length(col)==1)	col <- rep(col,length(unique(design$line)))
  if(is.null(spar)) spar <- switch(xtype, ordinal=0.33, sqrt=0.7, log10=0.7)
  library(MASS)
  if(!is.null(cluster)){
    e <- e[intersect(row.names(e),names(geneClusters)[which(geneClusters==cluster)]),]
  }
  if(yuse=="z" & is.null(correctFor)){
    o <- order(design$concentration)
    e <- e[,o]
    design <- design[o,]
    z <- t(scale(t(log(e+1))))
  }else{
    o <- order(design[[correctFor]],design$concentration)
    e <- e[,o]
    design <- design[o,]
    z <- do.call("cbind",lapply(unique(design[[correctFor]]),yuse=yuse,e=e,bg=design[[correctFor]],conc=design$concentration,FUN=function(x,e,bg,conc,yuse){
      if(yuse %in% c("FC","absFC")){
        tmp <- e[,which(bg==x)]
        w <- which(bg==x & conc==0)
        if(length(w)>1){
          ctl <- rowMeans(e[,w])
        }else{
          ctl <- e[,w]
        }
        ret <- log2( (1+e[,which(bg==x)])/(1+ctl) )
        if(yuse=="absFC") ret <- abs(ret)
        return(ret)
      }else{
        return(t(scale(t(log(e[,which(bg==x)]+1)))))
      }
    }))
    o <- order(design$concentration)
    e <- e[,o]
    design <- design[o,]
    z <- z[,o]
  }
  z <- z[apply(z,1,FUN=function(y){ all(!is.na(y))}),]
  ztop <- aggregate(apply(z,2,FUN=function(x){ as.numeric(quantile(x,0.75))}),by=list(design$concentration),FUN=mean); row.names(ztop) <- ztop[,1]
  zbot <- aggregate(apply(z,2,FUN=function(x){ as.numeric(quantile(x,0.25))}),by=list(design$concentration),FUN=mean); row.names(zbot) <- zbot[,1]
  z <- apply(z,2,na.rm=T,FUN=median)
  ylab <- switch(yuse, FC="log2(foldchange)", absFC="abs(log2(foldchange))", z="z-score")
  xlab <- switch(xtype, ordinal="Concentration rank", log10="log10(1+concentration)", sqrt="sqrt(concentration)")
  x <- as.character(unique(design$concentration))
  tmp <- order(unique(design$concentration))
  names(tmp) <- x
  xx <- switch(xtype,
               ordinal=tmp[as.character(design$concentration)],
               log10=log10(design$concentration+1),
               sqrt=sqrt(design$concentration))
  cols <- col[as.numeric(as.factor(as.character(design[[correctFor]])))]
  if(!showPoints) cols <- rgb(1,1,1,0)
  plot(xx,z,ylab=ylab,xlab=xlab,main=ifelse(is.null(main),paste("Overall pattern for cluster",cluster),main),pch=16,col=cols,bty="n")
  if(yuse %in% c("FC","absFC")) abline(h=0,lty="dashed")
  x2 <- unique(xx)
  ztop <- .plSmooth(x2,ztop[as.character(x),"x"],interpolate=interpolate,spar=0.2)
  zbot <- .plSmooth(x2,zbot[as.character(x),"x"],interpolate=interpolate,spar=0.2)
  polygon(c(ztop$x,rev(ztop$x),ztop$x[[1]]),c(ztop$y,rev(zbot$y),ztop$y[1]),col=maketrans(shadingCol),border=NA,xpd=T)
  if(!is.null(interpolate)){
    if(length(interpolate)>1){
      newX <- interpolate
    }else{
      newX <- c(x2[[1]])
      for(i in 1:(length(x2)-1)) newX <- c(newX,seq(x2[i],x2[i+1],length.out=(interpolate+2))[-1])
    }
  }else{
    newX <- x2
  }
  if(splitLines){
    j <- 0
    for(i in unique(design[[correctFor]])){
      j <- j+1
      w <- which(design[[correctFor]] == i)
      sm <- smooth.spline(xx[w],z[w],spar=spar)
      lines(predict(sm,newX), lwd=2, col=col[j])
    }
  }
  sm <- smooth.spline(xx,z,spar=spar)
  if(!splitLines){
    lines(predict(sm,newX), lwd=2, col=col[[1]])
  }
  z <- predict(sm,unique(xx))$y
  z <- z/(max(z)-min(z))
  names(z) <- unique(design$concentration)
  return(z)
}

.plSmooth <- function(x,y,spar=0.1,interpolate=0,force0=FALSE){
  if(interpolate>0){
    xx <- seq(from=min(x),to=max(x),length.out=length(x)*interpolate)
  }else{
    xx <- x
  }
  d <- data.frame(x=xx, y=predict(smooth.spline(x,y,spar=spar),xx)$y)
  if(force0) d$y[order(xx)[1]] <- 0
  return(d)
}

#' getFoldchangeMatrix
#'
#' Get a matrix of Foldchanges relative to the DMSO sample from the counts matrix of a given set of genes
#' 
#' @param x The matrix with the normalized gene counts.
#' @param design The design matrix with the info about each sample
#' @param by The line or replicate by which to calculate FC (e.g. must include a ctrl related to each treated)
#' @param ctrls Which samples are controls; if omitted, will attempt to fetch from the design matrix.
#' @param maxLog2FC If given, flattens the logFC to the specific value.
#' @param is.log Logical indicating whether `x` is already log-transformed (default FALSE)
#' 
#' 
#' @return A matrix of foldchanges
#'
#' @export
getFoldchangeMatrix <- function(x,design,by="line",ctrls=NULL,maxLog2FC=NULL,is.log=FALSE){
  if(is.null(ctrls)){
    if("concentration" %in% colnames(design)){
      cond <- design$concentration
    }else{
      if("EXPO2" %in% colnames(design)){
        cond <- design$EXPO2
      }else{
        cond <- design$Condition
      }
    }
    ctrls <- which(cond %in% c("DMSO","CNT","0"))
  }
  fc <- as.matrix(x)
  ep <- as.character(design[[by]])
  for(i in unique(ep)){
    w <- which(ep==i)
    wc <- intersect(ctrls,w)
    cm <- rowMeans(x[,wc,drop=F])
    if(is.log){
      fc[,w] <- fc[,w]-cm
    }else{
      for(j in w) fc[,j] <- apply(cbind(x[,j],cm),1,groups=c("A","B"),FUN=log2FC)
    }
  }
  if(!is.null(maxLog2FC)){
    w <- which(abs(fc)>maxLog2FC)
    fc[w] <- sign(fc[w])*maxLog2FC
  }
  fc
}

#' plotFoldchangeMatrix
#'
#' Plot heatmap of Foldchanges of each sample relative to the DMSO for a given set of genes
#' 
#' @param x The matrix with the foldchanges of the gene set.
#' @param design The design matrix with the info about each sample
#' @param geneClusters An optional vector of cluster assignment for all genes (values are cluster numbers, names are gene names), to add to the heatmap annotation.
#' @param show_rownames Logical; whether to show gene names (default FALSE).
#' @param col Color palette.
#' @param breaks Numeric vector indicating the values at which to change colors.
#' 
#' @export
plotFoldchangeMatrix <- function(x,design, geneClusters=NULL,show_rownames=F,col=colorRampPalette(c("blue", "black", "yellow"))(29),breaks=c(-6,seq(-2.5,2.5,length.out=28),6)){
  o <- order(design$line, design$concentration)
  design <- design[o,]
  x <- x[,o]
  cbreaks <- which(!duplicated(design$line))-1
  cbreaks <- cbreaks[which(cbreaks>0)]
  if(is.null(geneClusters)){
    anno <- NULL
  }else{
    anno <- data.frame(row.names=names(geneClusters), cluster=geneClusters)
  }
  pheatmap(getFoldchangeMatrix(x,design),cluster_col=F,show_rownames=show_rownames,annotation_row=anno,gaps_col=cbreaks,color=col, breaks=breaks, border_color=NA)
}


.getEnrichment <- function(sfor,sin,universe){
  sfor <- intersect(sfor,universe)
  sin <- intersect(sin, universe)
  expected <- length(sin)*(length(sfor)/length(universe))
  obs <- length(intersect(sfor,sin))
  return(obs/expected)
}

.getProp <- function(sof,sin){ length(intersect(sof,sin))/length(sin) }



#' plotListClusters
#'
#' Plots gene clusters of each intersection between `degs` and `categories`
#'
#' @param degs a vector of DEGs, or a named list of such vectors for each system
#' @param categories a named list of genesets
#' @param e a gene count matrix
#' @param design a design data-frame
#' @param splitBy the variable used for splitting samples into systems (default `ex2` or `paste(type,mix)`)
#' @param plotEnrichment logical, whether to show enrichment barplots (default TRUE)
#' @param showOnlyTopCluster whether to plot only the top cluster (default FALSE)
#' @param absFC whether to use absolute foldchanges (default FALSE)
#' @param useSpecificDEGs whether to use DEGs specific for each system (default TRUE)
#' @param doCluster whether to cluster genes or plot all genes (default TRUE, i.e. clusters...)
#' @param forceClusterNb number of clusters
#' @param ... further arguments passed to the `plotGenesClusters` function.
#'
#' @return nothing; creates a plot
#'
#' @export
plotListClusters <- function(degs, categories, e, design, splitBy=NULL, plotEnrichment=TRUE, showOnlyTopCluster=FALSE, absFC=FALSE, useSpecificDEGs=TRUE, doCluster=TRUE, forceClusterNb=NULL, ...){
  
  if(!all(row.names(design)==colnames(e))){
    warning("WARNING: the row.names of `design` doesn't match the colnames of `e`! We are assuming that they correspond, but double-check!")
  }
  
  if(is.null(splitBy)){
    if(!("ex2" %in% colnames(design))){
      design$ex2 <- paste(design$type, design$mix)
    }
  }else{
    design[["ex2"]] <- design[[splitBy]]
  }
  
  if(!is.list(degs)) stop("`degs` should be a list of vectors.")
  
  udegs <- unique(unlist(degs))
  
  if(length(unique(design$ex2))!=length(degs) || !all(sort(unique(design$ex2))==sort(names(degs)))){
    if(all(names(degs) %in% unique(design$ex2))){
      message("Some systems do not have attached lists of DEGs... these will be ignored.")
      w <- which(design$ex2 %in% names(degs))
      e <- e[,w]
      design <- design[w,]
    }else{
      stop(paste0("`degs` should either be a character vector or a list; if a list, the names of the elements should match the unique values of `splitBy`, in this case:\n",paste(unique(design$ex2),collapse=", ")))
    }
    
  }
  
  en <- donorm(filterGenes(e[,row.names(design)]))
  
  tested <- row.names(en)
  
  categories <- lapply(categories,y=tested, FUN=intersect)	
  
  nbsys <- length(unique(design$ex2))
  nbcats <- length(categories)
  exps <- unique(design$ex2)
  
  if(plotEnrichment){
    laym <- matrix(0,nrow=nbcats,ncol=2*nbsys)
    j <- 1
    for(i in 1:nbsys){
      laym[,i*2-1] <- j
      laym[,i*2] <- j+1:nbcats
      j <- j+1+nbcats
    }
    layout(laym)
  }else{
    layout(matrix(1:(nbsys*nbcats),nrow=nbcats,ncol=nbsys))
  }
  
  
  probs <- lapply(degs, categories=categories, universe=tested, FUN=function(x, categories, universe){
    sapply(categories,set2=x,universe=universe,FUN=overlap.prob)
  })
  enrichs <- lapply(degs, categories=categories, universe=tested, FUN=function(x, categories, universe){
    ee <- sapply(categories,sin=x,universe=universe,FUN=.getEnrichment)
    ee[which(ee==0)] <- NA
    return(log2(ee))
  })
  
  xlim <- c(min(log10(unlist(probs))),max(unlist(enrichs),na.rm=T))
  
  for(j in 1:nbsys){
    par(mar=c(3,10,3,2))
    if(plotEnrichment) enrichmentBarplot(enrichs[[j]],probs[[j]],names(categories),xlim=xlim)
    tmpd <- design[which(design$ex2==exps[[j]]),]
    par(mar=c(2,4,2,1))
    for(i in 1:nbcats){
      if(useSpecificDEGs){
        gg <- intersect(degs[[j]],categories[[i]])
      }else{
        gg <- intersect(udegs,categories[[i]])
      }
      tmpe <- en[gg,which(design$ex2==exps[[j]]),drop=F]
      fc <- getFoldchangeMatrix(tmpe, tmpd)
      if(absFC){
        fc <- abs(fc)
        ylab <- "abs(log2(FC))"
      }else{
        ylab="log2(foldchange)"
      }
      if(doCluster){
        if(nrow(tmpe)>3){
          if(!is.null(forceClusterNb)){
            co2 <- kmeans(fc,forceClusterNb)$cluster
          }else{
            co2 <- try(getConsClust(fc,6,plot=FALSE),silent=T)
            if(is(co2,"try-error")){
              dev.off()
              co2 <- kmeans(fc,min(5,max(2,floor(sqrt(nrow(fc))))))$cluster
            }
          }
        }else{
          co2 <- 1:nrow(fc)
          names(co2) <- row.names(tmpe)
        }
        plotGenesClusters(fc, tmpd, co2, ...)
        #plotGenesClusters(fc, tmpd, co2, main=paste(names(categories)[i],"-",exps[[j]]), ...)
        #plotClusterAggregation(fc, tmpd, co2, main=paste(names(categories)[i],"-",exps[[j]]), showOnlyTopCluster=showOnlyTopCluster,ylab=ylab)
      }else{
        plotGenesPatterns(fc, tmpd, main=paste(names(categories)[i],"-",exps[[j]]),ylab=ylab)
      }
    }
  }
}


#' enrichmentBarplot
#'
#' Barplot of fold-enrichment and pvalues.
#'
#' @param enr numeric vector of fold-enrichments
#' @param pvals numeric vector of p-values
#' @param categories character vector indicating the names of the sets
#' @param xlim xlim (default auto)
#'
#' @export
enrichmentBarplot <- function(enr, pvals, categories=NULL, xlim=NULL){
  dat <- data.frame(enr=enr, pvals=log10(pvals))
  dat$enr[which(is.infinite(dat$enr))] <- NA
  if(is.null(xlim)) xlim <- c(min(dat$pval),max(dat$enr,na.rm=T))
  barplot(t(dat),horiz=T,beside=T,xaxt="n",yaxt="n",xlab="Fold enrichment", xlim=xlim)
  abline(v=0)
  x1 <- seq(0,max(c(xlim,ceiling(dat$enr)),na.rm=T))
  axis(side=1,at=x1,labels=2^x1)
  x2 <- seq(min(c(xlim,floor(dat$pvals)),na.rm=T),0)
  axis(side=3,at=x2,labels=format(10^x2,digits=1))
  if(!is.null(categories)) axis(side=2,at=1:length(categories)*3-1,labels=categories,las=2,lty=0)
}


#' crossEnrichmentBarplot
#'
#' Plots Overlap, enrichment and significance between 2 lists of character vectors
#'
#' @param ll1 a list with DEGs vectors
#' @param ll2 a list with gene vectors
#' @param universe the gene vector to be used as the universe for computing enrichment and significance
#' @param noNegEnrichment Logical; whether to suppress negative enrichment (default FALSE)
#'
#' @return a list of matrices with overlap, enrichment and significance values
#'
#' @export
crossEnrichmentBarplot <- function(ll1, ll2, universe, noNegEnrichment=FALSE){
  ll1 <- lapply(ll1,y=universe,FUN=intersect)
  ll2 <- lapply(ll2,y=universe,FUN=intersect)
  enr <- matrix(0,nrow=length(ll1),ncol=length(ll2))
  pvals <- matrix(1,nrow=length(ll1),ncol=length(ll2))
  oo <- matrix(0,nrow=length(ll1),ncol=length(ll2))
  for(i in 1:length(ll1)){
    for(j in 1:length(ll2)){
      oo[i,j] <- length(intersect(ll1[[i]],ll2[[j]]))
      ee <- log2(.getEnrichment(ll1[[i]],ll2[[j]],universe=universe))
      if(ee<0){
        pvals[i,j] <- overlap.prob(ll1[[i]],ll2[[j]],universe,lower=T)
      }else{
        pvals[i,j] <- overlap.prob(ll1[[i]],ll2[[j]],universe)
      }
      if(is.infinite(ee)) ee <- 0
      enr[i,j] <- ee
    }
  }
  if(noNegEnrichment) enr[which(enr<0)] <- 0
  pvals <- -log10(pvals)
  row.names(enr) <- names(ll1)
  colnames(enr) <- names(ll2)
  row.names(pvals) <- names(ll1)
  colnames(pvals) <- names(ll2)
  cols <- getQualitativePalette(length(ll1))
  m <- list(overlaps=oo, enr=enr,pvals=pvals)
  autoLayout(3)
  barplot(m[[1]],beside=T,las=3,col=cols,ylab="Overlap")
  legend("topleft",bty="n",fill=cols, legend=names(ll1),xpd=T, cex=0.8)
  barplot(m[[2]],beside=T,las=3,col=cols,ylab="log2(enrichment)")
  barplot(m[[3]],beside=T,las=3,col=cols,ylab="-log10(p-value)")
  abline(h=-log10(0.05),lty="dashed")
  return(m)
}

#' plotGenesClusters
#'
#' Plots gene clusters dose-responses across EDC concentrations
#'
#' @param fcmat a vector of DEGs, or a named list of such vectors for each system
#' @param design a design data-frame
#' @param classification a named vector containing the cluster value for each gene (genes as names)
#' @param labels x-axis labels
#' @param showEachGene Logical; whether to plot each gene as thin transparent lines (default FALSE)
#' @param nbClusters Integer indicating the number of clusters. If omitted, will be estimated from the data.
#' @param minClusterSize Integer indicating the minimum cluster size (default 1).
#' @param spar Smoothing parameter (default 0.3).
#' @param interpolate Integer >= 0. Number of points to interpolate between actual points (default 3).
#' @param alpha alpha value for transparency (default 25)
#' @param q quantile to be plotted (default `c(0.125,0.875)`)
#' @param cols Colors.
#' @param ylim y-axis limits.
#' @param xlab x-axis label (default 'Concentration').
#' @param ylab y-axis label (default 'log2(foldchange)').
#' @param agFun Gene aggregation function (default `median`).
#' @param showNumber Logical; whether to show the number of genes in each cluster (default FALSE)
#' @param ... further arguments passed to the `plot` function.
#'
#' @export
plotGenesClusters <- function(fcmat, design, classification, labels=NULL, showEachGene=FALSE, nbClusters=NULL, minClusterSize=1, spar=0.3, interpolate=3, alpha=25, q=c(0.125,0.875), cols=NULL,ylim=NULL,xlab="Concentration",ylab="log2(foldchange)", agFun=median, showNumber=FALSE, ...){
  if(nrow(fcmat)==0){
    plot(1,1,col="white",xaxt="n",yaxt="n",xlab="",ylab="",bty="n",...)
    text(1,1,labels="No overlap")
    return()
  }
  fcmat <- aggregate(t(fcmat),by=list(x=design$concentration),FUN=mean)
  x <- fcmat[,1]
  row.names(fcmat) <- as.character(fcmat[,1])
  fcmat <- as.data.frame(t(fcmat[,-1,drop=FALSE]))
  nbConcs <- length(x)
  x <- 1:nbConcs
  
  cs1 <- table(classification)
  w <- which(cs1>=minClusterSize)
  cs <- 1+(10*cs1/sum(cs1))
  cs[which(cs1==1)] <- 1
  cs <- cs[w]
  cs1 <- cs1[w]
  co <- order(cs1, decreasing=T)
  cs <- cs[co]
  if(!is.null(nbClusters))	cs <- cs[1:min(length(cs,nbClusters))]
  cs <- rev(cs)
  classification <- classification[which(as.character(classification) %in% names(cs))]
  xx <- seq(from=1,to=nbConcs,length.out=nbConcs*interpolate)
  
  fcmat <- as.data.frame(t(apply(fcmat[names(classification),,drop=FALSE],1,x=x,spar=spar,interpolate=interpolate,force0=TRUE,FUN=function(y,x,xx,spar,interpolate,force0){
    .plSmooth(x,y,spar=spar,interpolate=interpolate,force0=force0)$y
  })))
  
  if(is.null(cols))	cols <- getQualitativePalette(length(cs))
  names(cols) <- names(cs)
  
  if(is.null(ylim)) ylim <- .getClustersYlim(fcmat, classification, q=q)
  plot(1,0,col="white",xlim=c(1,nbConcs), xlab=xlab, ylim=ylim,ylab=ylab,xaxt=ifelse(is.null(labels),"s","n"),bty="n",...)
  
  if(showEachGene){
    for(i in 1:nrow(fcmat))	lines(smooth.spline(xx,fcmat[i,],spar=0.1),col=maketrans(cols[classification[row.names(fcmat)[i]]],round(alpha/2)))
  }
  
  abline(h=0,lty="dashed")
  
  for(cc in names(cs)){
    e <- fcmat[names(classification)[which(classification==cc)],,drop=FALSE]
    if(!is.null(q) & nrow(e)>1){
      d <- .getClusterDisp(e, q=q)
      polygon(c(xx,rev(xx),xx[1]), c(d[,2],rev(d[,1]),d[1,2]), col=maketrans(cols[cc],alpha),border=NA,xpd=T)
    }
    #		lines(xx,apply(e,2,FUN=agFun),lwd=cs[cc],col=cols[cc],xpd=T)
    if(nrow(e)>1){
      lines(.plSmooth(xx,apply(e,2,FUN=agFun),interpolate=5,spar=spar),lwd=cs[cc],col=cols[cc],xpd=T)
    }else{
      lines(.plSmooth(xx,as.numeric(e),interpolate=5,spar=spar),lwd=cs[cc],col=cols[cc],xpd=T)
    }
    if(showNumber){
      nn <- nrow(e)
      if(nn==1)	nn <- row.names(e)
      text(max(xx),agFun(e[,ncol(e)]),labels=nn,pos=4,col=cols[cc],xpd=T,cex = 0.5)
    }
  }
  
  if(!is.null(labels)){
    axis(1, 1:nbConcs, labels)
  }
}

.getClusterDisp <- function(fcmat, q=c(0.125,0.875)){
  if(is.null(q))	return(t(apply(fcmat,2,FUN=median)))
  if(is.numeric(q) & length(q)==2) return(t(apply(fcmat,2,probs=q,FUN=quantile)))
  if(all(q=="sd")){
    m <- apply(fcmat,2,FUN=median)
    sds <- apply(fcmat,2,FUN=sd)
    return(cbind(m-sds,m+sds))
  }
  stop("Invalid `q`")
}

.getClustersYlim <- function(fcmat, classification, q=c(0.125,0.875)){
  range(sapply(unique(classification),fcmat=fcmat,cc=classification,q=q,FUN=function(x, fcmat, cc, q){
    range(as.matrix(.getClusterDisp(fcmat[names(cc)[which(cc==x)],], q=q)))
  }))
}


#' getConsClust
#'
#' kmeans clustering with consensus nb of clusters
#'
#' @export
getConsClust <- function(dat, maxClusters=10, method="kmeans", plot=TRUE){
  if(maxClusters==2) nb <- 2
  else{
    library(NbClust)
    dat <- dat[,which(apply(dat,2,FUN=function(x){ !all(x==0) }))]
    if(!plot) pdf(file=NULL)
    nb <- NbClust(dat, min.nc=2, max.nc=maxClusters, method=method)
    if(!plot) dev.off();
    tc <- table(nb$Best.nc[1,])
    nb <- max(as.numeric(names(tc)[which(tc==max(tc))]))
  }
  if(!(nb>1)){
    cc <- rep("0",nrow(dat))
    names(cc) <- row.names(dat)
    return(cc)
  }
  dat <- as.matrix(dat)
  dat[which(is.na(dat) | is.nan(dat))] <- 0
  cc <- kmeans(dat, nb)
  cc <- as.factor(cc$cluster)
  levels(cc) <- as.character(rank(-1*table(cc),ties.method="random"))
  cc <- as.character(cc)
  names(cc) <- row.names(dat)	
  return(cc)
}

#' filterGenes
#'
#' returns a filtered version of an expression matrix for the purpose of DEA.
#'
#' @param counts the expression matrix, with gene symbols as row names.
#' @param minCount the minimum number of reads (in at least `minSamples`) for a gene to be included.
#' @param minSamples the minimum number of samples with count > `minCount` in order for a gene to be included.
#' @param removeERCC logical, whether to remove ERCC spike-ins (default TRUE).
#' @param removeSmall logical, whether to remove SNORNAs and genes below 200nt (default TRUE).
#' @param removeRRNA logical, whether to remove ribosomal RNAs (default TRUE).
#' @param removeFusion logical, whether to remove gene fusions (default TRUE).
#' @param geneLengths a named vector of gene lengths, necessary if `removeSmall=T` and if the row names are not human gene symbols.
#'
#' @return a filtered matrix or data.frame (depending on `counts`).
#'
#' @export
filterGenes <- function(counts, minCount=20, minSamples=2, removeERCC=T, removeSmall=T, removeRRNA=T, removeFusion=T, geneLengths=NULL){
  if(is(counts, "SummarizedExperiment")){
    SE <- counts
    counts <- assay(SE)
  }else{
    SE <- NULL
  }
  nb1 <- nrow(counts)
  counts <- counts[which(apply(counts,1,minCount=minCount,minSamples=minSamples,function(x,minCount,minSamples){ sum(x>minCount) })>minSamples),]
  if(removeERCC)	counts <- counts[grep("^ERCC-",row.names(counts),invert=T),]
  if(removeSmall){
    if(is.null(geneLengths)){
      rm(geneLengths)
      #load("../Data/geneLengths.RData")
      if(!(length(intersect(row.names(counts),names(geneLengths))) > 10))	stop("The names of geneLengths do not correspond to the row names of `counts`. Use `geneLengths` to specific the genes' lengths for your species.")
    }else{
      if(class(geneLengths)=="character") load(geneLengths)
    }
    sg <- names(geneLengths)[which(geneLengths<200)]
    counts <- counts[which(!(row.names(counts) %in% sg)),]
    counts <- counts[grep("^SNOR",row.names(counts),invert=T),]
    counts <- counts[grep("^RNU",row.names(counts),invert=T),]
  }
  if(removeRRNA){
    counts <- counts[grep("^RNA5",row.names(counts),invert=T),]
    counts <- counts[grep("^RNA18",row.names(counts),invert=T),]
    counts <- counts[grep("^RNA28",row.names(counts),invert=T),]
    counts <- counts[grep("^RNA45",row.names(counts),invert=T),]
    counts <- counts[which(!(row.names(counts) %in% c("MT-RNR1","MT-RNR2"))),]
  }
  if(removeFusion) counts <- counts[grep("-",row.names(counts),invert=T,fixed=T),]
  message(paste(nb1-nrow(counts),"rows out of",nb1,"were removed."))
  if(!is.null(SE)) return(SE[row.names(counts),])
  counts
}



#' plPCA
#'
#' Runs PCA and plots the given components
#'
#' @param x The data.frame or matrix
#' @param normalize.genes Logical; whether to normalize rows, i.e. genes (default FALSE)
#' @param plot.components Numeric vector of length 2 or 3 indicating which components to plot. Default c(1,2)
#' @param plot.labels Logical; whether to plot labels (default TRUE)
#' @param tsplit Null, or a vector of length=ncol(x) indicating at most two groups on the basis of which to try to split the PCA. Ignored if length(plot.components) > 2
#' @param points.color Color of the points. Default "white" (not visible) if plot.labels is TRUE, or "black" otherwise. This can either be a single color, or a vector of length=ncol(x).
#' @param labels.color Color of the labels (default "black"). This can either be a single color, or a vector of length=ncol(x).
#' @param ... Any other argument to be passed to the plot function
#'
#' @return plots a figure
#'
#' @examples
#' # create random data
#' m <- matrix(rnorm(1000,10),nrow=100)
#' colnames(m) <- paste0("s",1:10)
#' # plot PCA
#' plPCA(m, plot.labels=T, col="white")
#'
#' @export
plPCA <- function(x, normalize.genes=F, plot.components=c(1,2), plot.labels=T, tsplit=NULL, doclusplot=0, points.color=NULL, labels.color="black", ...){
  if(!(length(plot.components) %in% c(2,3)))	stop("Plot components should be of size 2 or 3")
  if(is.null(points.color)) points.color <- ifelse(plot.labels,"white","black")
  x <- t(x[which(!(apply(x, 1, var)==0)),])
  if(normalize.genes)	x <- scale(x)
  pca <- prcomp(x)
  xlab <- paste("PC ", plot.components[1], " (", round(pca$sdev[plot.components[1]]/sum(pca$sdev)*100,0), "%)", sep="")
  ylab <- paste("PC ", plot.components[2], " (", round(pca$sdev[plot.components[2]]/sum(pca$sdev)*100,0), "%)", sep="")
  if(length(plot.components)==2){
    if(doclusplot>1){
      library(cluster)
      clusplot(pam(-1*pca$x[,plot.components], doclusplot), shade=TRUE, sub="", xlab=xlab, ylab=ylab, ...)
    }else{
      plot(pca$x[,plot.components[1]], pca$x[,plot.components[2]], xlab=xlab, ylab=ylab, col=points.color, ...)
    }
    if(plot.labels) text(pca$x[,plot.components[1]], pca$x[,plot.components[2]], labels=row.names(pca$x), cex=.8, font=2, col=labels.color)
    if(!is.null(tsplit)){
      if(length(tsplit)==nrow(x) & length(unique(tsplit))==2){
        divline <- supervisedPlotDivision(pca$x[,plot.components], tsplit)
        abline(a=divline$intercept, b=divline$slope, lty="dashed", lwd=2, col="grey")
      }else{
        warning("tsplit is not a vector of length=ncol(x) factorizable to two levels, and will therefore be ignored.")
      }
    }
  }else{
    library(scatterplot3d)
    zlab=paste("PC ", plot.components[3], " (", round(pca$sdev[plot.components[3]]/sum(pca$sdev)*100,0), "%)", sep="")
    p <- scatterplot3d(pca$x[,plot.components],type="h",lty.hplot="dashed", xlab=xlab, ylab=ylab, zlab=zlab, ...)
    if(plot.labels){
      label.coord <- p$xyz.convert(pca$x[,plot.components[1]], pca$x[,plot.components[2]], pca$x[,plot.components[3]])
      text(label.coord$x, label.coord$y, labels=row.names(pca$x),pos=4, cex=.7)
    }
  }
}

#' normalizes a dataset
#'
#' Wrapper to normalize a dataset, using the given method.
#'
#' @param dataset is a data.frame
#' @param method is a character string of either 'linear', 'housekeeping', 'quantile', or any method supported by edgeR's \code{\link[edgeR]{calcNormFactors}}.
#'
#' @return The normalized data.frame.
#'
#' @examples
#' donorm(matrix(1:12,nrow=4),"linear")
#'
#' @export
donorm <- function(dataset, method="TMM", returnCPM=FALSE){
  method <- match.arg(method, c("TMM","RLE","upperquartile","quantile","linear","geometric","none"))
  if(method=="none")	return(dataset)
  if(method=="geometric"){
    en <- exp(donorm(log(dataset+1),"linear"))
    return(en - min(1,min(en)))
  }	
  if(method=="linear"){
    nf <- getLinearNormalizers(dataset)
    return(t(t(dataset)*nf))
  }
  if(method=="quantile"){
    return(preprocessCore::normalize.quantiles(dataset))
  }
  # call edgeR
  library(edgeR)
  d <- DGEList(counts=dataset,group=colnames(dataset))
  d <- calcNormFactors(d, method=method)
  e <- cpm(d, normalized.lib.sizes=T)
  if(!returnCPM) e <- e*mean(d$samples$lib.size)/1000000
  return(as.data.frame(e))
}

#' getQualitativePalette
#'
#' Returns a qualitative color palette of given size
#'
#' @param nbcolors number of colors (from 1 to 22)
#'
#' @return A vector of colors
#'
#' @export
getQualitativePalette <- function(nbcolors){
  # based on Paul Tol's colors
  switch(as.character(nbcolors),
         "1"=c("#4477AA"),
         "2"=c("#4477AA", "#CC6677"),
         "3"=c("#4477AA", "#DDCC77", "#CC6677"),
         "4"=c("#4477AA", "#117733", "#DDCC77", "#CC6677"),
         "5"=c("#332288", "#88CCEE", "#117733", "#DDCC77", "#CC6677"),
         "6"=c("#332288", "#88CCEE", "#117733", "#DDCC77", "#CC6677","#AA4499"),
         "7"=c("#332288", "#88CCEE", "#44AA99", "#117733", "#DDCC77", "#CC6677","#AA4499"),
         "8"=c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#CC6677","#AA4499"),
         "9"=c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#CC6677", "#882255", "#AA4499"),
         "10"=c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#661100", "#CC6677", "#882255", "#AA4499"),
         "11"=c("#332288", "#6699CC", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#661100", "#CC6677", "#882255", "#AA4499"),
         "12"=c("#332288", "#6699CC", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#661100", "#CC6677", "#AA4466", "#882255", "#AA4499"),
         "13"=c("#882E72", "#B178A6", "#1965B0", "#5289C7", "#7BAFDE", "#4EB265", "#90C987", "#CAE0AB", "#F7EE55", "#F6C141", "#F1932D", "#E8601C", "#DC050C"),
         "14"=c("#882E72", "#B178A6", "#D6C1DE", "#1965B0", "#5289C7", "#7BAFDE", "#4EB265", "#90C987", "#CAE0AB", "#F7EE55", "#F6C141", "#F1932D", "#E8601C", "#DC050C"),
         "15"=c("#114477", "#4477AA", "#77AADD", "#117755", "#44AA88", "#99CCBB", "#777711", "#AAAA44", "#DDDD77", "#771111", "#AA4444", "#DD7777", "#771144", "#AA4477", "#DD77AA"),
         "16"=c("#114477", "#4477AA", "#77AADD", "#117755", "#44AA88", "#99CCBB", "#777711", "#AAAA44", "#DDDD77", "#771111", "#AA4444", "#DD7777", "#771144", "#AA4477", "#DD77AA", "black"),
         "17"=c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788"),
         "18"=c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788"),
         "19"=c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788", "black"),
         "20"= c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788"),
         "21"= c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788"),
         "22"= c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788", "black"),
         stop("Max 22 colors")
  )
}

#' maketrans
#'
#' Makes a color transparent
#'
#' @param tcol color
#' @param alpha alpha level (1 to 255, default 100)
#'
#' @return transparent color code
#'
#' @export
maketrans <- function (tcol, alpha = 100){
  c <- col2rgb(tcol)
  rgb(c["red", 1][[1]], c["green", 1][[1]], c["blue", 1][[1]], 
      alpha, maxColorValue = 255)
}

#' autoLayout
#'
#' Creates a layout for the given number of panels
#'
#' @param nb number of panels
#'
#' @return nothing, but creates a layout
#'
#' @export
autoLayout <- function(nb,byrow=FALSE){
  nc <- ceiling(sqrt(nb))
  nr <- ceiling(nb/nc)
  layout(matrix(1:(nr*nc),nrow=nr,byrow=byrow))
}


#' overlap.prob
#'
#' Calculates the probability of the observed overlap between two sets on the basis of the hypergeometric distribution.
#'
#' @param set1 a character vector
#' @param set2 a character vector
#' @param universe either a character vector containing the universe/background (recommended), or an integer indicating the size of the universe
#' @param lower logical; whether the lower tail should be considered (default FALSE)
#'
#' @return the probability of the overlap
#'
#' @examples
#' overlap.prob( c("A","B","C","D"), c("B","D","G"), LETTERS )
#'
#' @export
overlap.prob <- function(set1,set2,universe,lower=F){
  set1 <- as.character(set1)
  set2 <- as.character(set2)
  if(class(universe)=="character"){
    set1 <- intersect(set1,universe)
    set2 <- intersect(set2,universe)
    universe <- length(unique(universe))
  }
  set1 <- unique(set1)
  set2 <- unique(set2)
  ov <- sum(set1 %in% set2)
  phyper(max(0,ov-1), length(set1), universe-length(set1), length(set2), lower.tail=lower)
}


#' counts2fpkm
#'
#' Returns fpkm values from count matrix
#'
#' @param counts a count matrix (with samples as columns and genes as rows)
#' @param mean.frag.size mean fragment size (default 220) for effective size estimate
#'
#' @return a fpkm matrix
#'
#' @export
counts2fpkm <- function(counts, mean.frag.size = 220){
  data("geneLengths")
  ii <- intersect(names(geneLengths), row.names(counts))
  fpkm <- apply(counts[ii,], 2, gl=geneLengths[ii], mean.frag.size=mean.frag.size, FUN=function(x, gl, mean.frag.size){
    x <- as.numeric(x)
    gl <- gl - mean.frag.size
    x[which(gl < 0)] <- 0
    gl[which(gl < mean.frag.size)] <- mean.frag.size
    10^9 * x/(gl * sum(x))
  })
  row.names(fpkm) <- ii
  colnames(fpkm) <- colnames(counts)
  return(fpkm)
}


#' sdStrips
#'
#' Plots stripchart and standard deviations
#'
#' @param a A matrix of expression values, with samples as columns and genes as rows.
#' @param sampletypes A vector indicating the type of each samples in a.
#' @param colors A vector of colors, with names corresponding to the levels of `sampletypes`.
#' @param legendLocation Passed to `legend`, indicates the location of the legend (default 'topleft').
#' @param colerror The color of the error bars.
#' @param colwidth The width of the whiskers.
#' @param ... arguments passed to the stripchart function
#'
#' @export
sdStrips <- function(a, sampletypes, colors=c(Fetal="#D6D6D6", Ngn2="#BAA5E0", Org="#A8C0C2", iPSC="lightsalmon"), legendLocation="topleft", colerror=rgb(0,0,0,1), errwidth=0.5, ...){
  a <- t(a)
  dat <- as.numeric(as.matrix(a))
  st <- as.factor(rep(sampletypes, ncol(a)))
  gene <- as.factor(rep(colnames(a),each=nrow(a)))
  stripchart(dat~st+gene, vertical=T, pch=16, method="jitter", jitter=0.3, cex=0.6, col=colors[as.character(unique(st))], xaxt="n", ylab="log(FPKM)", bty="n")
  axis(1, at=(1:ncol(a))*length(unique(st))-length(unique(st))/2, labels=colnames(a), cex.axis=0.8)
  
  m <- aggregate(a,by=list(sampletypes),FUN=mean)
  row.names(m) <- m[,1]
  m <- as.numeric(as.matrix(m[levels(st),levels(gene)]))
  s <- aggregate(a,by=list(sampletypes),FUN=sd)
  row.names(s) <- s[,1]
  s <- as.numeric(as.matrix(s[levels(st),levels(gene)]))
  
  x <- 1:(ncol(a)*length(unique(st)))
  segments(x0=x-errwidth/2,x1=x+errwidth/2,y0=m,y1=m, lty=1, col=colerror, lwd=2)
  arrows(x0=x,x1=x,y0=m-s,y1=m+s, length=0.05, lwd=1, lty="83", code=3, col=colerror, angle=90)
  
  legend(legendLocation, fill=colors, legend=names(colors), bty="n")
  
}



#' estimateTFactivity
#'
#' @param en The normalized expression matrix (preferentially log-transformed, or logFCs)
#' @param dea The DEA results table
#' @param net The network (list of tables with target genes and scores for each TF)
#' @param CD The colData of the experiment
#' @param form The formula for differential TF activity (`dta`). Terms must be columns of 
#' CD. If omitted, dta will not be performed, but TF activity will nevertheless be 
#' calculated and camera analysis will be performed.
#' @param test.term The coefficient(s) to test for dta, can be either terms of the formula
#' or columns of the model matrix, but not a combination of the two.
#' @param as.is Logical; whether to use the expression values as-is; default FALSE, which 
#' means that gene z-scores are calculated. Set to TRUE if `en` are log-foldchanges.
#'
#' @return A list.
estimateTFactivity <- function(en, dea, net, CD=NULL, form=NULL, test.term=NULL, as.is=FALSE){
  if(!as.is) en <- t(scale(t(en)))
  
  gsets <- lapply(net, FUN=function(x){ as.character(x$Target_Gene) })
  cam <- cameraWrapper(dea, gsets)
  row.names(cam) <- sapply(strsplit(as.character(cam$term),":"),FUN=function(x) x[3])
  
  tfa <- t(sapply(net, FUN=function(x){
    row.names(x) <- x$Target_Gene
    i <- intersect(row.names(x), row.names(en))
    colSums(en[i,,drop=FALSE]*sqrt(as.numeric(x[i,"Edge_Weight"])))
  }))
  if(is.null(CD) || is.null(form) || is.null(test.term)){
    return(list( tfa=tfa, camera=cam ))
  }
  
  df <- as.data.frame(CD)[,labels(terms(form)),drop=FALSE]
  df$y <- 1
  form <- reformulate(labels(terms(form)), response="y" )
  isTestFullterm <- all(test.term %in% colnames(df))
  if(!isTestFullterm && !all(test.term %in% colnames(model.matrix(form, data=df)))){
    stop("`test.term` should either be (a) term(s) of the formula or column(s) of the model matrix.")
  }
  a <- t(apply(tfa, 1, FUN=function(x){
    df$y <- x
    tryCatch({
      mod <- lm(form, data=df)
      if(isTestFullterm){
        v <- c(beta="Estimate")
      }else{
        v <- c(beta="Estimate",pval="Pr(>|t|)")
      }
      tt2 <- unlist(sapply(test.term, FUN=function(x){
        paste0(x,levels(as.factor(df[[x]]))[-1])
      }))
      co <- coef(summary(mod))[tt2, v, drop=F]
      co2 <- as.numeric(t(co))
      names(co2) <- paste(rep(row.names(co),each=length(v)), names(v), sep=".")
      if(isTestFullterm){
        co2 <- c(co2, PValue=drop1(mod, test.term, test="F")[test.term,'Pr(>F)'])
      }
      co2
    }, error=function(e){ 
      warning(e);
      if(isTestFullterm) return(rep(NA, length(test.term)+1))
      rep(NA, length(test.term)*2)
    })
  }))
  a <- as.data.frame(a)
  if(isTestFullterm){
    a <- a[order(a$PValue),]
    a$FDR <- p.adjust(a$PValue)
  }else{
    wp <- grep("pval",colnames(a))
    a <- a[order(rowMeans(a[,wp,drop=FALSE])),]
    if(length(wp)==1){
      a$FDR <- p.adjust(a[,wp])
    }
  }
  return( list( tfa=tfa, dtfa=a, camera=cam ) )
}

vipsummary  <- function(mrs, dea, n=50){
  mrs2 <- summary(mrs,n)
  mrs2 <- cbind(mrs2,dea[row.names(mrs2),])
  DT::datatable(mrs2)
  DT::datatable(mrs2[which(!is.na(rowMeans(mrs2[,grep("logFC\\.",colnames(mrs2))]))),-(ncol(mrs2)-1)])
}




#' svacor
#'
#' A wrapper around SVA-based correction
#'
#' @param SE An object of class `SummarizedExperiment`. Alternatively, a matrix can be 
#' used, but many options will not be supported.
#' @param form The formula of the differential expression model
#' @param form0 An optional formula for the null model
#' @param mm If `form=NULL`, the model.matrix
#' @param mm0 An optional null model.matrix.
#' @param regressOutNull Logical; whether to regress out the variables of `form0` (default TRUE)
#' @param seqb Whether to use the `svaseq` (default FALSE; uses vst+sva)
#' @param ... Any other param passed to the sva command.
#'
#' @return A list with the slots: 
#' * `sv`: a table of the surrogate variables
#' * `cor`: the corrected data (for plotting)
#' * `mm`: the model.matrix containing, in addition to the specified experimental 
#' variables, all detected surrogate variables.
#' 
#' @export
svacor <- function(SE, form=NULL, form0=~1, mm=NULL, mm0=NULL, regressOutNull=TRUE, seqb=FALSE, ...){
  library(sva)
  library(SummarizedExperiment)
  if( (is.null(form) && is.null(mm)) ||
      (!is.null(form) && !is.null(mm)) ) stop("Only one of `form` or `mm` should be given.")
  if(is.null(mm)){
    if(is(SE,"SummarizedExperiment")){
      CD <- as.data.frame(colData(SE))
      mm <- model.matrix(form, data=CD)
    }else{
      stop("If `form` is used, `SE` should be a SummarizedExperiment.")
    }
  }
  if(is.null(mm0)){
    if(is.null(form0)){
      mm0 <- mm[,1,drop=F]
    }else{
      mm0 <- model.matrix(form0, data=CD)
    }
  }
  if(is(SE,"SummarizedExperiment")){
    if(!is.null(form) && !seqb){
      message("Using variance-stabilizing transformation")
      suppressPackageStartupMessages(library(DESeq2))
      dds <- DESeqDataSetFromMatrix(round(assay(SE)), as.data.frame(colData(SE)), form)
      dds <- estimateSizeFactors(dds)
      en <- as.matrix(assay(vst(dds, blind=FALSE)))
    }else{
      en <- as.matrix(donorm(assay(SE)))
      if(!seqb) en <- log1p(en)
    }
  }else{
    en <- as.matrix(SE)
  }
  if(seqb){
    sv <- svaseq(en, mm, mm0, ...)
  }else{
    sv <- sva(en, mm, mm0, ...)
  }
  if(sv$n.sv==0) return(NULL)
  colnames(sv$sv) <- paste0("SV",1:ncol(sv$sv))
  X <- cbind(mm, sv$sv)
  H <- solve(t(X)%*%X)%*%t(X)
  b <- (H%*%t(en))
  if(regressOutNull){
    cn <- setdiff(colnames(X),setdiff(colnames(mm), colnames(mm0)))  
  }else{
    cn <- setdiff(colnames(X),colnames(mm))
  }
  cn <- setdiff(cn, "(Intercept)")
  encor <- en - t(as.matrix(X[,cn]) %*% b[cn,])
  sv <- sv$sv
  mm2 <- cbind(mm[,1,drop=F],sv,mm[,-1,drop=F])
  return(list(sv=sv, cor=encor, mm=mm2))
}

getMsigSets <- function(){
  library(msigdbr)
  m <- msigdbr(species="Homo sapiens")
  m <- m[which(m$gs_cat %in% c("H","C2","C5")),]
  m$name2 <- paste0(m$gs_cat,":",m$gs_subcat, ":", m$gs_name)
  split(m$gene_symbol, m$name2)
}

ancols <- list(EXPO2=c("CNT"="#0000FF", "DMSO"="#2900D5", "0.1X"="#5500AA", "1X"="#7E0080", "10X"="#AA0054", "100X"="#D3002B", "1000X"="#FF0000", "BPA0.04X"="#117733", "BPA1X"="#999933", "Triclosan3nM"="#117744", "Triclosan100nM"="#999944", "VITC"="lightgrey", "T3"="yellow",  "T3MixN"="orange"))

options("SEtools_def_anno_colors"=ancols)


#to be updated
goseq.enrichment <- function (allGenes, deGenes, gotype = c("GO:BP", "GO:MF", "GO:CC"), 
                              cutoff = 0.1, cutoff.onFDR = TRUE, org = "hg19", minCatSize = 10, 
                              maxCatSize = 1000, maxResults = 200) 
{
  suppressPackageStartupMessages(library("goseq"))
  suppressPackageStartupMessages(library("GO.db"))
  suppressPackageStartupMessages(library("org.Hs.eg.db"))
  suppressPackageStartupMessages(library("AnnotationDbi"))
  gotype <- match.arg(gotype, c("GO:CC", "GO:MF", "GO:BP"), 
                      several.ok = T)
  org <- match.arg(org, c("hg19", "mm9"))
  emptyRes <- data.frame(term = vector(mode = "character", 
                                       length = 0), enrichment = vector(mode = "numeric", length = 0), 
                         PValue = vector(mode = "numeric", length = 0), FDR = vector(mode = "numeric", 
                                                                                     length = 0), genes = vector(mode = "character", length = 0))
  if (length(deGenes) < 3) 
    return(emptyRes)
  if (all(is.na(allGenes))) {
    message("Using all (annotated) genes as a background.")
    if (org == "hg19") {
      data("hg38.refseq")
    }
    else {
      data("mm10.refseq")
    }
    allGenes <- unique(as.character(gtf$gene_id))
  }
  genes <- as.integer(allGenes %in% as.character(deGenes))
  names(genes) <- allGenes
  pwf = nullp(genes, org, "geneSymbol", plot.fit = F)
  go.all <- suppressMessages(goseq(pwf, org, "geneSymbol", 
                                   test.cats = gotype))
  go.all <- go.all[which(go.all$numInCat >= minCatSize & go.all$numInCat < 
                           maxCatSize), ]
  go.all$FDR <- p.adjust(go.all$over_represented_pvalue, method = "BH")
  if (cutoff.onFDR) {
    go.en <- go.all[go.all$FDR < cutoff, ]
  }
  else {
    go.en <- go.all[go.all$over_represented_pvalue < cutoff, 
                    ]
  }
  go.en <- go.en[order(go.en$over_represented_pvalue), ]
  go.en$Term <- sapply(go.en$category, FUN = function(x) {
    if (is.null(x) | is.na(x) | x == "") {
      return("")
    }
    else {
      return(Term(x))
    }
  })
  go.en$Expected <- length(deGenes) * (go.en$numInCat/length(allGenes))
  go.en$Enrichment <- go.en$numDEInCat/go.en$Expected
  go.en <- go.en[, c("category", "Term", "numInCat", "numDEInCat", 
                     "Enrichment", "over_represented_pvalue", "FDR")]
  names(go.en) <- c("GO.ID", "Term", "Annotated", "Significant", 
                    "Enrichment", "PValue", "FDR")
  if (!(nrow(go.en) >= 0)) 
    return(emptyRes)
  if (nrow(go.en) > maxResults) 
    go.en <- go.en[1:maxResults, ]
  species <- ifelse(substr(org, 0, 2) == "hg", "Hs", "Mm")
  db <- paste0("org.", species, ".eg")
  library(package = paste0(db, ".db"), character.only = T)
  eg <- AnnotationDbi::mget(as.character(go.en$GO.ID), get(paste0(db, 
                                                                  "GO2ALLEGS")), ifnotfound = NA)
  go.en$genes <- sapply(eg, db = db, sig = as.character(deGenes), 
                        FUN = function(x, db, sig) {
                          x <- x[which(!is.na(x))]
                          if (length(x) == 0) 
                            return("")
                          x <- unique(as.character(unlist(AnnotationDbi::mget(as.character(x), 
                                                                              get(paste0(db, "SYMBOL"))))))
                          x <- intersect(x, sig)
                          paste(sort(x), collapse = ", ")
                        })
  return(go.en)
}



#to be updated
cameraWrapper <- function(dea, gsets=NULL, addSets=NULL, addDEgenes=TRUE, dea.thres=0.05, minG=5, reportMax=500){
  library(limma)
  if(is.null(gsets)){
    gsets <- getMsigSets()
  }
  if(!is.null(addSets)) gsets <- c(gsets, addSets)
  dea <- homogenizeDEAresults(dea)
  dea <- dea[which(!is.na(dea$PValue)),]
  gsets <- lapply(gsets, y=row.names(dea), FUN=intersect)
  gsets <- gsets[which(sapply(gsets,length)>=minG)]
  if("stat" %in% colnames(dea)){
    gs <- dea$stat
  }else{
    gs <- sign(dea$logFC)*-log10(dea$PValue)
  }
  names(gs) <- row.names(dea)
  gs <- gs[which(!is.na(gs))]
  Cres <- cameraPR(gs, gsets)
  w <- which(Cres$PValue < 0.01)
  if(length(w)==0) return(NULL)
  w <- w[1:min(reportMax, length(w))]
  Cres <- Cres[w,,drop=F]
  if(nrow(Cres)==0) return(NULL)
  dg <- t(sapply(gsets[row.names(Cres)], sig=row.names(dea)[which(dea$FDR<dea.thres)], FUN=function(x,sig){
    c(sum(x %in% sig), paste(intersect(x,sig),collapse=", "))
  }))
  Cres <- cbind(nDE=as.numeric(dg[,1]), Cres, genes=dg[,2])
  df <- t(sapply(strsplit(row.names(Cres),":",fixed=T), FUN=function(x){
    if(length(x)==2) return(c("custom",x[[2]],paste(x[3:length(x)],collapse=":")))
    if(length(x)==1) return(c("custom",NA,paste(x[3:length(x)],collapse=":")))
    c(x[[1]],x[[2]],paste(x[3:length(x)],collapse=":"))
  }))
  colnames(df) <- c("collection","subcat","term")
  Cres <- as.data.frame(cbind(df,Cres),stringsAsFactors=F)
  row.names(Cres) <- NULL
  Cres
}



#to be updated
homogenizeDEAresults <- function(x){
  if(is(x,"list")) return(lapply(x,homogenizeDEAresults))
  x <- as.data.frame(x)
  colnames(x)[which(colnames(x) %in% c("FDR","padj","adj.P.Val"))] <- "FDR"
  colnames(x)[which(colnames(x) %in% c("P.Value","pvalue","PValue"))] <- "PValue"
  colnames(x)[which(colnames(x) %in% c("log2FoldChange","logFC"))] <- "logFC"
  return(x[order(x$FDR),])
}

byheatmap <- function(x, scale = "none", ...){                                                                                                                                                                                               
  require("pheatmap")
  scale <- match.arg(scale, c("none", "row", "column"))
  x <- x[which(apply(x, 1, FUN = function(y){
    !all(is.na(y))
  })), which(apply(x, 2, FUN = function(y) {
    !all(is.na(y))
  }))]
  pheatmap(x, scale = scale, color = colorRampPalette(c("blue", "black", "yellow"))(29), border_color = NA, ...)                                                                                                                                        
}

