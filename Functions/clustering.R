#' prepLFC4clust
#'
#' Prepares a logFC matrix for features clustering by capping foldchanges, 
#' penalizing sign changes and optionally aggregating by condition
#' 
#' @param lfc A matrix of log-foldchanges, with samples as columns
#' @param groups An optional vector of groups of length `ncol(lfc)`
#' @param signIncr An increment added to penalize lfcs in opposite directions
#' @param minLFC Minimum abs(logFC) to be considered
#' @param maxLFC Maximum abs(logFC) cap
#'
#' @return A modified logFC matrix
#' @export
prepLFC4clust <- function( lfc, groups=NULL, signIncr=0.25, minLFC=0.15,
                        maxLFC=quantile(abs(as.numeric(lfc)),c(0.95),na.rm=T),
                        scale=TRUE, center=FALSE ){
  rn <- row.names(lfc)
  if(!is.null(maxLFC)){
    w <- abs(lfc)>maxLFC
    lfc[w] <- sign(lfc[w])*maxLFC
  }
  lfc <- sign(lfc)*sqrt(abs(lfc))
  if(scale || center) lfc <- t(base::scale(t(lfc), center=center, scale=scale))
  if(!is.null(groups)) lfc <- sapply( split(seq_len(ncol(lfc)), groups), 
                                      FUN=function(x)
                                            rowMedians(lfc[,x,drop=FALSE])
                                      )
  if(!is.null(minLFC)) lfc[abs(lfc)<minLFC] <- 0
  if(sum(rowSums(lfc!=0)==0))
    warning("Some rows are entirely zeros after capping.")
  if(!is.null(signIncr)) lfc <- lfc+sign(lfc)*signIncr
  row.names(lfc) <- rn
  lfc
}

#' prepLFC4clust2
#'
#' Prepares a logFC matrix for features clustering by capping foldchanges, 
#' penalizing sign changes and optionally aggregating by condition
#' 
#' @param lfc A matrix of log-foldchanges, with samples as columns
#' @param groups An optional vector of groups of length `ncol(lfc)`
#' @param signIncr An increment added to penalize lfcs in opposite directions
#' @param minLFC Minimum abs(logFC) to be considered
#' @param maxLFC Maximum abs(logFC) cap
#'
#' @return A modified logFC matrix
#' @export
prepLFC4clust2 <- function( lfc, groups=NULL, minLFC=0.15,
                           maxLFC=quantile(abs(as.numeric(lfc)),c(0.98),na.rm=T) ){
  rn <- row.names(lfc)
  if(!is.null(maxLFC)){
    w <- abs(lfc)>maxLFC
    lfc[w] <- sign(lfc[w])*maxLFC
  }
  if(!is.null(groups)) lfc <- sapply( split(seq_len(ncol(lfc)), groups), 
                                      FUN=function(x)
                                        rowMedians(lfc[,x,drop=FALSE])
  )
  z <- t(scale(t(lfc)))
  if(!is.null(minLFC)) lfc[abs(lfc)<minLFC] <- 0
  if(sum(rowSums(lfc!=0)==0))
    warning("Some rows are entirely zeros after capping.")
  row.names(lfc) <- rn
  lfcs <- sign(lfc)
  lfcs[abs(lfc)<minLFC] <- 0
  lfcr <- t(apply(lfc,1,FUN=rank))/ncol(lfc)
  colnames(lfcs) <- paste0("sign_",colnames(lfcs))
  colnames(lfcr) <- paste0("rank_",colnames(lfcr))
  colnames(z) <- paste0("zscore_",colnames(z))
  cbind(lfc,lfcs,lfcr,z)
}


#' clusterWrapper
#'
#' A wrapper for pam/kmeans clustering across values of k
#'
#' @param d Data to cluster (passed directly to the clustering function)
#' @param k.test Range of number of clusters to test (default 2:10)
#' @param method either 'kmeans' (default) or 'pam'
#' @param ... Passed to the clustering function.
#'
#' @return Plots diagnostic graph and returns a list of clustering results.
#' @export
clusterWrapper <- function(d, k.test=2:10, method=c("kmeans","pam"), ...){
  library(cluster)
  library(stats)
  switch( match.arg(method),
          kmeans=.kmeansWrapper(d, k.test=k.test, ...),
          pam=.pamWrapper(d, k.test=k.test, ...) )
}

# see clusterWrapper()
.pamWrapper <- function(d, k.test=2:10, method=c("pam","kmeans"), ...){
  res <- lapply(k.test, FUN=function(x){ pam(d, x, ...) })
  names(res) <- k.test
  sw <- sapply(res, FUN=function(x){
    c( sil.avg.width=x$silinfo$avg.width,
       separation=mean(x$clusinfo[,"separation"]),
       max_diss=max(x$clusinfo[,"max_diss"],na.rm=TRUE)
       )
  })
  csizes <- lapply(res, FUN=function(x) table(x$clustering))
  csizes <- sapply(csizes, y=rep(NA_integer_, max(sapply(csizes,length))),
                   FUN=function(x,y){ y[seq_along(x)] <- x; y } )
  layout(matrix(1:4,nrow=2))
  boxplot(as.data.frame(csizes), ylab="Cluster size", xlab="Number of clusters")
  plot(k.test, sw[1,], xlab="Number of clusters", ylab="Average cluster silhouette width", main="Silhouette", type="b", lwd=2, pch=16, col="red")
  plot(k.test, sw[2,], xlab="Number of clusters", ylab="Separation", main="Separation", type="b", lwd=2, pch=16, col="blue")
  plot(k.test, sw[3,], xlab="Number of clusters", ylab="Max distance to medoid", main="Max dissimilarity", type="b", lwd=2, pch=16, col="red")
  res
}

# see clusterWrapper()
.kmeansWrapper <- function(d, k.test=2:10, ...){
  ds <- dist(d)
  res <- lapply(k.test, FUN=function(x){ kmeans(d, x, ...) })
  names(res) <- k.test
  sw <- lapply( res, FUN=function(x) cluster::silhouette(x$cluster, ds))
  ds <- as.matrix(ds)
  sw <- list( sil.avg.width=sapply( sw, FUN=function(x) summary(x)$avg.width ),
              separation=sapply(res, ds=ds, FUN=.getSeparation),
              var.explained=sapply( res, FUN=function(x) x$betweenss/x$totss ) )
  csizes <- lapply(res, FUN=function(x) table(x$cluster))
  csizes <- sapply(csizes, y=rep(NA_integer_, max(sapply(csizes,length))),
                   FUN=function(x,y){ y[seq_along(x)] <- x; y } )
  layout(matrix(1:4,nrow=2))
  boxplot(as.data.frame(csizes), ylab="Cluster size", xlab="Number of clusters")
  plot(k.test, sw$sil.avg.width, xlab="Number of clusters", ylab="Average cluster silhouette width", main="Silhouette", type="b", lwd=2, pch=16, col="red")
  plot(k.test, sw$separation, xlab="Number of clusters", ylab="Separation", main="Separation", type="b", lwd=2, pch=16, col="blue")
  plot(k.test, sw$var.explained, xlab="Number of clusters", ylab="within.ss/tot.ss", main="Variance explained", type="b", lwd=2, pch=16, col="blue")
  res
}

# returns separation (minimum distance between non-siblings)
.getSeparation <- function( cl, # vector of clusters
                            ds  # distance matrix
                          ){
  if(is(cl,"kmeans")) cl <- cl$cluster
  ds <- as.matrix(ds)
  cl <- cl[row.names(ds)]
  sameCl <- sapply(seq_along(cl), FUN=function(i) cl[i]==cl)
  ds[sameCl] <- NA
  min(ds,na.rm=TRUE)
}

#' plotClusterLFCs
#' 
#' Plots per-cluster log-foldchanges/z-scores curves
#'
#' @param lfc The matrix of values to plot
#' @param clust A vector of cluster identities
#' @param summary The cluster summary to use
#' @param return.df Whether to return the data.frame rather than plotting (default FALSe)
#'
#' @return A ggplot, or a data.frame if `return.df==TRUE`
#' @export
plotClusterLFCs <- function(lfc, clust, summary=c("median","mean","smooth"), return.df=FALSE){
  summary <- match.arg(summary)
  lfc <- lfc[,grep("^zscore_|^sign_|^rank_", colnames(lfc),invert=TRUE)]
  d <- data.frame(cluster=rep(clust, ncol(lfc)), 
                  feature=factor(rep(row.names(lfc),ncol(lfc))),
                  class=factor(rep(colnames(lfc),each=nrow(lfc)), colnames(lfc)),
                  val=as.numeric(lfc))
  if(return.df) return(d)
  require(ggplot2)
  p <- ggplot(d, aes(class, val, group=feature)) + geom_line(alpha=0.3) + 
    facet_wrap(~cluster, scale="free_y") + ylab("") + xlab("")
  if(summary=="smooth"){
    p <- p + geom_smooth(formula=y~x, method="loess", aes(group=cluster), 
                         se=FALSE, size=1.5)
  }else{
    p <- p + stat_summary(fun.y="mean",geom="line",size=1.5, colour="blue", 
                          aes(group=cluster))
  }
  p
}