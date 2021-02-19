#' clusterEnrichment
#'
#' @param clusters A named list of sets of genes of interest
#' @param sets A named list of reference genesets (must use the same gene
#' identifiers as `clusters`)
#' @param universe The optional universe (by default, the union of clusters is
#' used)
#' @param minSetSize The minimum set size (in universe) to consider
#' @param threshold The p-value threshold
#' @param family The model family for deviance calculation.
#'
#' @return A data.frame
clusterEnrichment <- function(clusters, sets, universe=NULL, minSetSize=10, threshold=0.05, family=c("binomial","poisson")){
    family <- match.arg(family)
    if(is.null(universe)) universe <- unlist(clusters)
    sets <- lapply(sets, y=universe, FUN=intersect)
    clusters <- lapply(clusters, y=universe, FUN=intersect)
    sets <- sets[sapply(sets,length)>=minSetSize]
    sl <- sapply(sets, length)
    cs <- sapply(clusters, length)
    universe <- length(universe)
    b <- sapply(clusters, FUN=function(x){
        sapply(sets, FUN=function(y){
            length(intersect(x,y))
        })
    })
    dev <- sapply(1:nrow(b),sz=cs,mod=family,FUN=function(g, sz, mod){
        x <- b[g,]
        expected <- sz*sl[g]/universe
        enr <- log1p(x)-log1p(expected)
        p=sum(x)/sum(sz)
        # calculation of
        if(mod=="binomial"){
            term1<-sum(x*log(x/(sz*p)), na.rm=TRUE)
            nx<-sz-x
            term2<-sum(nx*log(nx/(sz*(1-p))), na.rm=TRUE)
            dev <- 2*(term1+term2)
        }else{
            dev <- 2*sum(x*log(x/(sz*p)),na.rm=TRUE)-2*sum(x-sz*p)
        }
        pval <- pchisq(dev, length(x)-1, lower.tail=FALSE)
        c(deviance=dev, p.value=pval, FDR=NA_real_, enr)
    })
    dev <- as.data.frame(t(dev))
    dev$FDR <- p.adjust(dev$p.value, method="fdr")
    colnames(dev)[4:ncol(dev)] <- paste("enrichment",colnames(dev)[4:ncol(dev)],sep=".")
    row.names(dev) <- row.names(b)
    dev <- dev[dev$p.value<threshold,]
    dev[order(dev$p.value),]
}

plotClusterEnrichment <- function(dev, k=5, top=NULL,
                                  color=colorRampPalette(c("blue","black", "yellow"))(50),
                                  ..., returnMatrix=FALSE){
    dev <- dev[order(dev$p.value),]
    if(is.null(top)){
        d <- 1-cor(t(dev[,grep("enrichment",colnames(dev))]))
        bg <- cluster_louvain(knn.graph(d,k))
        dev$cluster <- NA_integer_
        dev[bg$names,"cluster"] <- bg$membership
        s1 <- row.names(dev)[unique(apply(dev[,grep("enrichment",colnames(dev))],2,which.max))]
        s2 <- row.names(dev)[!duplicated(dev$cluster)]
        s <- union(s1,s2)
    }else{
        s <- row.names(dev)[seq_len(min(top,nrow(dev)))]
    }
    deve <- dev[s,grep("enrichment",colnames(dev))]
    colnames(deve) <- gsub("enrichment\\.","",colnames(deve))
    if(returnMatrix) return(deve)
    pheatmap(deve, color=color, border_color = NA, ...)
}




#' goclassify
#'
#' @param cl A named vector indicating the cluster for each gene, or a named list
#' containing the genes for each cluster
#' @param go A named list of genes for each geneset
#' @param minSize The minimum size for a geneset to be considered
#' @param max_depth passed to xgboost
#' @param N the number of sets to select
#' @param nrounds passed to xgboost
#' @param by how to select the features
#' @param onlyAnnotated Logical; whether to restrict elements of `cl` to those that are
#' present in some set of `go`
#'
#' @return a list
#'
#' @examples
#' go <- msigdbr::msigdbr("Mus musculus", "C5")
#' go <- go[which(go$gs_subcat=="BP"),]
#' go <- split(go$gene_symbol, go$gs_name)
#' cl <- lapply(c(c1=50, c2=100), FUN=function(x) sample(unique(unlist(go)),x))
#' res <- goclassify(cl, go)
goclassify <- function( cl, go, minSize=5, max_depth=5, N=30, nrounds=50,
                        by=c("Frequency","Gain","Coverage"), onlyAnnotated=TRUE,
                        onlyInSets=TRUE){
    if(!is.list(cl)) cl <- split(names(cl),cl)
    if(onlyAnnotated) cl <- lapply(cl, y=unique(unlist(go)), FUN=intersect)
    if(onlyInSets) go <- lapply(go, y=unique(unlist(cl)), FUN=intersect)
    by <- match.arg(by)
    bm <- sapply(go, y=unlist(cl), FUN=function(x,y) as.numeric(y %in% x))
    row.names(bm) <- unlist(cl)
    bm <- bm[,colSums(bm)>=minSize]
    # library(glmnet)
    # fits = cv.glmnet( bm2, as.factor(labs), family = "multinomial",
    #                   type.multinomial = "grouped", standardize=FALSE )
    # co <- coef(fits, fits$lambda.min)
    # co <- row.names(co)[co[,1]!=0][-1]
    #
    library(xgboost)
    labs <- rep(as.integer(names(cl)), sapply(cl,length))-1
    fit <- xgboost(bm, labs, params=list( booster="gbtree", max_depth=max_depth,
                                          subsample=0.75, colsample_bytree=1,
                                          objective="multi:softprob",
                                          eval_metric="mlogloss",
                                          min_child_weight=3,
                                          num_class=length(cl) ),
                   nrounds = nrounds, verbose=0)
    vi <- xgb.importance(model=fit)
    co <- as.character(vi$Feature[order(vi$Gain, decreasing=TRUE)[seq_len(N)]])
    fns <- list( overlap=function(x,y) length(intersect(x,y)),
                 proportion=function(x,y) length(intersect(x,y))/length(x),
                 enrichment=function(x,y){
                     expected <- length(y)*length(x)/nrow(bm)
                     length(intersect(x,y))/expected
                 },
                 jaccard=function(x,y) length(intersect(x,y))/length(union(x,y)),
                 ocoef=function(x,y) length(intersect(x,y))/min(length(x),length(y))
    )
    m <- lapply(fns, FUN=function(fn){
        t(sapply(go[co],FUN=function(y){
            sapply(cl, y=y, FUN=fn)
        }))
    })
    m$lengths <- sapply(cl,length)
    m
}

plot.goclassify <- function(m, n=max(6, ncol(m$overlap)), annotation_legend_param=list(),
                            what=c("proportion","log2enr","enrichment","jaccard"), transpose=FALSE, ...){
    n <- min(n, nrow(m$overlap))
    an <- data.frame(row.names=names(m$lengths), size=as.numeric(m$lengths))
    an <- HeatmapAnnotation(df=an,annotation_legend_param=annotation_legend_param,
                            which=ifelse(transpose,"row","column"))
    what <- match.arg(what)
    if(what=="log2enr"){
        m$log2enr <- log2(m$enrichment+0.1)
    }
    row.names(m[[what]]) <- gsub("^GO_","",row.names(m[[what]]))
    if(!is.null(breakStrings))
        row.names(m[[what]]) <- breakStrings(gsub("_"," ",row.names(m[[what]])))
    library(ComplexHeatmap)
    if(transpose){
        return(Heatmap( t(m[[what]][seq_len(n),]), right_annotation=an,
                 name=what, cell_fun=function(j, i, x, y, width, height, fill){
                     grid.text(m$overlap[j,i], x, y, gp=gpar(fontsize=10))
                 }, ...))
    }
    h <- Heatmap( m[[what]][seq_len(n),], top_annotation=an,
                  name=what, cell_fun=function(j, i, x, y, width, height, fill){
                      grid.text(m$overlap[i,j], x, y, gp=gpar(fontsize=10))
                  }, ...)
    draw(h, heatmap_legend_side="bottom", merge_legends = T)
}
