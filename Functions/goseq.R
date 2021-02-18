#' goseq.enrichment
#'
#' Gets GO terms enriched in a given set of genes using the `goseq` package, and thereby 
#' correcting for the transcript length bias of RNA-seq differential expression data.
#'
#' @param allGenes A character vector containing all the gene symbols in the background (i.e. all tested genes). Set it to NA to use all genes.
#' @param deGenes A character vector containing the dysregulated genes.
#' @param gotype A character vector indicating the GO ontologies to fetch. Can be any combination of 'GO:BP', 'GO:MF', 'GO:CC'.
#' @param cutoff A number between 0 and 1, indicating the FDR threshold below which terms will be retained. Default 0.1.
#' @param cutoff.onFDR Whether to apply the cutoff on FDR (default), otherwise it will be applied on PValue.
#' @param org The organism, e.g. "hg19" (human) or "mm9" (mouse), although other genomes versions should work too.
#' @param minCatSize The minimum size of categories (nb of genes with the annotation) to be considered. Default 10.
#' @param maxCatSize The maximum size of categories (nb of genes with the annotation) to be considered. Default 1000.
#' @param maxResults The maximum number of categories to return. Default 200.
#'
#' @return A data.frame.
#'
#' @export
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



#' goTreemap
#'
#' Produce a treemap summarizing GO enrichment results.
#'
#' @param go A data.frame of class GOresults, as produced by the `go.enrichment` or `goseq.enrichment` functions.
#' @param gotype The GO ontology of the categories: 'BP', 'MF' or 'CC'.
#' @param flat logical; whether to disable hierarchical colouring (default F, i.e. activated)
#' @param removeParents logical; whether or not to remove parent categories with enriched children. Default TRUE.
#' @param maxCats integer; the maximum number of categories to plot (default 25)
#' @param maxHeavyTry integer; max size of combinations to try the heavy way -- don't try more than 5!!
#' @param minParents integer; min number of parents for hierarchical colouring. Default 3.
#' @param inflate logical; whether to inflate labels to the size of the box
#' @param tolerated integer; max overlaping subterms between colouring parent terms. Default 0.
#' @param enrichedAncestorsOny logical; whether hierarchical colouring will select parents only from the enriched ones. Default FALSE.
#' @param returndf logical; whether to return the resulting data.frame. Default FALSE.
#'
#' @return Produces a figure, and (if returndf=TRUE) returns a data.frame.
#'
#' @export
goTreemap <- function(	go,
                       gotype="BP",
                       flat=F,
                       removeParents=T,
                       maxCats=25, 
                       maxHeavyTry=NULL,
                       minParents=NULL,
                       inflate=F,
                       tolerated=0,
                       enrichedAncestorsOnly=F,
                       returndf=F 
){
    require(GO.db)
    require(treemap)
    gotype <- gsub("GO:","",gotype,fixed=T)
    gotype <- match.arg(gotype, c("CC","MF","BP"), several.ok=T)
    go[,1] <- as.character(go[,1])
    go$logp <- -log10(as.numeric(go$FDR))
    allgos <- go[,1]
    
    message("Preparing ancestors")
    if(gotype=="BP"){
        allancestors <- as.character(unlist(AnnotationDbi::as.list(GOBPANCESTOR)[allgos]))
    }else if(gotype=="CC"){
        allancestors <- as.character(unlist(AnnotationDbi::as.list(GOCCANCESTOR[allgos])))
    }else if(gotype=="MF"){
        allancestors <- as.character(unlist(AnnotationDbi::as.list(GOMFANCESTOR[allgos])))
    }else{
        error("Invalid gotype argument")
    }
    if(removeParents){
        message("Removing categories that have enriched children")
        go <- go[!(go[,1] %in% allancestors),]
    }else{
        flat=TRUE
    }
    go$term2 <- apply(go[,c("Term","FDR")], 1, FUN=function(x){ paste(x[[1]]," (",format(as.numeric(x[[2]]), scientific=T, digits=2),")",collapse="",sep="") })
    if(flat==F){
        if(enrichedAncestorsOnly){
            enrichedancestors <- sapply(allancestors[allancestors %in% allgos],FUN=function(x){ as.character(Term(x))})
        }else{
            enrichedancestors <- sapply(allancestors,FUN=function(x){ as.character(Term(x))})
        }
    }
    
    if(removeParents==F | !is.null(maxCats)){
        if(is.null(maxCats))	maxCats <- 20
        go <- go[1:min(nrow(go),maxCats),]
    }
    allgos <- go[,1]
    
    if(flat==F){
        parents <- getGroupingParents(allgos, gotype="BP", maxHeavyTry=maxHeavyTry, minParents=minParents, inSubset=enrichedancestors, tolerated=tolerated)
        if(is.null(parents)){
            message("No non-trivial enriched ancestor; will plot without hierarchy.")
            flat <- T
        }
    }
    if(flat){
        # prepare dataframe for treemap
        go2 <- go[,c("GO.ID","term2","logp")]
        treemap(go2, c("term2"), vSize="logp", lowerbound.cex.labels = 0, bg.labels = "#CCCCCC00", inflate=inflate, position.legend="none", fontsize.title=0)
        return(go2)
    }
    
    # prepare dataframe for treemap
    go2 <- go[,c("GO.ID","term2","logp")]
    go2$parent <- sapply(go2[,1],parents=parents, allancestors=allancestors,FUN=function(x,parents,allancestors){
        pp <- try(as.character(allancestors[[as.character(x)]]), silent=T)
        if(class(pp)=="try-error")	return(NA)
        if(!any(pp %in% parents))	return(NA)
        return(pp[which(pp %in% parents)[1]])
    })
    go2 <- go2[order(go2$parent),]
    go2$parent[is.na(go2$parent)] <- "other"
    go2$parent <- as.factor(go2$parent)
    go2$sortid <- as.numeric(go2$parent)
    
    treemap(go2, c("term2"), vSize="logp", sortID="sortid", type="categorical", vColor="parent", lowerbound.cex.labels = 0, bg.labels = "#CCCCCC00", inflate=inflate, position.legend="bottom", fontsize.title=0)
    if(returndf) return(go2)
}

.countGenesInCat <- function(x,normGenes=T){
    if(is.null(x) | !is.matrix(x)){
        if(length(x)>0)	return(sum(x))
        return(NA)
    }
    x <- x[,which(colSums(x)>0)]
    if(is.null(x) | !is.matrix(x)){
        if(length(x)>0)	return(sum(x))
        return(NA)
    }	
    if(normGenes) x <- t(t(x)/colSums(x))
    rowSums(x)
}


#' getGroupingParents
#'
#' Get maximally-including, non-overlapping parent categories from a given set of categories
#'
#' @param goid A vector of GO category IDs
#' @param gotype The GO ontology of the categories: 'BP', 'MF' or 'CC'.
#' @param maxHeavyTry integer; max size of combinations to try the heavy way -- don't try more than 5!!
#' @param minParents integer; min number of parents for hierarchical colouring. Default 3.
#' @param tolerated integer; max overlaping subterms between colouring parent terms. Default 0.
#'
#' @return Either NULL or a character vector of parent categories
#'
#' @export
getGroupingParents <- function(goids, gotype="BP", maxHeavyTry=NULL, minParents=NULL, inSubset=NULL, tolerated=0){
    gotype <- gsub("GO:","",gotype,fixed=T)
    message("Preparing ancestor intersections")
    # create a dataframe of each categories' ancestors
    if(gotype=="BP"){
        anc <- AnnotationDbi::as.list(GOBPANCESTOR[goids])
    }else if(gotype=="CC"){
        anc <- AnnotationDbi::as.list(GOCCANCESTOR[goids])
    }else if(gotype=="MF"){
        anc <- AnnotationDbi::as.list(GOMFANCESTOR[goids])
    }
    anc <- sapply(anc,FUN=function(x){ as.character(Term(x))})
    tmpl <- sapply(anc,FUN=length)
    df2 <- data.frame(go=rep(names(tmpl),as.numeric(tmpl)), ancestor=as.character(unlist(anc)))
    if(!is.null(inSubset)) df2 <- df2[df2$ancestor %in% inSubset,]
    tt <- table(df2$ancestor)
    df3 <- df2
    df3$ancestor <- as.character(df3$ancestor)
    
    # remove ancestors that have less than XX enriched children
    df3 <- df3[!(df3$ancestor %in% names(tt)[tt<2]),]
    if(nrow(df3) > 70)	df3 <- df3[!(df3$ancestor %in% names(tt)[tt<3]),]
    # remove ubiquitous ancestors
    df3 <- df3[!(df3$ancestor %in% names(tt[tt>=(length(goids)-1)]) ),]
    if(nrow(df3)==0) return(NULL)
    
    # prepare a matrix of intersections between the different ancestor categories
    test <- aggregate(df3$ancestor,by=list(go=df3$go),FUN=as.character)
    test2 <- as.list(test$x)
    names(test2) <- test$go
    m <- matrix(0,nrow=length(unique(df3$go)),ncol=length(unique(df3$ancestor)))
    colnames(m) <- unique(df3$ancestor)
    row.names(m) <- test$go
    m <- t(mapply(test2,FUN=function(x){ as.numeric(colnames(m) %in% x)}))
    m2 <- matrix(0,nrow=length(unique(df3$ancestor)),ncol=length(unique(df3$ancestor)))
    colnames(m2) <- unique(df3$ancestor)
    row.names(m2) <- unique(df3$ancestor)
    for(i in 1:nrow(m2)){
        m2[i,] <- as.numeric(apply(m,2,FUN=function(x){ sum((m[,i]+x)==2)}))
    }
    if(ncol(m2) < 3)	return(NULL)
    
    i <- 1:ncol(m2)
    j <- 2
    while(length(i) > 35){
        i <- (1:ncol(m2))[rowSums(m2)>j]
        j <- j+1
    }
    message("Data prepared. Now checking for non-overlapping combinations of ",length(i)," parent terms.")
    
    # checking for non-overlapping combinations
    checkCombn <- function(x,tolerated=0){
        ov <- sum(m2[x,x])-sum(diag(m2[x,x]))
        if(ov<=tolerated){
            return(list(ov,sum(rowSums(m[,x])>0)/length(goids),x))
        }else{
            return(list(ov,NA,x))
        }
    }
    aco <- data.frame(V1=vector(mode="numeric",length=0), V2=vector(mode="numeric",length=0), V3=vector(mode="list",length=0))
    if(is.null(maxHeavyTry)){
        if(ncol(m2) < 22){
            maxHeavyTry <- min(5,floor(ncol(m2)/2))
        }else{
            maxHeavyTry <- 4
        }
    }
    if(!is.null(minParents)){
        if(minParents > maxHeavyTry){
            minParents <- maxHeavyTry
        }
    }else{
        minParents <- 3
    }
    for(nbc in min(length(i),max(2,minParents)):min(length(i),maxHeavyTry) ){
        print(paste("Checking with",nbc,"terms..."))
        co <- combn(i,nbc, tolerated=tolerated, FUN=checkCombn)
        if(is.null(dim(co))){
            co <- matrix(co,nrow=1)
        }else{
            co <- t(co)
        }
        co <- co[!is.na(co[,2]) & as.numeric(co[,2]) > 0.2,]
        aco <- rbind(aco,as.data.frame(co))
    }
    if(nrow(aco)==0)	return(NULL)
    aco[,2] <- as.numeric(aco[,2])
    aco$t <- lapply(aco[,3],FUN=function(x){ as.character(colnames(m2)[as.numeric(unlist(x))])})
    aco$l <- sapply(aco$t,FUN=length)
    aco <- aco[aco[,2]==max(aco[,2]),]
    if(nrow(aco)>1){
        # more than one combination with same coverage
        aco$nbi <- sapply(aco$t, FUN=function(x){
            return(sum(rowSums(m2[as.character(unlist(x)),])))
        })
        aco$nbc <- sapply(aco$t, FUN=function(x){ nchar(paste(x,collapse="",sep=""))})
        aco <- aco[order(-aco$l, aco$nbi, aco$nbc),]
    }
    parents <- as.character(unlist(aco$t[1]))
    
    if(length(parents) == maxHeavyTry & aco[1,2] < 1){
        message("Found an imperfect combination; trying to increment it with additional terms")
        
        # attempt to add a member to the combination
        parents <- .tryAddParent(parents,m2,co)
        parents <- .tryAddParent(parents,m2,co)
    }
    
    return(parents)
}

.tryAddParent <- function(parents,m2,co){
    tmp <- colnames(m2[,colSums(m2[parents,])==0])
    if(!is.null(tmp)){
        co <- t(sapply(tmp, FUN=function(x){ checkCombn(c(parents,x)) }))
        if(nrow(co[!is.na(co[,2]),])>0){
            co <- as.data.frame(co[!is.na(co[,2]),])
            co[,2] <- as.numeric(co[,2])
            co <- co[order(as.numeric(co[,2]),decreasing=T),]
            return(as.character(unlist(co[1,3])))
        }
    }
    return(parents)	
}

