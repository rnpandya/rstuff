options(max.print=100)

source("http://bioconductor.org/biocLite.R")
library(AnnotationHub)
library(biomaRt)
library(clusterProfiler)
library(doSNOW)
library(GenomicRanges)
library(KEGG.db)
library(lattice)
library(mgsa)
library(reshape)
library(Rsamtools)
library(rtracklayer)
library(rtracklayer)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(VariantAnnotation)

# parallel setup

cl<-makeCluster(4)
registerDoSNOW(cl)

# shutdown
# stopCluster(cl)

# functions

get_kegg_genes <- function (pattern, sym = "ensembl_gene_id") {
    kegg <- toTable(KEGGPATHNAME2ID)
    kegg2 <- toTable(KEGGPATHID2EXTID)
    kegg_cancer_id <- kegg$path_id[grep(pattern,kegg$path_name,ignore.case = TRUE)]
    kegg$path_name[kegg$path_id %in% kegg_cancer_id]
    kegg_cancer_id <- paste("hsa", kegg_cancer_id, sep="")
    kegg_cancer_gene <- as.integer(unique(kegg2$gene_or_orf_id[kegg2$pathway_id %in% kegg_cancer_id]))
    ensembl_cancer_id <- getBM(attributes=c(sym),filter="entrezgene",values=kegg_cancer_gene,mart=ensembl)
    sort(unique(ensembl_cancer_id[[1]]))
}

inject <- function (vals, fn) {
    res <- NULL
    for (v in vals) {
        if (is.null(res)) {
            res <- v
        } else {
            res <- fn(res, v)
        }
    }
    res
}

adjacent <- function (vals, fn3, endval) {
    n <- length(vals)
    sapply(1:n, function (i) {
        fn3(if (i>1) vals[i-1] else endval, vals[i], if (i<n) vals[i+1] else endval)
    })
}

# should give sorted input!!
dupcount <- function (a) {
    d <- duplicated(a)
    len <- length(a)
    sapply(1:len, function (i) {
        if (i == len) {
            !d[i]
        } else if (d[i]) {
            0
        } else {
            n <- 1
            for (j in (i + 1):len) {
                if (! d[j]) break
                n <- n + 1
            }
            n
        }
    })
}
countMutations <- function (muts, regions) {
    c <- countOverlaps(regions, muts)
    r <- regions[c>0]
    r$mutations <- c[c>0]
    r
}

lapply_par <- function (lis, fn) {
    res <- foreach (i=lis) %dopar% fn(i);
    names(res) <- names(lis);
    res
}

lapply2_par <- function (lis1, lis2, fn) {
    res <- foreach (i=lis1) %:% foreach (j=lis2) %dopar% fn(i, j);
    res <- lapply(res, function (x) {names(x) <- names(lis2); x})
    names(res) <- names(lis1);
    res
}

overlap_p <- function (samples, regions) {
    samples_p <- sapply(samples, length) / 3.1e9
    regions_w <- width(regions)
    overlaps <- lapply_par(samples, function (s) { require(GenomicRanges); findOverlaps(s, regions)})
    result <- array(0, dim=length(regions))
    for (s in 1:50) {
        u <- unique(subjectHits(overlaps[[s]]))
        result[u] <- result[u] - log(1-(1-samples_p[s])^regions_w[u])
    }
    result
}
