source("http://bioconductor.org/biocLite.R")

biocLite("TxDb.Hsapiens.UCSC.hg19.knownGene")
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

biocLite("VariantAnnotation")
library(VariantAnnotation)
library(GenomicRanges)
library(Rsamtools)

biocLite("AnnotationHub")
library(AnnotationHub)
ah <- AnnotationHub()
filters(ah) <- list(Species="Homo sapiens")

allvcf <- sapply(grep(".vcf.gz$", list.files("t:/ravip/tcga-laml-wgs/full/"),value=TRUE), function (f) {
        v <- readVcf(paste("t:/ravip/tcga-laml-wgs/full/",f,sep=""),"hg19")
        ok <- geno(v)$SSC[,2]>=20
        ok[is.na(ok)] <- FALSE
        v[ok]
    })

allvcf2 <- sapply(allvcf, rowData)

allvcf2c <- do.call(c,allvcf2)
allvcfov <- countOverlaps(encode,allvcf2c)
table(allvcfov)
plot(table(allvcfov)[5:40])
encode[allvcfov >= 30,]
library(rtracklayer)
enctss <- import.gff("c:/ravip/encode/gencode.v19.TSS.notlow.gff")
start(enctss) <- start(enctss) - 2000
end(enctss) <- end(enctss) + 200
sum(width(enctss))
tssov <- countOverlaps(enctss,allvcf2c)

table(tssov)
enctss[tssov > 2]
tssov[tssov>2]

library("biomaRt")
ensembl=useMart("ensembl")
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)

tssgrouplev <- levels(enctss$group)
enctssgrp <- unlist(lapply(enctss[tssov > 2]$group,
    function (lvl) strsplit(tssgrouplev[lvl], " ")), recursive=TRUE)
enctssgrp <- unlist(lapply(enctssgrp, strsplit, ","))
enctssgrp <- gsub("\\..*", "", enctssgrp)
lapply(c("ensembl_gene_id","ensembl_transcript_id"),
    function(id) getBM(attributes=c(id,'description'), values=enctssgrp, filters=id, mart=ensembl))

# dnase data - similar for tfbs
df <- read.table("c:/ravip/encode/wgEncodeRegDnaseClusteredV3.bed.gz",sep="\t")
dnase <- makeGRangesFromDataFrame(df, seqnames.field = "V1", start.field = "V2", end.field = "V3", ignore.strand = TRUE, keep.extra.columns = TRUE)
sum(width(reduce(dnase)))
sum(countOverlaps(reduce(dnase),allvcf2c))

enctss$ensembl_gene <- sapply(tssgrouplev[enctss$group], function(g) gsub("\\..*", "", strsplit(g, " ")[[1]][2]))

tssgo <- getBM(attributes = c("ensembl_gene_id", "go_id"), filters = "ensembl_gene_id", values = enctss$ensembl_gene, mart = ensembl)

ens69 <- ah$ensembl.release.69.gtf.homo_sapiens.Homo_sapiens.GRCh37.69.gtf_0.0.1.RData
cosmic_ens <- findOverlaps(ens69,rowData(cosmic_67))
cosmic_ens_names <- unique(ens69$gene_id[unique(queryHits(cosmic_ens))])
tss_cosmic <- enctss[sapply(enctss$ensembl_gene, function(g) g %in% cosmic_ens_names)]

sum(width(tss_cosmic))
sum(countOverlaps(tss_cosmic, allvcf2c))

go_hsc <- c(GOBPOFFSPRING[["GO:0071425"]], "GO:0071425")
go_hsc_gene <- tssgo$ensembl_gene_id[tssgo$go_id %in% go_hsc]
tss_hsc <- reduce(enctss[enctss$ensembl_gene %in% go_hsc_gene])
sum(width(tss_hsc))
sum(countOverlaps(tss_hsc, allvcf2c))

top20 <- c("JAK2", "NPM1", "FLT3", "KIT", "MYD88", "ABL1", "CALR", "NOTCH1", "TET2", "DNMT3A", "ASXL1", "TP53", "BRAF", "SF3B1", "IKZF1", "NRAS", "CDKN2A", "SRSF2", "ATM", "GATA1")
# translated manually with http://www.genenames.org/cgi-bin/symbol_checker
top20_hgncid <- c("HGNC:6192", "HGNC:7910", "HGNC:3765", "HGNC:6342", "HGNC:7562", "HGNC:76", "HGNC:1455", "HGNC:7881", "HGNC:25941", "HGNC:2978", "HGNC:18318", "HGNC:11998", "HGNC:1097", "HGNC:10768", "HGNC:10769", "HGNC:13176", "HGNC:7989", "HGNC:4170")
top20_ensembl <- getBM(attributes = "ensembl_gene_id", filter="hgnc_id", values = top20_hgncid, mart=ensembl)

tfbs20 <- reduce(tfbs[tfbs$V4 %in% top20])
sum(width(tfbs20))
sum(countOverlaps(tfbs20, allvcf2c))

### 2/9/15 second pass
source("http://bioconductor.org/biocLite.R")

biocLite("TxDb.Hsapiens.UCSC.hg19.knownGene")
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

biocLite("VariantAnnotation")
library(VariantAnnotation)
library(GenomicRanges)
library(Rsamtools)
library(rtracklayer)

biocLite("AnnotationHub")
library(AnnotationHub)
ah <- AnnotationHub()
filters(ah) <- list(Species="Homo sapiens")

allvcf <- sapply(grep(".vcf.gz$", list.files("t:/ravip/tcga-laml-wgs/full/"),value=TRUE), function (f) {
        v <- readVcf(paste("t:/ravip/tcga-laml-wgs/full/",f,sep=""),"hg19")
        ok <- geno(v)$SSC[,2]>=20
        ok[is.na(ok)] <- FALSE
        v[ok]
    })
allvcf2 <- sapply(allvcf, rowData)

exons <- exons(txdb)
genes <- genes(txdb)
introns <- setdiff(genes,exons)

enctss <- import.gff("c:/ravip/encode/gencode.v19.TSS.notlow.gff")
start(enctss) <- start(enctss) - 2000
end(enctss) <- end(enctss) + 200
enctss_r <- reduce(enctss)
enctss_r_w <- sum(width(enctss_r))
enctss$ensembl_gene <- sapply(levels(enctss$group)[enctss$group], function(g) gsub("\\..*", "", strsplit(g, " ")[[1]][2]))

library("biomaRt")
ensembl <- useMart("ensembl")
ensembl <- useDataset("hsapiens_gene_ensembl",mart=ensembl)
# can't re-find this list
# could also try list from http://www.tumorportal.org/tumor_types?ttype=AML
top20 <- c("JAK2", "NPM1", "FLT3", "KIT", "MYD88", "ABL1", "CALR", "NOTCH1", "TET2", "DNMT3A", "ASXL1", "TP53", "BRAF", "SF3B1", "IKZF1", "NRAS", "CDKN2A", "SRSF2", "ATM", "GATA1")
# translated manually with http://www.genenames.org/cgi-bin/symbol_checker
top20_hgncid <- c("HGNC:6192", "HGNC:7910", "HGNC:3765", "HGNC:6342", "HGNC:7562", "HGNC:76", "HGNC:1455", "HGNC:7881", "HGNC:25941", "HGNC:2978", "HGNC:18318", "HGNC:11998", "HGNC:1097", "HGNC:10768", "HGNC:10769", "HGNC:13176", "HGNC:7989", "HGNC:4170")
top20_ensembl <- getBM(attributes = "ensembl_gene_id", filter="hgnc_id", values = top20_hgncid, mart=ensembl)

# dnase data - similar for tfbs
df <- read.table("c:/ravip/encode/wgEncodeRegDnaseClusteredV3.bed.gz",sep="\t")
dnase <- makeGRangesFromDataFrame(df, seqnames.field = "V1", start.field = "V2", end.field = "V3", ignore.strand = TRUE, keep.extra.columns = TRUE)
df <- read.table("c:/ravip/encode/wgEncodeRegTfbsClusteredV3.bed.gz",sep="\t")
tfbs <- makeGRangesFromDataFrame(df, seqnames.field = "V1", start.field = "V2", end.field = "V3", ignore.strand = TRUE, keep.extra.columns = TRUE)
rm(df)

genemap <- getBM(attributes=c('ensembl_gene_id', 'hgnc_symbol'),
 filters = c('hgnc_symbol','with_hgnc'), values=list(unique(tfbs$V4),TRUE), mart=ensembl)

# graphing

biocLite("reshape")
library(reshape)
library(lattice)

allvcf2n <- sapply(allvcf2, length)
regions <- list(genes=genes,introns=introns,exons=exons,tss=tss,tfbs=tfbs,dnase=dnase)
regions <- sapply(regions, reduce)
# long!!
counts <- sapply(regions, function(r) sapply(allvcf2, function (v) sum(countOverlaps(r,v))))
region_scale <- sapply(regions, function(r) 1e-6*sum(width(r)))
muts <- data.frame(counts/t(array(region_scale,dim=c(length(regions),50))))
muts_rel <- muts / (sapply(allvcf2, length) / 3.1e3)
muts$id <- 1:50
muts_rel$id <- 1:50
xyplot(value~id,groups=variable, data = melt(muts,id.vars=c('id')))
xyplot(value~id,groups=variable, data = melt(muts_rel,id.vars=c('id')),auto.key=TRUE)
summary(muts_rel)
bwplot(value~variable, data=melt(muts_rel,id.vars=c('id'),auto.key=TRUE))

grep("kegg",listMarts()$biomart,ignore.case = TRUE, value=TRUE)
grep("cosmic",names(ah),ignore.case = TRUE,value=TRUE)

cosmic <- ah$goldenpath.hg19.database.cosmic_0.0.1.RData

biocLite("KEGG.db")
library(KEGG.db)
kegg <- toTable(KEGGPATHNAME2ID)
kegg2 <- toTable(KEGGPATHID2EXTID)
kegg_cancer_id <- kegg$path_id[grep("(cancer)|(carcinoma)|(leukemia)|(oma$)",kegg$path_name,ignore.case = TRUE)]
kegg$path_name[kegg$path_id %in% kegg_cancer_id]
kegg_cancer_id <- paste("hsa", kegg_cancer_id, sep="")
kegg_cancer_gene <- as.integer(unique(kegg2$gene_or_orf_id[kegg2$pathway_id %in% kegg_cancer_id]))
ensembl_cancer_id <- getBM(attributes=c("entrezgene","ensembl_gene_id","hgnc_id","hgnc_symbol"),filter="entrezgene",values=kegg_cancer_gene,mart=ensembl)

tfbs_cancer <- tfbs[tfbs$V4 %in% ensembl_cancer_id$hgnc_symbol]
regions$tfbs_cancer <- reduce(tfbs_cancer)
counts$tfbs_cancer <- array(sapply(allvcf2, function (v) sum(countOverlaps(regions$tfbs_cancer,v))))
# screws up type of counts, need to do this differently!!!

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

plot(hclust(dist(t(muts_rel[1:9]))))

tfbs_leukemia_gene_muts <- tapply(X=tfbs_leukemia$count, INDEX=tfbs_leukemia$V4, FUN=sum)
tfbs_leukemia_gene_muts <- sort(tfbs_leukemia_gene_muts[!is.na(tfbs_leukemia_gene_muts)])
tfbs$count <- countOverlaps(tfbs,allmuts)
tfbs_gene_muts <- tapply(X=tfbs$count, INDEX=tfbs$V4, FUN=sum)
tfbs_gene_muts <- sort(tfbs_gene_muts[!is.na(tfbs_gene_muts)])

tfbs$wide <- width(tfbs)
tfbs_gene_width <- tapply(X=tfbs$wide, INDEX=tfbs$V4, FUN=sum)
tfbs_gene_rel <- sort(1e6 * tfbs_gene_muts / tfbs_gene_width)

biocLite("clusterProfiler")
library(clusterProfiler)
top_entrez <- getBM(attributes=c("entrezgene"),filter=c("hgnc_symbol"),value=names(tfbs_gene_rel)[120:161],mart=ensembl)
top_entrez <- sapply(top_entrez,as.character)
enrichGO(top_entrez)

gene_tfbs <- tapply(tfbs, INDEX=tfbs$V4, FUN=c)
gene_tfbs <- sapply(gene_tfbs, reduce)
sort(sapply(gene_tfbs, function(g) sum(width(g))))

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

gencode <- import.gff("t:/ravip/encode/gencode.v21.annotation.gff3.gz")
gencode_mut <- sapply(allvcf2, countMutations, gencode)

dnase_mut <- sapply(allvcf2, countMutations, dnase)
dnase_mut_all <- inject(dnase_mut, c)
dnase_mut_all$dups <- dupcount(dnase_mut_all)

df <- read.table("c:/ravip/encode/allGeneCorrelations100000.p05.txt",sep="\t",header=TRUE)
dhsgene <- makeGRangesFromDataFrame(df, seqnames.field = "dhs_chr", start.field = "dhs_start", end.field = "dhs_end", ignore.strand = TRUE, keep.extra.columns = TRUE)

dhsgene_mut <- sapply(allvcf2, countMutations, dhsgene)
# todo: sum muts, avg other cols, etc.
dhsgene_mut <- sapply(dhsgene_mut, function(r) r[!duplicated(r)])
dhsgene_mut_all <- sort(inject(dhsgene_mut, c))
dhsgene_mut_all$dups <- dupcount(dhsgene_mut_all)

dnase_blood <- dnase[sapply(dnase$V7,function(x) length(intersect(strsplit(levels(dnase$V7)[x],","),c("12","13","21","22","23","93","116","117","15","85","86","95")))>0)]

# dhs from blood tissue (normal+tumor)
dnase_blood_mut <- lapply(allvcf2, subsetByOverlaps, reduce(dnase_blood))
summary(sort(sapply(dnase_blood_mut,length)/sapply(allvcf2,length))/(sum(width(reduce(dnase_blood)))/3.1e9))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.3048  0.4575  0.5528  0.5533  0.6462  0.8706 

# all dhs
dnase_mut <- lapply(allvcf2, subsetByOverlaps, reduce(dnase))
summary(sort(sapply(dnase_mut,length)/sapply(allvcf2,length))/(sum(width(reduce(dnase)))/3.1e9))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.4287  0.5315  0.5703  0.5855  0.6389  0.8259 

# dhs with p<.05 corr w/nearby gene expr
dhsgene_mut2 <- lapply(allvcf2, subsetByOverlaps, reduce(dhsgene))
summary(sapply(dhsgene_mut2,length)/sapply(allvcf2,length)/(sum(width(reduce(dhsgene)))/3.1e9))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.3089  0.4183  0.4719  0.4831  0.5392  0.7461 

# dhs with p<?? corr w/nearby gene expr
dhsgene2_mut <- lapply(allvcf2, subsetByOverlaps, reduce(dhsgene2))
summary(sapply(dhsgene2_mut,length)/sapply(allvcf2,length)/(sum(width(reduce(dhsgene2)))/3.1e9))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.3489  0.4681  0.4932  0.5141  0.5679  0.7422 

# all tfbs
tfbs_mut <- lapply(allvcf2, subsetByOverlaps, reduce(tfbs))
summary(sort(sapply(tfbs_mut,length)/sapply(allvcf2,length))/(sum(width(reduce(tfbs)))/3.1e9))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.4793  0.6617  0.7251  0.7150  0.7649  0.8957 

# all genes
genes_mut <- lapply(allvcf2, subsetByOverlaps, reduce(genes))
summary(sort(sapply(genes_mut,length)/sapply(allvcf2,length))/(sum(width(reduce(genes)))/3.1e9))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.5022  0.7315  0.7847  0.7736  0.8446  0.9668 

# exons only
exons_mut <- lapply(allvcf2, subsetByOverlaps, reduce(exons))
summary(sort(sapply(exons_mut,length)/sapply(allvcf2,length))/(sum(width(reduce(exons)))/3.1e9))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.3290  0.4018  0.4358  0.4401  0.4841  0.5658 

#
# same with higher-quality mutations (conf >= 40)
#

allvcfhi <- lapply(allvcf, function(v) rowData(v[geno(v)$SSC[,2]>=40]))

dnase_blood_mut <- lapply(allvcfhi, subsetByOverlaps, reduce(dnase_blood))
summary(sort(sapply(dnase_blood_mut,length)/sapply(allvcfhi,length))/(sum(width(reduce(dnase_blood)))/3.1e9))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0000  0.3082  0.5221  0.4831  0.6222  1.1440 
dnase_mut <- lapply(allvcfhi, function (v) subsetByOverlaps(v, reduce(dnase)))
summary(sort(sapply(dnase_mut,length)/sapply(allvcfhi,length))/(sum(width(reduce(dnase)))/3.1e9))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.3400  0.4931  0.5446  0.5583  0.6127  0.8180 
dhsgene_mut2 <- lapply(allvcfhi, subsetByOverlaps, reduce(dhsgene))
summary(sapply(dhsgene_mut2,length)/sapply(allvcfhi,length)/(sum(width(reduce(dhsgene)))/3.1e9))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.1870  0.3470  0.4090  0.4284  0.4966  0.8280 
dhsgene2_mut <- lapply(allvcfhi, subsetByOverlaps, reduce(dhsgene2))
summary(sapply(dhsgene2_mut,length)/sapply(allvcfhi,length)/(sum(width(reduce(dhsgene2)))/3.1e9))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.2510  0.3922  0.4477  0.4667  0.5454  0.8886 
tfbs_mut <- lapply(allvcfhi, subsetByOverlaps, reduce(tfbs))
summary(sort(sapply(tfbs_mut,length)/sapply(allvcfhi,length))/(sum(width(reduce(tfbs)))/3.1e9))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.2901  0.5138  0.5653  0.5896  0.6743  0.9389 
genes_mut <- lapply(allvcfhi, subsetByOverlaps, reduce(genes))
summary(sort(sapply(genes_mut,length)/sapply(allvcfhi,length))/(sum(width(reduce(genes)))/3.1e9))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.2845  0.5070  0.5699  0.5806  0.6646  0.9156 
exons_mut <- lapply(allvcfhi, subsetByOverlaps, reduce(exons))
summary(sort(sapply(exons_mut,length)/sapply(allvcfhi,length))/(sum(width(reduce(exons)))/3.1e9))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.1887  0.3244  0.3722  0.3718  0.4073  0.5486 

dhsgene_blood <- subsetByOverlaps(dhsgene, dnase_blood)
dhsgene_blood_mut <- lapply(allvcfhi, subsetByOverlaps, reduce(dhsgene_blood))
summary(sapply(dhsgene_blood_mut,length)/sapply(allvcfhi,length)/(sum(width(reduce(dhsgene_blood)))/3.1e9))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0000  0.0000  0.3309  0.3219  0.4877  1.1130 

# ChromHMM segmentation

df <- read.table("c:/ravip/encode/wgEncodeBroadHmmK562HMM.bed.gz",sep="\t")
chromhmm_k562 <- makeGRangesFromDataFrame(df, seqnames.field = "V1", start.field = "V2", end.field = "V3", ignore.strand = TRUE, keep.extra.columns = TRUE)
hmmreg <- tapply(chromhmm_k562, chromhmm_k562$V4, reduce)
chromhmm_k562_ov <- lapply(allvcfhi, findOverlaps, chromhmm_k562, ignore.strand=TRUE, select="all", type="any")
m <- lapply(chromhmm_k562_ov, function (h) { h2 <- chromhmm_k562[subjectHits(h)]; tapply(h2, h2$V4, length) })
m <- t(data.frame(m))
chromhmm_k562_len <- tapply(chromhmm_k562, chromhmm_k562$V4, function (r) sum(width(reduce(r))))
m <- m/t(array(chromhmm_k562_len,dim=c(15,50)))
m <- m/(sapply(allvcfhi, length)/3.1e9)
sort(unlist(lapply(data.frame(m),mean)))
order(unlist(lapply(data.frame(m),mean)))
m <- cbind(m,id=1:50)
bwplot(value~X2, data=melt(m[,c(1:5,8:15)],id.vars=c('id'),auto.key=TRUE))

#
# 2/25/15
#
# now trying things in parallel
#

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

allvcfhi <- lapply_par(grep(".vcf.gz$", list.files("t:/ravip/tcga-laml-wgs/full/"),value=TRUE), function (f) {
    require(VariantAnnotation)
    v <- readVcf(paste("t:/ravip/tcga-laml-wgs/full/",f,sep=""),"hg19")
    ok <- geno(v)$SSC[,2]>=40
    ok[is.na(ok)] <- FALSE
    rowData(v[ok])
})

chromhmm_k562_ov <- lapply2_par(hmmreg, allvcfhi, function (r, v) {
    require(VariantAnnotation)
    findOverlaps(v, r, ignore.strand=TRUE, select="all", type="any")
})