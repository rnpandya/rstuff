# reload .RData to start with

# public databases

ensembl <- useMart("ensembl")
ensembl <- useDataset("hsapiens_gene_ensembl",mart=ensembl)

ah <- AnnotationHub()
filters(ah) <- list(Species="Homo sapiens")
ens69 <- ah$ensembl.release.69.gtf.homo_sapiens.Homo_sapiens.GRCh37.69.gtf_0.0.1.RData

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
exons <- exons(txdb)
genes <- genes(txdb)
introns <- setdiff(genes,exons)

genemap <- getBM(attributes=c('ensembl_gene_id', 'hgnc_symbol'),
 filters = c('hgnc_symbol','with_hgnc'), values=list(unique(tfbs$V4),TRUE), mart=ensembl)

kegg <- toTable(KEGGPATHNAME2ID)
kegg2 <- toTable(KEGGPATHID2EXTID)
kegg_cancer_id <- kegg$path_id[grep("(cancer)|(carcinoma)|(leukemia)|(oma$)",kegg$path_name,ignore.case = TRUE)]
kegg_cancer_id <- paste("hsa", kegg_cancer_id, sep="")
kegg_cancer_gene <- as.integer(unique(kegg2$gene_or_orf_id[kegg2$pathway_id %in% kegg_cancer_id]))
ensembl_cancer_id <- getBM(attributes=c("entrezgene","ensembl_gene_id","hgnc_id","hgnc_symbol"),filter="entrezgene",values=kegg_cancer_gene,mart=ensembl)

gencode <- import.gff("t:/ravip/encode/gencode.v21.annotation.gff3.gz")

df <- read.table("t:/ravip/encode/allGeneCorrelations100000.p2.txt.gz",sep="\t",header=TRUE)
dhsgene <- makeGRangesFromDataFrame(df, seqnames.field = "dhs_chr", start.field = "dhs_start", end.field = "dhs_end", ignore.strand = TRUE, keep.extra.columns = TRUE)

df <- read.table("t:/ravip/encode/wgEncodeBroadHmmK562HMM.bed.gz",sep="\t")
chromhmm_k562 <- makeGRangesFromDataFrame(df, seqnames.field = "V1", start.field = "V2", end.field = "V3", ignore.strand = TRUE, keep.extra.columns = TRUE)
hmmreg <- tapply(chromhmm_k562, chromhmm_k562$V4, reduce)

enctss <- import.gff("t:/ravip/encode/gencode.v19.TSS.notlow.gff")
start(enctss) <- start(enctss) - 2000
end(enctss) <- end(enctss) + 200

df <- read.table("t:/ravip/encode/wgEncodeRegDnaseClusteredV3.bed.gz",sep="\t")
dnase <- makeGRangesFromDataFrame(df, seqnames.field = "V1", start.field = "V2", end.field = "V3", ignore.strand = TRUE, keep.extra.columns = TRUE)

df <- read.table("t:/ravip/encode/wgEncodeRegTfbsClusteredV3.bed.gz",sep="\t")
tfbs <- makeGRangesFromDataFrame(df, seqnames.field = "V1", start.field = "V2", end.field = "V3", ignore.strand = TRUE, keep.extra.columns = TRUE)

hsgo <- readGAF("c:/ravip/encode/gene_association.goa_human.gz")

# my working data

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
