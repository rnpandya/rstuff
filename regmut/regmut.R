#
# looking at RNA-seq
#
gx <- read.table('c:/ravip/pancan/scpilot/joint_fpkm_uq.tsv.gz',header=T)
colnames(gx) <- gsub('^X','', colnames(gx))
colnames(gx) <- gsub('\\.', '-', colnames(gx))
# filter to samples with LAML VCF & RNA-seq
gx2 <- gx[c(
  'feature',
  'bae3bd60-74cf-466b-b488-5918b8045ef7',
  '826d37cd-4da7-4edd-a317-4b8779a9005a',
  '53d81e42-430b-480f-a5a5-591e1e0c2d4c',
  '81929ba6-a393-4090-b480-d1e6ce1d1e59',
  '70bb0fbf-4e70-4555-b4c3-9b7446032a08',
  '1ba2510e-84ec-447c-b924-4e4c90984dec',
  'f7baaf60-149c-493c-b6f9-6ffc63b9b040',
  'f00a9542-9f5b-4847-836d-bfbf60f4302e',
  '76773f60-1e96-489b-ad14-f74e240fbcf4',
  '9526e905-5339-4e36-a3ef-2754ae9303f3',
  'c0655f95-468f-4e63-bfad-2c0e97f7dc72',
  '4d54f302-3b92-434c-8ff2-dde1ead2d089',
  '7ec9e97b-8bef-4fb0-8639-cec88e770dbd',
  '5a372c2d-ce65-47c3-9bdc-54d2dbbf7213',
  '4588e050-e15f-476f-9065-ba747bdd3736',
  '30ca50fc-bd27-4d7a-bfe1-d09b5f60b771',
  '99697721-c39e-4043-a42c-865d3df26da1',
  '808e9937-6986-41b4-8406-74701267514d',
  '60a9ce85-58a6-4a3a-a34e-5be41e6c1b9e',
  '5416ac01-d888-4ecf-814a-497cbdaf45e0',
  '1b88b184-bc65-4406-a286-0fc94d61eea4')]
scalenzlog <- function (x) log10(x+1.0)-mean(log10(x[x>0]))
for (i in 2:length(colnames(gx2))) { gx2[i] <- scalenzlog(gx2[i]) }
gx2$ensembl_gene_id <- gsub('\\.[0-9]+$','',levels(gx2$feature)[gx2$feature])

dnaRnaMap <- read.table('laml-dna-rna.tsv',sep='\t',header=T,stringsAsFactors=F)

dhsgenenames <- getBM(attributes=c('ensembl_gene_id','hgnc_symbol','entrezgene','description'),values=levels(dhsgene$gene_name),filters='hgnc_symbol',mart=ensembl)

# look at relative expression of genes regulated by mutated TFBS for each sample
# also generate random sample of TFBS for comparison

tfbs_gx <- lapply(1:(length(dnaRnaMap$rna)*2), function (i) {
    ix <- floor((i+1)/2)
    dna <- dnaRnaMap$dna[ix]
    rna <- dnaRnaMap$rna[ix]
    if (mod(i,2)==1) {
        v <- tfbs_count[dna][[1]]
    } else {
        v <- sample(tfbs, length(tfbs_count[dna][[1]]$ov),replace=F,prob=width(tfbs))
        v$ov <- 1
        v$tfgene <- v$gene
    }
    q <- findOverlaps(v, dhsgene)
    df <- data.frame(tfgene=v$tfgene[queryHits(q)],ov=v$ov[queryHits(q)],
        gene_name=dhsgene$gene_name[subjectHits(q)],cor=dhsgene$cor[subjectHits(q)])
    df2 <- df %>% group_by(gene_name) %>%
        summarise(mean_cor=mean(cor),sd_cor=sd(cor),iqr_cor=IQR(cor),sum_ov=sum(ov))
    df3 <- df2 %>% rename(hgnc_symbol=gene_name) %>% inner_join(dhsgenenames, by='hgnc_symbol')
    result <- inner_join(df3, gx2[c('ensembl_gene_id', rna)], by='ensembl_gene_id')
    rnax <- result[[rna]]
    rnax <- rnax[rnax>min(rnax)+.001]
    stats <- c(length(v$ov), length(unique(v$tfgene)), length(subjectHits(q)), length(unique(subjectHits(q))),
            length(df3$hgnc_symbol), mean(rnax), median(rnax))
    print(stats)
    list(df=result,stats=stats)
})
tfbs_gx_stats <- t(sapply(tfbs_gx, . %>% .$stats))
colnames(tfbs_gx_stats) <- c('tfbs_n', 'tfbs_uq', 'tfbs_dhs_n', 'tfbs_dhs_uq', 'tfbs_ensembl_n', 'gx_mean', 'gx_median')
tfbs_gx_stats <- as.data.frame(tfbs_gx_stats)
tfbs_gx_stats$rna <- sapply(1:length(tfbs_gx_stats$tfbs_n), function (i) dnaRnaMap$rna[1+floor((i-1)/2)])
tfbs_gx_stats$random <- sapply(1:length(tfbs_gx_stats$tfbs_n), function (i) mod(i,2)==0)
tfbs_gx <- lapply(tfbs_gx, . %>% .$df)
tfbs_gx_stats %>% group_by(random) %>% summarise(mult_tfbs=mean(tfbs_dhs_n/tfbs_dhs_uq),mean_gx=mean(gx_mean),mult_dhs=mean(tfbs_dhs_n/tfbs_n))

plots <- lapply(1:2*length(tfbs_gx), function(i) {
    dfx <- tfbs_gx[[floor((i+1)/2)]]
    rna <- dnaRnaMap$rna[floor((i+3)/4)]
    if (mod(i,2)==1) {
        ggplot(dfx, aes(x=dfx$mean_cor,y=dfx[[rna]])) +
            geom_point(aes(size=dfx$sum_ov)) + scale_size_area() +
            geom_text(aes(label=dfx$hgnc_symbol))
    } else {
        ggplot(dfx, aes(x=dfx[[rna]])) + geom_histogram(binwidth=0.2)
    }
})


pdf('plots-gx.pdf')
par(mfrow=c(2,2))
for (i in 2:length(gx2)) {
    hist(gx2[,i], main=colnames(gx2)[i])
}
dev.off()

#
# look at TCGA-505 data from ICGC for comparison
#
# has both GBM data and matched random mutations
#

gx <- read.table('expression.tsv.gz',sep='\t',stringsAsFactors=F,header=T)
colnames(gx) <- gsub('\\.', '-', colnames(gx))
dna <- read.table('c:/ravip/pancan/tcga505/mutations.tsv.gz', sep='\t', header=T, stringsAsFactors=F)
dna <- makeGRangesFromDataFrame(dna, keep.extra.columns = T, start.field = 'pos', end.field = 'pos')
seqlevelsStyle(dna) <- 'UCSC'
barcodes <- unique(dna$barcode)
scalenzlog <- function (x) log10(x+1.0)-mean(log10(x[x>0]+1))
for (i in 3:length(colnames(gx))) { gx[i] <- scalenzlog(gx[i]) }
tfbs_gx <- lapply(1:(length(barcodes)*2), function (i) {
    bc <- barcodes[floor((i+1)/2)]
    v <- dna[mcols(dna)$barcode==bc]
    if (mod(i,2)==1) {
        v <- data.frame(gene=tfbs$gene, ov=countOverlaps(tfbs, v),
                width=width(tfbs), chr=seqnames(tfbs), start=start(tfbs), end=end(tfbs)) %>%
            mutate(rate=1.0e6*ov/width, tfgene=gene) %>%
            select(-width) %>%
            makeGRangesFromDataFrame(keep.extra.columns=T)
    } else {
        v <- sample(tfbs, length(start(v)),replace=F,prob=width(tfbs))
        v$ov <- 1
        v$tfgene <- v$gene
    }
    q <- findOverlaps(v, dhsgene, maxgap=0)
    df <- data.frame(tfgene=mcols(v)$tfgene[queryHits(q)],ov=mcols(v)$ov[queryHits(q)],
        gene_name=mcols(dhsgene)$gene_name[subjectHits(q)],cor=mcols(dhsgene)$cor[subjectHits(q)])
    df2 <- df %>% group_by(gene_name) %>%
        summarise(mean_cor=mean(cor),sd_cor=sd(cor),iqr_cor=IQR(cor),sum_ov=sum(ov)) %>%
        rename(symbol=gene_name)
    result <- inner_join(df2, gx[c('symbol', bc)], by='symbol')
    rnax <- result[[bc]]
    rnax <- rnax[rnax>min(rnax)+.001]
    stats <- c(length(v$ov), length(unique(v$tfgene)), length(subjectHits(q)), length(unique(subjectHits(q))),
            length(df2$symbol), mean(rnax), median(rnax))
    print(stats)
    list(df=result,stats=stats)
})
