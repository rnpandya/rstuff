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

dnaRnaMap <- read.table('c:/ravip/beataml/r/laml-dna-rna.tsv',sep='\t',header=T,stringsAsFactors=F)

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
tfbs_gx_stats %>% group_by(random) %>% summarise(mult_tfbs=mean(tfbs_dhs_n/tfbs_dhs_uq),mean_gx=mean(gx_mean),mult_dhs=mean(tfbs_dhs_n/tfbs_n))
tfbs_gx <- lapply(tfbs_gx, . %>% .$df)

tfbs_gx_all <- do.call(rbind, lapply(tfbs_gx,
      function (x) x %>% gather('sample', 'expr', 9, na.rm=T))) %>%
  rename(sample=variable, expr=value)

tfbs_gx_stats <- tfbs_gx_all %>% group_by(hgnc_symbol) %>%
  summarise(sum_ov=sum(sum_ov),mean_ov=mean(sum_ov),median_ov=median(sum_ov),n_ov=n_distinct(sample)) %>% arrange(desc(mean_ov))

tfbs_gx_stats %>% filter(n_ov>2 & mean_ov>3) %>% .$hgnc_symbol %>% unique %>%
  write.table(col.names=F,row.names=F,quote=F,'c:/ravip/beataml/r/g3.txt')

# run g3.txt through GO enrichment at pantherdb.org
# 1st & 4th highest fold enrichment are myeloid differentiation
# genes from GO terms are
tfbs_myeloid <- c('GAB3', 'CBFA2T3', 'GATA2', 'FOXP1', 'SBNO2', 'IRF4', 'CCR7', 'TGFB1', 'TNFSF9',
    'EPHA2', 'LTBR', 'JAGN1', 'JMJD6', 'ERCC2', 'RUNX1', 'DNASE2', 'PTPN6', 'HOXB7')
# samples by #genes with TFBS muts
tfbs_gx_target <- tfbs_gx_all %>% filter(hgnc_symbol %in% tfbs_myeloid) %>% group_by(sample) %>%
    summarise(n_targets=n_distinct(hgnc_symbol)) %>% rename(rna=sample)

# aml genes
aml_genes <- c('NPM1', 'FLT3', 'DNMT3A', 'IDH1', 'IDH2', 'NRSA', 'KRAS', 'RUNX1', 'TET2', 'TP53', 'CEBPA',
      'WT1', 'PTPN11', 'KIT')
laml_entrez <- getBM(attributes = "entrezgene", filter="external_gene_name", values = aml_genes, mart=ensembl)
txdb_laml <- txdb_genes[txdb_genes$gene_id %in% sapply(laml_entrez, as.character),]
as.granges <- function (g) makeGRangesFromDataFrame(as.data.frame(g))
laml_ov <- sapply(vcfs, function (v) countOverlaps(as.granges(v), as.granges(txdb_laml)))

mut_targets <- data.frame(n_muts=laml_ov, dna=names(laml_ov)) %>%
    inner_join(dnaRnaMap) %>% inner_join(tfbs_gx_target)
ggplot(mut_targets, aes(x=n_muts, y=n_targets)) + geom_bin2d() + geom_density2d()

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

#setwd('c:/ravip/pancan/tcga505')
setwd('d:/sequence/pancan/tcga505')

gx <- read.table('expression.tsv.gz',sep='\t',stringsAsFactors=F,header=T)
colnames(gx) <- gsub('\\.', '-', colnames(gx))
scalenzlog <- function (x) {
  f <- x>0 & ! is.na(x)
  m <- mean(log10(x[f]+1))
  ifelse(f, log10(x+1)-m, m)
}
gx %<>% mutate_each(funs(scalenzlog), -ENSEMBL_ID, -symbol)
ROWLIMIT=10000 # -1 for all
dna <- read.table('d:/sequence/pancan/tcga505/mutations.tsv.gz', sep='\t', header=T, stringsAsFactors=T, nrows=ROWLIMIT) %>%
  makeGRangesFromDataFrame(keep.extra.columns = T, start.field = 'pos', end.field = 'pos')
seqlevelsStyle(dna) <- 'UCSC'

dna_rnd <- read.table('d:/sequence/pancan/tcga505/randomised_TCGA505_toshare.txt.gz', header=F, skip=1,
        col.names=c('barcode', 'cancer', 'chr', 'pos', 'ref_allele', 'var_allele'), nrows=ROWLIMIT) %>%
    makeGRangesFromDataFrame(keep.extra.columns = T, start.field = 'pos', end.field = 'pos')
seqlevelsStyle(dna_rnd) <- 'UCSC'

barcodes <- unique(dna$barcode)

dna_all <- dna_rnd %>% as.data.frame %>% mutate(var_freq=NA, varscan_p=NA, rnd=T) %>%
    rbind(dna %>% as.data.frame %>% mutate(rnd=F)) %>%
    makeGRangesFromDataFrame(keep.extra.columns=T)

#for (i in 3:length(colnames(gx))) { gx[i] <- scalenzlog(gx[i]) }
tfbs_gx <- lapply_par(1:(length(barcodes)*2), function (i) {
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

#
# build general map of TF and target genes within some range
#

df <- read.table("c:/ravip/encode/wgEncodeRegTfbsClusteredV3.bed.gz",sep="\t") %>%
    rename(chr=V1, start=V2, end=V3, gene=V4)
tfbs <- makeGRangesFromDataFrame(df, ignore.strand = T, keep.extra.columns = T)
df <- read.table("c:/ravip/encode/allGeneCorrelations100000.p05.txt.gz",sep="\t",header=TRUE)
dhsgene <- makeGRangesFromDataFrame(df, seqnames.field = "dhs_chr", start.field = "dhs_start", end.field = "dhs_end", ignore.strand = TRUE, keep.extra.columns = TRUE)
rm(df)
ov_ <- findOverlaps(tfbs, dhsgene, maxgap=0)
tfbs_ <- as.data.frame(tfbs[queryHits(ov_)])
dhs_ <- as.data.frame(dhsgene[subjectHits(ov_)])
tfbs_target <- makeGRangesFromDataFrame(data.frame(seqnames=tfbs_$seqnames, start=tfbs_$start, end=tfbs_$end,
    tf=tfbs_$gene, tf_width=tfbs_$width, target=dhs_$gene_name, cor=dhs_$cor, pval=dhs_$pval), keep.extra.columns=T)
rm(ov_, tfbs_, dhs_)

count_tf <- function (v) v %>% as.data.frame %>%
      select(target, tf, mut_start, tf_width) %>%
      unique %>%
      group_by(target) %>%
      summarise(ntf=length(tf), wtf=sum(tf_width), rtf=length(tf)/sum(tf_width))

s_ <- dna %>% as.data.frame %>% split(., .$barcode)
dna_tf <- lapply(s_, function (d) {
    v <- makeGRangesFromDataFrame(d, keep.extra.columns=T)
    seqlevelsStyle(v) <- 'UCSC'
    match_regions(v, tfbs_target)
})
names(dna_tf) <- names(s_)
dna_tf_sum <- lapply(dna_tf, count_tf) %>% Reduce(rbind, .)
dna_tf_sum %>% filter(ntf<50) %>% .$ntf %>% hist(breaks=1:50)

s_ <- split(dna_rnd, dna_rnd$barcode)
dna_rnd_tf <- lapply(s_, function (d) {
    v <- makeGRangesFromDataFrame(d, keep.extra.columns=T)
    seqlevelsStyle(v) <- 'UCSC'
    match_regions(v, tfbs_target)
})
names(dna_rnd_tf) <- names(s_)
rm(s_)
dna_rnd_tf_sum <- lapply(dna_rnd_tf, count_tf) %>% Reduce(rbind, .)
dna_rnd_tf_sum %>% filter(ntf<50) %>% .$ntf %>% hist

dna_tf1 %>% as.data.frame %>%
  group_by(target) %>%
  summarise(ntf=length(tf), mean_cor=mean(cor)) %>%
  filter(ntf>3) %>%
  .$target %>%
  unique %>%
  write.table(file='c:/ravip/beataml/r/dna_tf1.txt', row.names=F, quote=F, col.names=F)

# count tf hits in each sample (.x=donor, .y=random)

# first pass, use inner join

dna_tf_sumall_inner <- names(dna_tf) %>% lapply(function (s)
        count_tf(dna_tf[[s]]) %>%
        inner_join(count_tf(dna_rnd_tf[[s]]), by='target') %>%
        mutate(sample=s, r_ntf=log10(.01+ntf.x/(.01+ntf.y)), rd_ntf=(ntf.x-ntf.y)/ntf.y)) %>%
    Reduce(rbind, .)
dna_tf_sumall2_inner <- dna_tf_sumall_inner %>%
    group_by(target) %>%
    summarise(n=length(sample), mean_rn=mean(r_ntf), sd_rn=sd(r_ntf), mean_rd=mean(rd_ntf), sd_rd=sd(rd_ntf))

dna_tf_sumall <- names(dna_tf) %>% lapply(function (s)
        count_tf(dna_tf[[s]]) %>%
        full_join(count_tf(dna_rnd_tf[[s]]), by='target') %>%
        mutate(ntf.x=ifelse(is.na(ntf.x),0,ntf.x), ntf.y=ifelse(is.na(ntf.y),0,ntf.y),
            rtf.x=ifelse(is.na(rtf.x),0,rtf.x), rtf.y=ifelse(is.na(rtf.y),0,rtf.y), sample=s) %>%
        mutate(r_ntf=logoddsratio(ntf.x, ntf.y), rd_ntf=(ntf.x-ntf.y)/ntf.y, mean_rd=mean(rd_ntf), sd_rd=sd(rd_ntf),
            r_rtf=logoddsratio(rtf.x, rtf.y))) %>%
    Reduce(rbind, .)
# plot shows bias to higher TF mutation counts
dna_tf_sumall2 <- dna_tf_sumall %>%
    group_by(target) %>%
    summarise(n=length(sample), mean_rn=mean(r_ntf), sd_rn=sd(r_ntf), mean_rr=mean(r_rtf), sd_rr=sd(r_rtf))

dna_tf_sumall2 %>%
    filter(n > 10) %>%
    ggplot(aes(x=mean_rn)) + geom_histogram(binwidth=.1)

# export to text, run through GO, download to tcga505-go.xlsx
# filtering for >5 in sample, top fold enrichment terms are all immune system
# but not not signif w/bonferroni correction - ??
dna_tf_sumall2 %>%
  filter(n>10&mean_nd>5) %>%
  .$target %>%
  write.table(file='dna_tf_sumall2.txt', row.names=F, quote=F, col.names=F)

# look at outer join, fraction of tf muts in real
# n>24 is top 3/4 of samples
dna_tf_sumall %<>% mutate(xf=ntf.x/(ntf.x+ntf.y))
dna_tf_sumall2a <- dna_tf_sumall %>% group_by(target) %>% summarise(n=length(sample), mean_xf=mean(xf), sd_xf=sd(xf), iqr_xf=IQR(xf))
dna_tf_sumall2a %>% filter(n>=24) %>% ggplot(aes(x=mean_xf,y=sd_xf)) + geom_point(aes(size=n))

# simply look at targets w/most TF muts
dna_tf_sumall %>%
    filter(ntf.x >= 9) %>% # top quartile
    group_by(target) %>%
    summarize(n=length(sample), mean_ntf=mean(ntf.x)) %>% arrange(desc(n)) %>%
    head(.01*length(.$n)) %>%
    .$target %>%
    write.table(file='dna_tf_sumall-freq.txt',row.names=F, quote=F, col.names=F)
# enriched for chromatin organization, cell cycle, apoptosis, cancer regulation, cellular response to stress
,.
#
# now trying predictive mapping from pathways
#

target_ids <- getBM(attributes=c('hgnc_symbol', 'ensembl_gene_id', 'entrezgene'), filters='hgnc_symbol', values=unique(tfbs_target$target), mart=ensembl) %>%
    as.data.frame %>%
    mutate(entrezgene=as.integer(entrezgene))

kegg_path <- read.table('c:/ravip/rstuff/regmut/cancer_kegg.csv',sep=',',header=T) %>%
    inner_join(KEGGPATHID2EXTID %>% as.data.frame) %>%
    rename(entrezgene=gene_or_orf_id) %>%
    mutate(entrezgene=as.integer(entrezgene))

target_path <- target_ids %>% inner_join(kegg_path)
barcode_cancer <- dna %>% as.data.frame %>% select(barcode,cancer) %>% unique

dna_tf_sumpath <- names(dna_tf) %>%
    lapply(function (s) {
        cancer_code <- (barcode_cancer %>% filter(barcode == s) %>% .$cancer)[[1]]
        targets <- target_path %>% filter(cancer == cancer_code) %>% .$hgnc_symbol
        count_tf(dna_tf[[s]]) %>%
            mutate(in_pathway=target %in% targets) %>%   
            full_join(count_tf(dna_rnd_tf[[s]]), by='target') %>%
            mutate(ntf.x=ifelse(is.na(ntf.x),0,ntf.x), ntf.y=ifelse(is.na(ntf.y),0,ntf.y),
                rtf.x=ifelse(is.na(rtf.x),0,rtf.x), rtf.y=ifelse(is.na(rtf.y),0,rtf.y),
                sample=s, cancer=cancer_code) %>%
            mutate(r_ntf=logoddsratio(ntf.x, ntf.y), rd_ntf=(ntf.x-ntf.y)/ntf.y, mean_rd=mean(rd_ntf), sd_rd=sd(rd_ntf),
                r_rtf=logoddsratio(rtf.x, rtf.y))}) %>%
    Reduce(rbind, .)
# ok, this pretty much disproves my hypothesis - these graphs look very similar
dna_tf_sumpath %>%
    filter(abs(r_ntf)<6) %>%
    ggplot(aes(x=r_ntf)) + geom_histogram(binwidth=.1) + facet_grid(in_pathway ~ ., scales='free')

#
# try predicting random vs. cancer
# based on vector of ntf by target gene
#

# todo: gather should do this, couldn't figure out how
dna_tf_rfin <- dna_tf_sumall %>% select(target, sample, ntf.x) %>% rename(ntf=ntf.x) %>% mutate(rnd=F) %>%
    rbind(dna_tf_sumall %>% select(target, sample, ntf.y) %>% rename(ntf=ntf.y) %>% mutate(rnd=T)) %>%
    spread(target, ntf) %>%
    rename(barcode=sample) %>%
    join(barcode_cancer) %>%
    mutate(flag=1) %>%
    spread(cancer,flag) %>%
    select(-barcode)

#
# Ernst enh-prom files
#

setwd('c:/ravip/rstuff/regmut/data/ernst')
ernst <- list.files(, '*.txt') %>% lapply(function (f)
    read.table(f, sep='\t', quote=NULL, header=F, col.names=c('seqnames', 'start', 'end', 'refseqtx')) %>%
        mutate(cellline=strsplit(f, '_')[[1]][2], chromstate=paste(strsplit(f, '_')[[1]][3], strsplit(f, '_')[[1]][4], sep='_'))) %>%
    Reduce(rbind, .) %>%
    makeGRangesFromDataFrame(keep.extra.columns=T)
ernst$enh_width <- width(ernst)

#
# He enh-prom files
#

he <- list.files('c:/ravip/rstuff/regmut/data/he', '*.csv', full.names=T) %>% lapply(function (f)
        read.table(f, skip=1,header=T,quote='"',sep=',') %>%
            rename(seqnames=Chr, start=Start, end=End, tss.seqnames=Chr.1, tss.start=TSS, tss.strand=Strand, promoter=ID) %>%
            mutate(cellline=strsplit(f, '\\.')[[1]][1], enhancer=paste(seqnames, ':', start, '-', end, sep=''), promoter.loc=paste(tss.seqnames, tss.start, sep=':'))) %>%
    Reduce(rbind, .) %>%
    makeGRangesFromDataFrame(keep.extra.columns=T)
he$enh_width <- width(he)

# match mutations against prom-enh pairs, remove dups from diff cell lines, multiple muts in region (OK??)
dna_he <- match_regions(dna_all, he) # %>% as.data.frame %>% filter(., ! duplicated(.$enhancer, .$promoter.loc, .$barcode))

count_enh <- function (v) v %>% as.data.frame %>%
      select(promoter, enhancer, mut_start, enh_width) %>%
      unique %>%
      group_by(promoter) %>%
      summarise(nenh=length(enhancer), wenh=sum(enh_width), renh=length(enhancer)/sum(enh_width))

he_prom <- he %>% as.data.frame %>%
    group_by(promoter) %>%
    summarise(n_enh=length(enhancer), w_enh=sum(enh_width))

dna_he_sum <- dna_he %>%
    as.data.frame %>%
    group_by(promoter, barcode, rnd) %>%
    summarise(n_enh=length(enhancer), w_enh=sum(enh_width), r_enh=length(enhancer)/sum(enh_width))

dna_he_sum2 <- dna_he_sum %>%
    select(-w_enh, -r_enh) %>%
    spread(rnd, n_enh, fill=0) %>%
    group_by(promoter, barcode) %>%
    summarise(ne.sample=sum(`FALSE`), ne.random=sum(`TRUE`)) %>%
    mutate(nef=ne.sample/(ne.sample+ne.random), both=ne.sample>0 & ne.random>0)

dna_he_sum2 %>% filter(both) %>% ggplot(aes(x=nef)) + geom_histogram(binwidth=0.05)
dna_he_sum2 %>% group_by(promoter) %>% summarise(n=length(barcode), fmean=mean(nef)) %>% ggplot(aes(x=fmean)) + geom_histogram(binwidth=0.05)
dna_he_sum2 %>% group_by(promoter) %>% summarise(n=length(barcode), fmean=mean(nef)) %>% ggplot(aes(x=fmean, y=log10(n))) + geom_density2d()

# most mutated in donor vs. random
dna_he_prom %>% filter(n>10) %>% arrange(desc(fmean)) %>% head(200) %>% .$promoter %>% write.table(file='dna_he_prom.txt',row.names=F, col.names=F,quote=F)
# GO terms w/lowest p-values
#type B pancreatic cell apoptotic process
#positive regulation of epithelial cell apoptotic process
#positive regulation of type B pancreatic cell apoptotic process
#regulation of type B pancreatic cell apoptotic process
#regulation of epithelial cell apoptotic process
#epithelial cell apoptotic process
dna_he_mx <- dna_he_sum %>%
    select(-w_enh, -r_enh) %>%
    rename(ne_sample=n_enh) %>%
    left_join(he_prom %>% select(-w_enh) %>% rename(ne_prom=n_enh)) %>% ==
    mutate(fne=ne_sample/ne_prom) %>%
    select(-ne_sample,-ne_prom) %>%
    filter(!is.na(promoter)) %>%
    spread(promoter, fne, fill=0)

dna_he_mx  %>%
    mutate(rnd=ifelse(rnd,1,0)) %>%
    select(-barcode) %>%
    write.table('dna_he_mx.csv',row.names=F,sep=',',quote=F)

#
# trying motif PWM scoring
#

query(MotifDb,'ZBTB33')[[1]] -> pwm
subsetByOverlaps(dna, tfbs[1,]) -> muts
hs <- getBSgenome("BSgenome.Hsapiens.UCSC.hg19")

motifs <- levels(mcols(tfbs)$gene) %>% sapply(function (tf) MotifDb::query(MotifDb, tf)) %>%
  sapply(. %>% subset(organism=='Hsapiens'))
motifs <- motifs[sapply(motifs, length) > 0] %>% sapply(. %>% .[[1]])

names(motifs) # 94 of 161


best_score <- function (pwm, string, min.score='50%') {
  pfm <- round(pwm*100)
  d <- (nchar(string) - 2*dim(pwm)[2] - 1) / 2
  if (d > 0) {
    string <- substring(string, d + 1, nchar(string) - d)
  }
  max(mcols(matchPWM(pfm, string, with.score=T, min.score=min.score))$score/maxScore(pfm), 0)
}

best_motif <- function (motifs, string, min.score='50%') sapply(motifs,
  . %>% best_score(string, min.score=min.score)) %>% sort(decreasing=T) %>% head(1)

mutation_delta <- function (motifs, mut, min.score='50%') {
  refdna <- mut$near[[1]]
  bestref <- best_motif(motifs, refdna, min.score=min.score)
  refdna[floor(nchar(refdna)/2)+1] <- as(levels(mut$var_allele)[mut$var_allele], 'DNAString')
  bestmut <- best_score(motifs[[names(bestref)]], refdna, min.score=min.score)
  bestmut - bestref
}

#
# trying simpler approach with single TF match per
#
motifs <- subset(MotifDb, organism=='Hsapiens')
motif2 <- sapply(levels(tfbs$gene)[unique(tfbs$gene)], function(tf) head(MotifDb::query(motifs, tf),1))
motif2 <- motif2[sapply(motif2, length)>0] %>% sapply(function (m) round(m[[1]]*100))
tfbs2 <- tfbs[tfbs$gene %in% names(motif2)]
dna_tf2 <- match_regions(tfbs2, makeGRangesFromDataFrame(dna_all, keep.extra.columns = T))
dna_tf2$near <- as.character(getSeq(hs, dna_tf2+30))
dna_tf2$var_allele <- sub('-','',dna_tf2$var_allele) # get rid of - symbols
dna_tf2$var_allele <- sub('+','',dna_tf2$var_allele) # get rid of + symbols

best_score2 <- function (pfm, string) {
  d <- (nchar(string) - 2*dim(pfm)[2] - 1) / 2
  if (d > 0) {
    string <- substring(string, d + 1, nchar(string) - d)
  }
  max(mcols(matchPWM(pfm, string, with.score=T, min.score='50%'))$score/maxScore(pfm), 0)
}

mutation_delta2 <- function(motifs, mut) {
  pfm <- motifs[[as.character(mut$gene)]]
  refdna <- mut$near
  refscore <- best_score2(pfm, refdna)
  n = floor(nchar(refdna)/2)
  mutdna <- paste(substring(refdna,1,n), as.character(mut$var_allele), substring(refdna, n+1, 2*n+1), sep='')
  mutscore <- best_score2(pfm, mutdna)
  mutscore - refscore  
}

# try with a small sample
d <- dna_tf2[sample(9390551, 1000)]
system.time(d$tfdelta <- unlist(foreach(i=1:1000) %do% mutation_delta2(motif2, d[i,])))
d %>% as.data.frame %>% ggplot(aes(x=tfdelta)) + geom_histogram(binwidth=0.05) + facet_grid(rnd ~ .)
# essentially same distribution for rnd & not
# itneresting peaks around +/- 0.5 - maybe some TF that has simple matrix?
# todo: include TF name as well

# run it all!!
system.time(dna_tf2_delta <-
  foreach(i=1:length(dna_tf2$near), .packages=c('BSgenome', 'magrittr')) %dopar% mutation_delta2(motif2, dna_tf2[i,]))
dna_tf2$tfdelta <- unlist(dna_tf2_delta)
