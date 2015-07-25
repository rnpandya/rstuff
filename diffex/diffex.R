# TP53

df<-read.table('c:/ravip/rstuff/diffex/tp53_prelim graph.csv',sep=',',quote='"',header=T,stringsAsFactors=F)
df$single<-1:length(df$TAF)>254	
ggplot(df,aes(x=TAF,y=Tumor.Expression)) + geom_point(aes(color=single))
ggplot(df,aes(x=Tumor.Expression)) + geom_histogram(aes(color=single)) + facet_grid(single ~ .)
scpilot <- read.table('d:/sequence/pancan/santa_cruz_pilot.v2.2015_0504.tsv',sep='\t',header=T,stringsAsFactors=F)
df2 <- df %>% rename(Tumour.WGS.submitter.sample.ID.s.=X.75) %>% 
	select(matches('^[^X].*')) %>%
	inner_join(scpilot, by='Tumour.WGS.submitter.sample.ID.s.')

dxvcf <- lapply(df2$`Tumour.WGS.aliquot.ID.s.`, function (id) {
    fname <- list.files('c:/ravip/rstuff/diffex/vcf', paste(id, '.*snv_mnv.vcf.gz', sep=''), full.names=T)
    vcf <- tryCatch(rowRanges(readVcf(fname, 'hg19')), error=function (e) { print(fname); GRanges() } )
    vcf[vcf$FILTER == 'PASS']
})

df <- read.table("c:/ravip/encode/allGeneCorrelations100000.p05.txt.gz",sep="\t",header=TRUE)
dhsgene <- makeGRangesFromDataFrame(df, seqnames.field = "dhs_chr", start.field = "dhs_start", end.field = "dhs_end", ignore.strand = TRUE, keep.extra.columns = TRUE)

tfp53 <- dhsgene[dhsgene$gene_name=='TP53',]
seqlevelsStyle(tfp53) <- 'NCBI'
dxp53mut <- lapply(dxvcf, function (v) subsetByOverlaps(v, tfp53, maxgap=150))
df2$p53tfmut <- sapply(dxp53mut, length)
df2 %<>% mutate(isp53tfmut=p53tfmut > 0)
# doesn't look like any connection with p53 TF binding site mutations
ggplot(df2,aes(x=TAF,y=Tumor.Expression)) +
	geom_point(aes(size=p53tfmut,color=single)) + scale_size_continuous(range=c(2,6))
ggplot(df2,aes(x=Tumor.Expression)) + geom_histogram() + facet_grid(single ~ isp53tfmut)

df <- read.table("c:/ravip/encode/wgEncodeRegTfbsClusteredV3.bed.gz",sep="\t") %>%
    rename(chr=V1, start=V2, end=V3, gene=V4)
tfbs <- makeGRangesFromDataFrame(df, ignore.strand = T, keep.extra.columns = T)

rm(df)

dxtfbsov <- lapply(1:length(dxvcf), function (i) {
	data.frame(gene=tfbs$gene, ov=countOverlaps(tfbs, dxvcf[[i]]), width=width(tfbs),
		chr=seqnames(tfbs), start=start(tfbs), end=end(tfbs)) %>%
	mutate(rate=1.0e6*ov/width, vcf=i)})

dhstfbsov <- sapply(dxtfbsov, function (d) tryCatch(
 	sum(countOverlaps(tfp53, makeGRangesFromDataFrame(d %>% filter(ov>0) %>% select(-width)))),
		error=function(e) 0))

ggplot(df2,aes(x=TAF,y=Tumor.Expression)) +
     geom_point(aes(color=single,size=isdhstfbsov)) + scale_size_discrete(range=c(2,4))
ggplot(df2,aes(x=Tumor.Expression)) + geom_histogram() + facet_grid(isdhstfbsov ~ single, labeller = label_both)

# CDKN2A

df <- read.table('c:/ravip/rstuff/diffex/cdkn2a.tsv',sep='\t',header=T,stringsAsFactors=F)
ggplot(df, aes(x=tumorDNAFraction,y=tumorRNAFraction)) + geom_point(aes(color=isSingle))
cdkn2a_dhs <- dhsgene[mcols(dhsgene)$gene_name=='CDKN2A']

vf <- list.files('g:/pancan2/vcf', '.*snv_mnv.vcf.gz')
vf <- gsub('\\.svcp.*','',vf)
vfc <- intersect(df$Tumor_Sample_UUID, vf)
seqlevelsStyle(tfbs) <- 'UCSC'
cdkn_tfbs <- subsetByOverlaps(tfbs, cdkn2a_dhs, maxgap=0)

tfbsc <- lapply(vfc, function (f) { tryCatch
		list.files('g:/pancan2/vcf', paste(f, '.*snv_mnv.vcf.gz', sep=''), full.names = T)[1] %>%
		readVcf('hg19') %>% rowRanges %>% function (v) v[v$FILTER=='PASS']),
	error = function (e) GRanges()) })

sapply(tfbsc, . %>% (. %>% subsetByOverlaps(cdkn_tfbs) %>% length))
# 3 1 1 1 1 1 0 0 1 0 0
