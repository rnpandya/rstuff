df<-read.table('c:/ravip/rstuff/diffex/tp53_prelim graph.csv',sep=',',quote='"',header=T,stringsAsFactors=F)
df$single<-1:length(df$TAF)>254	
ggplot(df,aes(x=TAF,y=Tumor.Expression))+geom_point(aes(color=single))
ggplot(df,aes(x=Tumor.Expression))+geom_histogram(aes(color=single))+facet_grid(single ~ .)
scpilot <- read.table('c:/ravip/pancan/santa_cruz_pilot.v2.2015_0504.tsv',sep='\t',header=T,stringsAsFactors=F)
df2 <- df %>% rename(Tumour.WGS.submitter.sample.ID.s.=X.75) %>% 
	select(matches('^[^X].*')) %>%
	inner_join(scpilot, by='Tumour.WGS.submitter.sample.ID.s.')

dxvcf <- lapply(df2$`Tumour.WGS.aliquot.ID.s.`, function (id) {
    fname <- list.files('g:/pancan2/vcf', paste(id, '.*snv_mnv.vcf.gz', sep=''), full.names=T)
    vcf <- tryCatch(rowRanges(readVcf(fname, 'hg19')), error=function (e) GRanges())
    vcf[vcf$FILTER == 'PASS']
})

df <- read.table("c:/ravip/encode/allGeneCorrelations100000.p05.txt.gz",sep="\t",header=TRUE)
dhsgene <- makeGRangesFromDataFrame(df, seqnames.field = "dhs_chr", start.field = "dhs_start", end.field = "dhs_end", ignore.strand = TRUE, keep.extra.columns = TRUE)

tfp53 <- dhsgene[dhsgene$gene_name=='TP53',]
seqlevelsStyle(tfp53) <- 'NCBI'
dxp53mut <- lapply(dxvcf, function (v) subsetByOverlaps(v, tfp53, maxgap=150))
df2$p53tfmut <- sapply(dxp53mut, length)
# doesn't look like any connection with p53 TF binding site mutations
ggplot(df2,aes(x=TAF,y=Tumor.Expression))+geom_point(aes(size=p53tfmut,color=single))+scale_size_continuous(range=c(2,6))
ggplot(df2,aes(x=Tumor.Expression)) + geom_histogram() + facet_grid(single ~ isp53tfmut)
