source("http://bioconductor.org/biocLite.R")

if (FALSE) {
	install.packages("bitops")
	install.packages(c("DBI","RSQLite"))
	install.packages(c('BatchJobs','BBmisc','XML','RCurl'))
	biocLite("VariantTools")
	biocLite('TxDb.Hsapiens.UCSC.hg19.knownGene')
	biocLite('org.Hs.eg.db')
	biocLite('LungCancerLines')
	biocLite('BSgenome.Hsapiens.NCBI.GRCh38')
}

library(VariantTools)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(BSgenome.Hsapiens.NCBI.GRCh38)

# from vignette
library(LungCancerLines)
p53 <- gmapR:::exonsOnTP53Genome("TP53")
bams <- LungCancerLines::LungCancerBamFiles()
bam <- bams$H1993
tally.param <- TallyVariantsParam(gmapR::TP53Genome(), high_base_quality = 23L, which = range(p53) + 5e4, indels = TRUE, read_length = 75L)
called.variants <- callVariants(bam, tally.param)

# setup package for full genome
library(gmapR)
library(BSgenome.Hsapiens.NCBI.GRCh38)
gmapGenome <- GmapGenome(genome=Hsapiens,
	directory = "/mnt/resource/data/Hsapiens",
	name = "GRCh38",
	create = TRUE)
makeGmapGenomePackage(gmapGenome=gmapGenome,
	version="1.0.0",
	maintainer="ravip@microsoft.com",
	author="Ravi Pandya",
	destDir="HsapiensGmapGenome",
	license="Artistic-2.0",
	pkgName="GmapGenome.Hsapiens.NCBI.GRCh38")

#
#
# run the demo
#
#

#config <- list(localDir='/mnt/resource/data', hdfsDir='/user/hdiuser/data', gwasFile='gwascatalog-keep.tsv', ntasks=39, profile=FALSE, limit=0, shrink=0, maxsamples=2)
#config <- list(localDir='/mnt/resource/data/data10', hdfsDir='/user/hdiuser/data10', gwasFile='gwascatalog-keep.tsv', ntasks=39, profile=FALSE, limit=0, shrink=0, maxsamples=10)
#samples <- c('HG00096', 'HG00097', 'HG00099', 'HG00100', 'HG00101', 'HG00102', 'HG00103', 'HG00105', 'HG00106', 'NA21116')

# 3x13
config <- list(localDir='/mnt/resource/data/data13', hdfsDir='/user/hdiuser/data13', gwasFile='gwascatalog-keep.tsv', ntasks=39, profile=FALSE, limit=0, shrink=0, maxsamples=3)
config$samples <- c('HG00101', 'HG00102', 'HG00105')
config$chunks <- c('aa', 'ab', 'ac', 'ad', 'ae', 'af', 'ag', 'ah', 'ai', 'aj', 'ak', 'al', 'am')

# 2x19
config <- list(localDir='/mnt/resource/data/data19b', hdfsDir='/user/hdiuser/data19b', gwasFile='gwascatalog-keep.tsv', ntasks=38, profile=FALSE, limit=0, shrink=0, maxsamples=2)
config$samples <- c('HG00101', 'HG00102')
config$chunks <- c('aa', 'ab', 'ac', 'ad', 'ae', 'af', 'ag', 'ah', 'ai', 'aj', 'ak', 'al', 'am', 'an',	'ao', 'ap', 'aq', 'ar', 'as')

# 7x5
config <- list(localDir='/mnt/resource/data/data5', hdfsDir='/user/hdiuser/data5', gwasFile='gwascatalog-keep.tsv', ntasks=35, profile=FALSE, limit=0, shrink=0, maxsamples=7)
config$samples <- c('NA12282', 'NA20274', 'HG03583', 'HG04144', 'NA19222', 'HG01861', 'HG04229')
config$chunks <- c('aa', 'ab', 'ac', 'ad', 'ae')

rxSetComputeContext(RxHadoopMR(consoleOutput=TRUE))

start <- Sys.time()
results <- execVariantCalls(config)
difftime(Sys.time(), start)



#
# cleaning up GWAS catalog
#
# first, must pass NHGRI file through Excel
# and save as CSV to clean up missing columns
#

df <- read.csv('c:/ravip/build/gwascatalog.csv',stringsAsFactors=FALSE,header=TRUE)
df <- df[grep('(', df$Disease.Trait, fixed=TRUE, invert=TRUE),]
df <- df[grep('(^\\[?(NR)?\\]? ?$)|(^\\[[^[]+\\] ?$)', df$X95..CI..text.),]
df$OR.or.beta <- as.numeric(df$OR.or.beta)
df$Risk.Allele.Frequency <- as.numeric(df$Risk.Allele.Frequency)
df <- subset(df, !is.na(Chr_id) & Chr_id != 23 & !is.na(OR.or.beta) & OR.or.beta > 0 & OR.or.beta < 3)
df <- df[order(df$Disease.Trait, df$SNPs, df$OR.or.beta),]
df <- df[!duplicated(subset(df, TRUE, c(Disease.Trait, SNPs))),]
df$ALT <- toupper(substr(df$Strongest.SNP.Risk.Allele, nchar(df$Strongest.SNP.Risk.Allele), nchar(df$Strongest.SNP.Risk.Allele)))
df <- df[df$ALT %in% c('A', 'C', 'T', 'G'),]

# merge allele frequencies from additional database
extra.af<-read.csv('c:/ravip/build/extra_af.csv')
df$Risk.Allele.Frequency[is.na(df$Risk.Allele.Frequency)] <-
	merge(df, extra.af, by.x='SNPs', by.y='rsid', all.x=TRUE)$af[is.na(df$Risk.Allele.Frequency)]

# add prior incidence column where available
priors <- read.table('c:/ravip/build/priors.tsv', sep='\t', stringsAsFactors=FALSE,header=TRUE)
priors$diseaseKey <- tolower(priors$disease)
df$diseaseKey <- tolower(df$Disease.Trait)
df$incidence <- merge(df, priors, by=c('diseaseKey'), all.x=TRUE)$incidence
df$diseaseKey <- NULL

#only keep those with incidence & allele frequency
df <- subset(df, !(is.na(incidence)|is.na(Risk.Allele.Frequency)))

# calculate odds ratio adjustment for variants
AF <- df$Risk.Allele.Frequency
OR <- df$OR.or.beta
# p(disease | homo non-risk allele), approximate
p00 <- df$incidence / ((1 - AF)^2 + OR * 2 * AF * (1- AF) + (OR * AF) ^ 2)
df$logUpdate00 <- log10(p00 / df$incidence)
df$logOR <- log10(df$OR.or.beta)

# write it out
write.table(df, 'c:/ravip/build/gwascatalog-keep2.tsv', sep='\t')

bed <- data.frame(seqname=df$Chr_id, start=pmax(df$Chr_pos-2500, 0), end=df$Chr_pos+2500)
bed <- bed[order(bed$seqname, bed$start),]
write.table(bed, 'c:/ravip/build/gwascatalog-keep2.bed', sep='\t', quote=FALSE, row.names=FALSE, col.names=FALSE)
