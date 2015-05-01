
#
# multi-sample version of scripts
#

try({
	library(VariantTools)
	.libPaths(c(.libPaths(), '/mnt/resource/data/libs'))
	library(GmapGenome.Hsapiens.NCBI.GRCh38)
}, silent=TRUE)

namemap <- data.frame(sample=c(''), full=c(''))
try({
	namemap <- read.csv('namemap.csv')
	namemap$name = paste(namemap$first, namemap$last)
	if (sum(duplicated(namemap$name)) > 0) {
		warning("duplicated names in namemap.csv: ", do.call(paste, namemap$name[duplicated(namemap$name)]))
	}
}, silent=TRUE)


updateSampleList <- function (config) {
	if (length(config$chunks) == 0) {
		beds <- rxHadoopListFiles(paste(config$hdfsDir, "bed", sep="/"), intern=TRUE)
		config$chunks <- gsub(".*gwascatalog-(..).bed.*", "\\1", beds[2:length(beds)])
	}
	if (length(config$samples) == 0) {
		bams <- rxHadoopListFiles(paste(config$hdfsDir, "bam", sep="/"), intern=TRUE)
		config$samples <- gsub(".*/([^/]+)-aa.bam", "\\1", bams[grep("-aa.bam$", bams)])
		# validate results
		allchunks <- as.vector(outer(config$samples, config$chunks, paste, sep="-"))
		bamfiles <- gsub(".*/", "", bams[grep(".bam$", bams)])
		allbam <- paste(allchunks, ".bam", sep="")
		missing <- allbam[! (allbam %in% bamfiles)]
		baifiles <- gsub(".*/", "", bams[grep(".bam.bai$", bams)])
		allbai <- paste(allchunks, ".bam.bai", sep="")
		missing <- c(missing, allbai[! (allbai %in% baifiles)])
		if (length(missing) > 0) {
			warning("Missing files: ", paste(missing))
		}
	}
	    if (config$maxsamples > 0 & length(config$samples) > config$maxsamples) {
            config$samples <- config$samples[1:config$maxsamples]
    }
	config
}

loadGwas <- function (config) {
	df <- read.table(syncHdfsToLocal(config$gwasFile, config),sep="\t",header=TRUE,stringsAsFactors=FALSE)
	makeGRangesFromDataFrame(df, seqnames.field='Chr_id', start.field='Chr_pos', end.field='Chr_pos',keep.extra.columns=TRUE)
}

syncHdfsToLocal <- function (filepath, config) {
	loc <- file.path(config$localDir, filepath)
	if (! file.exists(loc)) {
		#print(c('rxHadoopCopyToLocal', paste(config$hdfsDir, filepath, sep="/"), loc))
		s <- Sys.time()
		dir.create(dirname(loc), showWarnings=FALSE, recursive=TRUE)
		rxHadoopCopyToLocal(paste(config$hdfsDir, filepath, sep="/"), loc)
		if (! file.exists(loc)) {
			stop("syncHdfsToLocal failed to create ", loc, " from ", paste(config$hdfsDir, filepath, sep="/"))
		}
	}
	loc
}

bedCache <- list()
loadBed <- function (chunk, config) {
	if (length(bedCache[[chunk]]) == 0) {
		bed.df <- tryCatch(read.table(syncHdfsToLocal(file.path("bed", paste("gwascatalog-", chunk, ".bed", sep="")), config),
			sep="\t", colClasses=c('character','integer','integer')),
			error = function (e) { m <- paste("loadBed", chunk, "failed with", e); warning(m); stop(m) })
		bedCache[[chunk]] <- makeGRangesFromDataFrame(bed.df, seqnames.field='V1', start.field='V2', end.field='V3')
		seqlevelsStyle(bedCache[[chunk]]) <- 'NCBI'
		seqlevels(bedCache[[chunk]], force=TRUE) <- seqlevels(GmapGenome.Hsapiens.NCBI.GRCh38)
	}
	bedCache[[chunk]]
}

loadBam <- function (chunk, config) {
	syncHdfsToLocal(file.path("bam", paste(chunk, ".bam.bai", sep="")), config)
	syncHdfsToLocal(file.path("bam", paste(chunk, ".bam", sep="")), config)
}

# per-node code, chunks is list of strings in form "sample-chunk"
callVariantsChunk <- function(chunks, config) {
	library(VariantTools)
	.libPaths(c(.libPaths(), '/mnt/resource/data/libs'))
	library(GmapGenome.Hsapiens.NCBI.GRCh38)
	gwas <- loadGwas(config)
	hostname <- Sys.getenv('HOSTNAME')
	# each chunk is string samplename-chromnumber
	dochunk <- function (chunk) {
		ranges <- loadBed(strsplit(chunk, "-")[[1]][2], config)
		chunkgwas <- subsetByOverlaps(gwas, ranges, ignore.strand=FALSE)
		variantParams <- TallyVariantsParam(GmapGenome.Hsapiens.NCBI.GRCh38, high_base_quality = 23L, indels = FALSE, read_length = 100L, which = ranges)
		variants <- callVariants(BamFile(loadBam(chunk, config)), variantParams, read.count=1)
		ov <- suppressWarnings(findOverlaps(chunkgwas, variants))
		gwashits <- chunkgwas[queryHits(ov)]
		varhits <- variants[subjectHits(ov)]
		# simple check for hom/het variant - VariantTools callGenotype requires cov
		gwashits$variantAlt <- alt(varhits)
		gwashits$genotype <- paste(alt(varhits),
			ifelse(varhits$raw.count.ref < varhits$raw.count.total/4, alt(varhits), gwashits$genomeref), sep="")
		gwashits <- subset(gwashits, ! duplicated(gwashits))
		gwasmiss <- subset(chunkgwas, !(chunkgwas %in% gwashits))
		gwasmiss$variantAlt <- "N"
		gwasmiss$genotype <- paste(gwasmiss$genomeref, gwasmiss$genomeref, sep="")
		result <- rbind(as.data.frame(gwashits), as.data.frame(gwasmiss))
		result$zygosity <- (result$ALT == substr(result$genotype, 1, 1)) + (result$ALT == substr(result$genotype, 2, 2))
		result$logUpdate <- result$logUpdate00 + result$zygosity * result$logOR
		result$sample <- strsplit(chunk, "-")[[1]][1]
		result$hostname <- hostname
		result #[,c('PUBMEDID', 'Disease.Trait', 'Reported.Gene.s.', 'SNPs', 'logOR', 'logUpdate', 'sample', 'incidence', 'Risk.Allele.Frequency', 'logUpdate00', 'zygosity', 'variantAlt', 'ALT', 'genotype', 'Region', 'hostname', 'Chr_id', 'Chr_pos')]
	}
	Reduce('rbind', lapply(chunks, dochunk))
}

# master code
execVariantCalls <- function(config) {
	config <- updateSampleList(config)
	allchunks <- sample(as.vector(outer(config$samples, config$chunks, paste, sep="-"))) # cross product, then random permutation
	allchunks <- lapply(1:config$ntasks, function (i) allchunks[floor((i-1)*length(allchunks)/config$ntasks+1): floor(i*length(allchunks)/config$ntasks)])
	detail <- rxExec(callVariantsChunk, rxElemArg(allchunks), config,
		execObjects = c('loadGwas', 'syncHdfsToLocal', 'bedCache', 'loadBed', 'loadBam'))
	detail <- Reduce('rbind', detail)
	detail <- subset(detail, !duplicated(detail[c('Disease.Trait','sample','Region')]	))
	detail <- merge(detail, namemap[c('sample','name')], by='sample', all.x=TRUE)
	detail$name[is.na(detail$name)] <- detail$sample
	summary <- aggregate(detail$logUpdate, by=list(detail$name, detail$Disease.Trait), FUN=sum)
	colnames(summary) <- c('sample', 'disease', 'logOR') # rename to match old columns for visualization etc.
	summary$oddsRatio <- 10^summary$logOR
	list(detail=detail, summary=summary)
}

#
# individual demo
#

individualVariants <- function (snps, gwas, details = FALSE) {
	
	# merge snps & gwas, filter
	m <- merge(gwas, snps, by.x=c('SNPs'), by.y=c('rsid'))
	m <- m[order(m$Disease.Trait, m$SNPs, m$logOR),]
	m <- m[!duplicated(subset(m, TRUE, c(Disease.Trait, SNPs))),]
	m <- subset(m, !(is.na(Risk.Allele.Frequency) | is.na(OR.or.beta)))

	# find most common disease risks, adjust odds based on prior
	# 0, 1, or 2 risk alleles
	m$zygosity <- (m$ALT == substr(m$genotype, 1, 1)) + (m$ALT == substr(m$genotype, 2, 2))
	m$logUpdate <- m$logUpdate00 + m$zygosity * m$logOR

	# aggregate and compute odds
	ma <- aggregate(m$logUpdate, by=list(m$Disease.Trait), FUN=sum)
	colnames(ma) <- c('disease', 'logUpdate')
	ma$incidence <- merge(ma, aggregate(m$incidence, by=list(m$Disease.Trait), FUN=min), by.x='disease', by.y='Group.1')$x
	ma$myIncidence <- ma$incidence * 10 ^ ma$logUpdate
	ma <- subset(ma, !is.na(ma$myIncidence))
	ma <- ma[order(abs(ma$myIncidence-ma$incidence),decreasing = TRUE),]
	if (details) m else ma
}

individualFileVariants <- function (snpfile, gwasfile, details=FALSE) {
	snps <- read.table(snpfile,sep="\t",stringsAsFactors=FALSE)
	colnames(snps) <- c('rsid', 'chr', 'pos', 'genotype')
	gwas <- read.table(gwasfile,sep="\t",header=TRUE,stringsAsFactors=FALSE)
	individualVariants(snps, gwas,details)
}
