library(ensembldb)		# needs to be first to avoid S4 inconsistencies
library(Rsubread)		# for read mapping

library(tibble)			# general tools
library(tidyr)
library(stringr)
library(dplyr)
library(readr)
library(keypress)

library(DESeq2)			# RNAseq with DESeq2
library(IHW)			# for p-value adjustment with IHW
library(ggplot2)
library(ashr)

ALIGN <- FALSE		# do we start from the aligned data?
BOTH <- TRUE        # is it required to pair both ends in the aligment?
use.online.annotation <- TRUE
VERBOSE <- TRUE
INTERACTIVE <- FALSE

ORG <- "COTURNIX"
ORG <- "GALLUS"

if (ORG == "COTURNIX" ) {
    # we will move inside './data/.' to work, so all paths are relative to it.
    reference <- 'Cjaponica'               # coturnix
    release <- "Coturnix_japonica_2.0"
    target.organism <- 'Coturnix_japonica'
    ncbi.taxid <- '93934'
    ncbi.genus <- 'Coturnix'
    ncbi.species <- 'japonica'
    ncbi.version <- '2.0'
    alignment.dir <- 'Alignments_Coturnix'
    ens.db.pkg <- "EnsDb.Cjaponica"
    ens.version <- "105"
    mart.name <- 'cjaponica_gene_ensembl'
    rnaseq.out <- 'rnaseq-test'
    org.package <- "org.Cjaponica.eg.db"
    KEGG_org <- 'cjo'

} else if (ORG == "GALLUS") {
    reference <- "GRCg6a"			# gallus
    #reference <- 'Gg_GRGc7b'
    release <- reference
    target.organism <- 'Gallus_gallus'
    ncbi.taxid <- '9031'
    genus <- 'Gallus'
    species <- 'gallus'
    version <- '6a'
    alignment.dir <- 'gg-aln-Rsubread-6a'
    ens.db.pkg <- "EnsDb.Ggallus"
    ens.version <- '106'
    mart.name <- 'ggallus_gene_ensembl'
    rnaseq.out <- 'rnaseq-test'
    org.package <- "org.Ggallus.eg.db"
    org.package <- "org.Gg.eg.db"
    KEGG_org <- 'gga'
}

n.genes <- 1000		# number of top genes to revise

fastq.data <- 'fastq-qc-paired'
#data <- fastq.data

out <-  alignment.dir

if (BOTH == T) {
    requireBothEnds <- T
    folder <- paste(rnaseq.out, "both_ends", sep='/')	# whether both ends should match or just any end
    logfile <- paste(folder, "log", "RNAseq.log", sep='/')
} else {
    requireBothEnds <- F
    folder <- paste(rnaseq.out, "any_end", sep='/')
    logfile <- paste(folder, "log", "RNAseq.log", sep='/')
}

# create needed directory hierarchy
dir.create(rnaseq.out, showWarnings=FALSE)
dir.create(folder, showWarnings=FALSE)
dir.create(file.path(folder, "log"), showWarning=FALSE)
dir.create(file.path(folder, "img"), showWarning=FALSE)
dir.create(file.path(folder, "DESeq2/raw"), showWarning=FALSE)
dir.create(file.path(folder, "DESeq2/shrunk"), showWarning=FALSE)
dir.create(file.path(folder, "DESeq2/signif"), showWarning=FALSE)
dir.create(file.path(folder, "DESeq2"), showWarning=FALSE)
dir.create(file.path(folder, "DESeq2/go"), showWarning=FALSE)
dir.create(file.path(folder, "DESeq2/pfam"), showWarning=FALSE)
dir.create(file.path(folder, "DESeq2/img"), showWarning=FALSE)
dir.create(file.path(folder, "DESeq2/cluster"), showWarning=FALSE)
dir.create(file.path(folder, "DESeq2/cluster/kmeans"), showWarning=FALSE)
dir.create(file.path(folder, "DESeq2/cluster/pam"), showWarning=FALSE)
dir.create(file.path(folder, "DESeq2/cluster/dbscan"), showWarning=FALSE)

as.png <- function(PLOT=NULL, 
                   file='out.png', width=1024, height=1024, 
                   overwrite=FALSE) {
    
    if (is.null(file)) {
        tryCatch(print(PLOT))
    } else if ( overwrite | ! file.exists(file)  ) {
        if (VERBOSE)
            cat("as.png(): creating", file, '\n')
        tryCatch( {
                png(file, width=width, height=height)
                print(PLOT)
            },
            finally=dev.off()
        )
    }
    return ()
}


align.fastq <- function(path, reference, aln.out, save.dir) {
    # get the fastq file names
    R1.fastq.files <- list.files(path=fastq.data, pattern='R1', full.names=TRUE)
    R2.fastq.files <- list.files(path=fastq.data, pattern='R2', full.names=TRUE)
    print(R1.fastq.files)
    print(R2.fastq.files)
    
    # we'll check if the output files exist to avoid repeating
    # work already done
    setwd(aln.out)
    dbname <- reference
    reference.fasta <- paste("./ref/", reference, ".fna", sep='')
    
    if (! file.exists(paste(reference, '.0.b.tab', sep=''))) {
        cat('\nBUILDING INDEX\n')
        # build the reference index inside the 'aln.out' directory
        buildindex(basename=dbname,reference=reference.fasta)
        dir()
    }
    
    r11 <- basename(R1.fastq.files[1])
    if (! file.exists(paste(aln.out, '/', r11, '.subread.BAM', sep=''))) {
        cat('\nALIGNING\n')
        # align the reads
        # IMPORTANT NOTE: WE NEED TWO FILES LISTING ALL THE FASTQ FILES TO ALIGN
        #	R1.fastq.files and R2.fastq.files
        align(index=dbname, readfile1=R1.fastq.files, readfile2=R2.fastq.files)
        
        # align will generate the output in the data directory, we
        # do not want to mix datasets, so we move the alignment results
        # to the corresponding output directory
        system(paste("mv ", fastq.data, "/*.BAM .", sep=""))
        system(paste("mv ", fastq.data, "/*.vcf .", sep=""))
        system(paste("mv ", fastq.data, "/*.summary .", sep=""))
    }
    setwd('..')
    
    # get the bam file names and inspect them to count the number
    # of reads that map to each genome position
    if ( ! file.exists(paste(save.dir, 'bam_files_stats.txt', sep='/'))) {
        bam.files <- list.files(path=aln.out, pattern='.BAM$', full.names = TRUE)
        bam.files
        props <- propmapped(files=bam.files)
        props
        write.table(props, file=paste(save.dir, 'bam_files_stats.txt', sep='/'), 
    	    row.names=T, col.names=T, sep='\t')
    }

}


compute.feature.counts <- function(files, reference, requireBothEnds, save.dir) {

    reference.gtf   <- paste("./ref/", reference, ".gtf", sep='')
    reference.gff   <- paste("./ref/", reference, ".gff", sep='')

    fc <- featureCounts(files=bam.files, 
	        annot.ext=reference.gtf, 
            isGTFAnnotationFile=T,
	        isPairedEnd=T, 
            requireBothEndsMapped=requireBothEnds, 
            primaryOnly=T, 
            ignoreDup=T, 
            useMetaFeatures=T)
    #  this means: process all BAM files
    #	Use as annotation the external file reference.gtf which is GTF
    #	BAM files contain paired end reads and we will only consider those
    #	where both ends match against the genome
    #    -------------------------------------------------------------
    #          ----->        <-----
    #	We will filter matches so that if a read may match more than one
    #	place in the genome, we will only consider the primary match and
    #	ignore other, duplicate matches
    #	We use meta-features
    #	We match against any feature in the GTF file, not only genes
    #	(this implies we will need to check the annotation carefully later)
    #	We do not remove chimeric fragments, but do actually count them too

    # SAVE FEATURE COUNTS
    # -------------------
    # now the variable fc contains 4 columns: annotation, target, counts and stats
    # but they exist only in the RAM memory, they are not stored somewhere safe,
    # so, we save them and create new variables so as to make it easier for us to 
    # manipulate the data
    write.csv(fc$counts, file=paste(save.dir, 'featureCounts.csv', sep='/'),
	      row.names=T, col.names=T, quote=F)
    write_delim(fc$stat, file=paste(save.dir, 'featureCounts_stat.txt', sep='/'), 
	        delim='      ')
    # save all the contents of 'fc' in an RDS file
    saveRDS(fc, file=paste(save.dir, '/featureCounts.rds', sep=''))
    #	'fc' can later be recovered with: fc <- readRDS(file='featureCounts.rds'
    # and save as well as Rdata file
    save(fc, file=paste(save.dir, '/featureCounts.RData', sep=''))
    #	'fc' can later be recovered with: fc <- load(file='featureCounts.RData')
}


#################################################################################
#
# Get the alignments and feature counts
#
#################################################################################

# Note: since we will be looking at infected cells with aberrant RNA viral genomes
# which might present recombinations of virus-host reads, we will restrict ourselves
# to matches in both ends.

if ( ALIGN == TRUE) {
    ### 
    ### JR ### This needs to be generalized
    ###
    out <- 'gg-aln-Rsubread-6a'
   # out <- 'Alignments_Coturnix'
    align.fastq(path=path, reference=reference, aln.out=out, save.dir=folder) 
} else {
    out <- 'gg-aln-Rsubread-6a'
    #out <- 'Alignments_Coturnix'
}

#bam.files <- list.files(path = out, pattern = '.bam$', full.names = TRUE)
bam.files <- list.files(path = out, pattern = '.BAM$', full.names = TRUE)

# get the feature counts
#	we use EnsEMBL annotation
#	this will tae the position counts and match them against the 
#	genes annotated in the GFF3 file for the reference, obtaining
#	a list of genes and the number of reads that mapped to them
#	as an indicator of their expression level
cat('\nCOMPUTING FEATURE COUNTS\n')

if (! file.exists(paste(folder, '/featureCounts.rds', sep=''))) {
    compute.feature.counts(files=bam.files, 
                           reference=reference, 
                           requireBothEnds=requireBothEnds, 
                           save.dir=folder) 
} else {
    fc <- readRDS(file=paste(folder, '/featureCounts.rds', sep=''))
}

#################################################################
#
#################################################################
                                                                #
#################################################################
#
#################################################################
                                                                #
#################################################################

# ---------------------------------------------------------------
# A N A L Y S I S     W I T H     D E S E Q 2
# ---------------------------------------------------------------

countDataFile <- paste(folder, "/featureCounts.csv", sep='')

countData <- read.csv(countDataFile,
		row.names=1) %>%
		as.matrix()
countData <- countData[rowSums(countData)>1, ]
head(countData)

colData <- read.delim(paste(rnaseq.out, "SampleInfo.tab", sep='/'))

colData$PFU <- as.factor(colData$PFU)
colData$Src <- as.factor(colData$Src)
colData$viral.dose <- as.factor(colData$viral.dose) #coturnix

# before doing the analysis, we will calculate normalized counts and
# save them for our record
# edgeR-type TMM normalization depends on the comparison design
# try both and see if there are differences
# sort metadata
metadata <- rbind(colData[13:15,], colData[1:3,], colData[10:12,], colData[7:9,], colData[4:6,])
rownames(metadata) <- 1:15
# sort countdata
counts <- cbind(countData[,13:15], countData[,1:3], countData[,10:12], countData[,7:9], countData[,4:6])

dds.pfu <- DESeqDataSetFromMatrix(counts, metadata,  design=~PFU)
dds.pfu <- estimateSizeFactors(dds.pfu)	# add norm factors to dds.pfu
sizeFactors(dds.pfu)
normalized_counts <- counts(dds.pfu, normalized=TRUE)	# get normalized counts
write.table(normalized_counts, file="rnaseq-test/normalized_counts.tab",
	        sep='\t', row.names=T, col.names=T) # better if row.names=F
#
tc <- t(normalized_counts)
tc <- cbind(c(rep(0, 3), rep(0.1, 3), rep(1.0,3), rep(10.0,3), rep(100.0,3)), tc)
# Note: DESeq2 doesn't actually use normalized counts, rather it uses
# the raw counts and models the normalization inside the Generalized Linear
# Model (GLM). These normalized counts will be useful for downstream
# visualization of results, but cannot be used as input to DESeq2 or any other
# tools that peform differential expression analysis which use the negative
# binomial model.
#
# To use these counts directly for variable reduction, we should first 
# select significant genes to reduce the number of explanatory variables
# and we should make sure that in the design, wt is the first element in
# countData.


name <- paste(folder, "DESeq2/dds.DESeq2.PFU", sep='/')

if ( ! file.exists(paste(name, "rds", sep='.')) ) {

    dds.pfu <- DESeqDataSetFromMatrix(countData, colData,  design=~PFU, tidy=F)

    dds.pfu <- DESeq(dds.pfu)

    # SAVE DESEQ ANALYSIS
    # -------------------
    #name <- paste(folder, "/DESeq2/dds.DESeq2", sep='')
    save(dds.pfu, file=paste(name, "RData", sep='.'))
    #	'dds' can later be recovered with: 
    #		dds <- load(file=paste(folder, "/dds.DESeq2.RData", sep=''))
    saveRDS(dds.pfu,  file=paste(name, "rds", sep='.'))
    # read as follows:    
} else {
    dds.pfu <- readRDS(file=paste(name, "rds", sep='.'))
}

cat("DESeq2 analysis produced the following fit")
print(resultsNames(dds.pfu))


name <- paste(folder, "DESeq2/dds.DESeq2.viral.dose", sep='/')
if ( ! file.exists(paste(name, "rds", sep='.')) ) {

    dds.vd <- DESeqDataSetFromMatrix(countData, colData,  design=~viral.dose, tidy=F)

    dds.vd <- DESeq(dds.vd)

    # SAVE DESEQ ANALYSIS
    # -------------------
    #name <- paste(folder, "/DESeq2/dds.DESeq2", sep='')
    save(dds.vd, file=paste(name, "RData", sep='.'))
    #	'dds' can later be recovered with: 
    #		dds <- load(file=paste(folder, "/dds.DESeq2.RData", sep=''))
    saveRDS(dds.vd,  file=paste(name, "rds", sep='.'))
    # read as follows:    
} else {
    dds.vd <- readRDS(file=paste(name, "rds", sep='.'))
}

cat("DESeq2 analysis produced the following fit")
print(resultsNames(dds.vd))

# get ens.ann from above

ens.ann <- get_ens_ann(rownames(countData), ens.db, folder)$ens.ann

ens.ann.1 <- ens.ann[ ! duplicated(ens.ann$GENEID), ]

# Add annotation to dds to keep everything in one place
#dds$ens.annot.1 <- ens.ann.1
#dds$bm.annot.1 <- b .annot.1


# SAVE SELECTED COMPARISONS
# -------------------------
# We will concentrate on the results between adjacent experiments
# where wt - PC - P - 1.0 - 0.1

dds_compare_annotate_save <- function(dds, 
	column, 
        x, 
        y, 
        filterFun=ihw, 
        alpha=0.01, 
        ensembl.ann, 
        biomart.ann,
        outDir='.' ) {

    out.file.base <- paste("raw", column, x, "vs", y, sep='_')
    # obtain comparison results
    cmp <- results(dds, 
                   contrast=c(column, x, y),
                   filterFun=filterFun,
                   alpha=alpha)
    # convert to data frame and add row names (ensembl.gene.id) 
    # as an additional column
    cmp.df.a <- data.frame(ensembl.gene.id=rownames(cmp), cmp)
    # ensembl.ann matches all ENSMEBL-ID in the cmp results 
    # (but lacks description)
    cmp.df.a <- cbind(cmp.df.a, 
                      ensembl.ann[ match(cmp.df.a$ensembl.gene.id, ensembl.ann$GENEID), ])
    # bm.annot fails to annotate some entries (why?)
    cmp.df.a <- cbind(cmp.df.a, 
                      biomart.ann[ match(cmp.df.a$ensembl.gene.id, biomart.ann$ensembl_gene_id), ])

    # and now save
    # save summary
    sink(paste(outDir, '/DESeq2/', out.file.base, '_summary.txt', sep=''), split=T)
    summary(cmp)
    sink()
    # unnanotated results object
    saveRDS(cmp, paste(outDir, "/DESeq2/", out.file.base, ".rds", sep=""))
    # annotated results as data frame (table)
    write.table(cmp.df.a,
	    paste(outDir, "/DESeq2/", out.file.base, "_annotated.tab", sep=""),
	    row.names=T, col.names=T, sep='\t')
    # histogram plot
    out.png <- paste(folder, '/DESeq2/img/DESeq2_', out.file.base, '_hist.png', sep='')
    as.png( {
        margins <- par("mar")
        par(mar=c(5, 5, 5, 5))
        hist(cmp.df.a$pvalue, main=out.file.base, breaks=1/alpha, xlab="p-value")
        par(mar=margins)
    }, out.png )

    return (cmp.df.a)
}


dds <- dds.pfu		# comparisons by column "PFU"

for (a in levels(as.factor(target$PFU)) ) {
    for (b in levels(as.factor(target$PFU)) ) {
        print(paste(a, b))
        if (a == b) next
        dds_compare_annotate_save( 
                        dds,
                        column="PFU",
                        x=a, y=b,
                        filterFun=ihw, alpha=0.01,
                        ensembl.ann=ens.ann.1,
                        biomart.ann=bm.annot.1,
                        outDir=folder
                        )

    }
}



# SAVE GENES WITH SIGNIFICANT CHANGES
# -----------------------------------

plot_and_save <- function(dds, 
                          column, x, y, 
                          filterFun=ihw, alpha=0.01,
        		  ensembl.ann, biomart.ann,
                          outDir='.',
                          save=TRUE ) {
    # base output file name
    out.base <- paste(column, ':_', x, '_x_', y, sep='')
    
    # Do the comparison and get the raw results (with baseMean, log2FC, p, padjusted)
    raw <- results(dds, 
                   contrast=c(column, x, y),
                   filterFun=filterFun,
                   alpha=alpha)
    # raw contains the data after comparing x and y as a DESEq2 result object
    # convert to a data frame and add row names (ensembl.gene.id) 
    # as an additional column
    raw.df <- data.frame(ensembl.gene.id=rownames(raw), raw)
    # ensembl.ann matches all ENSMEBL-ID in the raw results 
    # (but lacks description)
    raw.df <- cbind(
                raw.df, 
                ensembl.ann[ match(raw.df$ensembl.gene.id, ensembl.ann$GENEID), ]
                )
    # bm.annot fails to annotate some entries (why?)
    raw.df <- cbind(
                raw.df, 
                biomart.ann[ 
                  match(raw.df$ensembl.gene.id, biomart.ann$ensembl_gene_id),
                   ]
                )

    # and now save
    # save summary
    sink(paste(outDir, '/DESeq2/raw/raw_', out.base, '_summary.txt', sep=''), split=T)
    summary(raw)
    sink()
    # unnanotated results object
    saveRDS(raw, paste(outDir, "/DESeq2/raw/raw_", out.base, ".rds", sep=""))
    # annotated results as data frame (table)
    write.table(raw.df,
	    paste(outDir, "/DESeq2/raw/raw_", out.base, "_annotated.tab", sep=""),
	    row.names=T, col.names=T, sep='\t')
    # histogram plot
    out.png <- paste(folder, '/DESeq2/img/DESeq2_raw2', out.base, '_hist.png', sep='')
    as.png( {
        margins <- par("mar")
        par(mar=c(5, 5, 5, 5))
        hist(raw.df$pvalue, main=out.base, breaks=1/alpha, xlab="p-value")
        par(mar=margins)
    }, out.png )
    
    
    # now we'll shrink the data to improve visualization and ranking
    shrunk.lfc <- lfcShrink(dds, contrast=c(column, x, y), type="ashr")
    if (save) {
        ofile <- paste(outDir, "/DESeq2/img/DESeq2_", out.base, "_raw+shrunk_MA.png", sep='')
        cat("    plotting", ofile, '\n')
    } else ofile <- NULL
    as.png( {
        dim <- par("mfrow")
        par(mfrow=c(1,2))
        # plotMA shows the log2 fold changes attributable to a given variable
        # over the mean of normalized counts for all the samples in the dataset
        # Points above alpha are colored, outliers are shown as directed triangles    
        plotMA(raw, alpha=alpha)
        # it is useful to look at the shrunk l2fc values, which removes the noise
        # associated with l2fc changes from low-count genes
        plotMA(shrunk.lfc, alpha=alpha)
        # after plotMA, one may identify interesting points interactively using
        # identify() and clicking on them:
        # idx <- identify(res$baseMean, res$log2FoldChange)
        # rownames(res[idx, ])
        #
        # Alternatively, looking at the plot and deciding which coordinates are
        # of interest, one can use, e.g. 
        # res_wt_vs_PC[ (res_wt_vs_PC$log2FoldChange) > 4) 
        #               & (res_wt_vs_PC$baseMean > 1000), ]
        # or
        # res_wt_vs_PC[  (abs(res_wt_vs_PC$log2FoldChange) > 4) 
        #              & (res_wt_vs_PC$baseMean > 1000), ]
        par(mfrow=dim)
    }, ofile )
    
    # prepare for ggplot:
    #	convert to data frame
    #	add rownames as two new columns GENEID and ensembl_gene_id
    #   add annotation from EnsDb using GENEID
    #   add biomaRt annotation using ensembl_gene_id
    #   rename some columns for plotting
    shrunk.lfc$ensembl_gene_id <- rownames(shrunk.lfc)
    ann.shrunk <- as.data.frame(shrunk.lfc) %>%
    rownames_to_column("GENEID") %>% 
    left_join(ensembl.ann, "GENEID") %>% 
    left_join(biomart.ann, "ensembl_gene_id") %>% 
    rename(logFC=log2FoldChange, FDR=padj)   
    
    if (save) {
        ofile <- paste(outDir, "/DESeq2/shrunk/shrunk_", out.base, ".rds", sep='')
        cat("    saving", ofile, '\n')
        saveRDS(shrunk.lfc, file=ofile)
        ofile <- paste(outDir, "/DESeq2/shrunk/shrunk_", out.base, "_annotated.tab", sep='')
        cat("    saving", ofile, '\n')
        write.table(ann.shrunk, 
                    file=ofile,
                    row.names=T, col.names=T, sep='\t')
    }
    
    # find statistically significant changes
    signif <- raw[ raw$padj < alpha, ]
    signif$abs_lfc <- abs(signif$log2FoldChange)
    if (save) {
        ofile <- paste(folder, '/DESeq2/signif/signif_', out.base, "_α<", alpha, ".tab", sep='')
        cat("    saving", ofile, '\n')
        write.table(signif, 
                file=ofile, 
                sep='\t', row.names=T, col.names=T)
    }
    
    # order the data, firstly by the decreasing absolute 
    # value of log2FoldChange and secondly by the increasing pvalue....
    srt <- signif[ order(signif$abs_lfc,
                        signif$padj,
                        decreasing=c(T, F)), ]
    if (save) {
        ofile <- paste(folder, "/DESeq2/signif/signif_sorted_", out.base, ".tab", sep='')
        cat("    saving", ofile,'\n')
        write.table(srt, 
                file=ofile, 
	        sep='\t', row.names=T, col.names=T)
    }
    # annotate the sorted significant data
    srt.df <- data.frame(ensembl.gene.id=rownames(srt), srt)
    srt.df <- cbind(srt.df, 
                ensembl.ann[ match(srt.df$ensembl.gene.id, ensembl.ann$GENEID), ])
    srt.df <- cbind(srt.df, 
                biomart.ann[ match(srt.df$ensembl.gene.id, biomart.ann$ensembl_gene_id), ])
    if (save) {
        ofile <- paste(folder, "/DESeq2/signif/signif_sorted_", out.base, "_annotated.tab", sep='')
        cat("    saving", ofile,'\n')
        write.table(srt.df, 
                file=ofile, 
	        sep='\t', row.names=T, col.names=T)
    }

    
    return(list(result=raw, 
                shrunk=shrunk.lfc, 
                shrunk.annot=ann.shrunk, 
                signif=srt, 
                signif.annot=srt.df))
}

# Generating this takes long. We should save it and load from file
# when we run the script
ds.data <- list()
#x <- 0
contr <- "PFU"
grps <- levels(colData[ , contr ])
for (a in grps) {
    for (b in grps) {
        print(paste(a, b))
        if (a == b) next
        res <- plot_and_save( 
                        dds, contr,
                        x=a, y=b,
                        filterFun=ihw, alpha=0.01,
                        ensembl.ann=ens.ann.1,
                        biomart.ann=bm.annot.1,
                        outDir=folder,
						#save=TRUE
                        )
        print(names(res))
	name <- paste(contr, "_", a, "_", b, sep='')
        #x <- x + 1
        #ds.data[[x]] <- res
        ds.data[[name]] <- res
        #names(ds.data)[x] <- name
        #stop()
    }
}
print(names(ds.data))