library(ensembldb)		# needs to be first to avoid S4 inconsistencies
library(AnnotationHub)		# to seek annotation packages
library(AnnotationForge)	# to build our annotation package
library(GO.db)
library(PFAM.db)

library("biomaRt")		# an alternate approach to retrieve annotation

##############################################################################
#
#  PREPARE ANNOTATION SO IT IS READY WHEN WE GET THE RESULTS
#
##############################################################################

# Now is time to prepare to connect all the results we get with the existing
# information  we know from the literature, so we will retrieve infos from
# ENSEMBL and connect with the genes  we have kept. 

# we need to load the package that enables us to connect a file as an
# external database ensembldb. After that we create a variable for the DB
# storing it in the memory simple process: we extract the files from ENSEMBL
# and store the data to MySQL account and then with the script
# 'generate-EnsDBs.R' we create a sqlite file that has all the information
# needed to continue

cat("\n================================================\n")
cat(  "     B U I L D I N G   A N N O T A T I O N\n")
cat(  "================================================\n")

ah <- AnnotationHub()
# get ENSEMBL annotations for our target
#qr <- query(ah, c("EnsDb", target.organism))
# choose the most recent (last) one
#last.ref <- names(qr)[length(names(qr))]
#net.edb <- qr[[last.ref]]
#net.edb <- qr[["AH97993"]]	# Cjaponica v105

# with this we do not need to build it from GTF/GFF/mysql	

# SELECT ENSEMBL ANNOTATION DATABASE TO USE
if (use.online.annotation == TRUE) {
    dir.create('net.EnsDb', showWarnings=FALSE)
    if ( ! file.exists(paste('net.EnsDb', '/net.ens.db.RData', sep='')) ) {
        # get ENSEMBL annotations for our target
        qr <- query(ah, c("EnsDb", target.organism))
        # choose the most recent (last) one
        last.ref <- names(qr)[length(names(qr))]
        net.edb <- qr[[last.ref]]
        ens.db <- net.edb

        # save for future reference (downloading is too slow)
        saveRDS(ens.db, file=paste('net.EnsDb', '/net.ens.db.rds', sep=''))
        #	'ens.db' can later be recovered with: fc <- readRDS(file='net.ens.db.rds'
        # and save as well as Rdata file
        save(ens.db, file=paste('net.EnsDb', '/net.ens.db.RData', sep=''))
        #	'ens.db' can later be recovered with: fc <- load(file='net.ens.db.RData')
    } else {
        ens.db <- readRDS(file=paste('net.EnsDb', '/net.ens.db.rds', sep=''))
    }
} else {
    ens.db <- EnsDb(DBFile)
}

# check it
print(ens.db)
head(keys(ens.db, 'GENEID'))
columns(ens.db)


get_ens_ann <- function(gene_names, ens.db, folder='.') {
    out.file <- paste(folder, '/ensembl.annotation.txt', sep='')
    out.file.1 <- paste(folder, '/ensembl.annotation.1st.txt', sep='')
    if ( file.exists(out.file)) {
        ens.ann <- read.table(file=out.file, 
	                       sep='\t', header=T)
        ens.ann.1 <- read.table(file=out.file.1, 
	                       sep='\t', header=T)

    } else {
        ens.ann <- ensembldb::select(ens.db, 
                      column='GENID', keytype='GENEID', keys=gene_names, 
                      columns= c('SEQNAME', 'SYMBOL', 'DESCRIPTION', # no longer available
                                 'GENENAME', 'GENEID', 'ENTREZID', # empty
                                 'TXID', 'TXBIOTYPE', # these make the call fail
                                 'PROTEINID', 'UNIPROTID' # no longer available
                                )) 
        # SAVE ANNOTATION
        # ---------------
        write.table(ens.ann , file=out.file, 
	        sep='\t', row.names=T, col.names=T)

        # check if the amount of genes we have is the same as the number of 
        # the annotations that we have extracted
        if ( ! table(ens.ann$GENEID==rownames(fit)) ) {
            cat("Annotation does not match genes\n")
            cat("Using only one entry per gene\n")
            ens.ann.1 <- ens.ann[ ! duplicated(ens.ann$GENEID), ]
        } else {
            ens.ann.1 <- ens.ann
        }
        write.table(ens.ann.1, file=paste(folder, '/ensembl.annotation.1st.txt', sep=''), 
	        sep='\t', row.names=T, col.names=T)
    }
    return( list(ens.ann=ens.ann, ens.ann.1=ens.ann.1) )

}

# these should both give similar results
get_ensembl_ann <- function(gene_names) {
    return( ensembldb::select(ens.db, 
              keytype= 'GENEID', keys=rownames(fit), 
              columns= c('SEQNAME', 'SYMBOL', 'DESCRIPTION',
                         'GENENAME', 'GENEID', 'ENTREZID', 
                          'TXNAME', 'TXBIOTYPE', 
                          'PROTEINID', 'UNIPROTID')) )

} 



# ---------------------------------------------------------------
# O B T A I N   B I O M A R T   A N N O T A T I O N
# ---------------------------------------------------------------
if ( ! file.exists(paste(folder, 'biomaRt.annotation.1st.txt', sep='/')) ) {
    if ( ! file.exists(paste(folder, 'biomaRt.annotation.txt', sep='/')) ) {

        marts <- listMarts()
        head(marts)
        datasets <- listDatasets(useMart("ensembl"))
        head(datasets)

        ## set up connection to ensembl database
        ensembl=useMart("ENSEMBL_MART_ENSEMBL")

        # list the available datasets (species)
        listDatasets(ensembl) %>%  filter(str_detect(description, release))

        # already defined above
        #mart.name <- 'cjaponica_gene_ensembl'
        #mart.name <- 'ggallus_gene_ensembl'

        mart.db <- useMart("ensembl", mart.name)
        attributes <- listAttributes(mart.db)
        head(attributes)
        filters <- listFilters(mart.db)
        head(filters)

        # check the available "filters" - things you can filter for
        listFilters(mart.db) %>% 
            filter(str_detect(name, "ensembl"))

        # check the available "attributes" - things you can retreive
        attr <- listAttributes(mart.db)[,1:2]	# the first 200 are general
        listAttributes(mart.db)[,1:2] %>% 
            head(20)

        # we cannot get all the annotation at once because it times out
        #full.annot <- getBM(attributes=
        #                       c("ensembl_gene_id", "ensembl_transcript_id", 
        #			   "start_position", "end_position", 
        #                          "chromosome_name", "gene_biotype", 
        #                          "description", 
        #                          "entrezgene_id", "entrezgene_accession", "entrezgene_description", 
        #                          "go_id", "name_1006", "definition_1006", "go_linkage_type", "namespace_1003", 
        #                          "goslim_goa_accession", "goslim_goa_description", 
        #                          "pdb", 
        #                          "reactome", "uniprotswissprot"), 
        #                       mart=mart.db)

        # so we will retrieve the data in pieces, including ensembl_gene_id in
        # each piece so we can use it as key for merging the annotation
        if ( ! file.exists( paste(folder, '/biomart.ensembl.tab', sep='')) ) {
            bm.ensembl.annot <- getBM(attributes=c(
                                   "ensembl_gene_id", 
                                   "ensembl_transcript_id", 
                                   "start_position", "end_position", 
                                   "chromosome_name", "gene_biotype", 
                                   "description"),
                                mart=mart.db)
            write.table(bm.ensembl.annot, paste(folder, '/biomart.ensembl.tab', sep=''), 
	            row.names=T, col.names=T, sep='\t')
        } else {
            bm.ensembl.annot <- read.table(paste(folder, '/biomart.ensembl.tab', sep=''), 
	            header=T, sep='\t')
        }

        if ( ! file.exists(paste(folder, '/biomart.entrez.tab', sep='')) ) {
            bm.entrez.annot <- getBM(attributes=c(
                                   "ensembl_gene_id", 
                                    "entrezgene_id", "entrezgene_accession", "entrezgene_description"),
                                mart=mart.db)
            write.table(bm.entrez.annot, paste(folder, '/biomart.entrez.tab', sep=''), 
	            row.names=T, col.names=T, sep='\t')
        } else {
            bm.entrez.annot <- read.table(paste(folder, '/biomart.entrez.tab', sep=''), 
	            header=T, sep='\t')
        }

        if ( ! file.exists(paste(folder, '/biomart.go.tab', sep='')) ) {
            bm.go.annot <- getBM(attributes=c(
                                   "ensembl_gene_id", 
                                    "go_id", "name_1006", "definition_1006", "go_linkage_type", "namespace_1003"), 
                               mart=mart.db)
            write.table(bm.go.annot, paste(folder, '/biomart.go.tab', sep=''), 
	            row.names=T, col.names=T, sep='\t')
        } else {
            bm.go.annot <- read.table(paste(folder, '/biomart.go.tab', sep=''), 
	            header=T, sep='\t')
        }

        if ( ! file.exists(paste(folder, '/biomart.goslim.tab', sep='')) ) {
            bm.goslim.annot <- getBM(attributes=c(
                                   "ensembl_gene_id", 
                                    "goslim_goa_accession", "goslim_goa_description"),
                               mart=mart.db)
            write.table(bm.goslim.annot, paste(folder, '/biomart.goslim.tab', sep=''), 
	            row.names=T, col.names=T, sep='\t')
        } else {
            bm.goslim.annot <- read.table(paste(folder, '/biomart.goslim.tab', sep=''), 
	            header=T, sep='\t')
        }

        if ( ! file.exists(paste(folder, '/biomart.fam.tab', sep='')) ) {
            bm.fam.annot <- getBM(attributes=c(
                                   "ensembl_gene_id",
                                   "pfam",
                                   "pirsf",                                   "prints",
                                   "tigrfam"
                                   ), 
                               mart=mart.db)
            write.table(bm.fam.annot, paste(folder, '/biomart.fam.tab', sep=''), 
	            row.names=T, col.names=T, sep='\t')
        } else {
            bm.fam.annot <- read.table(paste(folder, '/biomart.fam.tab', sep=''), 
	            header=T, sep='\t')
        }
        if ( ! file.exists(paste(folder, '/biomart.prosite.tab', sep='')) ) {
            bm.prosite.annot <- getBM(attributes=c(
                                   "ensembl_gene_id",
                                   "scanprosite",
                                   "pfscan" 
                                   ), 
                               mart=mart.db)
            write.table(bm.prosite.annot, paste(folder, '/biomart.prosite.tab', sep=''), 
	            row.names=T, col.names=T, sep='\t')
        } else {
            bm.prosite.annot <- read.table(paste(folder, '/biomart.prosite.tab', sep=''), 
	            header=T, sep='\t')
        }
        
        if ( ! file.exists(paste(folder, '/biomart.sfam.tab', sep='')) ) {
            bm.sfam.annot <- getBM(attributes=c(
                                   "ensembl_gene_id",
                                   "superfamily"
                                   ), 
                               mart=mart.db)
            write.table(bm.sfam.annot, paste(folder, '/biomart.sfam.tab', sep=''), 
	            row.names=T, col.names=T, sep='\t')
        } else {
            bm.sfam.annot <- read.table(paste(folder, '/biomart.sfam.tab', sep=''), 
	            header=T, sep='\t')
        }

        if ( ! file.exists(paste(folder, '/biomart.extra.tab', sep='')) ) {
            bm.extra.annot <- getBM(attributes=c(
                                   "ensembl_gene_id", 
                                   "pdb",
                                   #"reactome", 
                                   "uniprotswissprot"), 
                               mart=mart.db)
            write.table(bm.extra.annot, paste(folder, '/biomart.extra.tab', sep=''), 
	            row.names=T, col.names=T, sep='\t')
        } else {
            bm.extra.annot <- read.table(paste(folder, '/biomart.extra.tab', sep=''), 
	            header=T, sep='\t')
        }

        # Now that we have all the pieces, merge them all together
        # into a single annotation variable
        #bm.annot <- merge(bm.ensembl.annot, bm.entrez.annot, by="ensembl_gene_id")
        #bm.annot <- merge(bm.annot, bm.go.annot, by="ensembl_gene_id")
        #bm.annot <- merge(bm.annot, bm.goslim.annot, by="ensembl_gene_id")
        #bm.annot <- merge(bm.annot, bm.extra.annot,  by="ensembl_gene_id")

        # SAVE ANNOTATION
        # ---------------
        write.table(bm.annot, file=paste(folder, '/biomaRt.annotation.txt', sep=''), 
	        sep='\t', row.names=T, col.names=T)
    } else {
        bm.ensembl.annot <- read.table(paste(folder, '/biomart.ensembl.tab', sep=''), 
	        header=T, sep='\t')
        bm.entrez.annot <- read.table(paste(folder, '/biomart.entrez.tab', sep=''), 
	        header=T, sep='\t')
        bm.go.annot <- read.table(paste(folder, '/biomart.go.tab', sep=''), 
	        header=T, sep='\t')
        bm.goslim.annot <- read.table(paste(folder, '/biomart.goslim.tab', sep=''), 
	        header=T, sep='\t')
        bm.fam.annot <- read.table(paste(folder, '/biomart.fam.tab', sep=''), 
	        header=T, sep='\t')
        bm.prosite.annot <- read.table(paste(folder, '/biomart.prosite.tab', sep=''), 
	        header=T, sep='\t')
        bm.sfam.annot <- read.table(paste(folder, '/biomart.sfam.tab', sep=''), 
	        header=T, sep='\t')
        bm.extra.annot <- read.table(paste(folder, '/biomart.extra.tab', sep=''), 
	        header=T, sep='\t')

        # this is too big to use
        #bm.annot <- read.table(file=paste(folder, '/biomaRt.annotation.txt', sep=''), 
	#    sep='\t', header=T)
    }

    # One possible way to do it would be to filter the queries above
    #	to retrieve the annotation matching ensembl_ids
    #
    # We can set a field to use to filter the output data
    # Set the filter type and values
    #ourFilterType <- "ensembl_gene_id"
    # and the values to select from that field
    #filterValues <- rownames(fit)
    #
    # and then obtain the specified annotation from records that match the values
    # specified in the filter field
    #fit.bm.extra.annot <- getBM(attributes=c(
    #                       "ensembl_gene_id", 
    #                       "pdb",
    #                       "reactome", 
    #                       "uniprotswissprot"), 
    #                   mart=mart.db,
    #                   filters=ourFilterType,
    #                   values=filterValues)
    #                   
    # deduplicate selecting the first annotation
    #fit.bm.extra.annot.1 <- fit.bm.extra.annot[ ! duplicated(fit.bm.extra.annot$ensembl_gene_id), ]
    # then we would repeat this for each annotation subset and merge all of them
    # at the end...

    # or we could deduplicate everything first and match aftwerards
    bm.ensembl.annot.1 <- bm.ensembl.annot[ ! duplicated(bm.ensembl.annot$ensembl_gene_id), ]
    bm.entrez.annot.1 <- bm.entrez.annot[ ! duplicated(bm.entrez.annot$ensembl_gene_id), ]
    bm.go.annot.1 <- bm.go.annot[ ! duplicated(bm.go.annot$ensembl_gene_id), ]
    bm.goslim.annot.1 <- bm.goslim.annot[ ! duplicated(bm.goslim.annot$ensembl_gene_id), ]
    bm.extra.annot.1 <- bm.extra.annot[ ! duplicated(bm.extra.annot$ensembl_gene_id), ]

    bm.annot.1 <- merge(bm.ensembl.annot.1, bm.entrez.annot.1, by="ensembl_gene_id")
    bm.annot.1 <- merge(bm.annot.1, bm.go.annot.1, by="ensembl_gene_id")
    bm.annot.1 <- merge(bm.annot.1, bm.goslim.annot.1, by="ensembl_gene_id")
    bm.annot.1 <- merge(bm.annot.1, bm.extra.annot.1,  by="ensembl_gene_id")

    write.table(bm.annot.1, 
                file=paste(folder, '/biomaRt.annotation.1st.txt', sep=''), 
	        sep='\t', row.names=T, col.names=T)
    # we save it to avoid repeating this in the future

    # or even do it all at once?
    ##bm.annot.1 <- bm.annot[ ! duplicated(bm.annot$ensembl_gene_id), ]
    ##write.table(bm.annot.1, 
    ##            file=paste(folder, '/biomart.annotation.1st.txt', sep=''), 
    ##            sep='\t', row.names=T, col.names=T)
} else {
    bm.annot.1 <- read.table(
    		file=paste(folder, '/biomaRt.annotation.1st.txt', sep=''), 
	        sep='\t', 
                header=T)
}
# as a byproduct we know at this point annotation files have been saved