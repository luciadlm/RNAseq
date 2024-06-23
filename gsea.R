library(GSEABase)
library(fgsea)
library(reactome.db)
library(pathview)

library(gprofiler2)
library(clusterProfiler)
require(DOSE)
library(enrichplot)
library(europepmc)

#' as.png
#'
#' run specified code to generate a graphic and save as PNG file
#' 
#' @param
#' 
#' @return
#'
#' @usage
#' 
#' @examples
#'
#' @author	(C) José R. Valverde, CNB-CSIC, 2019
#'
#' @license	EU-GPL
#'
#' @export
#
as.png <- function(PLOT=NULL, 
               file='out.png', width=1024, height=1024, 
               overwrite=TRUE, verbose=T) {

    if (is.null(file)) {
        tryCatch(print(PLOT))
    } else if ( overwrite || ! file.exists(file) ) {
        if (verbose){
            cat("as.png(): creating", file, "\n")
        }
	    tryCatch( {
                png(file, width=width, height=height)
                print(PLOT)
            },
            finally = dev.off()
        )
    }
    return()
}



###################################
# ------------------------------- #
# Do Gene Set Enrichment Analysis #
# ------------------------------- #
###################################

GO_fgsea <- function (ann.shrunk.lfc, 
		      max.size=250,
                      out.dir=paste(folder, 'go_fgsea', sep='/'), 
                      out.name='GO_fgsea',
                      use.description=TRUE,
                      top.n=20,
                      top.biblio=5,
                      verbose=FALSE) {
    # Do GSEA on GO terms using fgsea

    # Rank all genes on their fold change.
    #	Here we exclude genes for which we have no EntrezID and
    #	use shrunk LFC values
#    gseaDat <- filter(ann.shrunk.lfc, !is.na(ENTREZID))
    gseaDat <- filter(ann.shrunk.lfc, !is.na(GENEID))

    ranks <- gseaDat$lfc
    #names(ranks) <- gseaDat$ENTREZID
    names(ranks) <- gseaDat$GENEID
    head(ranks)
    #uranks <- ranks[!duplicated(sort(names(ranks)))]

    # plot all the ranked fold changes
    out.png <- paste(out.dir, '/', out.name, '.barplot.png', sep='')
    as.png(barplot(sort(ranks, decreasing=T)), out.png)
    
    # load pathways
    #pathways <- ann.shrunk.lfc$ENTREZID
    #pathways <- ann.shrunk.lfc$go_id
    #names(pathways) <- gseaDat$ENTREZID
    #names(pathways) <- gseaDat$GENEID
    #upathways <- pathways[!duplicated(sort(names(pathway)))]
    
    # create a list of go_id terms where each term contains a vector
    # of emsembl_gene_id in that term
    #   first recover the annotation (in case we are re-run and do not
    #   have it
    if ( ! exists(substitute(bm.go.annot)) ) {
        bm.go.annot <- read.table(paste(folder, '/biomart.go.tab', sep=''), 
	            header=T, sep='\t')
    }
    #if ( ! exists(substitute(bm.goslim.annot)) ) {
    #    bm.goslim.annot <- read.table(paste(folder, '/biomart.goslim.tab', sep=''), 
    #	            header=T, sep='\t')
    #}

    if (use.description == FALSE) {
        # split() will divide the gene-ids by their go-id
        pathways <- split(bm.go.annot$ensembl_gene_id, bm.go.annot$go_id)
    } else {
        # Do the same but with large gene ontology names
        # split() will divide the gene-ids by their go-name
        pathways.go <- split(bm.go.annot$ensembl_gene_id, bm.go.annot$name_1006)
    }
    # do analysis
    # the resulting table contains enrichment scores and p-values
    out.file <- sprintf("%s/go_fgsea_10-%d.RData", out.dir, max.size)
    if (file.exists(out.file)) {
        load(file=out.file)
    } else {
        fgseaRes <- fgsea(pathways=pathways.go, 
                          stats=ranks, 
		          minSize=10, 
		          maxSize=max.size, 
                          nPermSimple=100000 
                          )
        save(fgseaRes, file=out.file)
    }
    if (verbose == TRUE)
        head(fgseaRes[order(padj, -abs(NES)), ], n=10)

    if (verbose == TRUE) {
        # plot enrichment score
        sorted.fgsea.res <- fgseaRes[order(padj, -abs(NES)), ]
        sfr.names <- sorted.fgsea.res$pathway
        for (i in 1:top.n) {
            if (use.description == FALSE)
                descr <- bm.go.annot[bm.go.annot$go_id == sfr.names[i], "name_1006"]
            else
                descr <- sfr.names[i]
            print(
                plotEnrichment(pathways.go[[ sfr.names[i] ]], ranks) +
                    labs(title=descr)
                )
            ans <- readline("Press RETURN to continue: ")
            if (ans == "q") break
        }
    }
    # gsea table plot of top.n gene families
    #   top_n() is now deprecated
    topUp <- fgseaRes %>%
        filter(ES > 0) %>%
        top_n(top.n, wt=-padj)

    topDown <- fgseaRes %>%
        filter(ES < 0) %>%
        top_n(top.n, wt=-padj)

    topPathways <- bind_rows(topUp, topDown) %>%
        arrange(-ES)
    # last resort when "pos" is used    
    #topPathways <- sorted.fgsea.res[1:top.n, ]

    # do the plots and save descriptions
    out.file <- paste(out.dir, "/topUp.", top.n, '.txt', sep='')
    for (i in 1:top.n) {
        name <- topUp[i]$pathway
        if (use.description == FALSE)
            descr <- unique(bm.go.annot[bm.go.annot$go_id == name, "name_1006"])
        else
            descr <- name
        cat(i, descr, file=out.file, '\n', sep='\t', append=TRUE)
        #print(name)
        #out.png <- paste(out.dir, '/', out.name, '.topUp', i, '.enrichment.png', sep='')
        out.png <- sprintf("%s/%s.topUp.%03d.enrichment.png", 
        	out.dir, out.name, i)
        out.plt <-
            plotEnrichment(pathways.go[[ name ]] , ranks) +
            labs(title=descr)
        ggsave(out.png, out.plt)
    }
    out.file <- paste(out.dir, "/topDn.", top.n, '.txt', sep='')
    for (i in 1:top.n) {
        name <- topDown[i]$pathway
        if (use.description == FALSE)
            descr <- unique(bm.go.annot[bm.go.annot$go_id == name, "name_1006"])
        else
            descr <- name
        cat(i, descr, file=out.file, '\n', sep='\t', append=TRUE)
        #print(name)
        #out.png <- paste(out.dir, '/', out.name, '.topDown', i, '.enrichment.png', sep='')
        out.png <- sprintf("%s/%s.topDn.%03d.enrichment.png", 
        	out.dir, out.name, i)
        out.plt <-
            plotEnrichment(pathways.go[[ name ]] , ranks) +
            labs(title=descr)
        ggsave(out.png, out.plt)
    }

    out.png <- paste(out.dir, '/', out.name, '.GSEAtable.png', sep='')
    as.png(
        plotGseaTable(pathways.go[topPathways$pathway], 
                      ranks, 
                      fgseaRes, 
                      gseaParam = 0.5)
    , out.png, width=1024, height=100*top.n, overwrite=TRUE)
    
    
    # and now do a plot of the interest in citations during
    # the last ten years for the top 5 sets
    out.png <- paste(out.dir, '/', out.name, '.EUPMC.png', sep='')
    cur.year <- as.integer(format(Sys.Date(), "%Y"))
    terms <- topDown[1:top.biblio]$pathway
    as.png(
        pmcplot(terms, (cur.year-10):cur.year, proportion=FALSE),
        out.png)

    
    # finally, return fgseaRes
    return(fgseaRes)
}



GO_KEGG_clusterProfiler <- function(ann.shrunk.lfc, 
		      max.size=250,
                      out.dir=paste(folder, 'go_fgsea', sep='/'), 
                      out.name='GO_fgsea',
                      use.description=TRUE,
                      OrgDb = org.Cjaponica.eg.db,
                      kegg_organism = "cjo",	# (cjo = coturnix japonica)
                                             # (gga = gallus gallus)
                      top.n=10,
                      top.biblio=5,
                      verbose=FALSE) {
        
    gseaDat <- filter(ann.shrunk.lfc, !is.na(ENTREZID))
    ranks <- gseaDat$logFC
    names(ranks) <- gseaDat$ENTREZID
    ranks<-na.omit(ranks)

    s.ranks <- sort(ranks, decreasing=T)
    
    ### G O annotation
    #
    gse <- gseGO(geneList=s.ranks, 
                 ont ="ALL", 
                 keyType = "ENTREZID", 
                 minGSSize = 3, 
                 maxGSSize = 800, 
                 pvalueCutoff = 0.05, 
                 verbose = TRUE, 
                 OrgDb = OrgDb, 
                 pAdjustMethod = "fdr")
                 #pAdjustMethod = "none")
    if (dim(gse)[1] == 0) {
        # Try without correction issuing a warning
        gse <- gseGO(geneList=s.ranks, 
                 ont ="ALL", 
                 keyType = "ENTREZID", 
                 minGSSize = 3, 
                 maxGSSize = 800, 
                 pvalueCutoff = 0.05, 
                 verbose = TRUE, 
                 OrgDb = OrgDb, 
                 #pAdjustMethod = "fdr")
                 pAdjustMethod = "none")
        out.name <- paste(out.name, '.raw_p', sep='')
    }
    if (dim(gse)[1] == 0) {
        # give up
        cat('', file=paste(out.dir, '/',
                out.name, 'NO_ENRICHED_GO', sep=''))
#        return()
    } else {
        out.png <- paste(out.dir, '/', out.name, '.GOdotplot.png', sep='')
        as.png(
            dotplot(gse, showCategory=20, split=".sign") + facet_grid(.~.sign)
            , out.png)

        out.png <- paste(out.dir, '/', out.name, '.GOemmaplot.png', sep='')
        as.png(
            emapplot(pairwise_termsim(gse), showCategory = top.n)
            , out.png)

        out.png <- paste(out.dir, '/', out.name, '.GOridgeplot.png', sep='')
        as.png(
            ridgeplot(gse) + labs(x = "enrichment distribution")
            , out.png)

        out.png <- paste(out.dir, '/', out.name, '.EUPMC.png', sep='')
        cur.year <- as.integer(format(Sys.Date(), "%Y"))
        terms <- gse$Description[1:top.biblio]
        pmcplot(terms, (cur.year-10):cur.year, proportion=FALSE)

        for (i in 1:dim(gse)[1] ) {
        out.png <- paste(out.dir, '/', 
                out.name, '.GOcnetplot.', i, '.', gse$ID[i], '.png', sep='')
        # categorySize can be either 'pvalue' or 'geneNum'
        as.png(
            cnetplot(gse, categorySize="pvalue", foldChange=s.ranks, 
                     showCategory=i)
            , out.png)

            # Use the `Gene Set` param for the index in the title, and as the value for geneSetId
            out.png <- paste(out.dir, '/', 
                    out.name, '.GOgseaplot.', i, '.', gse$ID[i], '.png', sep='')
            gseaplot(gse, by = "all", title = gse$Description[i], geneSetID = i)
        }
        out.file <- paste(out.dir, '/', out.name, '.topGO.tab', sep='')
        write.table(gse, file=out.file, row.names=T, col.names=T, sep='\t')
    }
    
    
    ### K E G G annotation
    # let's try with KEGG (from the ENTREZID which is the same a ncbi-genid)
    kse <- gseKEGG(geneList     = s.ranks,
               organism     = kegg_organism,
               #nPerm        = 10000,
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "fdr",
               keyType       = "ncbi-geneid",
               nPermSimple = 100000)
               
    if (dim(kse)[1] == 0) {
        # Try without correction issuing a warning
        kse <- gseKEGG(geneList     = s.ranks,
                   organism     = kegg_organism,
                   #nPerm        = 10000,
                   minGSSize    = 3,
                   maxGSSize    = 800,
                   pvalueCutoff = 0.05,
                   pAdjustMethod = "none",
                   keyType       = "ncbi-geneid",
                   nPermSimple = 100000)
        out.name <- paste(out.name, 'raw_p', sep='')
    }
    if (dim(kse)[1] == 0) {
        # give up
        cat('', file=paste(out.dir, '/',
                out.name, 'NO_ENRICHED_KEGG', sep=''))
#        return()
    } else {
        out.png <- paste(out.dir, '/', out.name, '.KEGGdotplot.png', sep='')
        as.png( 
            dotplot(kse, showCategory = 10, title = "Enriched Pathways" , 
                    split=".sign") + facet_grid(.~.sign)
            , out.png)

        out.png <- paste(out.dir, '/', out.name, '.KEGGemmaplot.png', sep='')
        as.png(
            emapplot(pairwise_termsim(kse))
            , out.png)


        # categorySize can be either 'pvalue' or 'geneNum'
        out.png <- paste(out.dir, '/', out.name, '.KEGGcnetplot.png', sep='')
        as.png(
            cnetplot(kse, categorySize="pvalue", foldChange=s.ranks)
            , out.png)

        out.png <- paste(out.dir, '/', out.name, '.KEGGridgeplot.png', sep='')
        as.png(
            ridgeplot(kse) + labs(x = "enrichment distribution")
            , out.png)

        cur.dir <- getwd()
        setwd(out.dir)
        for (i in 1:dim(kse)[1]) {
            # for each of the pathways in kse

            # Use the `Gene Set` param for the index in the title, and as the value for geneSetId
            out.png <- paste(out.dir, '/', 
                out.name, '.GSEAplot.', i, '.', kse$ID[i], '.png', sep='')
            gseaplot(kse, by = "all", title = kse$Description[i], geneSetID = i)

            # Produce the native KEGG plot (PNG)
            dme <- pathview(gene.data=s.ranks, pathway.id=kse$ID[i], species = kegg_organism)

            # Produce a different plot (PDF) (different from the previous one)
            dme <- pathview(gene.data=s.ranks, pathway.id=kse$ID[i], species = kegg_organism, kegg.native = F)
        }
        setwd(cur.dir)

        out.file <- paste(out.dir, '/', out.name, '.topKEGG.tab', sep='')
        write.table(kse, file=out.file, row.names=T, col.names=T, sep='\t')
    }
}


# we do not know what is the best upper limit, so we'll try several
ms <- 500 # default value in GSEA, should work as well as the others
for (ms in c(50, 100, 250, 500)) {
    for (cmp.name in names(ds.data)) {
        # for each comparison name
        cmp.data <- ds.data[[ cmp.name ]]
        cat('Doing fGSEA on', cmp.name, '\n') 
        ann.shrunk.lfc <- cmp.data[[ "shrunk.annot" ]]
        #out.dir <- paste(folder, "DESeq2/GO_fgsea", cmp.name, sep='/')
        out.dir <- sprintf("%s/DESeq2/GO_fgsea/max.size=%03d/%s",
                folder, ms, cmp.name)
        dir.create(out.dir, showWarning=FALSE, recursive=TRUE)
        gogsea <- GO_fgsea(ann.shrunk.lfc, 
                   max.size=ms,
                   out.dir=out.dir, 
                   out.name='fgsea',
                   use.description=TRUE)
        # try dotplot(gogsea), emapplot(gogsea)... etc from clusterProfiler
        # we should add fgsea results to the cmp.data list
        ds.data[[cmp.name]][['go.fgsea']] <- gogsea

        out.dir <- sprintf("%s/DESeq2/GO+KEGG_cProf/max.size=%03d/%s",
                folder, ms, cmp.name)
        dir.create(out.dir, showWarning=FALSE, recursive=TRUE)
        GO_KEGG_clusterProfiler(
                      ann.shrunk.lfc, 
		      max.size=ms,
                      out.dir=out.dir, 
                      out.name='cProf',
                      use.description=TRUE,
                      verbose=FALSE)
    }
}
