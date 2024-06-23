    library(CA)
    library(vegan)
    library(ggvegan)

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
    
    ##########################################################################################
	#
    # Correspondence Analysis
	#
	##########################################################################################
	#
	# CA
	#
	# PCA maximizes the variance explained by each axis, while CA 
	# maximizes the correspondence between species and sample scores
	# (in our case the correspondence between cell types and genes).
	#
	# it should be as simple as issuing
	# data.ca <- CA(data)
	#
	# problem here is we need data to be a contingency table, i.e., the
	# sample names should be columns, and we cannot have duplicate names
	# so this might be better applied to log2FC data, although we can
	# make a contingency matrix by estimating the population counts from
	# the average of the three replicas of each sample.
	#
	### for CA using factominer
	#
	mtc <- aggregate( tc[ , 2:dim(tc)[2] ], by=list(tc$sample), FUN=mean )
	rownames(mtc) <- mtc$Group.1
	mtc$Group.1 <- NULL 
	#
	# CA will produce a plot, we will save it
	as.png(
		CA(mtc),
		'ca/mtc.ca.plot.png')
    mtc.ca <- CA(mtc),
	print(mtc.ca)
	# gives five dimensions (we have six levels)
	n.dim <- dim(mtc.ca$eig)[[1]]
	
	row <- get_ca_row(mtc.ca)
	as.png(corrplot(row$contrib, is.corr = FALSE), 
	       'impgenes/mtc_ca_contrib.png')
    
	# get eigenvalues / variances
	eig.val <- get_eigenvalue(mtc.ca)
	print(eig.val)
	as.png(
        fviz_screeplot(mtc.ca, addlabels = TRUE, ylim = c(0, 50)),
		'ca/mtc.ca.scree.png')
	as.png(
		fviz_screeplot(mtc.ca) + 
		geom_hline(yintercept = mean(eig.val), linetype = 2, color = "red"),
		'ca/mtc.ca.scree.ave.png')

	# biplot
	as.png(
    	fviz_ca_biplot(mtc.ca, repel = TRUE),
		'ca/mtc.ca.biplot.png')
	
	# graph of row variables
	# 	rows
	row <- get_ca_row(mtc.ca)
	row
    # 	row coordinates
	head(row$coord)
	#	quality
	head(row$cos2)
	#	contribution
	head(row$contrib)
	# 	plot
	as.png(
		fviz_ca_row(mtc.ca, repel = TRUE),
		'ca/mtc.ca.rows.png')
	as.png(
		fviz_ca_row(mtc.ca, col.row = "steelblue", shape.row = 15),
		'ca/mtc.ca.rows.tn.png')
	
	# quallity of representation of rows
	as.png(
		fviz_ca_row(mtc.ca, col.row = "cos2", 
					gradient.cols = c("#00AFBB", "#E7B800","#FC4E07"), 
					repel = TRUE),
		'ca/mtc.ca.rows.qc.png')
    # using transparency
	as.png(
		fviz_ca_row(mtc.ca, alpha.row = "cos2"),
		'ca/mtc.ca.rows.qc.tn.png')
	# visualize rows as big dots
	as.png(
		corrplot(row$cos2, is.corr = FALSE),
		'ca/mtc.ca.dot.rows.dot.png')
	# see as boxplot
	as.png(
		fviz_cos2(mtc.ca, choice='row', axes=1:2),
		'ca/mtc.ca.dot.rows.dot.tn.png')
	
	# contributions of rows to the dimensions
    as.png(
		corrplot(row$contrib, is.corr = FALSE),
		'ca/mtc.ca.dot.rows.contrib.png')
	# contribution of rows to dimension i
	# we shouldn't use top
	for (i in 1:n.dim) {
		as.png(
			fviz_contrib(mtc.ca, choice = "row", axes = i, top = 10),
			paste('ca/mtc.ca.dot.rows.contrib.PC', i, '.png', sep=''))
	}
	#fviz_contrib(mtc.ca, choice = "row", axes = 1, top = 10)	# DF1.I
	#fviz_contrib(mtc.ca, choice = "row", axes = 2, top = 10)	# DF1.PC DF1.PC.I
	#fviz_contrib(mtc.ca, choice = "row", axes = 3, top = 10)	# DF1
	#fviz_contrib(mtc.ca, choice = "row", axes = 4, top = 10)	# DF1.PC DF1.PC.I
	#fviz_contrib(mtc.ca, choice = "row", axes = 5, top = 10)	# DF1.P DF1.P.I
	#
	# PC2 and PC4 identify PC and PC.I cells
	# PC5 identifies P and P.I cells
	as.png(
	    fviz_contrib(mtc.ca, choice = "row", axes = c(2,4), top = 10),
	    'ca/mtc.ca.rows.contrib.dim24.png')
	# but PC2 and PC4 together also identify P cells in addition to PC and PC.I
	#
	# highlight most contributing rows
	as.png(
		fviz_ca_row(mtc.ca, col.row = "contrib", 
			gradient.cols = c("#00AFBB", "#E7B800","#FC4E07"), repel = TRUE),
		'ca/mtc.ca.rows.contrib.top.png')
	# with transparency
	as.png(
		fviz_ca_row(mtc.ca, alpha.row = "contrib", repel = TRUE),
		'ca/mtc.ca.rows.contrib.top.tn.png')
	
	# graphs of column variables
	#	NOTE:	we have too many genes to see anything clearly
	#  	As soon as we include columns (genes), we should consider
	# filtering plots by quality (adding 'select.col=list(cos2=0.8)'
	# or higher to display only most significant genes)
	#	NOTE: select.col=list(cos2=0.995|0.999) seems to allow some
	# seeing
	#
	col <- get_ca_col(mtc.ca)
	col
	head(col$coord)
	head(col$cos2)
	head(col$contrib)
	as.png(
		fviz_ca_col(mtc.ca),
		'ca/mtc.ca.cols.png')
	as.png(
		fviz_ca_col(mtc.ca, col.col = "cos2", 
			gradient.cols = c("#00AFBB", "#E7B800","#FC4E07"), repel = TRUE),
		'ca/mtc.ca.cols.qc.png')
	as.png(
		fviz_ca_col(mtc.ca, col.col = "cos2", 
		    select.col=list(cos2=0.995),
			gradient.cols = c("#00AFBB", "#E7B800","#FC4E07"), repel = TRUE),
		'ca/mtc.ca.cols.qc.dim12.0.995.png')
	as.png(
		fviz_ca_col(mtc.ca, col.col = "cos2", 
		    select.col=list(cos2=0.99), axes=c(2,4),
			gradient.cols = c("#00AFBB", "#E7B800","#FC4E07"), repel = TRUE),
		'ca/mtc.ca.cols.qc.dim24.0.99.png')
	as.png(
		fviz_cos2(mtc.ca, choice = "col", axes = 1:2),
		'ca/mtc.ca.cols.qc.tn.png')
	as.png(
		fviz_contrib(mtc.ca, choice = "col", axes = 1:2),
		'ca/mtc.ca.cols.contrib.png')

	# symmetric biplot
	as.png(
		fviz_ca_biplot(mtc.ca, repel = TRUE),
		'ca/mtc.ca.biplot.symm.png')
	# asymmetric biplot (using arrows)
	as.png(
		fviz_ca_biplot(mtc.ca, map = "rowprincipal", 
					   arrows = c(TRUE, TRUE), repel = TRUE),
		'ca/mtc.ca.biplot.asymm.arrows.png')
	# contribution biplot
	as.png(
		fviz_ca_biplot(mtc.ca, map = "colgreen", arrows = c(TRUE, FALSE), 
	            	   repel = TRUE),
		'ca/mtc.ca.biplot.contrib.png')
	
	# dimension description
	mtc.desc <- dimdesc(mtc.ca, axes = c(1,2))
	# dimension 1 on rows (cells)
	head(mtc.desc[[1]]$row, 4)
	# dimension 1 on columns (genes)
	head(mtc.desc[[1]]$col, 4)
	# etc for the rest.
	
	# prepare a minimal plot
	as.png( {
		plot(col$coord)
		text(row$coord, labels=rownames(row$coord), 
			col='red', cex=2, vfont=c("serif", "bold"))
		points(row$coord, col='yellow')
	}, file='ca/ca_dim_1+2.png')
	
	# having row$coord and col$coord we can save genes by cells/dims
	row.coord=row$coord
	col.coord=col$coord
	#ca_assign_cols_to_rows(col.coord[c('LOC107049418', 'IFNG'),], row.coord, 1)
	#gw <- ca_assign_cols_to_rows(col.coord, row.coord, 1)
	#for (i in names(gw)) cat(i, dim(gw[[i]]), '\n')
	for (i in 1:5) {
	    # get groups in that dimension
		cat('DIM', i, '\n')
		gw <- ca_assign_cols_to_rows(col.coord, row.coord, i)
		for (j in names(gw)) {
			cat('    GRP', j, dim(gw[[j]]), '\n')
			grp_genes <- gw[[j]]
			colnames(grp_genes) <- c('gene', 'distance')
			grp_genes <- enrich_genes(grp_genes, bio.ann, 
									'gene', 'entrezgene_accession')
			out <- paste('ca/CA_DIM_', i, '_', j, '.tab', sep='')
			write.table(grp_genes, file=out, row.names=F)
		}    
	}
	gw <- ca_assign_cols_to_rows(col.coord, row.coord, 1:2)
		for (j in names(gw)) {
			cat('    GRP', j, dim(gw[[j]]), '\n')
			grp_genes <- gw[[j]]
			colnames(grp_genes) <- c('gene', 'distance')
			grp_genes <- enrich_genes(grp_genes, bio.ann, 
									'gene', 'entrezgene_accession')
			out <- paste('ca/CA_DIM_1+2_', j, '.tab', sep='')
			write.table(grp_genes, file=out, row.names=F)
		}    
	gw <- ca_assign_cols_to_rows(col.coord, row.coord, c(2, 4))
		for (j in names(gw)) {
			cat('    GRP', j, dim(gw[[j]]), '\n')
			grp_genes <- gw[[j]]
			colnames(grp_genes) <- c('gene', 'distance')
			grp_genes <- enrich_genes(grp_genes, bio.ann, 
									'gene', 'entrezgene_accession')
			out <- paste('ca/CA_DIM_2+4_', j, '.tab', sep='')
			write.table(grp_genes, file=out, row.names=F)
		}    
	#
    #
	#
	### for CA using vegan
	#
	# with vegan we can use the same function as for CCA specifying
	# only the X matrix (if no Y matrix then CA is performed).
	mtc <- aggregate( tc[ , 2:dim(tc)[2] ], by=list(tc$sample), FUN=mean )
	rownames(mtc) <- mtc$Group.1
	mtc$Group.1 <- NULL 
	
	mtc.vca <- cca(X=mtc)
	# vegan is convenient because it facilitates access to results
	summary(mtc.vca)
	eigenvals(mtc.vca)
	summary(eigenvals(mtc.vca))
	# site (group) scores
	summary(mtc.vca)$sites 	# %>% round(3) %>% head()
	# species (genes) scores
	#	Note that the summary only returns the first 6 scores for each gene
	#   and site/sample, they can be extracted directly from mtc.vca
	#	but in our case that is not needed
	datos <- summary(mtc.vca)$species %>% round(3)
	datos_df <- as.data.frame(datos)
	df_sorted <- lapply(datos_df, function(x) head(sort(x, decreasing = TRUE), 30))
	
	# plot
	#	Note that as with CA the sheer amount of genes hides the
	#	BLUE dots for classes/samples
	as.png(
		plot(mtc.vca),				# same as for CA
		'ca/mtc.vca.plot.png')
	#
	# Gavin Simpson has written functions to organize the objects created 
	# by vegan in a consistent manner for graphing via ggplot2.
	# install.packages('remotes')
	# remotes::install_github("gavinsimpson/ggvegan")
	as.png(
		autoplot(mtc.vca),
		'ca/mtc.vca.aplot.png')
	
	# show sites names in the autoplot
	as.png(autoplot(mtc.vca, geom = "point") +
  		   geom_text(data = site_scores_df, 
		   aes(x = CA1, y = CA2, label = rownames(mtc)), 
		   vjust = -0.5, hjust = 0.7), 'impgenes/mtc.vca.aplot.rnam.png')

	# fortify() creates a dataframe with all the data needed to be
	# able to make ggplot2 plots
	mtc.vca.gg <- fortify(mtc.vca)
	as.png(
		ggplot(data = mtc.vca.gg,
    		   aes(x = CA1, y = CA2)) +
    		   geom_point(aes(shape = Score, colour = Score)) +
    		   facet_grid(facets = . ~ Score) +
    		   guides(shape = FALSE, colour = FALSE) +
    		   theme_bw(),
		'ca/mtc.vca.gplot.png')