 library(corrplot)
 library(factoextra)
 library(ctv)
 library(psych)
 library(psychTools)
 library(phyloseq)
 library(nFactors)
 
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
  # Now, we can proceed with the factor analysis:
  #  we expect to identify the most important groups of genes that determine the
  # outcome variable
  #
  ##########################################################################################
  #
  # EFA
  #
  # with EFA, we calculate which genes behave coordinately across
  # the full dataset
  #
  # We can do EFA on categorical variables if we can provide a
  # correlation matrix.
 #
  # Unfortunately this does not work with counts because of excessive
  # categories (seems to take all variables as categorical, possibly
  # bcause categories are defined as a column).
  # 
  # To do EFA we need a correlation matrix. We can generate one for
  # categorical variables using polychoric():
  # library(polycor)
  # pc.cor <- polychoric(dat)$rho
  # or
  # pc <- hetcor(data, ML=TRUE)
  # and then run EFA
  #
  # we can determine the number of factors with fa.parallel() or
  # directly:
  # fap <- fa.parallel(data, cor='poly')   # parallel analysis for dichotomous data
  # vss(pc$correlations, n.obs=N, rotate="varimax")  # very simple structure
  #
  # Then, with the polychoric correlation matrix we can use fa():
  # faPC <- fa(r=pc$correlations, nfactors=NF, n.obs=N, rotate="varimax")
  # faPC$loadings
  #
  # it was once also possible to use fa.poly() from 'psych')
  # (but it is now deprecated in favor of fa() )
  # faPCdirect <- fa.poly(data, nfactors=NF, rotate="varimax")
  # faPCdirect$fa$loadings
  #
  # and do plots:
  # factor.plot(faPCdirect$fa, cut=0.5)
  # fa.diagram(faPCdirect)
  #
  # For factor scores, look at package ltm which has a factor.scores() 
  # function specifically for polytomous outcome data.
  # see https://github.com/drizopoulos/ltm/blob/master/Examples/Scoring.R
  # library(ltm)
  # fs <- factor.scores(faPC)
  # plot(fs, include.items = TRUE, main = "EFA.poly")
  #
  # We can try to remove the 'sample' column and make it into the
  # row names, but we cannot have duplicate row.names, and if we
  # leave the three replicas, they would be considered as different
  # outcomes, not collated as a single one.
  #
  # In conclusion: in the case of Coturnix, EFA worked best with
  # count data, but in the case of Gallus, we cannot rely on raw
  # (or normalized, multi-sample) count data for either EFA or CA.
  #
  # A possible solution might be to take the average of all three counts
  # for each category, move the categories to row-names and try again
  # but in that case, we might as well try to use log2FC instead, which
  # (in theory) better summarizes our data, although it would reflect
  # differences between categories and not the categories themselves.
  #
  # To start with, we will try to approximate "population counts" from
  # our sample counts by taking the average of the normalized counts
  # using aggregate(), and then try to do EFA and CA with them. We
  # can do this because we are using normalized counts.
  mtc <- aggregate( tc[ , 2:dim(tc)[2] ], by=list(tc$sample), FUN=mean )
  rownames(mtc) <- mtc$Group.1
  mtc$Group.1 <- NULL 

  # and now we do EFA on the 'estimated population counts'
  #mtc.fac.5 <- fac(mtc,nfactors=5, rotate="varimax")
  #as.png(fa.plot(mtc.fac.5),
  #     'impgenes/mtc.fac.5.png')
 #as.png(fa.diagram(mtc.fac.5),
  #     'impgenes/mtc.fac.5.diagram.png')
  #mtc.fac.6 <- fac(mtc,nfactors=6, rotate="varimax")
  #as.png(fa.plot(mtc.fac.6),
  #     'impgenes/mtc.fac.6.png')
 #as.png(fa.diagram(mtc.fac.6),
  #	   'impgenes/mtc.fac.6.diagram.png')

WE_HAVE_LOTS_OF_TIME_AND_CPU <- FALSE
# NOTE we have 100 CPUs, 1TB memory, and still set this to FALSE!!!
if ( WE_HAVE_LOTS_OF_TIME_AND_CPU ) {
 # try to find out optimal number of hidden factors using fa.parallel()

    mtc.cor <- cor(mtc, use="pairwise.complete.obs")
    as.png(corrplot(mtc.cor, order='hclust', tl.col="black"), #the plot is wrong
        "efa/mtc_corplot.png", width = 800, height = 600)

    # this may take very long: it will run EFA many times, increasing the number
    # of factors to find an optimal number.
    # Maybe we should consider doing it with a subsample? Would that be representative?
    mtc.fa.paral <- fa.parallel(x=mtc.cor, fm="minres", fa="fa")
    print(mtc.fa.paral)
	} else {
	 # 800 is on the verge of tolerability for just finding out
		# a prediction for the number of factors/components
		# doing it with different size subsets will give us an 
		# idea of a potential number and its dependence on size
		# (if we assume genes are at random, size should have little
		# to no impact or at least tend to stabilize at some point)
    mtc100.fa.paral <- fa.parallel(x=mtc[ , 1:100], fm="minres", fa="fa")
    print(mtc100.fa.paral)
    mtc200.fa.paral <- fa.parallel(x=mtc[ , 1:200], fm="minres", fa="fa")
    print(mtc200.fa.paral)
    mtc400.fa.paral <- fa.parallel(x=mtc[ , 1:400], fm="minres", fa="fa")
    print(mtc400.fa.paral)
        mtc800.fa.paral <- fa.parallel(x=mtc[ , 1:800], fm="minres", fa="fa")
    print(mtc800.fa.paral)

	}
	# when run incrementally, there seems to be an agreement that
	# the optimal number of factors should be 2 according to fa.parallel()
	# using eigenvalue analysis with simulated data
	# and 6 by the elbow criterion as deduced from visual inspection of
	# the plots: 
	# we will try 6 and 5 (2 is too few).

  nFactors = 6
  st.seed(SEED)
  get_efa_factors(mtc, n.factors=6, annotation=bio.ann, 
      force.recalculation=F, 
          out.prefix='efa/mtc', 
          verbose=TRUE)
  # this works

  nFactors = 5
  st.seed(SEED)
  get_efa_factors(data=mtc, n.factors=5, annotation=bio.ann,
         force.recalculation=F, 
          out.prefix='efa/mtc', 
          verbose=TRUE)
  # this one raises an error:
  #  Error in e$vectors %*% diag(inv.sqrt.ev) : non-conformable arguments
  #
  # We can also try to do EFA with log2FC, this would tell us which
  # groups of genes are important to define each specific difference.
  #
  # Left for later.