##########################################################################################################
##########################################################################################################
##########################################################################################################
##########################################################################################################

# These are some convenience functions to help with printing messages, loading libraries and plotting

#' my.name()
#'
#' my.name() obtains the name of the function that called it
#'
#' This function will return the name of the calling function.
#' The CALLING function !!!
#' THE CALLING FUNCTION !!!!!!
#' It is intended as a way to obtain a function's own name for error
#' reporting. 
#  We could do it directly, within a function using 
#  my.name <- deparse(sys.call())
#  but it would be difficult to understand. Isolating this in a function
#  call makes it more readable.
#'
#' 
#' @return	the name of the calling function (one up in the call stack)
#'
#' @usage	me <- my.name()
#'
#' @examples
#'	f <- function() { print(my.name()) }
#'	f()
#'	##[1] "print(my.name())"
#'	## my.name() is being passed to print() as an argument and is thus
#'	## called by 'print', hence the output
#'	##
#'	f <- function() { me <- my.name() ; print(me) }
#'	f()
#'	##[1] "f()"
#'	## my.name() is called by f() to obtain its own name and then the name
#'	## is handled to print().
#'
#' @author	(C) CNB-CSIC
#'
#' @license	EU-GPL
#'
#' @export
#
my.name <- function(verbose=F) { 
	my.call.number <- sys.nframe()
	if (my.call.number == 1) {
	    # we are at the script/command line top level
		top.file <- sys.frame(1)$ofile
	    if (is.null(top.file))
		    caller <- "R"
		else {
		    caller <- top.file
            if (verbose == FALSE) caller <- basename(caller)
		}
	} else {
    	caller <- deparse(sys.call(-1)) 
		if (! verbose) caller <- gsub('\\(.*', '', as.character(caller))
		#who.called.me <- deparse( sys.calls()[[ my.parent.number ]] )
		#return( who.called.me )
    }
    return(caller)
}

# alternative (and simpler) implementation
me <- function(..., verbose=F) {
    if (sys.nframe() == 1) callingFun <- "R"
    else callingFun = as.list(sys.call(-1))[[1]]
	
    calledFun = as.list(sys.call())[[1]]
    if (verbose) message(paste(callingFun, " is calling ", calledFun, sep=""))
	return(callingFun)
}


my.caller <- function(verbose=F) {
	my.call.number <- sys.nframe()
	if (my.call.number <= 2) {
	    # we are at the script/command line top level
		# or colled by a function whose parent is the
		# top level
		top.file <- sys.frame(1)$ofile
	    if (is.null(top.file))
		    caller <- "R"
		else {
		    caller <- top.file
            if (verbose == FALSE) caller <- basename(caller)
		}
	} else {
    	caller <- deparse(sys.call(-2)) 
		if (! verbose) caller <- gsub('\\(.*', '', as.character(caller))
		#who.called.me <- deparse( sys.calls()[[ my.parent.number ]] )
		#return( who.called.me )
    }
    return(caller)
}


#' cat.err()
#'
#' err() prints a generic error message using cat()
#'
#' @param	abort	(boolean) whether the program should be stopped
#'			(defaults to FALSE)
#' @param	...	The message to print and cat options (see cat())
#'
#' @return	whatever cat() returns
#'
#' @usage	cat.err('some message\n', sep='')
#'
#' @examples
#'		f <- function(x) { if (x < 10) {print(x)} else {cat.err(x, "is too big\n")}}
#'		f(13)
#'		##ERROR in  f(13) : 13 is too big
#'
#' @author	(C) CNB-CSIC
#'
#' @license	EU-GPL
#'
#' @export
#
cat.err <- function(..., abort=FALSE, verbose=T) { 
    # error messages should be always printed
	# so we will ignore the verbose argument
	verbose <- TRUE

    # caller <- mycaller()
    #caller <- deparse( sys.calls()[[sys.nframe()-1]] )
    caller <- my.caller(verbose=verbose)
    cat('ERROR in ', caller, ":", ...)
    if (abort) {
        quit(save="no", status=1, runLast=FALSE)
    }
}


#' cat.warn()
#'
#' cat.warn() prints a generic warning message using cat()
#'
#' @param	...	The message to print and cat options (see cat())
#'
#' @return	whatever cat() returns
#'
#' @usage	cat.warn('some message\n', sep='')
#'
#' @examples
#'		f <- function(x) { if (x < 10) {print(x)} else {cat.warn(x, "is too big\n")}}
#'		f(13)
#'		##WARNING in  f(13) : 13 is too big
#'
#' @author	(C) CNB-CSIC
#'
#' @license	EU-GPL
#'
#' @export
#
cat.warn <- function(..., verbose=F) { 
    #caller <- deparse( sys.calls()[[sys.nframe()-1]] )
    #me <- as.list(sys.call())[[1]]
    #parent <- as.list(sys.call(-1))[[1]]
	caller <- my.caller(verbose=verbose)
	cat('WARNING in ', caller, ":", ...)
}


#' cat.info()
#'
#' cat.info() prints a generic information message using cat()
#'
#' @param	...	The message to print and cat options (see cat())
#'
#' @return	whatever cat() returns
#'
#' @usage	cat.warn('some message\n', sep='')
#'
#' @examples
#'		f <- function(x) { if (x < 10) {print(x)} else {cat.info(x, "is too big\n")}}
#'		f(13)
#'		##INFO  f(13) : 13 is too big
#'
#' @author	(C) CNB-CSIC
#'
#' @license	EU-GPL
#'
#' @export
#

cat.info <- function(..., verbose=F) { 
    #caller <- deparse( sys.calls()[[sys.nframe()-1]] )
    #me <- as.list(sys.call())[[1]]
    #parent <- as.list(sys.call(-1))[[1]]
    caller <- my.caller(verbose=verbose)
    cat('INFO ', caller, ":", ...)
}


# to be used instead of library(): this function ensures
# that the package is installed if not present in the system
#' use.package
#'
#' checks if a package is installed and loads it, if it is not then it 
#' will install it and then load it
#' 
#' @param	p	a single package name to load
#'
#' @param	silent	whether we want messages produced during the load displayed
#' 
#' @return	the package is guaranteed to be installed and loaded
#'
#' @usage	use.package(p, silent=T)
#' 
#' @examples	use.package('stringr')
#'
#' @author	(C) CNB-CSIC
#'
#' @license	EU-GPL
#'
#' @export
#
use.package <- function(p, silent = F, verbose = T) {
    ### NOTE: we should likely try CRAN first, Bioconductor afterwards,
    ### install.packages offline, devtools::install and install_github
    ### for now we will kepp it simple
    if (!is.element(p, installed.packages()[,1])) {
        # we prefer BiocManager because it can install both Bioconductor
        # and CRAN packages
        #
        #install.packages(p, dep = TRUE)
        if (!is.element('BiocManager', installed.packages()[,1])) {
            install.packages('BiocManager', dependencies=TRUE)
        }
        BiocManager::install(p, dep = TRUE)

    }
    # this will fail for 'banneR' which is neither on CRAN nor on
    # Bioconductor, so we'll check for it as a special case for now
    if (p == 'banneR')
       library(banneR)		# if not available we will do without it
    if (verbose)
        require(p, character.only = TRUE)
	else
        suppressPackageStartupMessages(require(p, character.only = TRUE))
}



#' use.packages
#'
#' given a vector with one or more package names, it first verifies if all
#' are installed, if not, then missing packages are installed, and finally
#' all the packages are loaded
#' 
#' @param	p	vector of package names to load
#' 
#' @return	the packages are guaranteed to be installed and loaded
#'
#' @usage	use.packages(pkgs)
#' 
#' @examples	use.packages(c('stringr', 'tidyr'))
#'
#' @author	(C) CNB-CSIC
#'
#' @license	EU-GPL
#'
#' @export
#
use.packages <- function(pkgs){
    new.pkgs <- pkgs[!(pkgs %in% installed.packages()[, "Package"])]
    if (length(new.pkgs))
        # if there are uninstalled packages
        install.packages(new.pkgs, dependencies = TRUE)

    sapply(pkgs, require, character.only = TRUE)
}


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


#' get_annotation
#'
get_annotation <- function(ensembl.db=NULL, biomart.db=NULL, verbose=F) {
    if (verbose) cat.info('loading annotation\n')
	
    ensembl.annot <- NULL
    bm.annot.1 <- NULL
    if (! is.null(biomart.db) && file.exists(biomart.db)) {
	    if (verbose) cat.info('loading', biomart.db, '\n')
        # get file extension and convert to lower case to simplify checks
		ext <- tolower(file_ext(biomart.db))
        if (ext == 'tab' || ext == 'tsv' || ext == 'txt') {
            # from .../signif: '../../biomaRt.annotation.1st.txt'
            bm.annot.1 <- read.table(file=biomart.db, header=T)
        } else if (ext == 'rds') {
            bm.annot.1 <- readRDS(biomart.db)
        } else {
            cat.warn('extension of biomart database must be txt, tsv, tab or rds\n')
        }
    }
    if (! is.null(ensembl.db) && file.exists(ensembl.db)) {
	    if (verbose) cat.info("Loading ENSMBL annotation from", ensembl.db, '\n')
        # get file extension and convert to lower case to simplify checks
        ext <- tolower(file_ext(ensembl.db))
        if (ext == 'tab' || ext == 'tsv' || ext == 'txt') {
            # from .../signif: '../../../../net.EnsDb/net.ens.db.rds'
            ensembl.annot <- readRDS(file=ensembl.db)
        } else if (ext == 'rds') {
            ensembl.annot <- readRDS(ensembl.db)
        } else {
            cat.warn('extension of ensembl database must be txt, tsv, tab or rds\n')
        }
    }
    # we could also load an endsb.org.db package
    # and then (but this currently doesn't work for some reason)
    #ens.annot <- ensembldb::select(ens.db, column='GENEID', keytype='GENEID', keys=rownames(lfc), columns='SEQNAME')
    # ensembldb::supportedFilters(ens.db)
    # columns(edb)
    # listColumns(edb)
    # keytypes(edb)
    # gids <- keys(edb, keytype = "GENEID")
    # length(gids)
    #ens.ann <- AnnotationDbi::select(ens.db, 
    #              column='GENID', keytype= 'GENEID', keys=rownames(lfc), 
    #              columns= c('SEQNAME', 'SYMBOL', 'DESCRIPTION', # no longer available
    #                         'GENENAME', 'GENEID', 'ENTREZID', # empty
    #                         'TXID', 'TXBIOTYPE', # these make the call fail
    #                         'PROTEINID', 'UNIPROTID' # no longer available
    #                          )
    #              )

    return( list(ensembl=ensembl.annot, biomart=bm.annot.1) )
}

#' enrich_genes
#'
enrich_genes <- function(data, ann, data.genes='ensembl_gene_id', ann.genes=data.genes) {
    # add the annotation as additional columns to data, using as the
    # column data.genes from data as index into ann.genes to find the
    # corresponding annotation.
    #
    # use as
    #	annotated.data <- enrich_genes(data, "ensembl.gene.id", ensembl.ann, "GENEID")
    #	annotated.data <- enrich_genes(data, "ensembl.gene.id", biomart.ann, "ensembl_gene_id")
    #   annotated.data <- enrich
	# for ENTREZGENE data
	# annotated.data <- enrich_genes(data, "gene.id", biomart.ann, "ensembl_gene_id")
	
    enriched_genes <- cbind(data, ann[match(data[ , data.genes], ann[ ,ann.genes]), ]) 

    return(enriched_genes)
}


load_coturnix_l2fc_data <- function(file="all_log2FC_transpose.txt") {
	# return a list with the datasets that we will use in the subsequent analysis


    # Load log2FC of significantly altered genes:
    #df <- read.table("all_log2FC.txt",header=T)
    # unused    df <- read.table("all_log2FC.txt",header=T)
    fd <- read.table(file, header=T)

    # this ensures that rows get the proper PFU irrespective
    # of order
    #	1. remove prefix
    pfu <- sub("wt.*pfu", "", fd$gene)
    #	2. remove suffix
    pfu <- sub("_l2fc", "", pfu)


    # set name of first column
    #names(fd)[1] <- 'wt.vs.PFU'
    #fd[,1] <- c(0.1, 1.0, 10.0, 100.0)

    # substitute column1 (file names) by column "$pfu"
    fd <- cbind(pfu=as.numeric(pfu), fd[ , -1])

    # define working data sets
    # 1. add a row of all zeros for the wild.type
    #    (inf.pfu/wt.vs.pfu = 0, l2fc=0 for all genes)
    lfc  <- as.data.frame(rbind(rep(0,dim(fd)[2]), fd))
    # 2. remove columns with NAs
    lfc.noNA <- lfc[ , colSums(is.na(lfc))==0]			# remove NAs
    clfc <- lfc.noNA
    # we now have
    #	fd 	-> the log2FC transposed data for infected cells (4x8594)
    #	lfc -> the log2FC transposed data for wt and infected cells (5x8594)
    #   clfc -> the log2FC transposed data for wt and infected cells without NA-containing columns (5x1182)
    return( list( signif.l2fc.all.data=lfc, signif.l2fc.common.data=clfc ) )
}


load_coturnix_norm_count_data <- function(file="normalized_counts.tab", select=NULL) {
    # load normalized counts of all genes
    nc <- read.table(file=file, row.names=1)
    tnc <- t(nc)
    # prepend a column with infection levels
    # upon read columns are not preserved
    #tnc <- cbind(c(rep(0, 3), rep(0.1, 3), rep(1.0,3), rep(10.0,3), rep(100.0,3)), tnc)
    # so we will fix the labels
    pfu <- rownames(tnc)
    pfu <- sub('^X', '', pfu)
    pfu <- sub('pfu.*', '', pfu)
    pfu <- sub('wt.*', '0.0', pfu)
    pfu <- as.numeric(pfu)
    
    # amend nc and tnc to include a 'pfu' row/column
    nc <- rbind(pfu=pfu, nc)
    tnc <- cbind(pfu=pfu, tnc)
    tnc <- as.data.frame( tnc )		# just in case

    if ( ! is.null(select)) {
        # select only counts of sigfigicantly altered genes:
        sig.nc <- rbind(pfu=nc[1, ], nc[rownames(nc) %in% select, ])
        sig.tc <- cbind(pfu=tnc[ ,1], tnc[ , colnames(tnc) %in% select, ] )
    #tc <- sig.tc
    } else {
        sig.nc <- nc
        sig.tc <- tnc
    }

    # tc now contains the significant transposed normalized counts	(15x8954)
    # nc contains the significant normalized counts	(8954x15)

    return(list( all_norm_counts=nc, transp_norm_counts=tnc,
                 signif_norm_counts=sig.nc,
                 signif_transp_norm_counts=sig.tc))
}

get_gallus_significant_genes <- function(pattern='signif.*.tab', verbose=T) {
    # one way is to load all saved significant comparisons and find
	# the union of the gene names
	files <- list.files(pattern='signif_sorted_sample.*_annotated.tab')
	genes <- c()
	data <- list()
	for (i in files) {
	    df <- read.table(i, header=T)
		# gene annotation in gallus was defined for all genes only in 'gene.name'
		genes <- union(genes, df$gene.name)
		data[[i]] <- df
		# safety check
		if (verbose) cat(i,sum(!is.na(df$gene)), 'OK', sum(is.na(df$gene)), 'KO', length(genes), 'tot\n')
	}
	data$sig.genes <- genes
	return(data)
}

# for gallus we will only use normalized count data
# but this function will load count data for all genes, irrespective
# of whether they are significant or not, that is why we need a select
# argument
load_gallus_norm_count_data <- function(file="normalized_counts.tab", select=NULL) { 
    # for gallus we have also annotated normalized counts that we can use
	nc <- read.table(file, header=T, row.names=1, sep='\t')
	
	# select significant genes
    if ( ! is.null(select)) {
        # select only counts of sigfigicantly altered genes:
        sig.nc <- nc[rownames(nc) %in% select, ]
    #tc <- sig.tc
    } else {
        sig.nc <- nc
    }
	return(list(total.counts=nc, significant.counts=sig.nc))
}

#
# -- Correlation Method
#
# Pearson Correlation
imp_cor <- function (formula, data, verbose=T) {
    # get the left hand side of the formula (the dependent variable) as a sub-formula
    fy <- f_lhs(formula)
    # convert to column names
    ny <- all.vars(fy)
    # get the right hand side (the predictor variables) as a sub-formula
    fx <- f_rhs(formula)
    # convert to column names
    nx <- all.vars(fx)
    if ( nx == '.' ) {
       nx <- colnames(data)[ colnames(data) != ny ]
    } 
    # get x and y
    y <- data[ , ny]
    x <- data[ , nx]

    if(verbose) cat.info("calculating correlation between variables\n")

    data_cor <- cor(x, y)	

    #if (verbose) cat.info("selecting variables with abs(cor) > 0.8\n")
    #highlyCorrelated <- rownames(data_cor[abs(data_cor[, 1]) > 0.8, , drop = FALSE])

    return(data_cor)
}



# -- relimp
#library(relaimpo)
#
# Apply a linear regression model and get the relative importance of 
# each variable using R² partitioned by averaging over orders. Observations
# with missing values are automatically excluded.
#
imp_relimp <- function(formula, data, tune=F, verbose=T) {
    if (tune == F) {
        if (verbose) cat.info("building linear model\n")
        # directly build a linear model
        regressor <- lm(formula, data=data) # fit a lm() model

    } else { 
        if (verbose) cat.info("building linear model with CARET\n")
         # tune with caret
        ### XXX WARNING!!!
        ### MAGIC NUMBERS!!!
        ### WE need to think if this is general or not
        control <- trainControl(method="repeatedcv", number=10, repeats=3)
        tuneGrid <- expand.grid(k=1:3, size=seq(5,20,by=1))	
        regressor <- train(formula, data=data, 
                           method="lm",  # Use "lm" for linear regression
                           preProcess="scale", 
                           trControl=control#,
                           #tuneGrid=tuneGrid
                           )
        # if verbose, plot the regressor
        # use try() to ignore errors at this step.
        if (verbose) try(plot(regressor))		
        if (verbose) cat.info("calculating importance\n")
        if (FALSE) {
            # this is variable importance using CARET
            #result <- calc.relimp(regressor$finalModel, type = "lmg", rela = TRUE) # check how to retrieve the lm from caret result
            result <- caret::varImp(regressor, conditional=TRUE) 	# conditional=True, adjusts for correlations 
    													 	        # between predictors
            return(result)
       } else {
            #to use relimp, we need to access directly the  model
            regressor <- regressor$finalModel
        }
    }
    if (verbose) try(summary(regressor))
    if (verbose) cat.info("calculating importance\n")
    relImportance <- calc.relimp(regressor, type = "lmg", rela = TRUE) # calculate relative importance scaled to 100

    return(sort(relImportance$lmg, decreasing=TRUE)) # relative importance
}

high_cor_pairs <- which(abs(cor_matrix) > 0.9, arr.ind = TRUE)

remove_high_cor <- function(data, threshold = 0.9) {
  cor_matrix <- cor(data)
  high_cor_pairs <- which(abs(cor_matrix) > threshold, arr.ind = TRUE)
  
  # Variables a eliminar (solo una de cada par)
  vars_to_remove <- unique(high_cor_pairs[high_cor_pairs[, 1] > high_cor_pairs[, 2], 1])
  
  return(data[, -vars_to_remove])
}

data_reduced <- remove_high_cor(data[, -1])  # Excluir la variable dependiente 
data_reduced <- cbind(data$y, data_reduced)
names(data_reduced)[1] <- "y"



# -- Step-wise Regression Method
# 
# apply a linear model and reduce automatically the number of variable using
# a stepwise approach: variables are sequentially entered or removed from
# the model, and their added informative contribution is used to decide
# whether to keep them or not. We use bidirectional elimination testing at
# each step for variables to include or exclude.
#
# Prone to overfitting (i.e. results may be too specific).
#
# If you have large number of predictors, consider splitting the Data in 
# chunks of 10 predictors with each chunk holding the responseVar.
#
imp_stepwise <- function(formula, data, direction="both", tune=F, verbose=T) {

    if (verbose) cat.info("building linear model\n")
    # for the base model we need a formula that looks as "dep.var ~ 1"
    # take LHS of the formula (the dependent variable) and convert (deparse) to string
    y=deparse(f_lhs(formula))
    base.formula <- as.formula(paste(y, "~ 1"))		# create the needed base formula
    base.mod <- lm(base.formula , data=data) 
    if (tune == F) 
        full.mod <- lm(formula, data=data) # full model with all predictors
    else  {
        control <- trainControl(method="repeatedcv", number=10, repeats=3)
        full.mod <- train(formula, data=tc, method='lm', trControl=control)$finalModel
    }

    if (verbose) cat.info("calculating importance\n")
    stepMod <- step(base.mod, 
		            scope = list(lower = base.mod, upper = full.mod), 
                    direction = direction, 
                    trace = 1, 
                    steps = 1000) 

    # according to str() we want to save $Coefficients
    return( stepMod$coefficients )
}
# The output might include levels within categorical variables, since
# 'stepwise' is a linear regression based technique.
# 
# If you have a large number of predictor variables, the above code may need to
# be placed in a loop that will run stepwise on sequential chunks of
# predictors. The shortlisted variables can be accumulated for further analysis
# towards the end of each iteration. This can be very effective method, if you
# want to
# 
# * Be highly selective about discarding valuable predictor variables. 
# 
# * Build multiple models on the response variable.




#
# -- Recursive Feature Elimination RFE Method
#
# Search for a subset of variables using a backward selection (removing
# variables) to find the optimal combination removing the ones with less
# importance according to some metric.
#
#	We tune RFE using caret, controlling the output sizes from 4 to 2^8
# an letting the algorithm select the optimal number. We use a custom
# random forests model.
#
#	Alternatives are linear models (lmFuncs), naive Bayes (nbFuncs), 
# bagged trees (treebagFuncs) and caret-available functions (caretFuncs)
#
#library(caret)

imp_rfe <- function (formula, data, verbose=T) {
    # WE HAVE A FORMULA    y ~ a + b + c ...  or y ~ .
    # get the left hand side of the formula (the dependent variable) as a sub-formula
    fy <- f_lhs(formula)
    # convert to column names
    ny <- all.vars(fy)
    # get the right hand side (the predictor variables) as a sub-formula
    fx <- f_rhs(formula)
    # convert to column names
    nx <- all.vars(fx)
    if ( nx == '.' ) {
       nx <- colnames(data)[ colnames(data) != ny ]
    } 
    # get x and y
    y <- data[ , ny]
    x <- data[ , nx]

    if (verbose) cat.info("building RFE model\n")
    control <- rfeControl(functions=rfFuncs, method="cv", number=10)
    # caution MAGIC NUMBERS HERE
    results <- rfe(x, y, sizes=2^(2:8), rfeControl=control)

    if (verbose) cat.info("calculating importance\n")
    predictors <- predictors(results)
    if (verbose) {
        # summarize the results
        # list the chosen features
        # plot the results
        plot(results, type=c("g", "o"))
    }


    return(predictors)
}



# -- MARS (earth package)
#
# apply a non-linear model using spline curves estimated automatically to
# model automatically non-linear dependencies and interactions between
# variables.
#
# The earth package implements variable importance based on Generalized cross
# validation (GCV), number of subset models the variable occurs (nsubsets) and
# residual sum of squares (RSS).
library(earth)

imp_mars <- function (formula, data, verbose=T) {
    if (verbose) cat.info("building MARS model\n")
    regressor <- earth(formula, data=data) 

    if (verbose) cat.info("calculating importance\n")
    ev <- evimp (regressor) # estimate variable importance

    return(ev)
}



# 
# -- Learning Vector Quantization (LVQ) Method for CATEGORICAL variables
#
# can be considered a special kind of neural network, precursor of SOM and
# related o kNN, that applies a winner takes it all learning approach. During
# learning it tends to select the most meaningful variables.
#
# LVQ uses dropout for training. 
#
#library(caret)

imp_lvq <- function(formula, data, verbose=T) {
    if (verbose) cat.info("building LVQ model\n")
    control <- trainControl(method="repeatedcv", number=10, repeats=3)
    # to avoid error 'subscript out of bounds'
    tuneGrid = expand.grid(k=1:3,size=seq(5,20,by=1))
    # train the model
    regressor <- train(formula, data=data, 
                      method="lvq", 			# error: wrong model for regression
                      preProcess="scale", 
                      trControl=control,
                      tuneGrid=tuneGrid)
    if (verbose) cat.info("calculating importance\n")
    # estimate variable importance
    importance <- caret::varImp(regressor, scale=FALSE)$importance

    return(importance)
}



#
# -- rpart (Recursive partitioning and regression trees)
#
# is an extension of decision trees from tables, and can be used for both
# classification, regression and survival like the popular and successful
# CART (Classification and Regression Trees, which builds a binary tree
# of binary or numeric data using the GINI index and cost-complexity pruning) 
# or C5.0 (which gives a binary or multibranch classification tree using
# information gain or entropy pruned by the binomial confidence limit)
#
#library(rpart)
#library(caret)

imp_rpart <- function(formula, data, tune=F, verbose=T) {
    if (tune == F) {
        if (verbose) cat.info("building rpart model\n")
        # this rcontrol was taken from iv_num in package 'riv'
        #	CAUTION: MAGIC NUMBERS
        minbucket <- nrow(data)/100
        if (minbucket < 1) minbucket <- 1
        rcontrol <- rpart.control(cp=0.0001, minbucket=minbucket)
        model <- rpart(data=tc, formula=formula, control=rcontrol)
	
	if (verbose) {
		png('figures/rpart_tree.png')
		plot(model)
		text(model, pretty=1)
		dev.off()
	}
        # we need now a way to calculate importance. woe/riv might
        # be a good one if it did work, so we will use the ones
        # rpart itself calculates
        # model$variable.importance is a named vector of weights with variables as names
        return(model$variable.importance)
    } else {	# (tune == TRUE)
        if (verbose) cat.info("building rpart model with CART\n")
        set.seed(SEED)
        # CAUTION: MAGIC NUMBERS
        trainControl= trainControl(method="repeatedcv", repeats=5)
        rPartMod <- train(formula, 
                          data=data, 
                          method="rpart",
                          importance=T,
                          trControl=trainControl)

        if (verbose) cat.info('calculating importance\n')
        rpartImp <- caret::varImp(rPartMod, scale=F)

        return(rpartImp)
    }
}





# -- random forest
#
# a random forest uses an ensemble learning method to combines the output 
# of multiple decision trees to reach a single result and works for both, 
# classification and regression (using a najority rule for classification
# or the average for regression).
#
# Random forest can be very effective to find a set of predictors that best
# explains the variance in the response variable.
# library(caret)
# library(randomForest)
# This function builds a Random Forest to find the most relevant variables
#
# A Random Forest is like 'rpart' but with many decision trees instead of only one
#
imp_random_forest <- function(formula=NULL, data, dep.var='', verbose=T, method='rangerRF', ...) {
    # create a RF model
    if (verbose) cat.info("building RF model\n")
    if (method == 'rangerRF') {
        # Ranger is a (not so) fast implementation of random forest particularly suited for high
        # dimensional data... 
        # ..and seems to work for <= 5 samples
        if (verbose) cat.info("calculating importance\n")	# we can do both together with ranger
        regressor <- ranger::ranger(formula, data=data, 
                                    importance='permutation', 
                                    local.importance=T,
                                    scale.permutation.importance=T)
        return(regressor$variable.importance)			# importance for each independent variable
        #return(regressor$variable.importance.local)	# importance for each variable and for each sample
 	} else if (method == 'rangerRF.dep.var' ) {
	    # formula is not a formula, but the name of the dependent variable
		# we need to use this when dimensionality is too high and the
		# formula interface fails
		regressor <- ranger::ranger(dependent.variable.name=dep.var, 
		                            data=tc,
									importance='permutation',
									local.importance=T,
									scale.permutation.importance=T)
        return(regressor$variable.importance)			# importance for each independent variable
   }
    # else (if method == 'RF')
    regressor <- randomForest(formula , data=data, importance=TRUE, ...)  # fit the random forest with 
                                                                          # provided parameters (if any)
    if (verbose) cat.info("calculating importance\n")
    # get variable importance, based on mean decrease in accuracy using
    # varImp in the caret package (we indicate it explicitly to avoid name clashes due to masking)
    if ( TRUE ) {
        result <- caret::varImp(regressor, conditional=TRUE) 	# conditional=True, adjusts for correlations 
    													 	    # between predictors
        return(result)
    } else {
        # unused for now, we will activate it eventually and add a function-call argument
        # to allow choosing
        varImp::varimpAUC(regressor)    # more robust towards class imbalance, returns a list
        return(result$importance)
    }
}



# -- Boruta Method
# 
# The 'Boruta' method can be used to decide if a variable is
# important or not based on RF. It is an all relevant method, not a 
# minimal optimal one (may include redundant variables).
#
# Random Ferns, originally proposed for vision can be considered a
# constrained decision tree ensemble, as reliable as randomized trees
# but faster, and compete with SVM in some tasks.
#
# We allow selecting each of the two methods as the base for Boruta.
#

# install.packages('Boruta')
#library(Boruta)

imp_boruta <- function (formula, data, getImp='getImpRfZ', out.png=NULL, verbose=2) {
    if (verbose) cat.info("building boruta model\n")

    boruta_output <- Boruta(formula, data, doTrace=verbose) # perform Boruta search

    # plot results
    if (verbose) 
	    as.png(plot(boruta_output,
		       cex.axis=.7, las=2, xlab="", main="Variable Importance"),
		       file=out.png)
    #as.png(
    #    plot(boruta_output, cex.axis=.7, las=2, xlab="", main="Variable Importance"),
    #    file="impgenes/imp_boruta_tc.png") 


    if (verbose) cat.info("calculating importance\n")
    # get confirmed and tentative factors
    # boruta_signif <- names(boruta_output$finalDecision)[boruta_output$finalDecision 
    # 									      %in% c("Confirmed", "Tentative")]
    # alternate (and more elegant) way
    boruta_signif <- getSelectedAttributes(boruta_output, withTentative=T)
    if (verbose) attStats(boruta_output)

    return(boruta_signif)
}





#
# -- Regularized Random Forest (RRF)
#
# regularized regression estimates a penalization function to produce
# more parsimonious models that have less variables with predictive power 
# similar to full models and robust to correlated variables, whereas RF 
# rely only on tree randomization. RRF can lead to substantially more
# accurate models.
#
imp_rrf <- function ( formula, data, verbose=T ) {
    if (verbose) cat.info("building RRF model\n")

    set.seed(SEED)

    if (verbose) cat.info("building rrf model\n")
    rrfMod <- train(formula, 
                data = data, 
                method = "RRF",
                importance = T) #this argument is required for varImp

    if (verbose) cat.info('calculating importance\n')
    rrfImp <- caret::varImp(rrfMod, scale=F)

    if (verbose) as.png(
				    plot(rrfImp, top = 20, main='Variable Importance'),
                    file="impgenes/imp_rrf_tc.png", verbose)

    if (verbose) plot(rrfImp$importance)

    return(rrfImp$importance)
}






#
#
# -- DALEX (Model Agnostic Language for Exploration and Explanation)
#
# DALEX is a model-agnostic framework to explore any model and 
# explain its behavior calculating the contribution of each 
# variable as changes in the expected model response given other 
# variables. 
#
# DALEX does not seem to work with ranger and categorical variables.
# RandomForest is limited by gene names, and takes much, much longer
# than ranger, but should work.
#
#library(randomForest)
#library(DALEX)
#
imp_dalex <- function(formula=NULL, data, dep.var='', method='rf', verbose=T) {
    if (verbose) cat.info("building DALEX model with", method, "\n")
    if (method == "rf") {
        regressor <- randomForest(formula , data=data, importance=TRUE) 
        # fit the random forest with default parameters
	} else if (method == "rangerRF") {
	    # DO NOT USE
        # Ranger is a (not so) fast implementation of random forest particularly suited for high
        # dimensional data... 
        # NOTE THAT DALEX DOES NOT WORK WITH RANGER
		regressor <- ranger::ranger(formula, data=data, 
                                    importance='permutation', 
                                               local.importance=T,
                                               scale.permutation.importance=T)
	} else if (method == 'rangerRF.dep.var' ) {
	    # DO NOT USE
	    # formula is not a formula, but the name of the dependent variable
		# we need to use this when dimensionality is too high and the
		# formula interface fails
        # NOTE THAT DALEX DOES NOT WORK WITH RANGER
		regressor <- ranger::ranger(dependent.variable.name=dep.var, 
		                            data=data,
									importance='permutation',
									local.importance=T,
									scale.permutation.importance=T)
    } else if (method == "bGLM") {
        regressor <- glm(formula , data=data, family='binomial')
    } else if (method == "mGLM") {
        cat.warn('multinomial models are unsupported\n')
        regressor <- glm(formula , data=data, family='multinomial')
    } else if (method == "gGLM") {
        regressor <- glm(formula , data=data, family='gaussian')
    } else {
        cat.err("unknown method", method, '\n')
    }

    if (verbose) cat.info("calculating importance\n")
    # Variable importance with DALEX
    # we need 'y', and to get it we need to parse the formula:
    #    all.vars returns all the variables that have been used in the formula
    #    so, the first element is the 'y' component (y ~ a ? b ? ...)
    if (dep.var == '') {
    	y <- all.vars(formula)[1]
    } else {
	    y = dep.var
	}
	# and now we use as y the column named y
    explained <- DALEX::explain(regressor, data=data, y=data[ , y], label=method)
    #x <- data[ , ! names(data) %in% c(y) ]
	#explained <- DALEX::explain(regressor, data=x, y=data[ , y], label=method)
    #
    # residuals <- model_performance(explained)
    # plot(residuals)
    # vip <- DALEX::variable_importance(explained) 
    # pdp <- DALEX::variable_response(explained, variable=*select one var*, type='pdp')

    DALEX=T
    if (DALEX) {
        if (verbose) cat.info('extracting variable importance (slow)\n')
        vip <- DALEX::variable_importance(explained)
		return( vip )
    } else {
        if (verbose) cat.info('extracting feature importance (slow)\n')
        # Get the variable importances
        #varimps <- variable_dropout(explained, type='raw')
        #print(varimps)
        #plot(varimps)
        vip <- DALEX::feature_importance(explained, type='variable_importance') # very slow
    }
	# vip contains the results of all the permutations
	# plot 20 most important
    if (verbose) {
		cat.info('plotting importance')
        nv<-dim(vip)[1] / length(levels(as.factor(vip$permutation)))
        as.png(plot(vip[c(1,(nv-20):nv),]), 'impgenes/top20_dalex.png')
    }
    # get mean loss and sort to get variable importance
    vip.mean.loss <- aggregate(vip$dropout_loss, list(vip$variable), mean)
    imp <- (vip.mean.loss[ order(vip.mean.loss$x, decreasing=F),])
    colnames(imp) <- c('variable', 'mean_dropout_loss')
    #print(imp)
    #head(imp)
    #head(sort(imp))
	# re-sort imp
	imp <- imp[order(imp$mean_dropout_loss, decreasing=T), ]
    return ( imp )
}




#
# -- VITA
#
#		PIPM implements the test approach of Altmann et 
# al. (2010) for the permutation variable importance measure 
# in a random forest for classification and regression by randomly
# shuffling the values of a single variable and observing the
# resulting degradation of the model score.
#
#library(vita)
#
imp_vita <- function (formula=NULL, data, dep.var="", method='RF', verbose=T) {
    if (verbose) cat.info("building RF model\n")

    if (method == 'rangerRF') {
        # Ranger is a (not so) fast implementation of random forest particularly suited for high
        # dimensional data... 
        # ..and seems to work for <= 5 samples
        regressor <- ranger::ranger(formula, data=data, 
                                    importance='permutation', 
                                               local.importance=T,
                                               scale.permutation.importance=T)
		# NOTE that PIPM implements the test approach of Altmann et 
		# al. (2010) for the permutation variable importance measure 
		# ‘VarImp’ in a random forest for classification and regression.
		# When we use ranger, we are already asking for the permutation
		# approach, so this is the same as using 'ramger' in imp_RF
		return(regressor$variable.importance)
    } else if (method == 'rangerRF.dep.var' ) {
	    # formula is not a formula, but the name of the dependent variable
		# we need to use this when dimensionality is too high and the
		# formula interface fails
        regressor <- ranger::ranger(dependent.variable.name=dep.var, 
		                            data=tc,
									importance='permutation',
									local.importance=T,
									scale.permutation.importance=T)
		# NOTE that PIPM implements the test approach of Altmann et 
		# al. (2010) for the permutation variable importance measure 
		# ‘VarImp’ in a random forest for classification and regression.
		# When we use ranger, we are already asking for the permutation
		# approach, so this is the same as using 'ranger' in imp_RF
		return(regressor$variable.importance)
	} else if (method == 'RF') {
        regressor <- randomForest(formula , data=data, importance=TRUE)  # fit the random forest with 
    }

    if (verbose) cat.info("calculating importance\n")
    y <- all.vars(formula)[1]
    pimp.varImp.reg <- PIMP(data, data[ , y], regressor,S =10, parallel=TRUE)
    return(pimp.varImp.reg$VarImp[order(pimp.varImp.reg$VarImp, decreasing=T),])
}




# -- LASSO regression (least absolute shrinkage and selection operator)
#
# selects variables and regularizes models to improve their precision and
# interpretability. Originally formulated for linear regression.
#
# Least Absolute Shrinkage and Selection Operator (LASSO) penalizes 
# with L1-norm.
#     It imposes a cost to having large weights (value of coefficients). 
#	  It is called L1 regularization, because the cost added, is
# proportional to the absolute value of weight coefficients.
# As a result, it reduces the coefficients of unwanted variables to zero leaving
# only the important ones.
#
#library(lasso)
#
# this works for both continuous and categorical variables
imp_lasso <- function (formula, data, verbose=T) {
    # get the left hand side of the formula (the dependent variable) as a sub-formula
    fy <- f_lhs(formula)
    # convert to column names
    ny <- all.vars(fy)
    # get the right hand side (the predictor variables) as a sub-formula
    fx <- f_rhs(formula)
    # convert to column names
    nx <- all.vars(fx)
    if ( nx == '.' ) {
       nx <- colnames(data)[ colnames(data) != ny ]
    } 
    # get x and y
    y <- data[ , ny]
    x <- data[ , nx]


    # make y binomial
    #y <- as.double(as.matrix(ifelse(y == 0.0, 0, 1))) 
    #y <- as.double(as.matrix(y))

    if(verbose) cat.info("building lasso model\n")
    set.seed(SEED)
    if (is.factor(y)) {
        if (length(levels(y)) == 2) {
            family <-'binomial'
        } else {
            family='multinomial'
        }
        if (verbose) cat.info("y is a factor, using a", family, "lasso model\n")
    } else {
        if (verbose) cat.info("assuming gaussian y\n")
        family='gaussian'
    }        
    cv.lasso <- cv.lasso(as.matrix(x), y, family=family, 
    		 		      alpha=1, 			# 0 -> ridge, 1-> LASSO
                          parallel=TRUE, 
                          standardize=TRUE #, 
                          #type.measure='auc'
                         )

    if (verbose) plot(cv.lasso)
    # plot with two vertical lines: first is the lambda with lowers mean squared error
    # and second the number of variables with highest deviance within 1 SD

    if (verbose) cat.info("calculating importance\n")
    # the best lambda value is stored in $lambda,min
    if (family != 'multinomial') {
        coef<- round(as.matrix(coef(cv.lasso, s=cv.lasso$lambda.min)), 2)
        coef <- coef[coef[ ,1] !=0, ]
    } else {
        coef <- list()
        for (i in levels(y)) {   # we get a set of coefficients for each level
            level.coef <- round(as.matrix(coef(cv.lasso, s=cv.lasso$lambda.min)[[ i ]]), 2)
            coef[[i]] <- level.coef[ level.coef[ ,1] != 0, ]
        }
    }

    # another way
    #coefList <- coef(cv.lasso, s='lambda.1se')
    #coefList <- data.frame(coefList@Dimnames[[1]][coefList@i+1],coefList@x)
    #names(coefList) <- c('var','val')


    # return coefficients with lambda value !=0
    return( coef )

}

library(glmnet)

# Preparar los datos para glmnet
x <- model.matrix(y ~ x1 + x2 + x3, data)[,-1]  # Matriz de dise?o sin el intercepto
y <- data$y  # Variable dependiente

# Ajustar el modelo LASSO
lasso_model <- cv.glmnet(x, y, alpha = 1)
plot(lasso_model)

# Obtener los coeficientes del mejor modelo LASSO
lasso_coefficients <- coef(lasso_model, s = "lambda.min")
print(lasso_coefficients)

library(Matrix)

# Extraer los nombres de las variables seleccionadas
selected_indices <- lasso_coefficients@i + 1  # +1 porque los ?ndices en R son 1-based
selected_variables <- lasso_coefficients@Dimnames[[1]][selected_indices]
selected_coefficients <- lasso_coefficients@x

# Crear un data frame con las variables seleccionadas y sus coeficientes
selected_model <- data.frame(Variable = selected_variables, Coefficient = selected_coefficients)

# Mostrar las variables seleccionadas y sus coeficientes
print(selected_model)

# -- xgboost (Extreme Gradient Boosting)
#
# Gradient boosting tries to find the function that best approximates
# the data by taking lots of simple functions and adding them together
# and training on the gradient of the error with respect to loss 
# predictions. In gradient boosting trees we consider as the simple model
# a decision tree. XGBoost looks at the distribution of variables across
# all data points in a leaf to reduce the search space, adds some
# regularizations and speed up the calculations.
#
# Gradient boosting trees can be more accurate than RF and RRF as they 
# can correct each other's errors.
#
#library(caret)
#library(xgboost)
#
imp_xgboost <- function(formula, data, verbose=T) {
    if (verbose) cat.info("building xgboost model\n")

    regressor <- train(formula, data, 
                       method = "xgbTree", 
                       trControl = trainControl("cv", number = 10))
                       # scale=TRUE)
    if (verbose) cat.info("calculating importance\n")
    result <- caret::varImp(regressor, conditional=T)

    # according to str() we want to save $importance from caret::varImp() results
    return(result$importance)
}




# -- Genetic Algorithm
#
# repeatedly uses a huge number of random selections of the variables and 
# calculates the RMSECV (RMS error of cross-validation) and select the best
# candidates discrading lower half combination on each iteration and "breeding"
# and mutating the best half until the ending criteria is met..
#
#library(caret)

imp_ga <- function (formula, data, verbose=T) {
    # get the left hand side of the formula (the dependent variable) as a sub-formula
    fy <- f_lhs(formula)
    # convert to column names
    ny <- all.vars(fy)
    # get the right hand side (the predictor variables) as a sub-formula
    fx <- f_rhs(formula)
    # convert to column names
    nx <- all.vars(fx)
    if ( nx == '.' ) {
       nx <- colnames(data)[ colnames(data) != ny ]
    } 
    # get x and y
    y <- data[ , ny]
    x <- data[ , nx]

    if (verbose) cat.info("building GA model\n")
    # Define control function
    ga_ctrl <- gafsControl(functions = rfGA,  method = "cv", repeats = 3)
    # another option is `caretGA`
    # Genetic Algorithm feature selection
    ga_obj <- gafs(x,
 	              y,
 	              iters = 100,  # normally much higher (100+)        
  	              gafsControl = ga_ctrl)

	if (verbose) {
    	save(ga_obj, file="impgenes/ga_obj.RData")
		# recover with load('impgenes/ga_obj.RData')
    }
	return(data.frame(variable=ga_obj$optVariables, 
					  importance=ga_obj$fit$importance#,
					  #importanceSD=ga_obj$fit$importanceSD
					  ))
}




# -- Simulated Annealing
#
# is a global search method that makes small perturbations (removing a
# small subset of model variable members) to an initial candidate solution 
# trying to improve it with an acceptance probability.
# A subotimal solution should eventually improve over interations.
#
#library(caret)
#
imp_sa <- function (formula, data, verbose=T) {
    # get the left hand side of the formula (the dependent variable) as a sub-formula
    fy <- f_lhs(formula)
    # convert to column names
    ny <- all.vars(fy)
    # get the right hand side (the predictor variables) as a sub-formula
    fx <- f_rhs(formula)
    # convert to column names
    nx <- all.vars(fx)
    if ( nx == '.' ) {
       nx <- colnames(data)[ colnames(data) != ny ]
    } 
    # get x and y
    y <- data[ , ny]
    x <- data[ , nx]

    if (verbose) cat.info("building SA model\n")
    # Define control function  
    sa_ctrl <- safsControl(functions = rfSA,
                           method = "repeatedcv",            
                           repeats = 3,            
                           improve = 5)

    set.seed(SEED)
    # Simulated Annealing feature selection
    sa_obj <- safs(x, 
    		       y,
                   safsControl = sa_ctrl)

    if (verbose) {
	    save(sa_obj, file="impgenes/sa_obj.RData")
		# recover with load('impgenes/sa_obj.RData')
	}
    return(data.frame(variable=sa_obj$optVariables, 
					  importance=sa_obj$fit$importance,
					  importanceSD=sa_obj$fit$importanceSD))
}