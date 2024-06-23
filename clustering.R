library(NbClust)
library(ggplot2)

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

#################################################
# H I E R A R C H I C A L   C L U S T E R I N G #
#################################################

hclust <- function(data, method = "complete", scale_data = TRUE) {
  if (escale_data) {
    data_scaled <- scale(data)
  } else {
    data_scaled <- data
  }
  
  dist_matrix <- dist(data_scaled)
  
  hclust_result <- hclust(dist_matrix, method)
  
  as.png(plot(hclust_result, main = "HClust Dendogram", xlab = "", sub = "", cex = 0.9), 
            "dendogram_hclust.png")
  
  as.png(heatmap(as.matrix(data_scaled), 
          Rowv = as.dendrogram(hclust_result), 
          symm = TRUE, 
          scale = "row", 
          margins = c(5, 10), 
          col = colorRampPalette(c("blue", "white", "red"))(100)), 
    "heatmap_hclust.png")
}

#####################################
# K M E A N S   C L U S T E R I N G #
#####################################

# find the optimal number of clusters for k-means
nclust <- function(data, dist = "euclidean", method = "complete", min_clusters = 2, max_clusters = 15) {
  data_scaled <- scale(data)
  
  nb <- NbClust(data_scaled, distance = dist, min.nc = min_clusters, max.nc = max_clusters, method = method)
  
  print(paste("El mejor número de clusters según el criterio mayoritario es:", nb$Best.nc[1]))
  print(nb$Best.nc)
  
  return(nb$Best.nc[1])
}

kmeansclust <- function(data, num_clusters = NULL, dist = "euclidean", method = "complete", min_clusters = 2, max_clusters = 15) {
  data_scaled <- scale(data)
  
  if (is.null(num_clusters)) {
    num_clusters <- nbclust(data, dist, method, min_clusters, max_clusters)
  }
  
  set.seed(290424)  
  kmeans_result <- kmeans(data_scaled, centers = num_clusters, nstart = 25)
  
  pca_result <- prcomp(data_scaled)
  data_pca <- as.data.frame(pca_result$x[, 1:2])
  data_pca$cluster <- as.factor(kmeans_result$cluster)
  
  as.png(ggplot(data_pca, aes(x = PC1, y = PC2, color = cluster)) +
    geom_point(size = 2) +
    labs(title = "K-Means Clustering Representado en 2D",
         x = "Componente Principal 1",
         y = "Componente Principal 2") +
    theme_minimal(), "kmeans_pca_clustering.png")
}





