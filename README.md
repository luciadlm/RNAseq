### Lucía de Lamadrid - Biotechnology (UPM)

This code was used to achieve the objectives of the Biotechnology Bachelor's thesis:
"ANÁLISIS COMPUTACIONAL DE LA EXPRESIÓN GÉNICA MEDIANTE RNA-SEQ EN LA RESPUESTA AVIAR A LA INFECCIÓN POR EL VIRUS IBDV: EXPLORACIÓN DE MÉTODOS DE APRENDIZAJE NO SUPERVISADOS" 
The methods were divided into six main stages: 
1. Data acquisition, including sample preparation and sequencing.
2. Alignment and quantification of reads with a reference genome.
3. Differential expression analysis using DESeq2.
4. Gene annotation through genomic databases.
5. Gene Set Enrichment Analysis (GSEA).
6. Application of machine learning methods for the selection and analysis of relevant variables.
These can be useful for interpreting differential expression data in a simplified manner, removing noise, and highlighting the most relevant genes.

Among the repository files, you will find:
- DGEanalysis.R: alignment, feature counting and differential gene expression analysis using DESEq2.
- annotation.R: subsequent gene annotation with ensembldb and biomaRt.
- gsea.R: GSEA with Go terms and KEGG pathways using both fgsea and clusterProfiler.
- clustering.R: performs clustering with hierarchical methods, including the visualization of dendrograms and heatmaps, and with k-means, estimating the optimal number of clusters.
- efa.R: performs exploratory factor analysis, several tests are conducted with different combinations of rotation and factorization methods. Parallel analysis is conducted to estimate the number of factors. 
- ca.R: performs correspondence analysis, this is used for categorical variables. It allows the visualization of various two-dimensional graphs.
- featureselection.R: includes functions to implement various feature selection methods (correlation, relaimpo, stepwise, lasso, mars, rpart, boruta, recursive feature elimination, dalex, vita, regularized random forest, xgboost, genetic algorithm, simulated annealing)
      
![image](https://github.com/luciadlm/RNAseq/assets/172217433/769218b5-67cc-4f48-8ab2-272cef165e48)



