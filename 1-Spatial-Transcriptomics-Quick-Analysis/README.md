<h1>Quick Spatial Transcriptomics Analysis</h1>
<p style="text-align: left;">
Hi everyone, in this video, we are going to perform a quick clustering for spatial transcriptomics data created using the 10X visium platform.
</p>
<h2>Load libraries</h2>
<pre>
  <code>
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
  </code>
</pre>
<h2>Perform the analysis</h2>
<pre>
  <code>
# download the `1-Spatial-Transcriptomics-Quick-Analysis` directory from the github repository 
#and copy/paste the absolute path in the following line. Hint: to get the absolute path, 
#you can navigate to the filtered_feature_bc_matrix folder that you have downloaded using
# the "Files" panel in Rstudio then you can press "More" and "copy file path to clipboard"
data_dir <- ('~/Desktop/BEC/1-Spatial-Transcriptomics-Quick-Analysis/Data')
slide1 <-Load10X_Spatial(data.dir=data_dir, slice='slide1')
# investigate the quality control parameters (RNA counts)
plot1 <- VlnPlot(slide1, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(slide1, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)
# perform normalization
slide1<-SCTransform(slide1, assay = "Spatial", verbose = FALSE)
saveRDS(slide1, "slide1.rds")
  </code>
</pre>
<h2>Perform clustering</h2>
<pre>
  <code>
slide1 <- readRDS("slide1.rds")
slide1 <- RunPCA(slide1, assay = "SCT", verbose = FALSE)
slide1 <- FindNeighbors(slide1, reduction = "pca", dims = 1:30)
slide1 <- FindClusters(slide1, verbose = FALSE)
slide1 <- RunUMAP(slide1, reduction = "pca", dims = 1:30)
p1 <- DimPlot(slide1, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(slide1, label = TRUE, label.size = 3)
p1 + p2
  </code>
</pre>
<p><img src="/plot1.png" alt="description" width="800" height="400" /></p>
<h2>Plot the expression of certain features</h2>
<pre>
  <code>
# let's plot Hpca which is a marker for the hippocampus
SpatialFeaturePlot(slide1, features = "Hpca", alpha = c(0.1, 1))
  </code>
</pre>
<p><img src="/plot2.png" alt="description" width="400" height="400" /></p>
<h2>Highlight cells that belong to certain clusters</h2>
<pre>
  <code>
SpatialDimPlot(slide1, cells.highlight = CellsByIdentities(object = slide1, idents = c(2, 1, 4)), facet.highlight = TRUE, ncol = 3)
  </code>
</pre>
<p><img src="/plot3.png" alt="description" width="800" height="400" /></p>
