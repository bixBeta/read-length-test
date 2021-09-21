# load integrated object 

source("scripts/packages.R")

cca.obj = readRDS("data/sobj.integrated.on.batchID.top.2000.features.Rds")

DefaultAssay(cca.obj) <- "integrated"

# Run the standard workflow for visualization and clustering
cca.obj <- ScaleData(cca.obj, verbose = T)
cca.obj <- RunPCA(cca.obj, npcs = 50, verbose = T)
ElbowPlot(cca.obj, ndims = 50)

cca.obj <- RunUMAP(cca.obj, reduction = "pca", dims = 1:50)
cca.obj <- FindNeighbors(cca.obj, reduction = "pca", dims = 1:50)
cca.obj <- FindClusters(cca.obj, resolution = c(0.4, 0.5, 0.6, 0.8))

cca.obj.meta = cca.obj@meta.data

save(cca.obj,cca.obj.meta, file = "data/02b_scale_pca_umap_neighbors_clusters_dim50_multiRes.Rdata")

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Idents(object = cca.obj) <- "integrated_snn_res.0.4"
umap.0.4 = DimPlot(cca.obj, reduction = "umap", label = TRUE, label.size = 6, raster = F)

Idents(object = cca.obj) <- "integrated_snn_res.0.5"
umap.0.5 = DimPlot(cca.obj, reduction = "umap", label = TRUE, label.size = 6, raster = F)

Idents(object = cca.obj) <- "integrated_snn_res.0.6"
umap.0.6 = DimPlot(cca.obj, reduction = "umap", label = TRUE, label.size = 6, raster = F)

Idents(object = cca.obj) <- "integrated_snn_res.0.8"
umap.0.8 = DimPlot(cca.obj, reduction = "umap", label = TRUE, label.size = 6, raster = F)

umap.0.4 + umap.0.5 + umap.0.6 + umap.0.8



