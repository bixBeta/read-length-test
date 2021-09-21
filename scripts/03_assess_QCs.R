load("data/02b_scale_pca_umap_neighbors_clusters_dim50_multiRes.Rdata")
source("scripts/packages.R")

# going with res 0.6 for now 
Idents(object = cca.obj) <- "integrated_snn_res.0.6"
cca.obj$seurat_clusters = cca.obj$integrated_snn_res.0.6
cca.obj.meta = cca.obj@meta.data

# get number of cells per cluster per sample
n_cells.cca.obj <- FetchData(cca.obj, 
                             vars = c("ident", "orig.ident")) %>%
  dplyr::count(ident, orig.ident) %>%
  tidyr::spread(ident, n)

write.csv(n_cells.cca.obj, "results/n_cells_per_cluster_per_sample.csv")
