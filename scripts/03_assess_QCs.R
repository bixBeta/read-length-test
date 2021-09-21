load("data/02b_scale_pca_umap_neighbors_clusters_dim50_multiRes.Rdata")
source("scripts/packages.R")

# going with res 0.6 for now 
Idents(object = cca.obj) <- "integrated_snn_res.0.6"
cca.obj$seurat_clusters = cca.obj$integrated_snn_res.0.6
cca.obj.meta = cca.obj@meta.data
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# get number of cells per cluster per sample
n_cells.cca.obj <- FetchData(cca.obj, 
                             vars = c("ident", "orig.ident")) %>%
  dplyr::count(ident, orig.ident) %>%
  tidyr::spread(ident, n)

write.csv(n_cells.cca.obj, "results/n_cells_per_cluster_per_sample.csv")

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

png(filename = "figures/cca.obj_umap_split_by_BatchID.png", width = 1920, height = 1080)
DimPlot(cca.obj, split.by = "batchID", raster = F, ncol = 3, label = T, label.size = 3)
dev.off()

png(filename = "figures/cca.obj_umap_split_by_batch.readlen.png", width = 1920, height = 1080/2)
DimPlot(cca.obj, split.by = "batch.readlen", raster = F, ncol = 3, label = T, label.size = 3)
dev.off()

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

metrics <-  c("nCount_RNA", "nFeature_RNA", "percent.mt")
DefaultAssay(cca.obj) <- "RNA"

FeaturePlot(cca.obj, 
            reduction = "umap", 
            features = metrics,
            pt.size = 0.4, 
            sort.cell = TRUE,
            min.cutoff = 'q10',
            label = TRUE, split.by = "batch.readlen", ncol = 3)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

rownames(n_cells_per_cluster_per_sample) <- n_cells_per_cluster_per_sample$orig.ident

n_cells.matrix = n_cells_per_cluster_per_sample %>% select(-1)
n_cells.matrix[is.na(n_cells.matrix)] <- 0

n_cells.matrix = n_cells.matrix / rowSums(n_cells.matrix) * 100
n_cells.matrix$sample = n_cells_per_cluster_per_sample$orig.ident

n_cells.matrix.joined = left_join(n_cells.matrix, cca.obj.meta, by = c("sample" = "orig.ident")) 

gg3.join = n_cells.matrix.joined %>% select(1:26,batchID, batch.readlen, sample, library)
gg3.melt = reshape2::melt(gg3.join)

colnames(gg3.melt) <- c("batchID", "batch.ReadLen", "sample", "library", "cluster.0.6", "value")

png("figures/distibution_of_cells_per_cluster_per_sample_per_cluster.png", width = 4000, height = 10000, res = 200)


ggplot(gg3.melt,
       aes(fill=batch.ReadLen, y=value, x=sample, label = batchID)) + 
  geom_bar(position="dodge", stat="identity") + theme(legend.position = "none", 
       axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_fill_viridis(discrete = T) +
  facet_grid(cluster.0.6 ~ library , scales = "free") + theme(axis.line=element_line())


dev.off()













#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++








