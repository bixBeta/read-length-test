library(Seurat)
library(dplyr)
library(tibble)
library(ggplot2)
library(UpSetR)
library(harmony)
library(patchwork)
library(tidyr)
library(viridis)
library(Matrix.utils)
library(SingleCellExperiment)
library(magrittr)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
cca.obj = readRDS("no__Filter__INTEGRATED_scaled_pca_umap_neighbors_clusters_dim50_multiRes.Rds")
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

DefaultAssay(cca.obj) <- "RNA"
cca.obj.meta = cca.obj@meta.data

Idents(object = cca.obj) <- "integrated_snn_res.0.4"
umap.0.4 = DimPlot(cca.obj, reduction = "umap", label = TRUE, label.size = 6, raster = F) + ggtitle("res 0.4")

Idents(object = cca.obj) <- "integrated_snn_res.0.5"
umap.0.5 = DimPlot(cca.obj, reduction = "umap", label = TRUE, label.size = 6, raster = F) + ggtitle("res 0.5")

Idents(object = cca.obj) <- "integrated_snn_res.0.6"
umap.0.6 = DimPlot(cca.obj, reduction = "umap", label = TRUE, label.size = 6, raster = F) + ggtitle("res 0.6")

Idents(object = cca.obj) <- "integrated_snn_res.0.8"
umap.0.8 = DimPlot(cca.obj, reduction = "umap", label = TRUE, label.size = 6, raster = F) + ggtitle("res 0.8")

png(filename = "figures__splitByBatch__noFilter/UMAP_res0.4_0.5_0.6_0.8.png", width = 1920, height = 1080)
umap.0.4 + umap.0.5 + umap.0.6 + umap.0.8
dev.off()

Idents(object = cca.obj) <- "integrated_snn_res.0.6"
cca.obj$seurat_clusters = cca.obj$integrated_snn_res.0.6

cca.obj.meta = cca.obj@meta.data

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

FeaturePlot(cca.obj, features = "percent.mt", label = T, split.by = "sex")
DimPlot(cca.obj, reduction = "umap", label = TRUE, label.size = 6, raster = F, split.by = "sex")


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# get number of cells per cluster per sample
n_cells_per_cluster_per_sample <- FetchData(cca.obj, 
                             vars = c("ident", "orig.ident")) %>%
  dplyr::count(ident, orig.ident) %>%
  tidyr::spread(ident, n)

write.csv(n_cells_per_cluster_per_sample, "results__splitByBatch__noFilter/n_cells_per_cluster_per_sample.csv")


rownames(n_cells_per_cluster_per_sample) <- n_cells_per_cluster_per_sample$orig.ident

n_cells.matrix = n_cells_per_cluster_per_sample %>% select(-1)
n_cells.matrix[is.na(n_cells.matrix)] <- 0

n_cells.matrix = n_cells.matrix / rowSums(n_cells.matrix) * 100
n_cells.matrix$sample = n_cells_per_cluster_per_sample$orig.ident

n_cells.matrix.joined = left_join(n_cells.matrix, cca.obj.meta, by = c("sample" = "orig.ident")) 

gg3.join = n_cells.matrix.joined %>% select(1:30,batchID, batch.readlen, sample, library, sex)
gg3.melt = reshape2::melt(gg3.join)

colnames(gg3.melt) <- c("batchID", "batch.ReadLen", "sample", "library", "sex", "cluster.0.6", "value")

png("figures__splitByBatch__noFilter/distibution_of_cells_per_cluster_per_sample_per_cluster.png", width = 4000, height = 10000, res = 200)


ggplot(gg3.melt,
       aes(fill=batch.ReadLen, y=value, x=sample, label = batchID)) + 
  geom_bar(position="dodge", stat="identity") + theme(legend.position = "none", 
                                                      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_fill_viridis(discrete = T) +
  facet_grid(cluster.0.6 ~ library , scales = "free") + theme(axis.line=element_line())


dev.off()

png("figures__splitByBatch__noFilter/SEX__distibution_of_cells_per_cluster_per_sample_per_cluster.png", width = 4000, height = 10000, res = 200)


ggplot(gg3.melt,
       aes(fill=batch.ReadLen, y=value, x=sample, label = batchID)) + 
  geom_bar(position="dodge", stat="identity") + theme(legend.position = "none", 
                                                      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_fill_viridis(discrete = T) +
  facet_grid(cluster.0.6 ~ sex , scales = "free") + theme(axis.line=element_line())


dev.off()


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# extract pseudoCounts
sce <- as.SingleCellExperiment(cca.obj)
assays(sce)

groups <- colData(sce)[, c("seurat_clusters", "orig.ident")]
pb <- aggregate.Matrix(t(counts(sce)), 
                       groupings = groups, fun = "sum") 

splitf <- sapply(stringr::str_split(rownames(pb), 
                                    pattern = "_",  
                                    n = 2),`[`, 1)

pb <- split.data.frame(pb, 
                       factor(splitf)) %>%
  lapply(function(u) 
    set_colnames(t(u), 
                 stringr::str_extract(rownames(u), "(?<=_)[:alnum:]+")))

pb.matrix = lapply(pb, function(x){as.matrix(x)})

save(pb,pb.matrix, file = "data__splitByBatch__noFilter/PseudoBulk_Counts.Rdata")
save.image(file = "data__splitByBatch__noFilter/noFilter__splitByBatch__SOBJ__scripts.Rdata")


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++




