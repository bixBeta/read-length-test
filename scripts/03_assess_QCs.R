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

vln.feat.list = list()
for (i in 1:26) {
  vln.feat.list[[i]] = 
    VlnPlot(cca.obj, features = "nFeature_RNA", split.by = "library", pt.size = 0, group.by = "orig.ident", idents = i-1)
  
  names(vln.feat.list)[[i]] <- paste0("cluster__",i-1)
}

wrap_plots(vln.feat.list,ncol = 1,widths = 1920, heights = 333)


vln.count.list = list()
for (i in 1:26) {
  vln.count.list[[i]] = 
    VlnPlot(cca.obj, features = "nCount_RNA", split.by = "library", pt.size = 0, group.by = "orig.ident", idents = i-1)
  
  names(vln.count.list)[[i]] <- paste0("cluster__",i-1)
}

wrap_plots(vln.count.list,ncol = 1,widths = 1920, heights = 333)



cca.obj.meta$phenotypeID[cca.obj.meta$phenotypeID == 1] <- "control"
cca.obj.meta$phenotypeID[cca.obj.meta$phenotypeID == 2] <- "case"

cca.obj <- AddMetaData(cca.obj, metadata = cca.obj.meta$phenotypeID, col.name = "phenotype")
cca.obj.meta = cca.obj@meta.data

cca.obj$rlt.Condition <- paste0(cca.obj$readlength, ".", cca.obj$phenotype, ".", cca.obj$day)
cca.obj.meta = cca.obj@meta.data


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



save.image("data/sep21.Rdata")

cca.obj.meta %>% 
  ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=percent.mt)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250) +
  facet_wrap(~orig.ident)

cca.obj$log10GenesPerUMI <- log10(cca.obj$nFeature_RNA) / log10(cca.obj$nCount_RNA)
cca.obj.meta = cca.obj@meta.data

# Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI
cca.obj.meta %>%
  ggplot(aes(x=log10GenesPerUMI, color = orig.ident, fill=orig.ident)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8)
