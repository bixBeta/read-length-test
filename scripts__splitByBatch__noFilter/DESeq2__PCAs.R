library(DESeq2)

preLim.colData = cca.obj.meta %>% rownames_to_column("temp") %>% 
  select(-temp, -library, -seurat_clusters, -percent.mt, -matches(c("integrate", "_RNA"))) %>%
  unique()


colData = preLim.colData[order(colData$orig.ident, colnames(pb.matrix[[1]])),]
rownames(colData) = colData$orig.ident

dds.list = list()
for (i in 1:length(pb.matrix)) {
  
  dds.list[[i]] <- DESeqDataSetFromMatrix(countData = pb.matrix[[i]],
                                           colData = colData[colData$orig.ident %in% colnames(pb.matrix[[i]]),], design = ~ day)
  names(dds.list)[[i]] <- paste0("dds__", names(pb.matrix)[[i]])
}

rld.list = 
  lapply(dds.list, function(x){
    rlog(x, blind = T)
  })


pca.list = list()
for (i in 1:length(rld.list)) {
  pca.list[[i]] <- DESeq2::plotPCA(rld.list[[i]], intgroup = c("batch.readlen","sex")) + ggtitle(names(rld.list)[[i]])
  names(pca.list)[[i]] <- paste0("pca__", names(rld.list)[[i]])
  
}

png(filename = "figures__splitByBatch__noFilter/pca__noFilter__logNorm.png", width = 4000, height = 2000)
wrap_plots(pca.list, ncol = 5, nrow = 6, heights = 800, widths = 800)
dev.off()


for (i in 1:length(pca.list)) {
  png(paste0("figures__splitByBatch__noFilter/PCA/", names(pca.list)[[i]] ,".png"),width = 500, height = 500)
  print(pca.list[[i]])
  Sys.sleep(1)
  dev.off()
} 


