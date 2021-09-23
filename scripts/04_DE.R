# all.pbmc.Markers <- FindAllMarkers(cca.obj, only.pos = TRUE, logfc.threshold = 0.25)
# write.csv(all.pbmcs.markers, "results/all.pbmcs.Markers.csv")
# 

source("scripts/packages.R")

DefaultAssay(cca.obj) <- "RNA"

long.v.short = list()
for (i in 1:length(levels(Idents(cca.obj)))) {
  long.v.short[[i]] <- FindMarkers(cca.obj, ident.1 = "long", 
                                            ident.2 = "short", group.by = "readlength",
                                            subset.ident = levels(Idents(cca.obj))[i], 
                                            only.pos = F, logfc.threshold = 0.1)
  
  names(long.v.short)[[i]] <- paste0("CLUSTER__", levels(Idents(cca.obj))[i])
  
}

filtered.long.v.short = 
  lapply(long.v.short, function(x){
    x %>% dplyr::filter(p_val_adj < 0.05)
  })

long.v.short.de.genes.list = unique(unname(unlist(sapply(filtered.long.v.short, function(x){rownames(x)}))))
DefaultAssay(cca.obj) <- "RNA"
all.genes = rownames(cca.obj)
avg.expression.cca.obj = 
  AverageExpression(cca.obj, assays = "RNA", group.by = "orig.ident", 
                    features = all.genes, 
                    return.seurat = T)

avg.scaled.matrix.all.genes = avg.expression.cca.obj@assays$RNA@scale.data
save.image("data/sep22.Rdata")


test = cca.obj.meta %>% select(1,5:11,17,18) %>% group_by(orig.ident)

l = seq(-2,2, .1)
heatmap.anno.long.v.short = unique(test) %>% select(1,2,4,5,8,10)

my_colour = list(
  rlt.Condition = 
    c(long.case.D1 = '#000075',
                    long.case.D2 = '#3cb44b',
                    long.control.D1 = '#ffe119',
                    long.control.D2 = '#4363d8',
                    short.case.D1 = '#f58231',
                    short.case.D2 = '#911eb4',
                    short.control.D1 = '#46f0f0',
                    short.control.D2 = '#f032e6'),
  
  batch.readlen = 
    c(early = '#bcf60c', 
    late.long = '#fabebe', 
    late.short = '#008080'), 

  
  sex = c(M = "blue", F = "gray")
)




png(filename = "figures/pheatmap.long.v.short.de.png", width = 4000, height = 20000, res = 100)
pheatmap::pheatmap(avg.scaled.matrix.all.genes[rownames(avg.scaled.matrix.all.genes) %in% long.v.short.de.genes.list,], scale = "none", 
         color=colorRampPalette(c("navy", "white", "red"))(41), breaks = l, 
         cluster_rows = T,
         cluster_cols = T,
         main = "long v short",
         annotation_col = heatmap.anno.long.v.short %>% column_to_rownames("orig.ident") %>% select(-library), 
         annotation_colors = my_colour,
         border_color = "grey60", 
         cellwidth = 20, cellheight = 20, 
         fontsize = 18, treeheight_col = 150)

dev.off()


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

long.case.D1.v.long.control.D1 = list()

for (i in 1:length(levels(Idents(cca.obj)))) {
  long.case.D1.v.long.control.D1[[i]] <- FindMarkers(cca.obj, ident.1 = "long.case.D1", 
                                   ident.2 = "long.control.D1", group.by = "rlt.Condition",
                                   subset.ident = levels(Idents(cca.obj))[i], 
                                   only.pos = F, logfc.threshold = 0.1)
  
  names(long.case.D1.v.long.control.D1)[[i]] <- paste0("CLUSTER__", levels(Idents(cca.obj))[i])
  
}



short.case.D1.v.short.control.D1 = list()

for (i in 1:length(levels(Idents(cca.obj)))) {
  short.case.D1.v.short.control.D1[[i]] <- FindMarkers(cca.obj, ident.1 = "short.case.D1", 
                                                     ident.2 = "short.control.D1", group.by = "rlt.Condition",
                                                     subset.ident = levels(Idents(cca.obj))[i], 
                                                     only.pos = F, logfc.threshold = 0.1)
  
  names(short.case.D1.v.short.control.D1)[[i]] <- paste0("CLUSTER__", levels(Idents(cca.obj))[i])
  
}

DefaultAssay(cca.obj) <- "RNA"
cca.obj = FindVariableFeatures(cca.obj, selection.method = "vst", nfeatures = 2000)
top100 <- head(VariableFeatures(cca.obj), 100)


