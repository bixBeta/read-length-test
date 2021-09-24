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

save(pb,pb.matrix, file = "data/PseudoBulk_Counts.Rdata")
save.image(file = "data/sep23.Rdata")
