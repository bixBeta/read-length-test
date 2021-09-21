library("Seurat")

sobj = readRDS("filtered_seurat_object_with_metaData_for_RLT.Rds")

sobj.list <- SplitObject(sobj, split.by = "batchID")

sobj.list <- lapply(X = sobj.list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = sobj.list)

sobj.anchors <- FindIntegrationAnchors(object.list = sobj.list, anchor.features = features)

sobj.integrated <- IntegrateData(anchorset = sobj.anchors)

saveRDS(sobj.integrated, "sobj.integrated.on.batchID.top.2000.features.Rds")
