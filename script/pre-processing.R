library(Seurat)
library(tidyverse)
library(harmony)

gc()
options(future.globals.maxSize = 30 * 1024^3)

# Load sample sheet
ss <- readxl::read_xlsx("data/ss.xlsx")

# Load Seurat objects
print("Loaded the objects from function.r")
seurat_list <- create_seurat_list(sample_sheet = ss,
                                  data_path = "data/raw")

# ============================================================
# FILTER: Calculate QC and subset EACH object
# ============================================================

seurat_list <- lapply(seurat_list, function(obj) {
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-", assay = "RNA")
  
  obj <- subset(
    obj, 
    subset = nFeature_RNA > 200 &
      nFeature_RNA < 9500 &
      percent.mt < 5
  )
  
  return(obj) 
})

# ============================================================
# NORMALIZE AND FIND VARIABLE FEATURES
# ============================================================

seurat_list <- lapply(X = seurat_list, FUN = function(x) {
  x <- NormalizeData(x)  # Normalize
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  return(x)
})

# ============================================================
# SEURAT INTEGRATION
# ============================================================
"features <- SelectIntegrationFeatures(object.list = seurat_list)
IntAnchors <- FindIntegrationAnchors(object.list = seurat_list, 
                                     anchor.features = features)
Int <- IntegrateData(anchorset = IntAnchors, k.weight = 50)"
print("Cancel CCA integration")
# ============================================================
# PRE-HARMONY: Check integration quality
# ============================================================


DefaultAssay(Int) <- "integrated"
Int <- ScaleData(Int)
Int <- RunPCA(Int, dims = 1:30)

p1 <- DimPlot(Int, reduction = "pca", group.by = "orig.ident")
p2 <- DimPlot(Int, reduction = "pca", group.by = "group")

png("analysis/pre_harmony_pca.png", res = 300,
    width = 10*300, height = 5*300)
p1 + p2
dev.off()

# ============================================================
# HARMONY CORRECTION
# ============================================================
Int <- RunHarmony(Int, group.by.vars = "group")

Int <- RunUMAP(Int, reduction = "harmony", dims = 1:30)

p1 <- DimPlot(Int, reduction = "harmony", group.by = "orig.ident")
p2 <- DimPlot(Int, reduction = "harmony", group.by = "group")

png("analysis/post_harmony_pca.png", res = 300,
    width = 10*300, height = 5*300)
p1 + p2
dev.off()

# ============================================================
# CLUSTERING
# ============================================================
Int <- FindNeighbors(Int, reduction = "harmony", dims = 1:30)
Int <- FindClusters(Int, resolution = 0.5)

# ============================================================
# UMAP VISUALIZATION
# ============================================================
Int <- RunUMAP(Int, reduction = "harmony", dims = 1:30)  # Already done above, can skip

p1 <- DimPlot(Int, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(Int, reduction = "umap", group.by = "group")
p3 <- DimPlot(Int, reduction = "umap", group.by = "seurat_clusters")

png("analysis/post_harmony_umap.png", res = 300,
    width = 15*300, height = 5*300)
p1 + p2 + p3
dev.off()

# ============================================================
# SAVE
# ============================================================
saveRDS(Int, file = "data/AllDNIntegration.rds")


sessionInfo()
