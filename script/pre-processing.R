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


#=============================================
# pre-integration Batch diagnosis
#===========================================
#===========================================
sc_merge_unint<- merge(x = seurat_list[[1]],
                       y = seurat_list[-1])
sc_merge_unint<-NormalizeData(sc_merge_unint)
sc_merge_unint<-FindVariableFeatures(sc_merge_unint,
                                     selection.method = "vst",
                                     nfeatures = 2000)
sc_merge_unint<-ScaleData(sc_merge_unint, vars.to.regress = "percent.mt")
sc_merge_unint<-RunPCA(sc_merge_unint,reduction.name ="unintegrated_pca")

pca_preint <- DimPlot(sc_merge_unint,
                      reduction = "unintegrated_pca",
                      group.by = "orig.ident")

pca_preint_group <- DimPlot(sc_merge_unint,
                            reduction = "unintegrated_pca",
                            group.by = "group")

png("analysis/pre_int_pca.png", res = 300,
    width = 10*300, height = 5*300)
pca_preint + pca_preint_group
dev.off()

# ============================================================
# SEURAT INTEGRATION
# ============================================================
features <- SelectIntegrationFeatures(object.list = seurat_list)
IntAnchors <- FindIntegrationAnchors(object.list = seurat_list, 
                                     anchor.features = features)
Int <- IntegrateData(anchorset = IntAnchors, k.weight = 50)

# CCA reducation 
DefaultAssay(Int) <- "integrated"

Int<-ScaleData(Int, vars.to.regress = "percent.mt")
Int<-RunPCA(Int, reduction.name = "cca_pca")
pca_int<-DimPlot(Int, reduction = "cca_pca", group.by = "orig.ident")
pca_int_group<-DimPlot(Int, reduction = "cca_pca", group.by = "group")

png("analysis/cca_pca.png", res = 300,
    width = 10*300, height = 5*300)
pca_int + pca_int_group
dev.off()

# find clusters 
Int<-RunUMAP(Int, reduction= "cca_pca", dims = 1:30)
Int <- FindNeighbors(Int, reduction = "cca_pca", dims = 1:30)
Int <- FindClusters(Int, resolution = 0.5)

p1 <- DimPlot(Int, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(Int, reduction = "umap", group.by = "group")
p3 <- DimPlot(Int, reduction = "umap", group.by = "seurat_clusters")

png("analysis/post_CCA_umap_main.png", res = 300,
    width = 15*300, height = 5*300)
p1 + p2 + p3
dev.off()


#save my cca int object
saveRDS(Int,"data/cca_int_object.rds")

print("CCA integration process done")


# ============================================================
# HARMONY CORRECTION
# ============================================================
sc.data <- RunHarmony(sc_merge_unint, group.by.vars = "orig.ident",
                      reduction.use = "unintegrated_pca",
                      dims.use = NULL,
                      reduction.save = "harmony"
                      )

pca_harmony <- DimPlot(sc.data, reduction = "harmony", group.by = "orig.ident")
pca_harmony_group <- DimPlot(sc.data, reduction = "harmony", group.by = "group")

png("analysis/post_harmony_pca.png", res = 300,
    width = 10*300, height = 5*300)
pca_harmony + pca_harmony_group
dev.off()



#======================get umap in harmony========================
sc.data <- RunUMAP(sc.data, reduction = "harmony", dims = 1:30)
sc.data <- FindNeighbors(sc.data, reduction = "harmony", dims = 1:30)
sc.data <- FindClusters(sc.data, resolution = 0.5)



p1 <- DimPlot(sc.data, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(sc.data, reduction = "umap", group.by = "group")
p3 <- DimPlot(sc.data, reduction = "umap", group.by = "seurat_clusters")

png("analysis/post_harmony_umap_main.png", res = 300,
    width = 15*300, height = 5*300)
p1 + p2 + p3
dev.off()




#save sc.data harmony integ 
saveRDS(sc.data,"data/sc.data_harmony_integ.rds")

print("both Seurat integration CCA and harmony similar PCA. but sample 
wheeze1 cluster in a cluster on its own in CCA,
      Will proceed with harmony.")
