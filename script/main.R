library(Seurat)
library(tidyverse)
library(harmony)

gc()
options(future.globals.maxSize = 30 * 1024^3)


prsc_merged("one sc_mergedegration analysis")
# Laod seurat list 

sc_merged <- merge(
  x = seurat_list[[1]],
  y = seurat_list[-1]
)
# ============================================================
# Normalize and Scale
# ============================================================


sc_merged <- ScaleData(sc_merged)
sc_merged <- RunPCA(sc_merged, dims = 1:30)

p1 <- DimPlot(sc_merged, reduction = "pca", group.by = "orig.ident")
p2 <- DimPlot(sc_merged, reduction = "pca", group.by = "group")

png("analysis/pre_harmony_pca_main.png", res = 300,
    width = 10*300, height = 5*300)
p1 + p2
dev.off()

sc_merged <- RunHarmony(sc_merged, group.by.vars = "orig.ident")

sc_merged <- RunUMAP(sc_merged, reduction = "harmony", dims = 1:30)

p1 <- DimPlot(sc_merged, reduction = "harmony", group.by = "orig.ident")
p2 <- DimPlot(sc_merged, reduction = "harmony", group.by = "group")

png("analysis/post_harmony_pca_main.png", res = 300,
    width = 10*300, height = 5*300)
p1 + p2
dev.off()

# ============================================================
# CLUSTERING
# ============================================================
sc_merged <- FindNeighbors(sc_merged, reduction = "harmony", dims = 1:30)
sc_merged <- FindClusters(sc_merged, resolution = 0.5)

# ============================================================
# UMAP VISUALIZATION
# ============================================================
sc_merged <- RunUMAP(sc_merged, reduction = "harmony", dims = 1:30)  # Already done above, can skip

p1 <- DimPlot(sc_merged, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(sc_merged, reduction = "umap", group.by = "group")
p3 <- DimPlot(sc_merged, reduction = "umap", group.by = "seurat_clusters")

png("analysis/post_harmony_umap_main.png", res = 300,
    width = 15*300, height = 5*300)
p1 + p2 + p3
dev.off()

#-================================
# Annotation 
#===================================

airway_markers <- list(
  
  #========================
  # Basal subsets
  #========================
  
  Basal_general = c(
    "KRT5", "TP63", "KRT14"
  ),
  
  Basal_1_less_mature = c(
    "TP63", "NPPC"
  ),
  
  Basal_2_more_mature = c(
    "KRT5", "TP63"
  ),
  
  Activated_Basal = c(
    "POSTN", "SERPINB2", "JAG1", "SPRY1"
  ),
  
  Cycling_Basal = c(
    "MKI67", "TOP2A", "KRT13"
  ),
  
  
  #========================
  # Secretory subsets
  #========================
  
  Club_cells = c(
    "SCGB1A1", "SCGB3A1", "MUC5B", "MUC5AC"
  ),
  
  Goblet_1 = c(
    "MUC5AC", "CEACAM5", "S100A4", "KRT4", "CD36"
  ),
  
  Goblet_2 = c(
    "MUC5AC", "CEACAM5", "IDO1", "NOS2",
    "IL19", "CSF3", "CXCL10"
  ),
  
  
  #========================
  # Ciliated subsets
  #========================
  
  Ciliated_general = c(
    "FOXJ1", "PIFO", "SNTN", "HYDIN", "DNAH12"
  ),
  
  Ciliated_1 = c(
    "PROS1", "FXYD1", "KCTD12"
  ),
  
  Ciliated_2 = c(
    "CCL20", "ATP12A", "AP2B1", "SYT5"
  ),
  
  Mucous_Ciliated = c(
    "FOXJ1", "PIFO", "MUC5AC", "CEACAM5"
  ),
  
  
  #========================
  # Specialized epithelial
  #========================
  
  Ionocytes = c(
    "FOXI1", "CFTR", "SCNN1B", "ASCL3"
  ),
  
  Neuroendocrine = c(
    "CHGA", "ASCL1", "INSM1", "HOXB5"
  ),
  
  Tuft_cells = c(
    "DCLK1", "ASCL2"
  ),
  
  Alveolar_1 = c(
    "AGER", "MGP", "CAV1"
  ),
  
  Alveolar_2 = c(
    "SFTPC", "SFTPB", "SFTPA2",
    "SERPINA1", "NKX2-1"
  ),
  
  
)


airway_markers_nodup <- list()
seen <- character()

for (nm in names(airway_markers)) {
  genes <- airway_markers[[nm]]
  genes <- genes[!genes %in% seen]
  genes <- genes[genes %in% rownames(sc_merged)]
  airway_markers_nodup[[nm]] <- genes
  seen <- c(seen, genes)
}

airway_markers_nodup <- airway_markers_nodup[lengths(airway_markers_nodup) > 0]

p <- DotPlot(sc_merged, features = airway_markers_nodup) + RotatedAxis()+ NoLegend()
png("analysis/dotplot.png", res = 300,
    width = 25*300, height = 5*300)
p
dev.off()




#==================================================
# Add cell cycle score
#=================================================
# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
sc_merged <- CellCycleScoring(sc_merged,
                           s.features = s.genes,
                           g2m.features = g2m.genes,
                           set.ident = TRUE)
DimPlot(sc_merged,group.by = "Phase")
saveRDS(sc_merged, file = "data/sc_merged.rds")



Idents(sc_merged)<- sc_merged$seurat_clusters
sc_merged<-JoinLayers(sc_merged)
markers <- FindAllMarkers(sc_merged, only.pos = T, 
                          min.pct = 0.25, min.diff.pct = 0.25)
top.20<- markers %>% dplyr::group_by(cluster)%>% slice_head(n=20)
write.csv(markers,"data/markers.csv")
#=========================================================
# See the group dist
#========================================================
group_cal <-list("Control"="red",
"Wheeze"="blue",
"RSV"="violet",
"WheezeRSV"="darkgreen")
group <- unique(sc_merged$group)


for (ch in group){
  print(ch)
  cells_highlight = colnames(sc_merged)[sc_merged$group==ch]
  p <- DimPlot(
    sc_merged,
    cells.highlight = cells_highlight,
    cols.highlight = group_cal[[ch]],
    cols = "grey85"
  ) +
    ggtitle(ch)
  png(paste0("analysis/group_feratureplot",ch,".png"), res = 300,   width = 6*300, height = 5*300)
      
  print(p)
  dev.off()
}

