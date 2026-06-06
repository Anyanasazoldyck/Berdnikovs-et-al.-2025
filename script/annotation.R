library(Seurat)
library(tidyverse)
library(harmony)
library(ggplot2)
# load theme

#-----------------------------
# Plot themes
#-----------------------------



general_theme <- theme(
  text = element_text(size = 12, family = "ArialMT")
)

umap_theme <- theme(
  plot.title = element_text(hjust = 0.5),
  legend.position = "none",
  axis.text = element_blank(),
  axis.ticks = element_blank(),
  axis.title = element_text(size = 18),
  text = element_text(size = 14, family = "ArialMT")
)




# set dir
setwd("D:/May_journal")


#load the data 
sc.data <- readRDS("D:/May_journal/data/sc.data_harmony_integ.rds")
# plot by group 
group_col <- c(
  "Control"   = "#4ECDC4",  # teal
  "Wheeze"    = "#FFA552",  # orange
  "RSV"       = "#7B6DCC",  # purple
  "WheezeRSV" = "#6DBE57"   # green
)
groups <- unique(sc.data$group)

for (ch in groups) {
  
  print(ch)
  
  cells_highlight <- colnames(sc.data)[sc.data$group == ch]
  
  p <- DimPlot(
    sc.data,
    cells.highlight = cells_highlight,
    cols.highlight = group_col[[ch]],
    cols = "grey85"
  ) &
    ggtitle(ch)&
    umap_theme
  
  png(
    paste0("analysis/group_featureplot_", ch, ".png"),
    res = 300,
    width = 5 * 300,
    height = 5 * 300
  )
  
  print(p)
  dev.off()
}

print(dim(sc.data))#36601 13242
#==================================================
# Add cell cycle score
#=================================================
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
sc.data <- CellCycleScoring(sc.data,
                              s.features = s.genes,
                              g2m.features = g2m.genes,
                              set.ident = TRUE)
p<-DimPlot(sc.data,group.by = "Phase")&umap_theme&theme(legend.position = "right")
ggsave("analysis/CellCyclePhase.png", p,dpi=300,
       width = 6, height = 5)

Idents(sc.data)<-sc.data$seurat_clusters
DimPlot(sc.data, reduction="umap")






#======================
#increase reslotion 
#=======================
DefaultAssay(sc.data) <- "RNA"

sc.data <- FindNeighbors(sc.data, reduction = "harmony", dims = 1:30)

sc.data <- FindClusters(
  sc.data,
  graph.name = "RNA_snn",
  resolution = c(0.8)
)
sc.data<- RunUMAP(sc.data,
                  reduction = "harmony",
                  dims = 1:30,
                  reduction.name = "umap")

DimPlot(sc.data, reduction = "umap", label = T)

# I need 17 Clustert
#==========================================
# Explore more sub clust 
#========================================

p <- DimPlot(sc.data, group.by = "seurat_clusters",
                     reduction = "umap", label=T, label.size = 4)&umap_theme&
  labs(x= "UMAP-1", y="UMAP-2",title = "Seurat Clusters")



sc.data<- FindSubCluster(sc.data,
                         cluster = c(5),
                         resolution = 0.4,
                         graph.name="RNA_snn")

main_umap <- DimPlot(sc.data, group.by = "sub.cluster",
                     reduction = "umap", label=T, label.size = 4)&umap_theme&
  labs(x= "UMAP-1", y="UMAP-2",title = "Seurat Clusters")

ggsave("analysis/umap_final.png", main_umap,dpi=300,
       width = 6, height = 5)


#=========================================
# Return to FindAllMarkers
#=========================================
Idents(sc.data)<- sc.data$sub.cluster
levels(sc.data)

sc.data<-JoinLayers(sc.data)
markers <- FindAllMarkers(sc.data, only.pos = T, 
                          min.pct = 0.25, min.diff.pct = 0.25)
top.20<- markers %>% dplyr::group_by(cluster)%>% slice_head(n=20)
write.csv(markers,"data/markers.csv")
write.csv(top.20,"data/markers_top.20.csv")
top.5<- markers %>% dplyr::group_by(cluster)%>% slice_head(n=5)
p<- DoHeatmap(sc.data,features = top.5$gene)
ggsave("analysis/dohm_top5.png", p,dpi=300,
       width = 10, height = 10)
#=======================
# Explor markers 
#========================
# Marker vectors from Vieira Braga et al. and markers from LungMAP-style epithelial annotation

cell.state.markers <- list(
  Basal1 = c("TP63", "NPPC"),
  Basal2 = c("KRT5", "KRT14"),
  Deuterosomal<- c("DEUP1", "FOXJ1"),
  club<-c("SCGB1A1" ,"SCGB3A1"),
  Ionocytes <- c("FOXI1", "CFTR"),
  
  Cycling <- c("MKI67", "TOP2A", "CENPF"),
  Goblet1 = c("CEACAM5", "S100A4", "MUC5AC", "KRT4", "CD36"),
  Goblet2 = c("CEACAM5", "S100A4", "MUC5AC", "IDO1", "NOS2", "IL19", "CSF3", "CXCL10"),
  Ciliated1 = c("PROS1", "FXYD1"),
  Ciliated2 = c("CCL20", "ATP12A", "COX7A1", "AP2B1", "SYT5")
)

vln.p <- VlnPlot(sc.data,
                 features = unique(unlist(cell.state.markers)),
                 stack = T, flip = T) + 
  theme(legend.position = "none", 
        axis.title.y = element_text(size = rel(1), angle = 0), 
        axis.text.y = element_text(size = rel(1))) 
vln.p
ggsave("analysis/markers_vln2.png", vln.p,dpi=300,
       width = 10, height = 10)

#--------------------------
# annotate 
#==========================
Idents(sc.data)<-sc.data$sub.cluster
levels(sc.data)
new.ident <- c(
 "0_0"= "Secretory_Precursor",  
 
 "0_1"="progenitors",# 0
  "2"="Ciliated",          # 2
  "3"="Ciliated",          # 3
  "4"="Club Progenrator",        # 4
  "5"="Club",              # 5
  "6_1"="Basal.2",             # 6_1
  "7"="Goblet",            # 7
  "8"="Deuterosomal",      # 8
  
 "9"="Ciliated",   # 9
  
 "10"="Goblet",            # 10
  "11"="Ciliated",          # 11
  "12"="Ciliated",          # 12
  "6_0"="Basal.1",             # 6_0
  "6_2"="Basal.2",             # 6_2
  "6_4"="Activated_Basal",   # 6_4
  "6_3"="Cycling_Basal"      # 6_3
)

sc.data <- RenameIdents(sc.data, new.ident)
sc.data <- AddMetaData(
  sc.data,
  metadata = Idents(sc.data),
  col.name = "cell_types"
)




#===================
# Plot final umap ----
#==================

p2<-dittoSeq::dittoBarPlot(sc.data, var="group",
                       group.by = "cell_types")

p<- DimPlot(sc.data,
            reduction = "umap",
            group.by = "cell_types",
            label = F)
p1<- DimPlot(sc.data,
             reduction = "umap",
             group.by = "sub.cluster",
             label = T)

ggsave("analysis/umap_final.png", p,dpi=300,
       width = 6, height = 5)

ggsave("analysis/umap_deurat_clust.png", p1,dpi=300,
       width = 6, height = 5)
ggsave("analysis/composition.png", p2,dpi=300,
       width = 10, height = 5)
saveRDS(sc.data,"data/final_harmony_seurat.rds")
#===========================
# Find top marker of each cell type
#============================

