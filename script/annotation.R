library(Seurat)
library(tidyverse)
library(harmony)


# set dir
setwd("D:/May_journal")


#load the data 
sc.data<- readRDS("data/sc.data.rds")


#check epithelial marker expression
p<- FeaturePlot(sc.data,features = "EPCAM", reduction = "umap", order =T)
p

ggsave("analysis/EPCAM_featureplot.png",p, dpi=300,
       width = 5, height = 5)
#-================================
# Felipe A. Vieira Braga 
#===================================

  airway_markers <- list(
    
    Ciliated_1 = c(
      "FOXJ1",
      "PIFO",
      "KCTD12"
    ),
    
    Type_2_alveolar = c(
      "SFTPC"
    ),
    
    Club = c(
      "SCGB1A1",
      "MUC5AC",
      "MUC5B"
    ),
    
    Basal_1 = c(
      "KRT5",
      "TP63",
      "NPPC"
    ),
    
    Basal_2 = c(
      "KRT5",
      "TP63"
    ),
    
    Goblet_1 = c(
      "MUC5AC"
    ),
    
    Goblet_2 = c(
      "MUC5AC",
      "CSF3"
    ),
    
    Ciliated_2 = c(
      "APOD"
    ),
    
    Ionocytes = c(
      "FOXI1",
      "CFTR"
    )) 

epi_markers<- unique(unlist(airway_markers))
  
  
plots <- FeaturePlot(
  sc.data,
  features = epi_markers,
  reduction = "umap",
  combine = FALSE,
  pt.size = 0.5,
  order=T,
  min.cutoff = 0,
  max.cutoff = 4,
  keep.scale = "all"
  
)

plots<- lapply(plots, function(p){
  p&theme_void()&
    theme(
      legend.position = "right",
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 12
      ))

  
})
  
combined_plot <- patchwork::wrap_plots(
  plots,
  ncol = 3,
  guides = "collect"
) &
  theme(
    legend.position = "right"
  )

ggsave("analysis/featureplot.png", combined_plot,dpi=300,
       width = 8, height = 5)

#==================================================
# Add cell cycle score
#=================================================
# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
sc.data <- CellCycleScoring(sc.data,
                              s.features = s.genes,
                              g2m.features = g2m.genes,
                              set.ident = TRUE)
DimPlot(sc.data,group.by = "Phase")
#saveRDS(sc.data, file = "data/sc.data.rds")

#=========================================
# FindAllMarkers
#=========================================

Idents(sc.data)<- sc.data$seurat_clusters
sc.data<-JoinLayers(sc.data)
markers <- FindAllMarkers(sc.data, only.pos = T, 
                          min.pct = 0.25, min.diff.pct = 0.25)
top.20<- markers %>% dplyr::group_by(cluster)%>% slice_head(n=20)
write.csv(markers,"data/markers.csv")


