#=================================================================================
# Analysis
#=================================================================================
sc.data<- readRDS("data/final_harmony_seurat.rds")
Idents(sc.data)<-sc.data$cell_types
table(sc.data$group)
"Control       RSV    Wheeze WheezeRSV 
3070      4004      1880      2450 "
dim(sc.data)
#===========
# Composition analysis
#============
# Basal ====
sc.basal <- subset (sc.data,
                    idents  = c("Basal.1","Basal.2","Basal.Cycling"))
# 1. Prepare data
df_pop <- Seurat::FetchData(sc.basal, vars = c("group", "cell_types")) %>%
  dplyr::count(group,cell_types) %>%
  group_by(group) %>%
  mutate(freq = n / sum(n))

# 2. Automatically generate enough colors if group_cols is broken
unique_groups <- unique(df_pop$group)
cell_cols <- c("darkblue","steelblue","lightblue","lightyellow")

# 3. Build the plot
p <- ggplot(
  df_pop,
  aes(
    y = group,
    x = freq * 100,
    fill = cell_types
  )
) +
  geom_col(
    color = "black",
    linewidth = 0.3
  ) +
  scale_fill_manual(values = cell_cols) +
  labs(
    x = "Composition (%)",
    y = NULL,
    fill = NULL
  ) + # <--- Fixed: The '+' must be on this line!
  theme_classic(base_family = "ArialMT") +
  theme(
    axis.text = element_text(size = 12, colour = "black"),
    axis.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    legend.position = "top"
  )

# Display the plot
print(p)


ggsave("analysis/basal composition.png", p,dpi=300,
       width = 8, height = 4)







#============
# Dotplot basal
#============
basal_markers= c("TP63", "NPPC","KRT5", 
                 "KRT14","MKI67", "TOP2A",
                 "STMN1","SPRY1","POSTN")
feature.cols <- c(
  "#08306B",  # dark blue
  "#2171B5",  # blue
  "#FFFFBF",  # yellow
  "#FDAE61",  # orange
  "#D73027"   # red
)
p<-FeaturePlot(sc.data,
           features = basal_markers,
           cols =  feature.cols,
           order=T)&umap_theme

p <- p + patchwork::plot_layout(guides = "collect")
ggsave("analysis/basal featureplot.png", p,dpi=300,
       width = 10, height = 10)

#========================================
# progen 
#========================================
levels(sc.data)
# Basal ====
sc.sec <- subset (sc.data,
                  idents  = c("Early Progintor","Club precursor","Suprabasal" ))
# 1. Prepare data
df_pop <- Seurat::FetchData(sc.sec, vars = c("group", "cell_types")) %>%
  dplyr::count(group,cell_types) %>%
  group_by(group) %>%
  mutate(freq = n / sum(n))

# 2. Automatically generate enough colors if group_cols is broken
unique_groups <- unique(df_pop$group)
cell_cols <- c("darkred","maroon","pink")

# 3. Build the plot
p <- ggplot(
  df_pop,
  aes(
    y = group,
    x = freq * 100,
    fill = cell_types
  )
) +
  geom_col(
    color = "black",
    linewidth = 0.3
  ) +
  scale_fill_manual(values = cell_cols) +
  labs(
    x = "Composition (%)",
    y = NULL,
    fill = NULL
  ) + # <--- Fixed: The '+' must be on this line!
  theme_classic(base_family = "ArialMT") +
  theme(
    axis.text = element_text(size = 12, colour = "black"),
    axis.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    legend.position = "top"
  )

# Display the plot
print(p)


ggsave("analysis/progen composition.png", p,dpi=300,
       width = 8, height = 4)


#============
# Dotplot basal
#============
Progen_markers= c("HES4", "NEAT1","KRT13", 
                 "SPRR2A", "DUSP5",  "SPRR3",
                 "CALML5","KRTDAP","KRT6A")
feature.cols <- c(
  "#08306B",  # dark blue
  "#2171B5",  # blue
  "#FFFFBF",  # yellow
  "#FDAE61",  # orange
  "#D73027"   # red
)
p<-FeaturePlot(sc.data,
               features = Progen_markers,
               cols =  feature.cols,
               order=T)&umap_theme

p <- p + patchwork::plot_layout(guides = "collect")
ggsave("analysis/progen featureplot.png", p,dpi=300,
       width = 10, height = 10)



progen<- DotPlot(sc.sec,
                 features =Progen_markers,
                 cols = c("lightyellow","darkred"))
ggsave("analysis/progen dotplot.png", progen,dpi=300,
       width = 10, height = 4)
#====================================================
# Club + Goblet
#==================================================

sc.clubgob <- subset (sc.data, idents  = c("Club"  ,"Goblet.1" ,"Goblet.2" ,
                                           "Squamous Secretory"))

df_pop <- Seurat::FetchData(sc.clubgob, vars = c("group", "cell_types")) %>%
  dplyr::count(group,cell_types) %>%
  group_by(group) %>%
  mutate(freq = n / sum(n))

# 2. Automatically generate enough colors if group_cols is broken
unique_groups <- unique(df_pop$group)
cell_cols <- c("brown","darkorange","yellow","lightyellow")

# 3. Build the plot
p <- ggplot(
  df_pop,
  aes(
    y = group,
    x = freq * 100,
    fill = cell_types
  )
) +
  geom_col(
    color = "black",
    linewidth = 0.3
  ) +
  scale_fill_manual(values = cell_cols) +
  labs(
    x = "Composition (%)",
    y = NULL,
    fill = NULL
  ) + # <--- Fixed: The '+' must be on this line!
  theme_classic(base_family = "ArialMT") +
  theme(
    axis.text = element_text(size = 12, colour = "black"),
    axis.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    legend.position = "top"
  )

# Display the plot
print(p)


ggsave("analysis/club_goblet composition.png", p,dpi=300,
       width = 8, height = 4)

#============
# Dotplot basal
#============
secretory.lineage.markers <- list(
  Goblet = c("CEACAM5", "S100A4", "MUC5AC"),
  Goblet1 = c("CEACAM5", "S100A4", "MUC5AC", "KRT4"),
  Goblet2 = c("CEACAM5", "S100A4", "MUC5AC", "HLA-B", "CXCL16", "C15orf48"),
  Club = c("SCGB1A1", "SCGB3A1", "MSMB", "BPIFB1", "MUC5B")
)

secretory.lineage.markers<- unique(unlist(secretory.lineage.markers))

feature.cols <- c(
  "#08306B",  # dark blue
  "#2171B5",  # blue
  "#FFFFBF",  # yellow
  "#FDAE61",  # orange
  "#D73027"   # red
)
p<-FeaturePlot(sc.data,
               features = secretory.lineage.markers,
               cols =  feature.cols,
               order=T)&umap_theme

p <- p + patchwork::plot_layout(guides = "collect")
ggsave("analysis/secretory.lineage.markers.featureplot.png", p,dpi=300,
       width = 12, height = 10)



progen<- DotPlot(sc.clubgob,
                 features =secretory.lineage.markers,
                 cols = c("lightyellow","darkred"))
ggsave("analysis/secretory.lineage.markers.dotplot.png", progen,dpi=300,
       width = 14, height = 4)
#=================
# ciliated 
#================
levels(sc.data)
# Basal ====
# Ciliated subset
sc.ciliated <- subset(
  sc.data,
  idents = c("Mature Ciliated", "Deuterosomal", "Early Ciliated")
)

# Verify subset
table(Idents(sc.ciliated))
levels(sc.ciliated)
df_pop <- Seurat::FetchData(sc.ciliated, vars = c("group", "cell_types")) %>%
  dplyr::count(group,cell_types) %>%
  group_by(group) %>%
  mutate(freq = n / sum(n))

# 2. Automatically generate enough colors if group_cols is broken
unique_groups <- unique(df_pop$group)
cell_cols <- c(
  "Early Ciliated" = "#5E3C99",  # purple
  "Deuterosomal"   = "#8073AC",  # lavender
  "Mature Ciliated"     = "#1B9E77"  # teal green
 
)

# 3. Build the plot
p <- ggplot(
  df_pop,
  aes(
    y = group,
    x = freq * 100,
    fill = cell_types
  )
) +
  geom_col(
    color = "black",
    linewidth = 0.3
  ) +
  scale_fill_manual(values = cell_cols) +
  labs(
    x = "Composition (%)",
    y = NULL,
    fill = NULL
  ) + # <--- Fixed: The '+' must be on this line!
  theme_classic(base_family = "ArialMT") +
  theme(
    axis.text = element_text(size = 12, colour = "black"),
    axis.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    legend.position = "top"
  )

# Display the plot
print(p)


ggsave("analysis/ciliated composition.png", p,dpi=300,
       width = 8, height = 4)

#============
# Dotplot basal
#============
ciliated.markers <- list(
  Mature_Ciliated = c(
    "SNTN"
  ),
  
  Deuterosomal<- c("DEUP1", "FOXJ1")
)

ciliated.markers<- unique(unlist(ciliated.markers))

feature.cols <- c(
  "#08306B",  # dark blue
  "#2171B5",  # blue
  "#FFFFBF",  # yellow
  "#FDAE61",  # orange
  "#D73027"   # red
)
p<-FeaturePlot(sc.data,
               features = ciliated.markers,
               cols =  feature.cols,
               order=T)&umap_theme

p <- p + patchwork::plot_layout(guides = "collect")
ggsave("analysis/ciliated.markers.featureplot.png", p,dpi=300,
       width = 12, height = 10)



progen<- DotPlot(sc.ciliated,
                 features =ciliated.markers,
                 cols = c("lightyellow","darkred"))
ggsave("analysis/ciliated.markers.dotplot.png", progen,dpi=300,
       width = 10, height = 4)

