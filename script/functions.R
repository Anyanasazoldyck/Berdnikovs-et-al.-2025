#====================================================
# Create Seurat objects from multiple .h5 files
#====================================================

create_seurat_list <- function(sample_sheet, data_path = "data/raw") {
  
  seurat_list <- list()
  
  for (sample_name in sample_sheet$Sample_Code) {
    
    #========================================
    # Match metadata
    #========================================
    
    meta_row <- sample_sheet[
      sample_sheet$Sample_Code == sample_name,
    ]
    
    original_id <- meta_row$Sample_Name
    
    
    #========================================
    # File path
    #========================================
    
    file_path <- file.path(
      data_path,
      paste0(sample_name, ".h5")
    )
    
    message("Reading: ", file_path)
    
    
    #========================================
    # Read 10X h5
    #========================================
    
    data <- Read10X_h5(
      file_path,
      use.names = TRUE,
      unique.features = TRUE
    )
    
    
    #========================================
    # Create Seurat object
    #========================================
    
    sc_obj <- CreateSeuratObject(
      counts = data,
      project = original_id,
      assay = "RNA"
    )
    
    
    #========================================
    # Add metadata
    #========================================
    
    sc_obj$sample_name <- original_id
    sc_obj$sample_code <- sample_name
    sc_obj$group <- meta_row$Group
    
    
    #========================================
    # Store in list
    #========================================
    
    seurat_list[[original_id]] <- sc_obj
  }
  
  return(seurat_list)
}
  
 

#add mitocount ####
#this function adds percent.mt to each seurat object in the seurat_list
add_metadata_to_seurat_objects = function(obj_list, mito_pattern ){
  
  message("Adding percent.mt to metadata")
  
  obj_list = lapply(obj_list, function(seurat_obj){
    seurat_obj[["percent.mt"]] <- PercentageFeatureSet(
      seurat_obj,
      pattern = mito_pattern
    )
    seurat_obj
  })
  
  return(obj_list)
}


#doublet finder #####


#---------------doublet finder per sample------------------------
#Source statsquid





# Functions ===================================================================
#----------------------------------------------------------#
# run_doubletfinder_custom Source : https://biostatsquid.com/doubletfinder-tutorial/
#----------------------------------------------------------#
# run_doubletfinder_custom runs Doublet_Finder() and returns a dataframe with the cell IDs and a column with either 'Singlet' or 'Doublet'
run_doubletfinder_custom <- function(seu_sample_subset, multiplet_rate = NULL){
  
  
  if(is.null(multiplet_rate)){
    print('multiplet_rate not provided....... estimating multiplet rate from cells in dataset')
    
    # 10X multiplet rates table
    #https://rpubs.com/kenneditodd/doublet_finder_example
    multiplet_rates_10x <- data.frame('Multiplet_rate'= c(0.004, 0.008, 0.0160, 0.023, 0.031, 0.039, 0.046, 0.054, 0.061, 0.069, 0.076),
                                      'Loaded_cells' = c(800, 1600, 3200, 4800, 6400, 8000, 9600, 11200, 12800, 14400, 16000),
                                      'Recovered_cells' = c(500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000))
    
    print(multiplet_rates_10x)
    
    multiplet_rate <- multiplet_rates_10x %>% dplyr::filter(Recovered_cells < nrow(seu_sample_subset@meta.data)) %>% 
      dplyr::slice(which.max(Recovered_cells)) %>% # select the min threshold depending on your number of samples
      dplyr::select(Multiplet_rate) %>% as.numeric(as.character()) # get the expected multiplet rate for that number of recovered cells
    
    print(paste('Setting multiplet rate to', multiplet_rate))
  }
  
  # Pre-process seurat object with standard seurat workflow --- 
  sample <- NormalizeData(seu_sample_subset)
  sample <- FindVariableFeatures(sample)
  sample <- ScaleData(sample)
  sample <- RunPCA(sample, nfeatures.print = 10)
  
  # Find significant PCs
  stdv <- sample[["pca"]]@stdev
  percent_stdv <- (stdv/sum(stdv)) * 100
  cumulative <- cumsum(percent_stdv)
  co1 <- which(cumulative > 90 & percent_stdv < 5)[1] 
  co2 <- sort(which((percent_stdv[1:length(percent_stdv) - 1] - 
                       percent_stdv[2:length(percent_stdv)]) > 0.1), 
              decreasing = T)[1] + 1
  min_pc <- min(co1, co2)
  
  # Finish pre-processing with min_pc
  sample <- RunUMAP(sample, dims = 1:min_pc)
  sample <- FindNeighbors(object = sample, dims = 1:min_pc)              
  sample <- FindClusters(object = sample, resolution = 0.1)
  
  # pK identification (no ground-truth) 
  #introduces artificial doublets in varying props, merges with real data set and 
  # preprocesses the data + calculates the prop of artficial neighrest neighbours, 
  # provides a list of the proportion of artificial nearest neighbours for varying
  # combinations of the pN and pK
  sweep_list <- paramSweep(sample, PCs = 1:min_pc, sct = FALSE)   
  sweep_stats <- summarizeSweep(sweep_list)
  bcmvn <- find.pK(sweep_stats) # computes a metric to find the optimal pK value (max mean variance normalised by modality coefficient)
  # Optimal pK is the max of the bimodality coefficient (BCmvn) distribution
  optimal.pk <- bcmvn %>% 
    dplyr::filter(BCmetric == max(BCmetric)) %>%
    dplyr::select(pK)
  optimal.pk <- as.numeric(as.character(optimal.pk[[1]]))
  
  ## Homotypic doublet proportion estimate
  annotations <- sample@meta.data$seurat_clusters # use the clusters as the user-defined cell types
  homotypic.prop <- modelHomotypic(annotations) # get proportions of homotypic doublets
  
  nExp.poi <- round(multiplet_rate * nrow(sample@meta.data)) # multiply by number of cells to get the number of expected multiplets
  nExp.poi.adj <- round(nExp.poi * (1 - homotypic.prop)) # expected number of doublets
  
  # run DoubletFinder
  sample <- doubletFinder(seu = sample, 
                          PCs = 1:min_pc, 
                          pK = optimal.pk, # the neighborhood size used to compute the number of artificial nearest neighbours
                          nExp = nExp.poi.adj) # number of expected real doublets
  # change name of metadata column with Singlet/Doublet information
  colnames(sample@meta.data)[grepl('DF.classifications.*', colnames(sample@meta.data))] <- "doublet_finder"
  
  # extracr meta_data
  double_finder_res <- data.frame(
    barcode = rownames(sample@meta.data),
    doublet_finder = sample@meta.data$doublet_finder,
    orig.ident = sample@meta.data$orig.ident,
    row.names = NULL
  )
  
}




# plot fgsea 

plot_fgsea <- function(fgsea_df) {
  fgsea_df$sig   <- ifelse(fgsea_df$padj < 0.05, "sig", "nonsignificant")
  fgsea_df$color <- ifelse(
    fgsea_df$sig == "nonsignificant", "nonsignificant",
    ifelse(fgsea_df$NES > 0, "upregulated", "downregulated")
  )
  
  ggplot(fgsea_df, aes(x = NES, y = reorder(pathway, NES), fill = color)) +
    geom_col(width = 0.3) +
    scale_x_continuous(breaks = seq(0,10, by = 0.5))+
    geom_vline(xintercept = 0, color = "black") +
    scale_fill_manual(values = c(
      "nonsignificant" = "gray",
      "upregulated"    = "salmon",
      "downregulated"  = "steelblue"
    ),labels=NULL) +
    labs(x = "NES", y = NULL) +
    theme_minimal() +
    theme(
      legend.position="none",
      axis.text.x.bottom = element_text(size = 10,  family = "ArialMT"),
      panel.grid   = element_blank(),
      legend.title = element_blank(),
      axis.text.y  = element_text(size = 15,  family = "ArialMT"),
      strip.text   = element_text(size = 20, family = "ArialMT",colour = "black")
    )
}
?scale_fill_manual
# Read back all per-cell-type fgsea results and combine
fgsea_all <- list()

for (ct in names(de_results_list)) {
  f <- file.path(output_dir_mast,
                 paste0("Top_fgsea_", gsub(" ", "_", ct), ".xlsx"))
  if (file.exists(f)) {
    tbl_gsea <- readxl::read_xlsx(f)
    tbl_gsea$cluster <- ct
    fgsea_all[[ct]] <- tbl_gsea
  }
}



# plot fgsea by color 

plot_fgsea_by_cluster <- function(fgsea_df) {
  
  fgsea_df <- fgsea_df %>%
    mutate(
      sig = ifelse(padj < 0.05, "sig", "nonsignificant"),
      color = ifelse(sig == "nonsignificant", "nonsignificant",
                     ifelse(NES > 0, "upregulated", "downregulated"))
    ) %>%

        group_by(cluster) %>%
    arrange(NES, .by_group = TRUE) %>%
    ungroup() %>%

    
        mutate(
      pathway_label = factor(pathway)
    )
  
  ggplot(fgsea_df, aes(x = NES, y = pathway_label, fill = cluster)) +
    geom_col(width = 0.7) +
    geom_vline(xintercept = 0, color = "black", linewidth = 0.8) +
    scale_x_continuous(breaks = seq(-4, 4, by = 0.5)) +
    labs(x = "NES", y = NULL) +
    theme_minimal() +
    theme(
      legend.position = "right",
      axis.text.x     = element_text(size = 10, family = "ArialMT"),
      panel.grid      = element_blank(),
      axis.text.y     = element_text(size = 10, family = "ArialMT"),
      strip.text      = element_text(size = 14, family = "ArialMT", colour = "black"),
      legend.title    = element_blank(),
      legend.text     = element_text(size = 10, family = "ArialMT")
    )
}
