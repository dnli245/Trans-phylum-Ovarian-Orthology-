library(Seurat)
library(cowplot)
library(ggplot2)
library(SeuratDisk)
library(dplyr)
library(openxlsx)
set.seed(100)
# load and merge datasets
human_ov_wk7 <- Read10X('GSE143380_RAW/GSM4257926_Female_Week_7_gene_matrix')
human_ov_wk9 <- Read10X('GSE143380_RAW/GSM4257927_Female_Week_9_gene_matrix')
human_ov_wk10 <- Read10X('GSE143380_RAW/GSM4257928_Female_Week_10_gene_matrix')
human_ov_wk13 <- Read10X('GSE143380_RAW/GSM4257929_Female_Week_13_gene_matrix')
human_ov_wk16 <- Read10X('GSE143380_RAW/GSM4257930_Female_Week_16_gene_matrix')

human_ov_07 <- CreateSeuratObject(counts = human_ov_wk7, project = "hmov_7", min.cells = 3, min.features = 200)
human_ov_09 <- CreateSeuratObject(counts = human_ov_wk9, project = "hmov_9", min.cells = 3, min.features = 200)
human_ov_10 <- CreateSeuratObject(counts = human_ov_wk10, project = "hmov_10", min.cells = 3, min.features = 200)
human_ov_13 <- CreateSeuratObject(counts = human_ov_wk13, project = "hmov_13", min.cells = 3, min.features = 200)
human_ov_16 <- CreateSeuratObject(counts = human_ov_wk16, project = "hmov_16", min.cells = 3, min.features = 200)


hm_combined <- merge(human_ov_07, c(human_ov_09, human_ov_10, human_ov_13, human_ov_16),
                     add.cell.ids = c("wk7","wk9","wk10","wk13","wk16"),
                     project = "embryonic_hmov_combined")
human_ov <- hm_combined



ovary.1_1.data <- Read10X_h5("GSE118127_RAW/GSM3319032_sample_1-1_filtered_gene_bc_matrices_h5.h5", use.names = TRUE, unique.features = TRUE)
ovary.1_2.data <- Read10X_h5("GSE118127_RAW/GSM3319033_sample_1-2_filtered_gene_bc_matrices_h5.h5", use.names = TRUE, unique.features = TRUE)
ovary.1_3.data <- Read10X_h5("GSE118127_RAW/GSM3319034_sample_1-3_filtered_gene_bc_matrices_h5.h5", use.names = TRUE, unique.features = TRUE)
ovary.1_4.data <- Read10X_h5("GSE118127_RAW/GSM3319035_sample_1-4_filtered_gene_bc_matrices_h5.h5", use.names = TRUE, unique.features = TRUE)
ovary.1_5.data <- Read10X_h5("GSE118127_RAW/GSM3319036_sample_1-5_filtered_gene_bc_matrices_h5.h5", use.names = TRUE, unique.features = TRUE)
ovary.1_6.data <- Read10X_h5("GSE118127_RAW/GSM3319037_sample_1-6_filtered_gene_bc_matrices_h5.h5", use.names = TRUE, unique.features = TRUE)
ovary.1_7.data <- Read10X_h5("GSE118127_RAW/GSM3319038_sample_1-7_filtered_gene_bc_matrices_h5.h5", use.names = TRUE, unique.features = TRUE)
ovary.1_8.data <- Read10X_h5("GSE118127_RAW/GSM3319039_sample_1-8_filtered_gene_bc_matrices_h5.h5", use.names = TRUE, unique.features = TRUE)
ovary.3_5.data <- Read10X_h5("GSE118127_RAW/GSM3319046_sample_3-5_filtered_gene_bc_matrices_h5.h5", use.names = TRUE, unique.features = TRUE)
ovary.3_6.data <- Read10X_h5("GSE118127_RAW/GSM3319047_sample_3-6_filtered_gene_bc_matrices_h5.h5", use.names = TRUE, unique.features = TRUE)

ovary.3_14.data <- Read10X_h5("GSE118127_RAW/GSM3319041_sample_3-14_filtered_gene_bc_matrices_h5.h5", use.names = TRUE, unique.features = TRUE)
ovary.3_15.data <- Read10X_h5("GSE118127_RAW/GSM3319042_sample_3-15_filtered_gene_bc_matrices_h5.h5", use.names = TRUE, unique.features = TRUE)
ovary.3_16.data <- Read10X_h5("GSE118127_RAW/GSM3319043_sample_3-16_filtered_gene_bc_matrices_h5.h5", use.names = TRUE, unique.features = TRUE)
ovary.3_17.data <- Read10X_h5("GSE118127_RAW/GSM3319044_sample_3-17_filtered_gene_bc_matrices_h5.h5", use.names = TRUE, unique.features = TRUE)
ovary.3_18.data <- Read10X_h5("GSE118127_RAW/GSM3319045_sample_3-18_filtered_gene_bc_matrices_h5.h5", use.names = TRUE, unique.features = TRUE)

ovary.1_1.data <- CreateSeuratObject(counts = ovary.1_1.data, project = "1.1", min.cells = 3, min.features = 200)
ovary.1_2.data <- CreateSeuratObject(counts = ovary.1_2.data, project = "1.2", min.cells = 3, min.features = 200)
ovary.1_3.data <- CreateSeuratObject(counts = ovary.1_3.data, project = "1.3", min.cells = 3, min.features = 200)
ovary.1_4.data <- CreateSeuratObject(counts = ovary.1_4.data, project = "1.4", min.cells = 3, min.features = 200)
ovary.1_5.data <- CreateSeuratObject(counts = ovary.1_5.data, project = "1.5", min.cells = 3, min.features = 200)
ovary.1_6.data <- CreateSeuratObject(counts = ovary.1_6.data, project = "1.6", min.cells = 3, min.features = 200)
ovary.1_7.data <- CreateSeuratObject(counts = ovary.1_7.data, project = "1.7", min.cells = 3, min.features = 200)
ovary.1_8.data <- CreateSeuratObject(counts = ovary.1_8.data, project = "1.8", min.cells = 3, min.features = 200)
ovary.3_5.data <- CreateSeuratObject(counts = ovary.3_5.data, project = "3.5", min.cells = 3, min.features = 200)
ovary.3_6.data <- CreateSeuratObject(counts = ovary.3_6.data, project = "3.6", min.cells = 3, min.features = 200)
ovary.3_14.data <- CreateSeuratObject(counts = ovary.3_14.data, project = "3.14", min.cells = 3, min.features = 200)
ovary.3_15.data <- CreateSeuratObject(counts = ovary.3_15.data, project = "3.15", min.cells = 3, min.features = 200)
ovary.3_16.data <- CreateSeuratObject(counts = ovary.3_16.data, project = "3.16", min.cells = 3, min.features = 200)
ovary.3_17.data <- CreateSeuratObject(counts = ovary.3_17.data, project = "3.17", min.cells = 3, min.features = 200)
ovary.3_18.data <- CreateSeuratObject(counts = ovary.3_18.data, project = "3.18", min.cells = 3, min.features = 200)

hm_adult_ov <- merge(ovary.1_1.data, c(ovary.1_2.data, ovary.1_3.data, 
                                       ovary.1_4.data, ovary.1_5.data, ovary.1_6.data, 
                                       ovary.1_7.data, ovary.1_8.data, ovary.3_5.data,
                                       ovary.3_6.data, ovary.3_14.data, ovary.3_15.data,
                                       ovary.3_16.data, ovary.3_17.data, ovary.3_18.data),
                     add.cell.ids = c("1.1","1.2","1.3","1.4","1.5","1.6", "1.7", "1.8",
                                      "3.5", "3.6", "3.14", "3.15", "3.16","3.17", "3.18"),
                     project = "adult_hmov_combined")

hm_combined <- merge(human_ov, hm_adult_ov,project = "combined_hmov")

obj <- hm_combined
# qc
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
obj<- subset(obj, subset = nFeature_RNA > 200 & percent.mt < 20)


obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj)
#
obj <- FindNeighbors(obj, dims = 1:30, reduction = "pca")
obj <- FindClusters(obj, resolution = 0.5, cluster.name = "unintegrated_clusters")

obj <- RunUMAP(obj, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")

DimPlot(obj, reduction = "umap.unintegrated", group.by = c("orig.ident", "unintegrated_clusters"))

# batch correction
obj <- RunHarmony(
  object = obj,
  group.by.vars = "orig.ident",
  dims = 1:30
)


obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30) #original 1:30
obj <- FindClusters(obj, resolution = 2.5, cluster.name = "harmony_clusters") # original 0.5 


obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")
p1 <- DimPlot(
  obj,
  reduction = "umap.harmony",
  group.by = c("orig.ident", "unintegrated_clusters", "harmony_clusters"),
  combine = TRUE, label.size = 2
)

pdf(file = "hm_harmony_analysis/cluster_plot.pdf",width=25,height=10)
p1
dev.off()


# stage categorization
obj$stage <- ifelse(obj$orig.ident %in% c("1.1","1.2","1.3","1.4","1.5","1.6", "1.7", "1.8",
                                          "3.5", "3.6", "3.14", "3.15", "3.16","3.17", "3.18"), "hm_adult", obj$orig.ident)
stage_plot <- DimPlot(obj, reduction = "umap.harmony", group.by = "stage")
pdf(file = "hm_harmony_analysis/stage_plot.pdf",width=12, height = 8)
stage_plot
dev.off()

# separate plots per timepoint 
for(stage in unique(obj$stage)) {
  subset <- obj[,obj$stage == stage]
  subset_plot <- DimPlot(subset, reduction = "umap.harmony", label = TRUE, label.size = 6, repel = TRUE)
  pdf(paste0("hm_harmony_analysis/",stage, ".pdf"), width = 12, height = 8)
  print(subset_plot)
  dev.off()
}
obj <- JoinLayers(obj)


# cluster annotation

cell_types <- c(
  "Stroma",          # 0
  "Stroma",          # 1
  "Stroma",          # 2
  "Pre-granulosa", # 3
  "Common Progenitor",          # 4
  "Endothelial",   # 5
  "Stroma",   # 6
  "Smooth muscle",   # 7
  "Stroma",     # 8
  "Stroma",      # 9
  "Endothelial",      # 10
  "Pre-granulosa",          # 11
  "Stroma",   # 12
  "Immune",   # 13
  "Endothelial",   # 14
  "Stroma",      # 15
  "Pre-granulosa",          # 16
  "Epithelial",   # 17
  "Pre-granulosa",   # 18
  "Stroma",      # 19
  "Epithelial",          # 20
  "Granulosa",          # 21
  "Granulosa",          # 22
  "Stroma",          # 23
  "Germline",       # 24
  "Smooth muscle",          # 25
  "Epithelial",          # 26
  "Immune",          # 27
  "Germline",          # 28
  "Stroma",          # 29
  "Pre-granulosa",          # 30
  "Epithelial",          # 31
  "Stroma",          # 32
  "Epithelial",          # 33
  "Endothelial",          # 34
  "Stroma",        # 35
  "Germline",        # 36
  "Stroma",        # 37
  "Immune",   # 38
  "Granulosa",    # 39
  "Perivascular",     # 40
  "Smooth muscle",     # 41
  "Immune"     # 42
)

cell_types_diff <- c(
  "Stroma 1",          # 0
  "Stroma 2",          # 1
  "Stroma 3",          # 2
  "Pre-granulosa 1",   # 3
  "Common Progenitor", # 4
  "Endothelial 1",     # 5
  "Stroma 4",          # 6
  "Smooth muscle 1",   # 7
  "Stroma 5",          # 8
  "Stroma 6",          # 9
  "Endothelial 2",     # 10
  "Pre-granulosa 2",   # 11
  "Stroma 7",          # 12
  "Immune 1",          # 13
  "Endothelial 3",     # 14
  "Stroma 8",          # 15
  "Pre-granulosa 3",   # 16
  "Epithelial 1",      # 17
  "Pre-granulosa 4",   # 18
  "Stroma 9",          # 19
  "Epithelial 2",      # 20
  "Granulosa 1",       # 21
  "Granulosa 2",       # 22
  "Stroma 10",         # 23
  "Germline 1",        # 24
  "Smooth muscle 2",   # 25
  "Epithelial 3",      # 26
  "Immune 2",          # 27
  "Germline 2",        # 28
  "Stroma 11",         # 29
  "Pre-granulosa 5",   # 30
  "Epithelial 4",      # 31
  "Stroma 12",         # 32
  "Epithelial 5",      # 33
  "Endothelial 4",     # 34
  "Stroma 13",         # 35
  "Germline 3",        # 36
  "Stroma 14",         # 37
  "Immune 3",          # 38
  "Granulosa 3",       # 39
  "Perivascular",      # 40
  "Smooth muscle 3",   # 41
  "Immune 4"           # 4
)
obj$annotated_celltype <- cell_types[as.numeric(as.character(obj$seurat_clusters)) + 1]
obj$cell_idents <- cell_types_diff[as.numeric(as.character(obj$seurat_clusters)) + 1]

pdf(file = "hm_harmony_analysis/annotated_harmony_umap",width=15,height=10)
DimPlot(obj, group.by = "annotated_celltype", label = TRUE, reduction = "umap.harmony")
dev.off()

pdf(file = "hm_harmony_analysis/harmony_celltypes_split_umap.pdf",width=12,height=8)
DimPlot(obj, group.by = "cell_idents", label = TRUE, reduction = "umap.harmony",label.size = 4.8, repel = TRUE)
dev.off()

features = c("NANOG","POU5F1", # PGC
             "SMC1B", "ZGLP1", "STRA8","SYCP1",
             "ZP3", "FIGLA", # oocyte from adult ov 
             "NR2F1", "GATA2",  # common progenitor
             "PDGFRA", "TCF21", # stroma
             "WNT6", "FOXL2", #pre-granulosa
             "AMH", # granulosa
             "KRT19","UPK3B", # epithelial
             
             "CD3D", "KLRB1","CD14", "CD68", # immune young 
             "CD69", "ITGB2", "CD2", "CD3G", "CD8A", # adult immune more if needed
             
             "VWF", "PECAM1", "CDH5", # endothelial
             "ACTA2", "RGS5", # smooth muscle
             "MCAM", "RERGL", "TAGLN" # perivascular 
             
)

# rewrite idents
Idents(obj) <- obj$cell_idents
levs <- levels(obj)

# 2. Re-order them alphabetically, case-insensitive
levs <- levs[order(tolower(levs))]

# 3. Write the new order back
levels(obj) <- levs
p_dot <- DotPlot(
  obj,
  features = features,
  cols      = c("lightgrey", "blue"),
  col.min   = -2.5,
  col.max   =  2.5,
  dot.min   = 0,
  dot.scale = 6
)

pdf(file = "hm_harmony_analysis/dotplot_grouped.pdf",width=25,height=10)
p_dot
dev.off()


saveRDS(obj, "hm_harmony_analysis/harmony_object.rds")
obj <- readRDS("hm_harmony_analysis/harmony_object.rds")

# remove stromal cells
obj_no_stroma <- obj[,obj$annotated_celltype != "Stroma"]



obj_no_stroma[["RNA3"]] <- as(obj_no_stroma[["RNA"]], Class = "Assay")
DefaultAssay(obj_no_stroma) <- "RNA3"
obj_no_stroma[["RNA"]] <- NULL
obj_no_stroma <- RenameAssays(obj_no_stroma, RNA3 = "RNA")

# 2. Drop the scale.data slot
#    We can use DietSeurat to remove scale.data but keep counts & data:
obj_no_stroma <- DietSeurat(
  obj_no_stroma,
  assays = "RNA",
  counts = TRUE,
  data = TRUE,
  scale.data = FALSE
)

# 3. Save and convert
SaveH5Seurat(obj_no_stroma, "hm_harmony_combined_no_stroma_8.14.h5seurat", overwrite = TRUE)
Convert("hm_harmony_combined_no_stroma_8.14.h5seurat", dest = "h5ad", overwrite = TRUE)

# generate markers
Idents(obj) <- obj$annotated_celltype

generate_markers <- function(seurat, species, folder_name, conversion_table) {
  for(i in unique(Idents(seurat))) {
    id = as.character(i)
    print(i)
    marker_table <- FindMarkers(seurat, ident.1 = id)
    marker_table <- marker_table %>%
      mutate(pct_ratio= pct.1 / pct.2)
    marker_table <- cbind(
      rownames(marker_table),
      marker_table)
    
    colnames(marker_table)[1] <- "names"
    marker_table <- marker_table %>%
      arrange(desc(pct_ratio))
    write.xlsx(marker_table, file = paste0(folder_name, "/", species, "_cluster_", id, ".xlsx"))
    print(i)
    
  }
}
generate_markers(obj, "hm","hm_harmony_analysis/hm_seurat_markers_w_germ",hm_conversion_table)

# write metadata 

write.csv(
  obj@meta.data,
  file = "hm_harmony_analysis/hm_metadata.csv",
  row.names = TRUE
)


# Germ cell subset 
germ_obj <- obj[, obj$annotated_celltype == "Germline"]

# recluster
germ_obj <- FindNeighbors(germ_obj, reduction = "harmony", dims = 1:30) #original 1:30
germ_obj <- FindClusters(germ_obj, resolution = 2.2, cluster.name = "germ_clusters") # original 0.5 


germ_obj <- RunUMAP(germ_obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")


features = c("NANOG","POU5F1",# PGC
             "ZGLP1", "STRA8", "SYCP1",# meiotic germ cell
             "ZP3", "FIGLA", "DAZL", "MAEL", "PRDM1", "DPPA3", "IFITM3", "DDX4", "SOX17","NANOS3",
             "KLF4", "SALL4", "TFAP2C")
f = FeaturePlot(germ_obj, features, reduction = "umap.harmony")
d = DimPlot(germ_obj, reduction = "umap.harmony", label = TRUE)
e = DimPlot(germ_obj, reduction = "umap.harmony", group.by = 'orig.ident',label = TRUE)


germ_cell_types <- c(
  "PGC", # 0
  "PGC", # 1
  "Oocyte", # 2
  "PGC", # 3
  "PGC", # 4
  "Oocyte", # 5 
  "PGC", # 6
  "PGC", # 7 
  "PGC", #8 
  "Oocyte", # 9 
  "Oocyte", # 10 
  "PGC", # 11
  "PGC", # 12
  "PGC", # 13
  "PGC" # 14
)

germ_cell_idents <- c(
  "PGC 1",     # 0
  "PGC 2",     # 1
  "Oocyte 1",  # 2
  "PGC 3",     # 3
  "PGC 4",     # 4
  "Oocyte 2",  # 5
  "PGC 5",     # 6
  "PGC 6",     # 7
  "PGC 7",     # 8
  "Oocyte 3",  # 9
  "Oocyte 4",  # 10
  "PGC 8",     # 11
  "PGC 9",     # 12
  "PGC 10",    # 13
  "PGC 11"     # 14
)
germ_obj$annotated_celltype <- germ_cell_types[as.numeric(as.character(germ_obj$seurat_clusters)) + 1]
germ_obj$cell_idents <- germ_cell_idents[as.numeric(as.character(germ_obj$seurat_clusters)) + 1]

# rewrite idents
Idents(germ_obj) <- germ_obj$cell_idents
levs <- levels(germ_obj)

# 2. Re-order them alphabetically, case-insensitive
levs <- levs[order(tolower(levs))]

# 3. Write the new order back
levels(germ_obj) <- levs
germ_obj_none <- germ_obj
germ_obj <- germ_obj[,germ_obj$annotated_celltype != "None"]

germ_dot <- DotPlot(
  germ_obj,
  features = features,
  cols  = c("lightgrey", "blue"),
  col.min   = -2.5,
  col.max   =  2.5,
  dot.min   = 0,
  dot.scale = 6
)
pdf(file = "hm_germ_cell_experiment_8.14/hm_dotplot_grouped_11.4.pdf",width=15,height=10)
germ_dot
dev.off()

# plot number of cells in each cluster over time 

library(reshape2)
plot_cluster_trajectory <- function(obj, cluster_category, stage_category) {
  cats <- unique(obj@meta.data[[cluster_category]])
  stages <- unique(obj@meta.data[[stage_category]])
  
  df <- matrix(0, nrow = length(stages), ncol = length(cats))
  df <- as.data.frame(df)
  colnames(df) <- cats
  rownames(df) <- stages
  
  for (cat in cats) {
    cat_subset <- obj[, obj@meta.data[[cluster_category]] == cat]
    for (stage in stages) {
      count <- sum(cat_subset@meta.data[[stage_category]] == stage)
      df[stage, cat] <- count
    }
  }
  
  df$Stage <- rownames(df)
  df_long <- reshape2::melt(df, id.vars = "Stage", variable.name = "Cluster", value.name = "Count")
  
  df_long$Stage <- factor(df_long$Stage, levels = stages)
  
  p <- ggplot(df_long, aes(x = Stage, y = Count, color = Cluster, group = Cluster)) +
    geom_line(size = 1.2) +
    geom_point(size = 2) +
    theme_minimal(base_size = 14) +
    labs(x = "Stage", y = "Cell count", title = "Cluster trajectory over developmental stages") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5)
    )
  
  return(p)
}
p <- plot_cluster_trajectory(germ_obj, "cell_idents", "stage")
p

pdf(file = "hm_germ_cell_experiment_8.14/hm_umap_specific.pdf",width=5,height=5)
DimPlot(germ_obj, reduction = "umap.harmony", group.by = c("cell_idents"), label = TRUE, label.size = 5, repel = TRUE)
dev.off()

pdf(file = "hm_germ_cell_experiment_8.14/hm_umap_broad.pdf",width=5,height=5)
DimPlot(germ_obj, reduction = "umap.harmony", group.by = c("annotated_celltype"), label = TRUE, label.size = 6, repel = TRUE)
dev.off()

write.csv(
  germ_obj_none@meta.data,
  file = "germ_cell_experiment/hm_germ_metadata.csv",
  row.names = TRUE
)

saveRDS(germ_obj, "hm_germ_cell_experiment_8.14/germ_obj.rds")

germ_obj <- readRDS("hm_germ_cell_experiment_8.14/germ_obj.rds")

germ_obj[["RNA3"]] <- as(germ_obj[["RNA"]], Class = "Assay")
DefaultAssay(germ_obj) <- "RNA3"
germ_obj[["RNA"]] <- NULL
germ_obj <- RenameAssays(germ_obj, RNA3 = "RNA")

# 2. Drop the scale.data slot
#    We can use DietSeurat to remove scale.data but keep counts & data:
germ_obj <- DietSeurat(
  germ_obj,
  assays = "RNA",
  counts = TRUE,
  data = TRUE,
  scale.data = FALSE
)

# 3. Save and convert
SaveH5Seurat(germ_obj, "hm_germ_8.14.h5seurat", overwrite = TRUE)
Convert("hm_germ_8.14.h5seurat", dest = "h5ad", overwrite = TRUE)


#ratio experiment:

hm7 <- germ_obj[, germ_obj$orig.ident == "hmov_7"]
hm7_oocyte <- hm7[, hm7$germ_celltype == "Oocyte"]
hm9 <- germ_obj[, germ_obj$orig.ident == "hmov_9"]
hm9_oocyte <- hm9[, hm9$germ_celltype == "Oocyte"]
hm10 <- germ_obj[, germ_obj$orig.ident == "hmov_10"]
hm10_oocyte <- hm10[, hm10$germ_celltype == "Oocyte"]

# insert new germ annotations back into original object
hm_new <- obj
hm_new <- AddMetaData(hm_new, metadata = germ_obj@meta.data)

hm_new$annotated_celltype[hm_new$annotated_celltype == "Germline"] <- 
  hm_new$germ_celltype[hm_new$annotated_celltype == "Germline"]
hm_new$annotated_celltype <- as.character(hm_new$annotated_celltype)
hm_new <- hm_new[,!is.na(hm_new$annotated_celltype)]

Idents(hm_new) <- hm_new@meta.data$annotated_celltype

hm_conversion_table <- read.csv('hm_conversion_table.txt',sep='\t', header=FALSE)

# generate new markers for new germ annotatins in original object
generate_markers(hm_new,"hm","hm_harmony_analysis/hm_seurat_markers_w_germ", hm_conversion_table)
generate_markers(germ_obj,"hm","hm_harmony_analysis/hm_germ_markers", hm_conversion_table)




