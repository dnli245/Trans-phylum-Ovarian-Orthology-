
library(reticulate)
library(Seurat)
library(anndata)
library(openxlsx)
library(presto)
library(dplyr)
library(ggplot2)
library(patchwork)
library(harmony)
set.seed(100)
sp_em_integ <- readRDS("SpInteg.rds") # embryonic dataset

sp_em_integ <- UpdateSeuratObject(sp_em_integ)
sp_data <- readRDS("sp_seurat_object.rds") # sp adult ovary dataset

# merge
sp_merge <- merge(sp_em_integ, sp_data, project = "combined_spov")

# QC
sp_merge<- subset(sp_merge, subset = nFeature_RNA > 200)

sp_merge <- NormalizeData(sp_merge, normalization.method = "LogNormalize", scale.factor = 10000)

sp_merge <- FindVariableFeatures(sp_merge, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(sp_merge)
sp_merge <- ScaleData(sp_merge, features = all.genes)

sp_merge <- RunPCA(sp_merge, features = VariableFeatures(object = sp_merge))

ElbowPlot(sp_merge)



sp_merge <- FindNeighbors(sp_merge, dims = 1:20)
sp_merge <- FindClusters(sp_merge, resolution = 1.7)

# batch correction
obj <- IntegrateLayers(
  object = sp_merge, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE,
  group.by = "orig.ident"
)

obj <- RunHarmony(
  object = obj,
  group.by.vars = "orig.ident",
  dims = 1:30
)


obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30) #original 1:30
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")
obj <- FindClusters(obj, resolution = 2.5, cluster.name = "harmony_clusters")
obj_adult <- subset(obj, subset = orig.ident == "0")
obj_em <- subset(obj, subset = orig.ident != "0")
obj@meta.data$"orig_annotation" <- c(obj_em$integrated_snn_res.0.5, obj_adult$cell_idents)

p1 <- DimPlot(
  obj,
  reduction = "umap.harmony",
  group.by = c("orig_annotation"),
  combine = TRUE, label.size = 2, label = TRUE
)

# insert specific germ cell annotations from germ cell subset analysis
germ_obj <- readRDS("germ_obj.rds")

obj_new <- obj

# Make sure we add metadata aligned by cell names only
meta <- germ_obj@meta.data
# keep only rows that exist in obj_new and preserve order
meta <- meta[colnames(obj_new), , drop = FALSE]
obj_new <- AddMetaData(obj_new, metadata = meta)

# make character vectors (avoid factor assignment issues)
obj_new$annotated_celltype <- as.character(obj_new$annotated_celltype)
obj_new$germ_celltype      <- as.character(obj_new$germ_celltype)

# build a clean index: no NAs allowed in the subscript
idx <- which(!is.na(obj_new$annotated_celltype) &
               obj_new$annotated_celltype == "Germline" &
               !is.na(obj_new$germ_celltype))

# assign only on valid positions
obj_new$orig_annotation[idx] <- obj_new$germ_celltype[idx]


p2 <- DimPlot(
  obj_new,
  reduction = "umap.harmony",
  group.by = c("orig_annotation"),
  combine = TRUE, label.size = 2, label = TRUE
)

obj_new$annotation_named <- as.character(obj_new$orig_annotation)

# Define embryonic cluster labels
embryo_labels <- c(
  "embr_Aboral ectoderm",         
  "embr_Oral ectoderm",           
  "embr_Ciliated cells",          
  "embr_Neural",                  
  "embr_Oral ectoderm",           
  "embr_Aboral ectoderm / neural",
  "embr_Endoderm",                
  "embr_Ciliated cells",          
  "embr_Endoderm",                
  "embr_Neural",                  
  "embr_Histone-enriched (early cells)",
  "embr_Pigment cells",           
  "embr_Oral ectoderm",           
  "embr_Oral ectoderm",           
  "embr_Endoderm",                
  "embr_Protease-enriched (late gastrula)",
  "embr_Skeleton",                
  "embr_Neural",                  
  "embr_Neural",                  
  "embr_Skeleton",                
  "embr_PGC",                
  "embr_germline?"
)

# Replace only numeric cluster entries (0–21)
num_idx <- suppressWarnings(!is.na(as.numeric(obj_new$orig_annotation)))
obj_new$annotation_named[num_idx] <- 
  embryo_labels[as.numeric(obj_new$orig_annotation[num_idx]) + 1]

# For any remaining textual labels (adult data), add "adult_" prefix
text_idx <- which(!num_idx)
obj_new$annotation_named[text_idx] <- paste0("adult_", obj_new$annotation_named[text_idx])

# Check results
table(obj_new$annotation_named, useNA = "ifany")

# cell-type annotated umap
p3 <- DimPlot(
  obj_new,
  reduction = "umap.harmony",
  group.by = c("annotation_named"), label.size = 3, label = TRUE, repel = TRUE
)


pdf(file = "sp_experiment/sp_integrated_umap.pdf",width=10,height=10)
p3
dev.off()


adult_and_em_germ <- subset(obj_new, subset = annotation_named == "embr_PGC" | orig.ident == "0" | annotation_named == "embr_germline?")

# germ cell plot
p4 <- DimPlot(
  adult_and_em_germ,
  reduction = "umap.harmony",
  group.by = c("annotation_named"), label.size = 3, label = TRUE, repel = TRUE
)

pdf(file = "sp_experiment/sp_integrated_just_embr_germ.pdf",width=10,height=10)
p4
dev.off()

# identify markers in overlapping PGC/germ cell cluster
# we identified overlapping PGC cluster as # 38 
generate_markers <- function(seurat, species, folder_name, conversion_table) {
  i = '38'
  id = as.character(i)
  print(i)
  marker_table <- FindMarkers(seurat, ident.1 = id)
  marker_table <- marker_table %>%
    mutate(pct_ratio= pct.1 / pct.2)
  marker_table <- cbind(
    rownames(marker_table),
    marker_table)
  
  colnames(marker_table)[1] <- "names"
  marker_table <- merge(marker_table, conversion_table)
  marker_table <- marker_table %>%
    arrange(desc(pct_ratio))
  write.xlsx(marker_table, file = paste0(folder_name, "/", species, "_cluster_", id, ".xlsx"))
  print(i)
  
  
}

generate_markers_conserved <- function(seurat, species, folder_name, conversion_table,
                                       cluster_id = 38) {
  # Split cells into two groups
  seurat$group <- ifelse(seurat$orig.ident == "0", "group0", "groupOther")
  
  # Split object by group
  split_objs <- SplitObject(seurat, split.by = "group")
  
  # Run FindMarkers on each group
  m0 <- FindMarkers(split_objs$group0, ident.1 = cluster_id)
  mO <- FindMarkers(split_objs$groupOther, ident.1 = cluster_id)
  
  # Merge results by gene
  merged <- merge(m0, mO, by = "row.names", suffixes = c(".0", ".other"))
  colnames(merged)[1] <- "gene"
  
  # Keep genes that are up or down in both
  merged <- merged %>%
    filter((avg_log2FC.0 > 0 & avg_log2FC.other > 0) |
             (avg_log2FC.0 < 0 & avg_log2FC.other < 0))
  
  # Add conversion info
  merged <- merge(merged, conversion_table, by = "gene", all.x = TRUE)
  
  # Sort by stronger average fold change
  merged <- merged %>%
    mutate(direction = ifelse(avg_log2FC.0 > 0, "Up", "Down")) %>%
    arrange(desc(abs((avg_log2FC.0 + avg_log2FC.other) / 2)))
  
  # Save to Excel
  write.xlsx(merged, file = paste0(folder_name, "/", species, "_cluster_", cluster_id, "_consistent.xlsx"))
  
  message("✅ Saved consistent markers for cluster ", cluster_id)
}


Idents(obj_new) <- obj_new@meta.data$harmony_clusters
obj_new <- JoinLayers(obj_new)

generate_markers(obj_new,"sp","sp_experiment/sp_seurat_integrated_harmony_markers", sp_conversion_table)
sp_conversion_table <- read.csv('sp_conversion_table.txt',sep='\t', header=FALSE)
names(sp_conversion_table)[3] <- "gene"
generate_markers_conserved(obj_new, "sp", "sp_experiment/sp_seurat_integrated_harmony_markers", sp_conversion_table)

# feature plots of markers in each timepoint
adult_obj <- subset(obj_new, subset = orig.ident == "0")
embr_obj <- subset(obj_new, subset =  orig.ident != "0")

plot_dual_feature <- function(feature) {
  folder <- paste0("sp_experiment/integrated_features/",feature)
  dir.create(folder)
  p_adult <- DimPlot(
    adult_obj,
    reduction = "umap.harmony",
    group.by = c("annotation_named"), label.size = 3, label = TRUE, repel = TRUE)
  p_embr <- DimPlot(
    embr_obj,
    reduction = "umap.harmony",
    group.by = c("annotation_named"), label.size = 3, label = TRUE, repel = TRUE)
  p_adult_feature <- FeaturePlot(adult_obj, c(feature))
  p_embr_feature <- FeaturePlot(embr_obj, c(feature))
  
  pdf(file = paste0(folder, "/", "adult_umap.pdf"),width=10,height=10)
  print(p_adult)
  dev.off()
  
  pdf(file = paste0(folder, "/", "embr_umap.pdf"),width=10,height=10)
  print(p_embr)
  dev.off()
  
  pdf(file = paste0(folder, "/", "adult_feature.pdf"),width=10,height=10)
  print(p_adult_feature)
  dev.off()
  
  pdf(file = paste0(folder, "/", "embr_feature.pdf"),width=10,height=10)
  print(p_embr_feature)
  dev.off()
  
  
  
}


plot_feature_integrated <- function(gene_list, seurat, folder) {
  for(gene in gene_list) {
    p <- FeaturePlot(seurat, gene)
    
    pdf(file = paste0(folder, "/", gene, ".pdf"),width=10,height=10)
    print(p)
    
    dev.off()
    
  }
}
gene_list <- c("Nanos2", "LOC100893220","LOC587419","LOC594212", "LOC100893206")
plot_feature_integrated(gene_list, obj_new, "sp_experiment/sp_integrated_feature_plots")
