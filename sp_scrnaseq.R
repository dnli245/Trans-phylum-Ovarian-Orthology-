library(reticulate)
library(Seurat)
library(anndata)
library(openxlsx)
library(presto)
library(dplyr)
library(ggplot2)
library(patchwork)


set.seed(100)
# set python path
use_python("your-path", required = TRUE)

# load sp adult ovary h5ad file as seurat object
sp_data <- read_h5ad("ourspecies/h5ad/sp_object.h5ad") 

sp_data <- CreateSeuratObject(counts = t(as.matrix(sp_data$X)), meta.data = sp_data$obs)

# QC
sp_data<- subset(sp_data, subset = nFeature_RNA > 200)

sp_data <- NormalizeData(sp_data, normalization.method = "LogNormalize", scale.factor = 10000)

sp_data <- FindVariableFeatures(sp_data, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(sp_data)
sp_data <- ScaleData(sp_data, features = all.genes)

sp_data <- RunPCA(sp_data, features = VariableFeatures(object = sp_data))

ElbowPlot(sp_data)

sp_data <- FindNeighbors(sp_data, dims = 1:20)
sp_data <- FindClusters(sp_data, resolution = 1.7)

sp_data <- RunUMAP(sp_data, dims = 1:20)

DimPlot(sp_data, reduction = "umap", label = TRUE)

# cluster annotation
features1 = c(
  "LOC578206", "LOC592537", 
  "vasa",
  "Nanos2", 
  "LOC577374",
  "LOC373434",
  "LOC586414" ,
  "LOC100892008",
  "LOC764071", 
  "LOC115924179",
  "LOC579828",
  "LOC105439489",
  "LOC105444039",
  "LOC100889874", "LOC100890503", "LOC581509", "LOC589794",
  "LOC100892701", "LOC584177", "LOC100893584",
  "LOC373488", "LOC582391","LOC582578","LOC587715",
  "LOC583660",
  "LOC100889113","LOC105442321","185/333",
  "LOC584097", "LOC586465",
  "LOC590620", "LOC587668", "LOC581904",
  "LOC100893082", "LOC581500", "LOC105441547", "LOC105443458","LOC576786",
  "LOC576448", "LOC105445992",
  "LOC579559",
  "LOC751846"
  
)

nicenames1 <- c(
  "LOC578206" = "neurotrypsin", "LOC592537" = "macoilin",
  
  "vasa",
  "Nanos2", 
  "LOC577374" = "Lefty",
  "LOC373434" = "seawi",
  
  "LOC764071" = "Synaptomenal complex 3",
  "LOC586414" = "meiotic nuclear division protein 1 homolog",
  "LOC100892008" = "meiosis-specific with OB domain-containing protein",
  "LOC115924179"="REC8 homolog",
  "LOC579828"="SPO11-like",
  "LOC105439489" = "MEI4-like",
  "LOC105444039" = "meiosis regulator and mRNA stability factor 1-like",
  "LOC100889874"= "steroid-17", "LOC100890503"="20 lyase", "LOC581509" = "SPARC", "LOC589794" = "collagen",
  "LOC100892701" = "oncoprotein-induced transcript 3 protein", "LOC584177" = "papilin isoform X1", "LOC100893584" = "papilin",
  "LOC373488"= "univin", "LOC582391"="arylsulfatase","LOC582578" = "frizzled-5","LOC587715" ="cilia and flagella 44",
  "LOC583660"="Vtgn1",
  "LOC100889113" = "echinoidin", "LOC105442321" = "toll-like R3", "185/333" = "185/333",
  "LOC584097" = "deleted in malignant brain tumoers 1 protein", "LOC586465" = "bactericidal permeability-increasing protein",
  "LOC590620" = "kappa-B", "LOC587668" = "titin", "LOC581904" = "neurexin",
  "LOC100893082"= "actin", "LOC581500" = "actin", "LOC105441547" = "voltage-dependent calcium channel", "LOC105443458" = "myoneurin-like", "LOC576786" = "myosin light chain kinase, smooth muscle",
  
  
  "LOC576448" = "muscle-specific protein 20", "LOC105445992" = "myosin light chain kinase",
  "LOC579559" = "60S ribosomal protein L12",
  "LOC751846" = "FOXL2"
)


# build dot plot 


cell_types <- c(
  "Uncharacterized", "Uncharacterized", "Uncharacterized", "Uncharacterized", "Epithelial",
  "Follicle", "Uncharacterized", "Neuronal", "Immune", "Epithelial",
  "Epithelial", "Epithelial", "Germline", "Germline", "Follicle", "Germline",
  "Germline", "Epithelial", "Epithelial", "Muscle",
  "Immune", "Epithelial", "Germline"
)
cell_types_diff <- c(
  "Uncharacterized", "Uncharacterized", "Uncharacterized", "Uncharacterized", "Epithelial 1",
  "Follicle 1", "Uncharacterized", "Neuronal", "Immune 1", "Epithelial 2",
  "Epithelial 3", "Epithelial 4", "Germline 1", "Germline 2", "Follicle 2", "Germline 3",
  "Germline 4", "Epithelial 4", "Epithelial 5", "Muscle 1",
  "Immune 2", "Epithelial 6", "Germline 5"
)
# assign new annotations
sp_data$annotated_celltype <- cell_types[as.numeric(as.character(sp_data$seurat_clusters)) + 1]
sp_data$cell_idents <- cell_types_diff[as.numeric(as.character(sp_data$seurat_clusters)) + 1]

# save umap plots
annotated_umap <- DimPlot(sp_data, label = TRUE, repel = TRUE, label.size = 5.5)
pdf(file = "sp_experiment/celltypes_split_umap",width=9,height=8)
annotated_umap
dev.off()

sp_data$orig_seurat_clusters <- Idents(sp_data) # save original clusters
sp_data$cell_idents <- factor(sp_data$cell_idents)
Idents(sp_data) <- "cell_idents"

saveRDS(sp_data, "sp_seurat_object.rds")


sp_data <- readRDS("sp_seurat_object.rds")

# dotplots

p2 <- DotPlot(
  sp_data,
  features = features1,
  cols      = c("lightgrey", "blue"),
  col.min   = -2.5,
  col.max   =  2.5,
  dot.min   = 0,
  dot.scale = 6
)
p2 <- p2 + 
  scale_x_discrete(labels = nicenames1[features1]) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
pdf(file = "sp_experiment/dotplot_grouped.pdf",width=20,height=15)
p2
dev.off()

# marker gene identification 
Idents(sp_data) <- sp_data$annotated_celltype # or sp_data$cell_idents
sp_conversion_table <- read.csv('sp_conversion_table.txt',sep='\t', header=FALSE)
names(sp_conversion_table)[3] <- "names"
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
    marker_table <- merge(marker_table, conversion_table)
    marker_table <- marker_table %>%
      arrange(desc(pct_ratio))
    write.xlsx(marker_table, file = paste0(folder_name, "/", species, "_cluster_", id, ".xlsx"))
    print(i)
    
  }
}

generate_markers(sp_data, "sp","sp_experiment/sp_seurat_markers_12.8",sp_conversion_table)

# write metadata
write.csv(
  sp_data@meta.data,
  file = "sp_experiment/sp_metadata.csv",
  row.names = TRUE
)

sp_data <- sp_data[,sp_data$annotated_celltype != "Uncharacterized"]

sp_data[["RNA3"]] <- as(sp_data[["RNA"]], Class = "Assay")
DefaultAssay(sp_data) <- "RNA3"
sp_data[["RNA"]] <- NULL
sp_data <- RenameAssays(sp_data, RNA3 = "RNA")

# 2. Drop the scale.data slot
#    We can use DietSeurat to remove scale.data but keep counts & data:
sp_data <- DietSeurat(
  sp_data,
  assays = "RNA",
  counts = TRUE,
  data = TRUE,
  scale.data = FALSE
)

# 3. Save and convert
SaveH5Seurat(sp_data, "sp_no_ribosomal_w_germ_8.15.h5seurat", overwrite = TRUE)
Convert("sp_no_ribosomal_w_germ_8.15.h5seurat", dest = "h5ad", overwrite = TRUE)


set.seed(100) # 12.8


# germ object subset 

germ_obj <- sp_data[,sp_data$annotated_celltype == "Germline"]
germ_obj <- FindNeighbors(germ_obj, dims = 1:20) #original 1:20
germ_obj <- FindClusters(germ_obj, resolution = 0.3, cluster.name = "germ_clusters") # original 1.5


germ_obj <- RunUMAP(germ_obj, dims = 1:30)


pdf(file = "germ_cell_experiment/sp_umap_specific.pdf", width=5, height=5)
DimPlot(germ_obj, label = TRUE, label.size = 4)
dev.off()


germ_cell_types <- c(
  "germline 1",
  "germline 2",
  "germline 3",
  "germline 4",
  "germline 5"
)
germ_obj$cell_idents <- germ_cell_types[as.numeric(as.character(germ_obj$seurat_clusters)) + 1]
germ_obj <- germ_obj[,germ_obj$cell_idents != "none"]
Idents(germ_obj) <- germ_obj$cell_idents
levs <- levels(germ_obj)
levs <- levs[order(tolower(levs))]
levels(germ_obj) <- levs


pdf(file = "germ_cell_experiment/sp_umap_specific_res_0.3_4_clusters.pdf", width=5, height=5)
DimPlot(germ_obj, label = TRUE, label.size = 4)
dev.off()


pdf(file = "germ_cell_experiment/sp_umap_broad.pdf",width=4,height=4)
DimPlot(germ_obj, group.by = "annotated_celltype")
dev.off()

generate_markers(germ_obj,"sp","germ_cell_experiment/sp_germ_cell_markers", sp_conversion_table)


# cluster annotation 
germ_features <- c(
  "vasa",
  "Nanos2", 
  "LOC577374",
  "LOC593520",
  "LOC373434",
  "LOC586414" ,
  "LOC100892008",
  "LOC764071", # Synaptomenal complex 3,
  "LOC115924179",
  "LOC579828",
  "LOC105439489",
  "LOC105444039",
  "LOC583022",
  "LOC100890845"
  
)
germ_nice <- c(
  "vasa",
  "Nanos2", 
  "LOC577374" = "Lefty",
  "LOC373434" = "seawi",
  "LOC593520" = "sox4",
  "LOC764071" = "Synaptomenal complex 3",
  "LOC583022" = "prdm9",
  "LOC100890845" = "prdm9 2",
  "LOC586414" = "meiotic nuclear division protein 1 homolog",
  "LOC100892008" = "meiosis-specific with OB domain-containing protein",
  "LOC115924179"="REC8 homolog",
  "LOC579828"="SPO11-like",
  "LOC105439489" = "MEI4-like",
  "LOC105444039" = "meiosis regulator and mRNA stability factor 1-like"
)

germ_dot <- DotPlot(
  germ_obj,
  features = germ_features,
  cols      = c("lightgrey", "blue"),
  col.min   = -2.5,
  col.max   =  2.5,
  dot.min   = 0,
  dot.scale = 6
)
germ_dot <- germ_dot + 
  scale_x_discrete(labels = germ_nice[germ_features]) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
pdf(file = "germ_cell_experiment/sp_germ_dotplot_8.15.pdf",width=15,height=10)
germ_dot
dev.off()


# dotplot of genes chosen for in situ hybridization
germ_dot_2 <- DotPlot(
  germ_obj,
  features = c("LOC577601","LOC589177", "LOC590208", "LOC115929968"),
  cols      = c("lightgrey", "blue"),
  col.min   = -2.5,
  col.max   =  2.5,
  dot.min   = 0,
  dot.scale = 6
)

pdf(file = "germ_cell_experiment/sp_germ_in_situ_dotplot_11_17.pdf",width=15,height=10)
germ_dot_2
dev.off()

germ_cell_types <- c(
  "germline 1",
  "germline 2",
  "germline 3",
  "germline 4",
  "germline 5",
  "germline 6",
  "germline 7",
  "germline 8",
  "germline 9",
  "germline 10",
  "germline 11"
)

germ_cell_types <- c(
  "germline 1",
  "germline 2",
  "germline 3",
  "germline 4",
  "germline 5",
  "none",
  "germline 6",
  "germline 7",
  "germline 8",
  "germline 9"
)

germ_obj$cell_idents <- germ_cell_types[as.numeric(as.character(germ_obj$seurat_clusters)) + 1]
germ_obj <- germ_obj[,germ_obj$cell_idents != "none"]
Idents(germ_obj) <- germ_obj$cell_idents

levs <- levels(germ_obj)
# 2. Re-order them alphabetically, case-insensitive
levs <- levs[order(tolower(levs))]

# 3. Write the new order back
levels(germ_obj) <- levs

write.csv(
  germ_obj@meta.data,
  file = "germ_cell_experiment/sp_germ_metadata.csv",
  row.names = TRUE
)

germ_umap = DimPlot(germ_obj, label = TRUE, label.size = 4) 
pdf(file = "germ_cell_experiment_real_germ_names/sp_umap_specific.pdf",width = 5,height=5)
germ_umap
dev.off()

saveRDS(germ_obj, "germ_obj.rds")
germ_obj <- readRDS("germ_obj.rds")

# germ cell marker exploration
features_germ_stem <- c(
  "LOC115926787", "LOC100889926", "LOC115919015",
  "LOC100889435", "LOC585801", "LOC115919242",
  "LOC592055", "LOC115917881", "LOC105445841",
  "LOC115925588"
)

nicenames_germ_stem <- c(
  "LOC115926787" = "musashi homolog 2",
  "LOC100889926" = "ubiquitin 3 mind bomb",
  "LOC115919015" = "activin",
  "LOC100889435" = "eva1",
  "LOC585801"    = "nova1",
  "LOC115919242" = "protogenin-a like",
  "LOC592055"    = "kiss-1",
  "LOC115917881" = "prospero",
  "LOC105445841" = "centriolin",
  "LOC115925588" = "harb1"
)


sp_germ_stem <- DotPlot(
  germ_obj,
  features = features_germ_stem,
  cols      = c("lightgrey", "blue"),
  col.min   = -2.5,
  col.max   =  2.5,
  dot.min   = 0,
  dot.scale = 6
)
sp_germ_stem <- sp_germ_stem + 
  scale_x_discrete(labels = nicenames_germ_stem[features_germ_stem]) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
pdf(file = "germ_cell_experiment/sp_stem_germ_dotplot_11.7.pdf",width=15,height=10)
sp_germ_stem
dev.off()

pdf(file = "germ_cell_experiment/sp_stem_germ_Featureplot_11.7.pdf",width=15,height=10)
FeaturePlot(germ_obj, features = features_germ_stem)
dev.off()

features_germ_oocyte <- c(
  "op",
  "mgb",
  # "LOC577801",
  #  "LOC575381",
  "LOC586390",
  "LOC115924179"
  #  "LOC105444710",
  # "LOC586403",
  #"LOC105442193",
  #"LOC579828"
  
)
nicenames_germ_oocyte <- c(
  "op" = "op",
  "mgb" = "rendezvin",
  #"LOC577801" = "Hif1a",
  #"LOC575381" = "Msh4",
  "LOC586390"= "Msh5",
  "LOC115924179" = "REC8"
  #"LOC105444710" = "sycp1",
  #"LOC586403" = "dmc1",
  #"LOC105442193" = "Meiob",
  #"LOC579828" = "spo11"
  
)

# oocyte markers
sp_germ_oocyte <- DotPlot(
  germ_obj,
  features = features_germ_oocyte,
  cols      = c("lightgrey", "blue"),
  col.min   = -2.5,
  col.max   =  2.5,
  dot.min   = 0,
  dot.scale = 6
)
sp_germ_oocyte <- sp_germ_oocyte + 
  scale_x_discrete(labels = nicenames_germ_oocyte[features_germ_oocyte]) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
pdf(file = "germ_cell_experiment/sp_oocyte_germ_dotplot_11.10.pdf",width=15,height=10)
sp_germ_oocyte
dev.off()

pdf(file = "germ_cell_experiment/sp_oocyte_germ_Featureplot_11.10.pdf",width=15,height=10)
FeaturePlot(germ_obj, features = features_germ_oocyte)
dev.off()

# save as h5ad
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
SaveH5Seurat(germ_obj, "sp_germ_8.15.h5seurat", overwrite = TRUE)
Convert("sp_germ_8.15.h5seurat", dest = "h5ad", overwrite = TRUE)
# insert germ annotations back to generate markers for each germline cluster
sp_new <- sp_data
sp_new <- AddMetaData(sp_new, metadata = germ_obj@meta.data)

sp_new$annotated_celltype[sp_new$annotated_celltype == "Germline"] <- 
  sp_new$germ_celltype[sp_new$annotated_celltype == "Germline"]
sp_new$annotated_celltype <- as.character(sp_new$annotated_celltype)
sp_new <- sp_new[,!is.na(sp_new$annotated_celltype)]

Idents(sp_new) <- sp_new@meta.data$cell_idents
generate_markers(sp_new,"sp","sp_experiment/sp_seurat_markers_w_germ_10.24", sp_conversion_table)

generate_markers <- function(seurat, species, folder_name, conversion_table) {
  id = "germ 7,8,9"
  marker_table <- FindMarkers(seurat, ident.1 = c("germ 7","germ 8", "germ 10"))
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
}
generate_markers(germ_obj,"sp","germ_cell_experiment/sp_germ_cell_markers", sp_conversion_table)

Features <- c("LOC115926787", # musashi homolog 2
              "LOC100889926",   # ubiquitin 3 mind bomb 
              "LOC115919015", # activin 
              "LOC100889435", # eva1
              "LOC585801", # nova1
              "LOC115919242", # protogenin-a like
              "LOC592055",# kiss-1
              "LOC115917881", # prospero
              "LOC105445841", # centriolin
              "LOC115925588" # harb1
)



# Create DotPlot
p2 <- DotPlot(
  sp_data,
  features  = features_wnt,
  cols      = c("lightgrey", "blue"),
  col.min   = -2.5,
  col.max   =  2.5,
  dot.min   = 0,
  dot.scale = 6
)

in_situ_list <- c("LOC373241", "LOC590208", "LOC580357", "LOC589177", 
                  "LOC577601", "LOC115919533", "LOC593520","LOC577374")

# Feature plot for FISH genes
feature_plot_all <- function(gene_list, full_data, germ_subset) {
  for(gene in gene_list) {
    f1 <- FeaturePlot(full_data, "vasa")
    f2 <- FeaturePlot(full_data, gene)
    g1 <- FeaturePlot(germ_subset, "vasa")
    g2 <- FeaturePlot(germ_subset, gene)
    
    
    pdf(paste0("sp_experiment/FISH_genes/", gene, ".pdf"), width = 40, height = 30)
    print(annotated_umap + f1 +f2 + germ_umap + g1 + g2)
    dev.off()
  }
}
feature_plot_all(in_situ_list, sp_data, germ_obj)



