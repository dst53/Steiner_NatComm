#' ---
#' title: "stRNA-seq analysis for "Identification of a gene expression 
#' signature of vascular invasion and recurrence in stage I lung adenocarcinoma 
#' via bulk and spatial transcriptomics"
#' output: github_document
#' ---

# expected run time is about 1-2 hours on a machine with 12 CPUs and 16gb RAM per core

#######
# SETUP
#######

# load R packages (all installed under R 4.2.1)

library(plyr)
library(dplyr)
library(Seurat)
library(future)
library(ggplot2)
library(RColorBrewer)
library(grid)
library(readbitmap)
library(hdf5r)
library(readxl)
library(SingleCellExperiment)
library(rmarkdown)
library(dittoSeq)
library(harmony)
library(patchwork)
library(scater)
library(scran)
library(proxy)
library(pheatmap)
library(GSVA)
library(msigdbr)
library(ggpubr)
library(mclust)
library(viridis)
library(reshape2)
library(ape)
library(stringr)
library(data.table)
library(cowplot)
library(edgeR)
library(sva)
library(BiocParallel)
library(HDF5Array)
library(scuttle)
library(parallel)
library(tidyr)
library(rstatix)
library(ggcorrplot)
library(nlme)
library(textshape)
library(scMerge)
library(pROC)
library(plotROC)
library(corrplot)
library(lme4)
library(ggfortify)
library(lmerTest)
library(ggnewscale)
library(forcats)
library(estimate)
library(CePa)
library(clusterProfiler)
library(SPATA2)
library(purrr)
library(spgwr)
library(sf)
library(tibble)
library(scMC)
library(enrichR)
library(GWmodel)
library(emmeans)
library(colorspace)
library(ggprism)
library(SeuratData)
library(SeuratDisk)
library(openxlsx)
library(here)
library(R.utils)
library(FNN)
library(psych)

source(here("utils.R"))

# change the current future plan to access parallelization (Seurat specific)

cores <- detectCores()

plan("multisession", workers = 4)

plan()

# define source data file for figures and tables

source_data <- createWorkbook()

##############
# DATA LOADING
##############

# read GEO metadata - placeholder until GEO series is public

# gse <- getGEO("GSE273528", GSEMatrix = TRUE)

geo_stRNAseq_metadata <- read_excel(here('data','Steiner_NatComm_stRNAseq_GEO_submission.xlsx'), 
                                      sheet = "Metadata", range = c('A41:AP57'))

##########
# TABLE S2
##########

tableS2 <- geo_stRNAseq_metadata %>% 
  select(`sample id`, `*library name`, `sequence batch`, age, sex, `smoking status`, 
         `capture frame VI`, `novel grade`, `TNM 8th edition stage`)

write.csv(tableS2, file = here('data','tableS2.csv'))



# define visium sample names

sample_names <- geo_stRNAseq_metadata$`*library name`

# define paths to samples

visium_path <-  here('data','GSE273378_RAW')

# unzip files and create new directory structure

image_files <- list.files(path = visium_path, 
                          pattern = "_tissue_lowres_image.png", 
                          full.names = TRUE)

image_folder <- here('data','GSE273378_RAW')

# Decompress the files

# Function to decompress and organize files

decompress_and_organize_files <- function(visium_path, file_name_pattern) {

  file_pattern <- paste0("_", file_name_pattern, ".gz")
  
  image_files <- list.files(path = visium_path, pattern = file_pattern, full.names = TRUE)
  
  decompressed_files <- sapply(image_files, function(x) {
    
    file_prefix <- sub(paste0(".*_(LM_.*)_", file_name_pattern, ".gz$"), "\\1", basename(x))
    
    subfolder_name <- file.path(here('data','GSE273378_RAW'), file_prefix, "outs", "spatial")
    
    if (!dir.exists(subfolder_name)) {
      dir.create(subfolder_name, recursive = TRUE)
    }
    
    decompressed_file <- file.path(subfolder_name, file_name_pattern)
    
    gunzip(x, destname = decompressed_file, overwrite = TRUE, remove = FALSE)
    
    return(decompressed_file)
  })
  
  return(decompressed_files)
}

decompressed_tissue_positions <- decompress_and_organize_files(visium_path, "tissue_positions_list.csv")

decompressed_lowres_files <- decompress_and_organize_files(visium_path, "tissue_lowres_image.png")

decompressed_hires_files <- decompress_and_organize_files(visium_path, "tissue_hires_image.png")

decompressed_scalefactors_files <- decompress_and_organize_files(visium_path, "scalefactors_json.json")

decompressed_detected_tissue_image_files <- decompress_and_organize_files(visium_path, "detected_tissue_image.jpg")

decompressed_aligned_fiducials_files <- decompress_and_organize_files(visium_path, "aligned_fiducials.jpg")

decompressed_pathology_files <- decompress_and_organize_files(visium_path, "pathology.csv")


# Function to decompress and organize filtered_feature_bc_matrix

organize_filtered_matrix <- function(visium_path, file_name_pattern) {
  
  file_pattern <- paste0("_", file_name_pattern)
  
  image_files <- list.files(path = visium_path, pattern = file_pattern, full.names = TRUE)
  
  decompressed_files <- sapply(image_files, function(x) {
    
    file_prefix <- sub(paste0(".*_(LM_.*)_", file_name_pattern), "\\1", basename(x))
    
    subfolder_name <- file.path(here('data','GSE273378_RAW'), file_prefix, "outs", "filtered_feature_bc_matrix")
    
    if (!dir.exists(subfolder_name)) {
      dir.create(subfolder_name, recursive = TRUE)
    }
    
    copied_file <- file.path(subfolder_name, paste0(file_name_pattern))
    
    file.copy(x, copied_file, overwrite = TRUE)
    
    return(copied_file)
  })
  
  return(decompressed_files)
}

organized_features <- organize_filtered_matrix(visium_path, "features.tsv.gz")

organized_barcodes <- organize_filtered_matrix(visium_path, "barcodes.tsv.gz")

organized_matrix <- organize_filtered_matrix(visium_path, "matrix.mtx.gz")


outs_paths <- list()

outs_paths <- lapply(X = list(sample_names), FUN = function(x) {
  x <- file.path(visium_path, x, "outs")
})


# paths to visium H&E images

visium_images <- lapply(X = outs_paths, FUN = function(x) {
  x <- file.path(x, "spatial")
})

# paths to 10x filtered spot-count matrix (post-10x tissue detection)

matrix_paths <- lapply(X = outs_paths, FUN = function(x) {
  x <- file.path(x, "filtered_feature_bc_matrix")
})

# paths to spot-level pathology annotations

pathology_paths <- lapply(X = outs_paths, FUN = function(x) {
  x <- file.path(x, "spatial/pathology.csv")
})


visium_lowres_images <- lapply(X = visium_images[[1]], FUN = function(x) {
  x <- Read10X_Image(x,
                     image.name = "tissue_lowres_image.png",
                     filter.matrix = TRUE)
})


visium_filtered_features <- lapply(X = matrix_paths[[1]], FUN = function(x) {
  x <- Read10X(x)
})

# Create visium seurat objects

visium_seurat_objects <- list()

for (i in 1:length(sample_names)) {
  
  visium_seurat_objects[[i]] <- CreateSeuratObject(counts = visium_filtered_features[[i]], 
                                                   assay="Spatial",
                                                   slice = 'slice1',
                                                   image = visium_lowres_images[[i]])
  
  Seurat::DefaultAssay(object = visium_lowres_images[[i]]) <- 'Spatial'  
  
  visium_seurat_objects[[i]][["slice1"]] <- visium_lowres_images[[i]]
  
  # add spot-level pathology annotations
  
  pathology <- read.csv(pathology_paths[[1]][i], row.names = 1)
  
  visium_seurat_objects[[i]]$pathology <- pathology
  
}



# reorder objects

names(visium_seurat_objects) <- sample_names

sample_names <- sample_names[order(as.numeric(sub('LM_SD_*','',sub('LM_SD_1216_*', '', sample_names))))]

visium_seurat_objects <- visium_seurat_objects[sample_names]


# add visium library IDs

for (i in 1:length(visium_seurat_objects)){
  
  visium_seurat_objects[[i]]@meta.data$orig.ident <- names(visium_seurat_objects)[i]
  
}


# read in nuclei count estimations per spot from stardist (except for sample 6)
# generated with stardist.py script on pre-filtered visium spots (raw data)

stardist_cell_counts <- list.files(path = here('data'),
                                   pattern = '_stardist_cell_counts\\.csv$')

stardist_cell_counts <- na.omit(stardist_cell_counts[match(sample_names, 
                                                           sub("_stardist_cell_counts.csv", 
                                                               "", stardist_cell_counts))])

stardist_cell_counts <- paste(here('data'), 
                              stardist_cell_counts, sep = '/')

data_list <- map(stardist_cell_counts, ~read.csv(.))

# remove sample 6

visium_seurat_objects_no6_unfiltered <- visium_seurat_objects[-6]

for (i in 1:length(visium_seurat_objects_no6_unfiltered)){
  
  visium_seurat_objects_no6_unfiltered[[i]]@meta.data$stardist_cell_counts <- data_list[[i]]$cell_count
  
}



#################
# Quality Control
#################

sum(sapply(visium_seurat_objects, ncol)) # 45,354 tissue spots in total dataset

# filter out spots with less than 250 genes

visium_seurat_objects <- lapply(X = visium_seurat_objects, FUN = function(x) {
  
  x <- subset(x, nFeature_Spatial > 250)
  
})

sum(sapply(visium_seurat_objects, ncol))  # (43,421 spots > 250 feature counts)

# save data checkpoint

# saveRDS(visium_seurat_objects, file = here('data','visium_seurat_objects_filtered.rds'))

visium_seurat_objects <- readRDS(here('data','visium_seurat_objects_filtered.rds'))


# merge visium samples into one seurat object

visium_seurat_merge <- merge(x = visium_seurat_objects[[1]], 
                             y = visium_seurat_objects[-1], 
                             add.cell.ids = names(visium_seurat_objects))


colorss <- colorRampPalette(brewer.pal(name="Paired", n = 12))(16)

# add additional metadata to seurat object

visium_seurat_merge@meta.data$column_name <- rownames(visium_seurat_merge@meta.data)

geo_stRNAseq_metadata <- geo_stRNAseq_metadata %>%
  rename(orig.ident = `*title`)

tmp <- merge(visium_seurat_merge@meta.data, geo_stRNAseq_metadata, by = "orig.ident")

tmp <- tmp[match(colnames(visium_seurat_merge), tmp$column_name),]

visium_seurat_merge@meta.data <- cbind(visium_seurat_merge@meta.data, tmp)





# stRNA-seq pathology annotations by VI status

######################
# EXTENDED DATA FIG 4A
######################

# use chi squared test

visium_seurat_merge_bub <- visium_seurat_merge


# remove sample LM_SD_6 for low mean counts due to tissue loss during workflow

Idents(visium_seurat_merge_bub) <- visium_seurat_merge_bub$orig.ident

visium_seurat_merge_bub <- subset(visium_seurat_merge_bub, idents = 'LM_SD_6', invert = TRUE)

dim(visium_seurat_merge_bub) 

# convert pathology annotations to factor, set normal lung as reference level

visium_seurat_merge_bub$pathology <- factor(visium_seurat_merge_bub$pathology, 
                                            levels = c("Normal Lung",unique(visium_seurat_merge_bub$pathology)[-which(unique(visium_seurat_merge_bub$pathology) == "Normal Lung")]))


# set unannotated areas to NA to drop from factor levels

levels(visium_seurat_merge_bub$pathology)[which(levels(visium_seurat_merge_bub$pathology) == "")] <- "Unannotated"

# remove path annotations with < 200 spots

visium_seurat_merge_bub$pathology <- visium_seurat_merge_bub$pathology[visium_seurat_merge_bub$pathology %!in% plyr::count(visium_seurat_merge_bub$pathology)$x[which(plyr::count(visium_seurat_merge_bub$pathology)$freq < 200)]]

visium_seurat_merge_bub$pathology <- droplevels(visium_seurat_merge_bub$pathology)

# downsample all remaining path annotatations to 200 spots

for(i in 1:length(levels(visium_seurat_merge_bub$pathology))){
  
  all_path <- rownames(visium_seurat_merge_bub@meta.data)[which(visium_seurat_merge_bub@meta.data$pathology == levels(visium_seurat_merge_bub$pathology)[i])]
  
  set.seed(123)
  
  downsampled_path <- sample(rownames(visium_seurat_merge_bub@meta.data)[which(visium_seurat_merge_bub@meta.data$pathology == levels(visium_seurat_merge_bub$pathology)[i])], size=200, replace=F)
  
  remove <- all_path[all_path %!in% downsampled_path]
  
  visium_seurat_merge_bub <- visium_seurat_merge_bub[,-which(colnames(visium_seurat_merge_bub) %in% remove)]
  
}


region_freq <- as.data.frame.matrix(table(visium_seurat_merge_bub$orig.ident, visium_seurat_merge_bub$pathology))

region_freq$orig.ident <- rownames(region_freq)

geo_stRNAseq_metadata$predominant_VI <- ifelse(geo_stRNAseq_metadata$`novel grade` == "VI",1,0)

region_VI_status <- geo_stRNAseq_metadata %>% dplyr::select(orig.ident, predominant_VI)

region_freq <- merge(region_freq, region_VI_status, by = "orig.ident")

# test at sample level not spot level


region_freq <- region_freq %>% `row.names<-`(.,NULL) %>% 
  column_to_rownames('orig.ident') %>% group_by(predominant_VI) %>% 
  summarise(across(everything(), ~sum(., na.rm = T))) %>% as.data.frame(.) %>%
  `row.names<-`(.,NULL) %>% column_to_rownames('predominant_VI')


model_stats <- list()

for(i in 1:length(levels(visium_seurat_merge_bub$pathology))){
  
  tmp <- data.frame('Other' = rowSums(region_freq[, -which(colnames(region_freq) == levels(visium_seurat_merge_bub$pathology)[i])]), 
                    'Region' = region_freq[,levels(visium_seurat_merge_bub$pathology)[i]])
  
  m <- chisq.test(tmp)
  
  model_stats[[i]] <- data.frame('region' = levels(visium_seurat_merge_bub$pathology)[i])
  
  model_stats[[i]][,'p_value'] <- m$p.value
  
  model_stats[[i]]['Direction'] <- ifelse(tmp$Region[1] < tmp$Region[2], "+","-")
  
}


model_stats_all <- bind_rows(model_stats, .id = "column_label")

# use this to make path labels in same order as bar plot

model_stats_all <- model_stats_all %>% dplyr::filter(region %!in% c('Collapsed Lung','Lymphoid Aggregate'))

model_stats_all$region <- factor(model_stats_all$region, 
                                 levels = c('Unannotated','Abnormal Pleura','Acinar',
                                            'Stroma','Cribriform','Desmoplastic Stroma',
                                            'Normal Vessel','Lepidic','Micropapillary','Normal Bronchus',
                                            'Normal Lung','Papillary','Solid','STAS','VI','VPI'))

levels(model_stats_all$region)[which(levels(model_stats_all$region) == 'Normal Vessel')] <- 'Vessel'

levels(model_stats_all$region)[which(levels(model_stats_all$region) == 'VI')] <- 'VI Focus'

levels(model_stats_all$region)[which(levels(model_stats_all$region) == 'VPI')] <- 'VPI Focus'


# bubble plot heatmap to chi squared test p values 

efig4A <- model_stats_all

# convert to -log10 p values

efig4A$FDR <- p.adjust(efig4A$p_value, n = length(efig4A$p_value))

efig4A$FDR <- -log10(efig4A$FDR)

efig4A$cluster <- as.factor(1)

addWorksheet(source_data, "E. Figure 4A")

writeData(source_data, sheet = "E. Figure 4A", 
          x = "E. Figure 4A", 
          startCol =  1, startRow = 1)

writeData(source_data, sheet = 'E. Figure 4A',
          x = efig4A, startCol = 1, startRow = 3)

efig4A_gg <- ggplot(data = efig4A, aes(x = cluster, y = region)) +
  geom_point(aes(fill = Direction, size = FDR), shape = 21) +
  scale_size(breaks=c(0,2,20,50,100),
             labels=c(0,2,20,50,100), 
             range = c(0,9)) +
  scale_fill_manual(values = c('-' = "#390099", '+' = "#FF0054")) +
  theme_classic() +
  theme(
    legend.position = 'right',
    panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
    plot.title = element_text(size=22, hjust = 0.5),
    axis.title.x = element_text(size=20),
    axis.title.y = element_text(size=20),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size=20),
    axis.text=element_text(size=20),
    axis.ticks=element_blank(),
    legend.title = element_text(size=20),
    legend.key.size = unit(0.8, "cm"),
    legend.text = element_text(size=20),
    text = element_text(family = 'Calibri')) +
  labs(title = 'StRNA-seq pathology \n annotations by VI status\n',
       x = element_blank(), 
       y = NULL, fill = element_blank()) +
  scale_y_discrete(limits = rev) + 
  scale_x_discrete(position = "top") + 
  guides(size=guide_legend(title="-log₁₀(p.adj)"),
         fill = guide_legend(title='Association with VI',
                             override.aes = list(size=8))) # +
guides(size='none', fill ='none') 

efig4A_gg

plot_dimensions <- dev.size("in")

aspect_ratio <- plot_dimensions[1] / plot_dimensions[2]

aspect_ratio <- 0.8363914

desired_height <- 7  

desired_width <- desired_height * aspect_ratio

TIFF_dpi <- 600

ggsave(here('figures','efig4A.tiff'),
       plot = efig4A_gg,
       device = 'tiff',
       width = desired_width,
       height = desired_height,
       units = 'in', dpi = TIFF_dpi)



######################
# EXTENDED DATA FIG 4B
######################

median(visium_seurat_merge$nFeature_Spatial)

efig4B <- visium_seurat_merge@meta.data %>%
  select('nFeature_Spatial', 'sample id')

addWorksheet(source_data, "E. Figure 4B")

writeData(source_data, sheet = "E. Figure 4B", 
          x = "E. Figure 4B", 
          startCol =  1, startRow = 1)

writeData(source_data, sheet = 'E. Figure 4B',
          x = efig4B, startCol = 1, startRow = 3)

efig4B_gg <- VlnPlot(visium_seurat_merge, features = c("nFeature_Spatial"), 
                     group.by = "sample id", ncol = 2, pt.size=0, cols = colorss)  + 
  labs(title=NULL, x = '\nSample', y = 'StRNA-seq genes/spot\n') +
  theme(
    plot.title = element_text(size=22, hjust = 0.5, face='plain'),
    axis.title.x = element_text(size=20),
    axis.title.y = element_text(size=20),
    axis.text.x = element_text(size=20, angle = 45),
    legend.position="none", 
    axis.text=element_text(size=20),
    legend.title = element_text(),
    text = element_text(family = 'Calibri'))

efig4B_gg

# Get the current plot dimensions in inches
plot_dimensions <- dev.size("in")

# Calculate aspect ratio
aspect_ratio <- plot_dimensions[1] / plot_dimensions[2]

aspect_ratio <- 2.22164

# Define a reasonable height (in inches)
desired_height <- 7 

# Calculate corresponding width to maintain aspect ratio
desired_width <- desired_height * aspect_ratio

# Define the resolution
TIFF_dpi <- 600

# Save the plot with adjusted dimensions
ggsave(here('figures','efig4B.tiff'),
       plot = efig4B_gg,
       device = 'tiff',
       width = desired_width,
       height = desired_height,
       units = 'in', dpi = TIFF_dpi)


######################
# EXTENDED DATA FIG 4C
######################

median(visium_seurat_merge$nCount_Spatial)

efig4C <- visium_seurat_merge@meta.data %>%
  select('nCount_Spatial', 'sample id')

addWorksheet(source_data, "E. Figure 4C")

writeData(source_data, sheet = "E. Figure 4C", 
          x = "E. Figure 4C", 
          startCol =  1, startRow = 1)

writeData(source_data, sheet = 'E. Figure 4C',
          x = efig4C, startCol = 1, startRow = 3)

efig4C_gg <- VlnPlot(visium_seurat_merge, features = c("nCount_Spatial"), 
                     group.by = "sample id", ncol = 2, pt.size=0, cols = colorss)  + 
  labs(title=NULL, x = '\nSample', y = 'StRNA-seq UMIs/spot\n') +
  theme(
    plot.title = element_text(size=22, hjust = 0.5, face='plain'),
    axis.title.x = element_text(size=20),
    axis.title.y = element_text(size=20),
    axis.text.x = element_text(size=20, angle = 45),
    legend.position="none", 
    axis.text=element_text(size=20),
    legend.title = element_text(),
    text = element_text(family = 'Calibri'))

efig4C_gg

# Save the plot with adjusted dimensions

ggsave(here('figures','efig4C.tiff'),
       plot = efig4C_gg,
       device = 'tiff',
       width = desired_width,
       height = desired_height,
       units = 'in', dpi = TIFF_dpi)



# remove sample LM_SD_6 for low mean counts due to tissue loss during workflow

Idents(visium_seurat_merge) <- visium_seurat_merge$orig.ident

visium_seurat_merge <- subset(visium_seurat_merge, idents = 'LM_SD_6', invert = TRUE)

dim(visium_seurat_merge)




#####################################
# GLOBAL NORMALIZATION OF ALL SAMPLES
#####################################


# increase max size of global objects passed to future package

options(future.globals.maxSize = 8000 * 1024^2)

# perform sctransform-based normalization on merged dataset

visium_seurat_merge <- SCTransform(visium_seurat_merge, assay = "Spatial",
                                   do.scale = FALSE, do.center = TRUE, 
                                   verbose = TRUE)

# save data checkpoint

# saveRDS(visium_seurat_merge, file = here('data','visium_seurat_merge_no6.rds'))

visium_seurat_merge <- readRDS(here('data','visium_seurat_merge_no6.rds'))


visium_seurat_merge <- ScaleData(visium_seurat_merge, assay = "SCT")

# dimensionality reduction and embedding

visium_seurat_merge <- RunPCA(visium_seurat_merge, assay = "SCT", verbose = TRUE,
                              npcs = 50)

visium_seurat_merge <- RunUMAP(visium_seurat_merge, dims = 1:30, 
                               assay = "SCT",
                               reduction = "pca")

# save data checkpoint

# saveRDS(visium_seurat_merge, file = here('data','visium_seurat_merge_sct_pca_scale_umap_no6.rds'))

visium_seurat_merge <- readRDS(here('data','visium_seurat_merge_sct_pca_scale_umap_no6.rds'))


# comparison with matched bulk RNAseq

######################
# EXTENDED DATA FIG 4D
######################

# create pseudobulk of samples (treat all sections as separate tumors)

visium_seurat_merge_sce <- as.SingleCellExperiment(visium_seurat_merge)

visium_seurat_merge_sce_pseudo <- aggregateAcrossCells(visium_seurat_merge_sce, 
                                                       id=colData(visium_seurat_merge_sce)[,c("block.id")])

dim(visium_seurat_merge_sce_pseudo)

# rename to compare with bulk sample IDs

colnames(visium_seurat_merge_sce_pseudo@assays@data$counts) <- sapply(strsplit(colnames(visium_seurat_merge_sce_pseudo@assays@data$counts), "_"), "[", 1)


# load bulk RNAseq data for same tumors from validation cohort

validation_dge <- readRDS(here('data','validation_dge.rds'))

gene_annotations <- readRDS(here('data','gene_annotations.rds'))

# convert ensemble gene names to hgnc

validation_dge <- ensemble_to_hgnc(validation_dge)


# subset validation dge object to samples that also had stRNA-seq data

colnames(validation_dge$counts) <- validation_dge$samples$`patient id`

validation_dge$counts <- validation_dge$counts[,which(colnames(validation_dge$counts) %in%
                                                        colnames(visium_seurat_merge_sce_pseudo@assays@data$counts))]

dim(validation_dge$counts) # 13 samples

# put sample names in same order (also duplicates bulk data for patient 21448 and 30868)

validation_dge$counts <- validation_dge$counts[,colnames(visium_seurat_merge_sce_pseudo@assays@data$counts)]

dim(validation_dge$counts)

# subset both assays to same genes

gene_intersect <- intersect(rownames(validation_dge$counts), rownames(visium_seurat_merge_sce_pseudo@assays@data$counts)) 

length(gene_intersect) 

validation_dge$counts <- validation_dge$counts[gene_intersect,]

visium_seurat_merge_sce_pseudo@assays@data$counts <- visium_seurat_merge_sce_pseudo@assays@data$counts[gene_intersect,]

# log transform and cpm normalize both count matrices

validation_dge$counts <- edgeR::cpm(as.matrix(validation_dge$counts +1), log=T)

visium_seurat_merge_sce_pseudo@assays@data$counts <- edgeR::cpm((visium_seurat_merge_sce_pseudo@assays@data$counts +1), log=T)

# calculate mean gene expression per assay per gene

avg_gene_bulk <- rowMeans(validation_dge$counts)

avg_gene_pseudo <- rowMeans(visium_seurat_merge_sce_pseudo@assays@data$counts)

# correlation plot of average expression of each gene (entire dataset)

efig4D <- data.frame('Avg_exp_bulk' = avg_gene_bulk,
                     'Avg_exp_visium' = avg_gene_pseudo)

addWorksheet(source_data, "E. Figure 4D")

writeData(source_data, sheet = "E. Figure 4D", 
          x = "E. Figure 4D", 
          startCol =  1, startRow = 1)

writeData(source_data, sheet = 'E. Figure 4D',
          x = efig4D, startCol = 1, startRow = 3)


# Compute Spearman correlation for each of the 15 sample pairs
efig4D_2 <- sapply(1:15, function(i) {
  cor(validation_dge$counts[, i],
      visium_seurat_merge_sce_pseudo@assays@data$counts[, i],
      method = "spearman")
})

# Label the results for clarity and print the correlations
names(efig4D_2) <- paste0("Sample_", 1:15)
print(efig4D_2)

addWorksheet(source_data, "E. Figure 4D_2")

writeData(source_data, sheet = "E. Figure 4D_2", 
          x = "E. Figure 4D_2", 
          startCol =  1, startRow = 1)

writeData(source_data, sheet = 'E. Figure 4D_2',
          x = efig4D_2, startCol = 1, startRow = 3)


efig4D_gg <- ggplot(efig4D, aes(x = Avg_exp_bulk, y = Avg_exp_visium)) +
  geom_point(color = "#E64B35FF", alpha = 0.75) +
  labs(title="Mean gene expression in validation cohort\n(n=15 matched samples)\n", 
       x = "\nBulk RNAseq log CPM", y = "Pseudobulk stRNA-seq log CPM\n") +
  theme_classic() +
  theme(
    plot.title = element_text(size=22, hjust = 0.5),
    axis.title.x = element_text(size=20),
    axis.title.y = element_text(size=20),
    legend.position="none", 
    axis.text=element_text(size=20),
    text = element_text(family = 'Calibri')) +
  geom_smooth(method = "lm", se = TRUE, formula = y ~ x, color = "black") + 
  stat_cor(method="spearman",
           aes(label = paste(after_stat(r.label),
                             sub("e","%.% 10^",after_stat(p.label)), 
                             sep = "~`,`~")), 
           size = 6, family = 'Calibri') 

efig4D_gg

plot_dimensions <- dev.size("in")

aspect_ratio <- plot_dimensions[1] / plot_dimensions[2]

aspect_ratio <- 0.9561934

desired_height <- 7  

desired_width <- desired_height * aspect_ratio

TIFF_dpi <- 600

ggsave(here('figures','efig4D.tiff'),
       plot = efig4D_gg,
       device = 'tiff',
       width = desired_width,
       height = desired_height,
       units = 'in', dpi = TIFF_dpi)



###################
# LOAD VI SIGNATURE
###################

VI_genesets <- list() 

# read in VI signature and clusters

VI_gene_cluster <- readRDS(here('data','VI_gene_clusters.rds'))

# read in VI predictor

VI_predictor_genes <- readRDS(here('data','VI_predictor_genes.rds'))

length(which(VI_predictor_genes %in% rownames(visium_seurat_merge))) # 40/48 genes


VI_predictor_genes_up <- readRDS(here('data','VI_predictor_genes_up.rds'))

VI_predictor_genes_dn <- readRDS(here('data','VI_predictor_genes_dn.rds'))

VI_genesets$VI_predictor_genes <- VI_predictor_genes


VI_genesets$VI_predictor_genes_up <- VI_predictor_genes_up

VI_genesets$VI_predictor_genes_dn <- VI_predictor_genes_dn


VI_genesets$VI_gene_cluster_1 <- rownames(VI_gene_cluster %>% filter(Cluster == 1))

VI_genesets$VI_gene_cluster_2 <- rownames(VI_gene_cluster %>% filter(Cluster == 2))

VI_genesets$VI_gene_cluster_3 <- rownames(VI_gene_cluster %>% filter(Cluster == 3))

VI_genesets$VI_gene_cluster_4 <- rownames(VI_gene_cluster %>% filter(Cluster == 4))



###################################################
# VI SIGNATURE CORRELATION WITH MATCHED BULK RNASEQ
###################################################


avg_VI_gene_bulk <- rowMeans(validation_dge$counts[which(rownames(validation_dge$counts) %in% rownames(VI_gene_cluster)),])

avg_VI_gene_pseudo <- rowMeans(visium_seurat_merge_sce_pseudo@assays@data$counts[which(rownames(visium_seurat_merge_sce_pseudo@assays@data$counts) %in% rownames(VI_gene_cluster)),])

# correlation plot of average expression of each gene (entire dataset)

efig4E <- data.frame('Avg_exp_bulk' = avg_VI_gene_bulk,
                   'Avg_exp_visium' = avg_VI_gene_pseudo)

addWorksheet(source_data, "E. Figure 4E")

writeData(source_data, sheet = "E. Figure 4E", 
          x = "E. Figure 4E", 
          startCol =  1, startRow = 1)

writeData(source_data, sheet = 'E. Figure 4E',
          x = efig4E, startCol = 1, startRow = 3)

######################
# EXTENDED DATA FIG 4E
######################


efig4E_gg <- ggplot(efig4E, aes(x = Avg_exp_bulk, y = Avg_exp_visium)) +
      geom_point(color = "#E64B35FF", alpha = 0.75) +
      labs(title="VI gene expression (426 genes) in validation\ncohort (n=15 matched samples)\n", 
           x = "\nBulk RNAseq log CPM", y = "Pseudobulk stRNA-seq log CPM\n") +
      theme_classic() +
      theme(
        plot.title = element_text(size=22, hjust = 0.5),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        legend.position="none", 
        axis.text=element_text(size=20),
        text = element_text(family = 'Calibri')) +
      geom_smooth(method = "lm", se = TRUE, formula = y ~ x, color = "black") + 
      stat_cor(method="spearman",
               aes(label = paste(after_stat(r.label),
                                 sub("e","%.% 10^",after_stat(p.label)), 
                                 sep = "~`,`~")), 
               size = 6, family = 'Calibri') 

plot_dimensions <- dev.size("in")

aspect_ratio <- plot_dimensions[1] / plot_dimensions[2]

aspect_ratio <- 0.9561934

desired_height <- 7  

desired_width <- desired_height * aspect_ratio

TIFF_dpi <- 600

ggsave(here('figures','efig4E.tiff'),
       plot = efig4E_gg,
       device = 'tiff',
       width = desired_width,
       height = desired_height,
       units = 'in', dpi = TIFF_dpi)





# spot scoring of merged seurat object with signatures

visium_seurat_merge <- AddModuleScore(visium_seurat_merge, features = VI_genesets,
                                      name = '_score')

names(visium_seurat_merge@meta.data)[grep("_score", names(visium_seurat_merge@meta.data))] <- names(VI_genesets)

# pathology annotations

###########
# FIGURE 2A
###########

tmp <- visium_seurat_merge

tmp@meta.data$predominant_VI <- ifelse(tmp@meta.data$`novel grade` == 'VI',1,0)

Idents(tmp) <- tmp@meta.data$predominant_VI

tmp$pathology <- factor(tmp$pathology)


for (i in 2:length(levels(tmp$pathology))){
  
  tmp@meta.data[,levels(tmp$pathology)[i]] <- ifelse(tmp$pathology == levels(tmp$pathology)[i],1,0)
  
}


# set unannotated areas to NA to drop from factor levels

levels(tmp$pathology)[which(levels(tmp$pathology) == "")] <- "Unannotated"

# collapse levels (could also create a 'Other' level)

tmp$pathology <- fct_collapse(tmp$pathology, Vessel = c("Hyal Vessel","Hyper Vessel","Normal Vessel"),
                              `Normal Lung` = c("Normal Lung","Normal Pleura","Pulmonary Macrophage"),
                              `Stroma` = c("Stroma","Cartilage","Lymphoid Aggregate","Lymphoid Aggregate (VI)","Collapsed Lung","Necrosis"))

levels(tmp$pathology)[which(levels(tmp$pathology) == "VI")] <- "VI Focus"

levels(tmp$pathology)[which(levels(tmp$pathology) == "VPI")] <- "VPI Focus"

fig2A <- tmp@meta.data

fig2A$predominant_VI <- as.factor(fig2A$predominant_VI)

colourCount = length(unique(fig2A$pathology))

getPalette = colorRampPalette(brewer.pal(9, "Set1"))

use_colors <- c(Unannotated = '#999999',
                `Abnormal Pleura` = '#7E6E85',
                Acinar = '#AC5782',
                Stroma = '#C66764',
                Cribriform = '#E1C62F',
                `Desmoplastic Stroma`= '#B16C29',
                Vessel = '#F17EB4',
                Lepidic = '#56A255',
                Micropapillary = '#E3712B',
                `Normal Bronchus` = '#AFEEEE',
                `Normal Lung` = '#3881B0',
                Papillary = '#FFA10D',
                Solid = '#FFE528',
                STAS = '#CB8CAD',
                `VI Focus` = '#E41A1C',
                `VPI Focus` = '#874F6F')

levels(fig2A$predominant_VI) <- c('VI-','VI+')

addWorksheet(source_data, "Figure 2A")

writeData(source_data, sheet = "Figure 2A", 
          x = "Figure 2A", 
          startCol =  1, startRow = 1)

writeData(source_data, sheet = 'Figure 2A',
          x = fig2A, startCol = 1, startRow = 3)

fig2A_gg <- ggplot(fig2A, aes(
  x = `sample id`
)) +
  geom_bar(aes(fill = pathology), position = "fill", width = 0.95) +
  scale_fill_manual(values = use_colors) +
  coord_cartesian(ylim = c(1.025,0)) +
  new_scale_fill() +
  geom_col(data = dplyr::distinct(fig2A, `sample id`, predominant_VI), 
           aes(y = -.05, fill = predominant_VI), width = 1) +
  scale_fill_manual(values = c("#390099","#FF0054")) +
  geom_col(data = dplyr::distinct(fig2A, `sample id`, predominant_VI), 
           aes(y = -.01), fill = "white", color = "white", width = 1) +
  facet_grid(. ~ predominant_VI,
             scales = "free_x", space = "free_x") +
  theme(
    strip.text = element_blank(),
    strip.background = element_blank(),
    plot.title = element_text(size=24,hjust = 0.5),
    axis.text.y = element_text(size=20),
    axis.text.x = element_text(size=20, vjust = 10),
    axis.ticks = element_blank(),
    panel.background = element_blank(),
    axis.title.x = element_text(size=24),
    axis.title.y = element_text(size=24),
    legend.key.size = unit(1.25, 'cm'),
    legend.text = element_text(size=20),
    legend.position = 'left',
    legend.title = element_blank(),
    text = element_text(family = 'Calibri')) +
  scale_y_continuous(name = NULL, 
                     sec.axis = sec_axis(~., name = "% stRNA-seq spots\n")) +
  guides(y = "none") +
  labs(x = "Sample") +
  ggtitle("Pathology annotations in stRNA-seq dataset\n")

fig2A_gg

# Get the current plot dimensions in inches

plot_dimensions <- dev.size("in")

# Calculate aspect ratio

aspect_ratio <- plot_dimensions[1] / plot_dimensions[2]

aspect_ratio <- 1.550499

# Define a reasonable height (in inches)

desired_height <- 10  

# Calculate corresponding width to maintain aspect ratio

desired_width <- desired_height * aspect_ratio

# Define the resolution

TIFF_dpi <- 600

# Save the plot with adjusted dimensions

ggsave(here('figures','fig2A.tiff'),
       plot = fig2A_gg,
       device = 'tiff',
       width = desired_width,
       height = desired_height,
       units = 'in', dpi = TIFF_dpi)



# plot mean VI cluster enrichment scores by pathology region

###########
# FIGURE 2C
###########

# using mixed linear effects model with sample ID as random effect

mixed_effects_stats <- list()

visium_seurat_merge_bub <- visium_seurat_merge

visium_seurat_merge_bub$predominant_VI <- ifelse(visium_seurat_merge_bub$`novel grade` == 'VI',1,0)

# convert pathology annotations to factor, set normal lung as reference level

visium_seurat_merge_bub$pathology <- factor(visium_seurat_merge_bub$pathology, levels = c("Normal Lung",unique(visium_seurat_merge_bub$pathology)[-which(unique(visium_seurat_merge_bub$pathology) == "Normal Lung")]))


# set unannotated areas to NA to drop from factor levels

levels(visium_seurat_merge_bub$pathology)[which(levels(visium_seurat_merge_bub$pathology) == "")] <- "Unannotated"

# remove path annotations with < 200 spots

visium_seurat_merge_bub$pathology <- visium_seurat_merge_bub$pathology[visium_seurat_merge_bub$pathology %!in% plyr::count(visium_seurat_merge_bub$pathology)$x[which(plyr::count(visium_seurat_merge_bub$pathology)$freq < 200)]]

visium_seurat_merge_bub$pathology <- droplevels(visium_seurat_merge_bub$pathology)

# downsample all remaining path annotatations to 200 spots

for(i in 1:length(levels(visium_seurat_merge_bub$pathology))){
  
  all_path <- rownames(visium_seurat_merge_bub@meta.data)[which(visium_seurat_merge_bub@meta.data$pathology == levels(visium_seurat_merge_bub$pathology)[i])]
  
  set.seed(123)
  
  downsampled_path <- sample(rownames(visium_seurat_merge_bub@meta.data)[which(visium_seurat_merge_bub@meta.data$pathology == levels(visium_seurat_merge_bub$pathology)[i])], size=200, replace=F)
  
  remove <- all_path[all_path %!in% downsampled_path]
  
  visium_seurat_merge_bub <- visium_seurat_merge_bub[,-which(colnames(visium_seurat_merge_bub) %in% remove)]
  
}

levels(visium_seurat_merge_bub@meta.data$pathology)[which(levels(visium_seurat_merge_bub@meta.data$pathology) == 'Unannotated')] <- NA

visium_seurat_merge_bub@meta.data$pathology <- droplevels(visium_seurat_merge_bub@meta.data$pathology)

clusters <- c('VI_gene_cluster_1', 'VI_gene_cluster_2',
              'VI_gene_cluster_3', 'VI_gene_cluster_4')

for(i in 1:length(clusters)){
  
  print(i)
  
  data <- visium_seurat_merge_bub@meta.data
  
  scores <- aggregate(data[,clusters[i]] ~ pathology, data = data, mean)
  
  colnames(scores)[2] <- clusters[i]
  
  data <- data.frame(y = data[,clusters[i]], x = data$pathology, sample = data$`sample id`, VI = data$predominant_VI)
  
  m <- lmer(y ~ 1 + x + (1 | sample), data = data)
  
  mixed_effects_stats[[i]] <- data.frame(summary(m)$coefficients)
  
  mixed_effects_stats[[i]]$pathology <- sapply(strsplit(rownames(mixed_effects_stats[[i]]),'x'),"[",2)
  
  mixed_effects_stats[[i]]$pathology[1] <- 'Normal Lung' # rename intercept
  
  mixed_effects_stats[[i]] <- merge(mixed_effects_stats[[i]], scores, by = 'pathology')
  
  mixed_effects_stats[[i]]$scores <- scores$`data[, clusters[i]]`
  
  rownames(mixed_effects_stats[[i]]) <- mixed_effects_stats[[i]]$pathology
  
  mixed_effects_stats[[i]] <- mixed_effects_stats[[i]] %>% dplyr::select(clusters[i], Pr...t..)
  
  colnames(mixed_effects_stats[[i]])[which(colnames(mixed_effects_stats[[i]]) == 'Pr...t..')] <- paste('p_value',clusters[i],sep = '_')
  
}



mixed_effects_stats_all <- bind_cols(mixed_effects_stats)


# bubble plot heatmap to show scores and p values together

mixed_effects_stats_all$region <- rownames(mixed_effects_stats_all)

tmp <- mixed_effects_stats_all %>% dplyr::select(region, VI_gene_cluster_1, VI_gene_cluster_2,
                                                 VI_gene_cluster_3, VI_gene_cluster_4)

fig2C <- gather(tmp, cluster, score, VI_gene_cluster_1:VI_gene_cluster_4, factor_key=TRUE)

tmp <- mixed_effects_stats_all %>% dplyr::select(p_value_VI_gene_cluster_1, p_value_VI_gene_cluster_2,
                                                 p_value_VI_gene_cluster_3, p_value_VI_gene_cluster_4)

data_long2 <- gather(tmp, cluster_p, p_value, p_value_VI_gene_cluster_1:p_value_VI_gene_cluster_4, factor_key=TRUE)

fig2C <- cbind(fig2C, data_long2)


# need to adjust p values for multiple hypothesis testing

fig2C$FDR <- p.adjust(fig2C$p_value, n = length(fig2C$p_value))

fig2C$FDR <- -log10(fig2C$FDR)

levels(fig2C$cluster) <- c('VI Cluster 1','VI Cluster 2','VI Cluster 3','VI Cluster 4')

fig2C$region <- factor(fig2C$region, levels = c('Normal Lung','Normal Vessel','Lymphoid Aggregate','Stroma','Desmoplastic Stroma','Normal Bronchus',
                                                'Collapsed Lung','Abnormal Pleura','Lepidic','Acinar','Papillary','Micropapillary','Cribriform','Solid',
                                                'VPI','VI'))

levels(fig2C$region)[which(levels(fig2C$region) == 'VI')] <- 'VI Focus'

levels(fig2C$region)[which(levels(fig2C$region) == 'VPI')] <- 'VPI Focus'


# scale by cluster 

fig2C <- fig2C %>% 
  group_by(cluster) %>% 
  mutate(score = (score - mean(score)) / sd(score))



# plot raw enrichment scores

# remove from plot enrichment scores > 0.01 FDR

fig2C <- fig2C %>%
  mutate(FDR_filtered = ifelse(FDR > 2, FDR, NA))

addWorksheet(source_data, "Figure 2C")

writeData(source_data, sheet = "Figure 2C", 
          x = "Figure 2C", 
          startCol =  1, startRow = 1)

writeData(source_data, sheet = 'Figure 2C',
          x = fig2C, startCol = 1, startRow = 3)

fig2C_gg <- ggplot(data = fig2C, aes(x = cluster, y = region)) +
  geom_point(aes(fill = score, size = FDR_filtered), shape = 21) +
  scale_size(breaks=c(3,20,50,100),labels=c(3,20,50,100), range = c(2,9)) +
  scale_fill_gradient2(low = "#007FFF", mid = "white", high = "red") +
  theme_classic() +
  theme(
    legend.position = 'right',
    panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
    plot.title = element_text(size=22, hjust = 0.5),
    axis.title.x = element_text(size=20),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size=20),
    axis.text=element_text(size=20),
    axis.ticks=element_blank(),
    legend.title = element_text(size=20),
    legend.key.size = unit(0.6, "cm"),
    legend.text = element_text(size=20),
    text = element_text(family = 'Calibri')) +
  labs(title = 'VI cluster expression by pathology annotation\n(n=15 samples)',
       x = element_blank(), y = 'Visium Region', 
       fill = 'Enrichment Score \n (Scaled)') +
  scale_y_discrete(limits = rev) + 
  scale_x_discrete(position = "top") +
  guides(size=guide_legend(title="-log₁₀(p.adj)"))  +
  annotate(
    xmin = c(-Inf,1.5,2.5,3.5), 
    xmax = c(1.5,2.5,3.5,Inf),
    ymin = 16.5,
    ymax = 17.25,
    geom = "rect",
    fill = c(unname(VI_cluster_colors))#[row_ord])
  ) +
  geom_hline(yintercept=c(16.5), color=c('black'), linewidth = 1)

fig2C_gg

# Get the current plot dimensions in inches

plot_dimensions <- dev.size("in")

# Calculate aspect ratio

aspect_ratio <- plot_dimensions[1] / plot_dimensions[2]

aspect_ratio <- 0.8267516

# Define a reasonable height (in inches)

desired_height <- 8  

# Calculate corresponding width to maintain aspect ratio

desired_width <- desired_height * aspect_ratio

# Define the resolution
TIFF_dpi <- 600

# Save the plot with adjusted dimensions

ggsave(here('figures','fig2C.tiff'),
       plot = fig2C_gg,
       device = 'tiff',
       width = desired_width,
       height = desired_height,
       units = 'in', dpi = TIFF_dpi)


# association of VI clusters with VI within high-grade annotated spots only

#############
# FIGURE 2D_1
#############

# convert pathology annotations to individual categories

tmp <- visium_seurat_merge

tmp@meta.data$predominant_VI <- ifelse(tmp@meta.data$`novel grade` == 'VI',1,0)

Idents(tmp) <- tmp@meta.data$predominant_VI

tmp$pathology <- factor(tmp$pathology)

tmp$predominant_VI <- factor(tmp$predominant_VI)

for (i in 2:length(levels(tmp$pathology))){
  
  tmp@meta.data[,levels(tmp$pathology)[i]] <- ifelse(tmp$pathology == levels(tmp$pathology)[i],1,0)
  
}


# set unannotated areas to NA to drop from factor levels

levels(tmp$pathology)[which(levels(tmp$pathology) == "")] <- "Unannotated"

tmp <- tmp@meta.data %>% dplyr::filter(pathology %in% c('Solid','Micropapillary',
                                                        'Cribriform'))

data_long <- gather(tmp, cluster, score, VI_gene_cluster_1:VI_gene_cluster_4, 
                    factor_key=TRUE)

names(data_long)[names(data_long) == "predominant_VI"] <- "VI Status"


# scale by cluster (for visualization purposes)

data_long <- data_long %>%
  group_by(cluster) %>%
  mutate(score = (score - mean(score)) / sd(score))

levels(data_long$`VI Status`) <- c('VI-','VI+')

use_colors_VI <- c(`VI-` = "#390099",
                   `VI+` = "#FF0054")


fig2D_1 <- data_long %>% filter(cluster == 'VI_gene_cluster_1')

addWorksheet(source_data, "Figure 2D-1")

writeData(source_data, sheet = "Figure 2D-1", 
          x = "Figure 2D-1", 
          startCol =  1, startRow = 1)

writeData(source_data, sheet = 'Figure 2D-1',
          x = fig2D_1, startCol = 1, startRow = 3)

data <- data.frame(y = fig2D_1$score, x = fig2D_1$`VI Status`, 
                   sample = fig2D_1$`sample id`, pathology = fig2D_1$pathology)

m <- lmer(y ~ 1 + x + (1 | sample) + pathology, data = data)


# type II anova 
Anova(m, type = 2)

mixed_effects_stats <- data.frame(Anova(m, type = 2, test.statistic = 'Chisq'))

mixed_effects_stats <- mixed_effects_stats %>% select(Pr..Chisq.)

p <- ggplot(data = fig2D_1, aes(x = `VI Status`, y = score, 
                                fill = `VI Status`)) +
  theme_classic() + 
  coord_cartesian(clip = "off") +
  scale_fill_manual(values = c(unname(use_colors_VI))) +
  geom_violin(lwd=0.5, width = 0.75) +
  geom_boxplot(color = "black", fill = 'white', width = .15, 
               outlier.alpha = 0, lwd=0.5) +
  theme(
    plot.title = element_text(size=20, hjust = 0.5),
    panel.background = element_blank(),
    axis.text.y=element_text(size=20),
    axis.text.x=element_text(size=20),
    axis.ticks.x=element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size=20),
    strip.text.x = element_text(size=0),
    legend.position = 'right',
    legend.title = element_blank(),
    legend.text = element_text(size=20),
    text = element_text(family = 'Calibri')) +
  labs(title = "High-grade pattern\nannotated spots only\n",
       x = "\nVI Status", y = "Enrichment score\n") +
  ylim(min(fig2D_1$score),max(fig2D_1$score)+2) +
  facet_grid(cols = vars(cluster),
             scales = "free_y", space = "free_y") +
  geom_bracket(
    xmin = 'VI-', xmax = "VI+", y.position = max(fig2D_1$score)+1,
    label = paste('Anova, p =',round(mixed_effects_stats$Pr..Chisq.[1], digits = 3)), 
    tip.length = c(0.02, 0.02), vjust = -0.25,
    size = 0.5,
    label.size = 6,
    inherit.aes = FALSE,
    family = 'Calibri') +
  scale_x_discrete(labels = c("VI-\n(n=8 samples)","VI+\n(n=7 samples)"))


# color individual  facets

g <- ggplot_gtable(ggplot_build(p))

strip_both <- which(grepl('strip-', g$layout$name))

fills <- c(unname(VI_cluster_colors))

k <- 1

for (i in strip_both) {
  
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  
  k <- k+1
  
}



plot_dimensions <- dev.size("in")

# Calculate aspect ratio
aspect_ratio <- plot_dimensions[1] / plot_dimensions[2]

aspect_ratio <- 1.2

# Define the resolution
dpi <- 600
width <- 5.75
height <- 5.75 / aspect_ratio

# Save the plot using tiff device with high resolution

tiff(filename = here('figures','fig2D_1.tiff'),
     width = width, height = height, units = "in", res = dpi)

grid.draw(g)

dev.off()

#############
# FIGURE 2D_2
#############

fig2D_2 <- data_long %>% filter(cluster == 'VI_gene_cluster_2')

addWorksheet(source_data, "Figure 2D-2")

writeData(source_data, sheet = "Figure 2D-2", 
          x = "Figure 2D-2", 
          startCol =  1, startRow = 1)

writeData(source_data, sheet = 'Figure 2D-2',
          x = fig2D_2, startCol = 1, startRow = 3)

data <- data.frame(y = fig2D_2$score, x = fig2D_2$`VI Status`, 
                   sample = fig2D_2$`sample id`, pathology = fig2D_2$pathology)

m <- lmer(y ~ 1 + x + (1 | sample) + pathology, data = data)


# type II anova 
Anova(m, type = 2)

mixed_effects_stats <- data.frame(Anova(m, type = 2, test.statistic = 'Chisq'))

mixed_effects_stats <- mixed_effects_stats %>% select(Pr..Chisq.)

p <- ggplot(data = fig2D_2, aes(x = `VI Status`, y = score, 
                                fill = `VI Status`)) +
  theme_classic() + 
  coord_cartesian(clip = "off") +
  scale_fill_manual(values = c(unname(use_colors_VI))) +
  geom_violin(lwd=0.5, width = 0.75) +
  geom_boxplot(color = "black", fill = 'white', width = .15, 
               outlier.alpha = 0, lwd=0.5) +
  theme(
    plot.title = element_text(size=20, hjust = 0.5),
    panel.background = element_blank(),
    axis.text.y=element_text(size=20),
    axis.text.x=element_text(size=20),
    axis.ticks.x=element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size=20),
    strip.text.x = element_text(size=0),
    legend.position = 'right',
    legend.title = element_blank(),
    legend.text = element_text(size=20),
    text = element_text(family = 'Calibri')) +
  labs(title = "High-grade pattern\nannotated spots only\n",
       x = "\nVI Status", y = "Enrichment score\n") +
  ylim(min(fig2D_2$score),max(fig2D_2$score)+2) +
  facet_grid(cols = vars(cluster),
             scales = "free_y", space = "free_y") +
  geom_bracket(
    xmin = 'VI-', xmax = "VI+", y.position = max(fig2D_2$score)+1,
    label = paste('Anova, p =',round(mixed_effects_stats$Pr..Chisq.[1], digits = 3)), 
    tip.length = c(0.02, 0.02), vjust = -0.25,
    size = 0.5,
    label.size = 6,
    inherit.aes = FALSE,
    family = 'Calibri') +
  scale_x_discrete(labels = c("VI-\n(n=8 samples)","VI+\n(n=7 samples)"))


# color individual  facets

g <- ggplot_gtable(ggplot_build(p))

strip_both <- which(grepl('strip-', g$layout$name))

fills <- c(unname(VI_cluster_colors))

k <- 2

for (i in strip_both) {
  
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  
  k <- k+1
  
}

# Save the plot using tiff device with high resolution

tiff(filename = here('figures','fig2D_2.tiff'),
     width = width, height = height, units = "in", res = dpi)

grid.draw(g)

dev.off()


#############
# FIGURE 2D_3
#############

fig2D_3 <- data_long %>% filter(cluster == 'VI_gene_cluster_3')

addWorksheet(source_data, "Figure 2D-3")

writeData(source_data, sheet = "Figure 2D-3", 
          x = "Figure 2D-3", 
          startCol =  1, startRow = 1)

writeData(source_data, sheet = 'Figure 2D-3',
          x = fig2D_3, startCol = 1, startRow = 3)

data <- data.frame(y = fig2D_3$score, x = fig2D_3$`VI Status`, 
                   sample = fig2D_3$`sample id`, pathology = fig2D_3$pathology)

m <- lmer(y ~ 1 + x + (1 | sample) + pathology, data = data)


# type II anova 
Anova(m, type = 2)

mixed_effects_stats <- data.frame(Anova(m, type = 2, test.statistic = 'Chisq'))

mixed_effects_stats <- mixed_effects_stats %>% select(Pr..Chisq.)

p <- ggplot(data = fig2D_3, aes(x = `VI Status`, y = score, 
                                fill = `VI Status`)) +
  theme_classic() + 
  coord_cartesian(clip = "off") +
  scale_fill_manual(values = c(unname(use_colors_VI))) +
  geom_violin(lwd=0.5, width = 0.75) +
  geom_boxplot(color = "black", fill = 'white', width = .15, 
               outlier.alpha = 0, lwd=0.5) +
  theme(
    plot.title = element_text(size=20, hjust = 0.5),
    panel.background = element_blank(),
    axis.text.y=element_text(size=20),
    axis.text.x=element_text(size=20),
    axis.ticks.x=element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size=20),
    strip.text.x = element_text(size=0),
    legend.position = 'right',
    legend.title = element_blank(),
    legend.text = element_text(size=20),
    text = element_text(family = 'Calibri')) +
  labs(title = "High-grade pattern\nannotated spots only\n",
       x = "\nVI Status", y = "Enrichment score\n") +
  ylim(min(fig2D_3$score),max(fig2D_3$score)+2) +
  facet_grid(cols = vars(cluster),
             scales = "free_y", space = "free_y") +
  geom_bracket(
    xmin = 'VI-', xmax = "VI+", y.position = max(fig2D_3$score)+1,
    label = paste('Anova, p =',round(mixed_effects_stats$Pr..Chisq.[1], digits = 3)), 
    tip.length = c(0.02, 0.02), vjust = -0.25,
    size = 0.5,
    label.size = 6,
    inherit.aes = FALSE,
    family = 'Calibri') +
  scale_x_discrete(labels = c("VI-\n(n=8 samples)","VI+\n(n=7 samples)"))


# color individual  facets

g <- ggplot_gtable(ggplot_build(p))

strip_both <- which(grepl('strip-', g$layout$name))

fills <- c(unname(VI_cluster_colors))

k <- 3

for (i in strip_both) {
  
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  
  k <- k+1
  
}

# Save the plot using tiff device with high resolution

tiff(filename = here('figures','fig2D_3.tiff'),
     width = width, height = height, units = "in", res = dpi)

grid.draw(g)

dev.off()


#############
# FIGURE 2D_4
#############

fig2D_4 <- data_long %>% filter(cluster == 'VI_gene_cluster_4')

addWorksheet(source_data, "Figure 2D-4")

writeData(source_data, sheet = "Figure 2D-4", 
          x = "Figure 2D-4", 
          startCol =  1, startRow = 1)

writeData(source_data, sheet = 'Figure 2D-4',
          x = fig2D_4, startCol = 1, startRow = 3)

data <- data.frame(y = fig2D_4$score, x = fig2D_4$`VI Status`, 
                   sample = fig2D_4$`sample id`, pathology = fig2D_4$pathology)

m <- lmer(y ~ 1 + x + (1 | sample) + pathology, data = data)


# type II anova 
Anova(m, type = 2)

mixed_effects_stats <- data.frame(Anova(m, type = 2, test.statistic = 'Chisq'))

mixed_effects_stats <- mixed_effects_stats %>% select(Pr..Chisq.)

pval <- grobTree(textGrob('Anova, p =    ',
                          x=0.4,  y=0.9, 
                          gp=gpar(fontsize=17, 
                                  fontfamily = 'Calibri')),
                 textGrob(eval(bquote(scientific_10(round(mixed_effects_stats$Pr..Chisq.[1], digits = 6)))),
                          x=0.635, y=0.91,
                          gp=gpar(fontsize=17, 
                                  fontfamily = 'Calibri')))

p <- ggplot(data = fig2D_4, aes(x = `VI Status`, y = score, 
                                fill = `VI Status`)) +
  theme_classic() + 
  coord_cartesian(clip = "off") +
  scale_fill_manual(values = c(unname(use_colors_VI))) +
  geom_violin(lwd=0.5, width = 0.75) +
  geom_boxplot(color = "black", fill = 'white', width = .15, 
               outlier.alpha = 0, lwd=0.5) +
  theme(
    plot.title = element_text(size=20, hjust = 0.5),
    panel.background = element_blank(),
    axis.text.y=element_text(size=20),
    axis.text.x=element_text(size=20),
    axis.ticks.x=element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size=20),
    strip.text.x = element_text(size=0),
    legend.position = 'right',
    legend.title = element_blank(),
    legend.text = element_text(size=20),
    text = element_text(family = 'Calibri')) +
  labs(title = "High-grade pattern\nannotated spots only\n",
       x = "\nVI Status", y = "Enrichment score\n") +
  ylim(min(fig2D_4$score),max(fig2D_4$score)+2) +
  facet_grid(cols = vars(cluster),
             scales = "free_y", space = "free_y") +
  geom_bracket(
    xmin = 'VI-', xmax = "VI+", y.position = max(fig2D_4$score)+1,
    label = paste(''), 
    tip.length = c(0.02, 0.02), vjust = -0.25,
    size = 0.5,
    label.size = 6,
    inherit.aes = FALSE,
    family = 'Calibri') +
  annotation_custom(pval) +
  scale_x_discrete(labels = c("VI-\n(n=8 samples)","VI+\n(n=7 samples)"))

# color individual  facets

g <- ggplot_gtable(ggplot_build(p))

strip_both <- which(grepl('strip-', g$layout$name))

fills <- c(unname(VI_cluster_colors))

k <- 4

for (i in strip_both) {
  
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  
  k <- k+1
  
}


plot_dimensions <- dev.size("in")

# Calculate aspect ratio

aspect_ratio <- plot_dimensions[1] / plot_dimensions[2]

aspect_ratio <- 1.2

# Define the resolution

dpi <- 600
width <- 5.75
height <- 5.75 / aspect_ratio


# Save the plot using tiff device with high resolution
tiff(filename = here('figures','fig2D_4.tiff'),
     width = width, height = height, units = "in", res = dpi)

grid.draw(g)

dev.off()





################################################
# INDIVIDUAL SAMPLE NORMALIZATION AND CLUSTERING
################################################

# Log-normal based normalization for spots

visium_seurat_objects <- lapply(X = visium_seurat_objects, FUN = function(x) {
  
  x <- NormalizeData(x, verbose = TRUE, assay = "Spatial")
  
  
})

# save data checkpoint 

# saveRDS(visium_seurat_objects, file = here('data','visium_seurat_objects_lognormal.rds'))

visium_seurat_objects <- readRDS(here('data','visium_seurat_objects_lognormal.rds'))


# Dimensionality reduction and clustering

visium_seurat_objects <- lapply(X = visium_seurat_objects, FUN = function(x) {
  
  x <- FindVariableFeatures(x)
  
  x <- ScaleData(x, assay = "Spatial", verbose = TRUE)
  
  x <- RunPCA(x, assay = "Spatial", verbose = TRUE)
  
  ElbowPlot(x, ndims = 30)
  
  x <- FindNeighbors(x, reduction = "pca", dims = 1:30)
  
  # scale cluster resolution parameter with number of spots
  
  x <- FindClusters(x, verbose = TRUE, 
                    resolution = round(((log(ncol(x))) / 10), digits = 1))
  
  x <- RunUMAP(x, reduction = "pca", dims = 1:30)
  
})

# save data checkpoint

# saveRDS(visium_seurat_objects, file = here('data','visium_seurat_objects_lognormal_clustered.rds'))

visium_seurat_objects <- readRDS(here('data','visium_seurat_objects_lognormal_clustered.rds'))



################################################################
# INDIVIDUAL SAMPLE SIGNATURE ENRICHMENT WITHIN SPATIAL CLUSTERS
################################################################

visium_seurat_objects <- lapply(X = visium_seurat_objects, FUN = function(x) {
  
  x <- AddModuleScore(x, features = VI_genesets, name = '_score')
  
})

for (i in 1:length(visium_seurat_objects)){
  
  names(visium_seurat_objects[[i]]@meta.data)[grep("_score", names(visium_seurat_objects[[i]]@meta.data))] <- names(VI_genesets)
  
}

# spatially weighted regression analysis of VI clusters

######################
# EXTENDED DATA FIG 4F
######################

for (i in 1:length(visium_seurat_objects)){
  
  # get image coordinates
  
  cat('Getting coordinates', i, 'of', length(visium_seurat_objects),'\n')
  
  visium_coordinates <- visium_seurat_objects[[i]]@images$slice1@coordinates
  
  # add geneset scores for analysis
  
  visium_coordinates <- visium_coordinates %>%
    bind_cols(FetchData(visium_seurat_objects[[i]], c('VI_gene_cluster_1', 
                                                      'VI_gene_cluster_2', 
                                                      'VI_gene_cluster_3', 
                                                      'VI_gene_cluster_4')))
  
  # convert to spatial points dataframe
  
  coordinates(visium_coordinates) <- c('row','col')
  
  # perform spatially weighted regression analysis
  
  gws_results <- gwss(visium_coordinates, vars = c(names(visium_coordinates)[grep("VI_gene_cluster", names(visium_coordinates))]), bw = 5)
  
  gws_results <- data.frame(gws_results$SDF)
  
  # extract spearman correlation summary stats
  
  gws_results <- gws_results %>% 
    select(starts_with('Spearman'))
  
  visium_seurat_objects[[i]]@meta.data <- cbind(visium_seurat_objects[[i]]@meta.data, gws_results)
  
  
}


visium_gw_cor_means_all <- list()

for (i in 1:length(visium_seurat_objects)){
  
  # take mean of each spatially weighted correlation pair
  
  visium_gw_cor_means <- visium_seurat_objects[[i]]@meta.data %>%
    select(starts_with("Spearman")) %>%
    dplyr::summarise(across(everything(), ~mean(., na.rm = TRUE)))
  
  names(visium_gw_cor_means) <- paste(unique(visium_seurat_objects[[i]]$orig.ident), names(visium_gw_cor_means), sep = '_')
  
  visium_gw_cor_means_all[[i]] <- data.frame(t(visium_gw_cor_means))
  
}

# plot global correlation

visium_gw_cor_means_all_df <- bind_rows(visium_gw_cor_means_all, .id = "orig.ident")

# remove sample 6

visium_gw_cor_means_all_df <- visium_gw_cor_means_all_df %>%
  filter(orig.ident != 6)

geo_stRNAseq_metadata$Library <- as.integer(sub(".*_", "", geo_stRNAseq_metadata$`*library name`))

visium_gw_cor_means_all_df <- visium_gw_cor_means_all_df %>% 
  rownames_to_column(var = 'corr') %>%
  mutate(clusters = sapply(strsplit(corr,'rho_'),"[",2)) %>%
  dplyr::rename(Library = orig.ident) %>%
  mutate(Library = as.integer(Library)) %>%
  inner_join(geo_stRNAseq_metadata, by = 'Library')


# calculate mean over all samples

visium_gw_cor_means_all_df_mean <- visium_gw_cor_means_all_df %>%
  group_by(clusters) %>%
  dplyr::summarise(mean_clusters = mean(t.visium_gw_cor_means.))


# Extract the gene cluster names

unique_clusters <- unique(unlist(strsplit(visium_gw_cor_means_all_df_mean$clusters,"\\.")))

# Create a square matrix to store the correlations

num_clusters <- length(unique_clusters)

cor_matrix <- matrix(0, nrow = num_clusters, ncol = num_clusters,
                     dimnames = list(unique_clusters, unique_clusters))

# Fill in the correlation matrix

for (i in seq_len(length(visium_gw_cor_means_all_df_mean$mean_clusters))) {
  
  row_cluster <- unique(unlist(strsplit(visium_gw_cor_means_all_df_mean$clusters[i],"\\.")))[1]
  
  col_cluster <- unique(unlist(strsplit(visium_gw_cor_means_all_df_mean$clusters[i],"\\.")))[2]
  
  cor_value <- visium_gw_cor_means_all_df_mean$mean_clusters[i]
  
  cor_matrix[row_cluster, col_cluster] <- cor_value
  
  cor_matrix[col_cluster, row_cluster] <- cor_value
  
  print(i)
  
}

cor_matrix[cor_matrix == 0] <- 1

# plot mean spatially weighted correlation across samples

bluered <- colorRampPalette(c("#007FFF","white","red"))(length(seq(-1, 1, by = 0.01)))

rownames(cor_matrix) <- gsub('_',' ',rownames(cor_matrix))

colnames(cor_matrix) <- gsub('_',' ',colnames(cor_matrix))

efig4F <- cor_matrix

addWorksheet(source_data, "E. Figure 4F")

writeData(source_data, sheet = "E. Figure 4F", 
          x = "E. Figure 4F", 
          startCol =  1, startRow = 1)

writeData(source_data, sheet = 'E. Figure 4F',
          x = efig4F, startCol = 1, startRow = 3)

efig4F_gg <- pheatmap::pheatmap(efig4F,
                                color = diverging_hcl(length(seq(-1, 1, by = 0.01)), 
                                                      palette = 'Blue-Red 3'), 
                                border_color = 'white', scale = "none", show_rownames = TRUE,
                                show_colnames = TRUE, cluster_rows = F, cluster_cols = F,
                                clustering_distance_rows = 'euclidean', 
                                clustering_distance_cols = 'euclidean',
                                legend = TRUE, annotation_legend = FALSE,
                                angle_col = 0,
                                fontsize = 16,
                                clustering_method = "ward.D2", 
                                main = "Mean Spatially Weighted\nCorrelation (n=15)\n",
                                breaks = seq(-1,1, by = 0.01),
                                fontfamily = 'Calibri')


dev.off()

efig4F_gg$gtable$grobs[[1]]$gp = gpar(fontface = 'plain', fontsize = 25)

efig4F_gg$gtable$grobs[[5]]$gp = gpar(fontface = 'plain')

efig4F_gg


plot_dimensions <- dev.size("in")

aspect_ratio <- plot_dimensions[1] / plot_dimensions[2]

aspect_ratio <- 1.123786

desired_height <- 7.75

desired_width <- desired_height * aspect_ratio

TIFF_dpi <- 600

ggsave(here('figures','efig4F.tiff'),
       plot = efig4F_gg,
       device = 'tiff',
       width = desired_width,
       height = desired_height,
       units = 'in', dpi = TIFF_dpi)



# convert seurat objects to spata objects for plotting pathology annotations
# remove unannotated spots for plotting purposes

visium_spata_objects_path <- list()

visium_spata_objects_path <- lapply(X = visium_seurat_objects, FUN = function(x) {
  
  Idents(x) <- x$pathology
  
  x <- subset(x, idents = "", invert = TRUE)
  
  x <- SPATA2::asSPATA2(
    object = x,
    sample_name = x$orig.ident[1],
    img_name = "slice1", 
    platform = "VisiumSmall",
    assay_name = "Spatial",
    assay_modality = 'gene'
  )
  
})


# supervised pathology labels only

for (i in 1:length(visium_spata_objects_path)){
  
  x <- names(visium_spata_objects_path)[i]
  
  visium_spata_objects_path[[i]]@meta_obs$pathology <- as.factor(visium_spata_objects_path[[i]]@meta_obs$pathology)
  
  # unannotated spots
  
  levels(visium_spata_objects_path[[i]]@meta_obs$pathology)[which(levels(visium_spata_objects_path[[i]]@meta_obs$pathology) == "")] <- 'Unannotated'
  
  # collapse levels (could also create a 'Other' level)
  
  visium_spata_objects_path[[i]]@meta_obs$pathology <- fct_collapse(visium_spata_objects_path[[i]]@meta_obs$pathology, Vessel = c("Hyal Vessel","Hyper Vessel","Normal Vessel"),
                                                                      `Normal Lung` = c("Normal Lung","Normal Pleura","Pulmonary Macrophage"),
                                                                      `Lymphoid Aggregate` = c("Lymphoid Aggregate","Lymphoid Aggregate (VI)"))
  
  levels(visium_spata_objects_path[[i]]@meta_obs$pathology)[which(levels(visium_spata_objects_path[[i]]@meta_obs$pathology) == "VI")] <- "VI Focus"
  
  levels(visium_spata_objects_path[[i]]@meta_obs$pathology)[which(levels(visium_spata_objects_path[[i]]@meta_obs$pathology) == "VPI")] <- "VPI Focus"
  
}

# color scheme for pathology annotations

use_colors <- c(`Abnormal Pleura` = '#64638E',
                `Lymphoid Aggregate` = '#94539E',
                Cartilage = '#B85D6F',
                `Collapsed Lung` = '#5C9A5C',
                Necrosis = '#FFB415',
                Acinar = '#D97187',
                Stroma = '#BB614F',
                Cribriform = '#D0A62D',
                `Desmoplastic Stroma`= '#AF6729',
                Vessel = '#F781BF',
                Lepidic = '#49A75B',
                Micropapillary = '#DE6F33',
                `Normal Bronchus` = '#419583',
                `Normal Lung` = '#3983AC',
                Papillary = '#FF8502',
                Solid = '#FFE428',
                STAS = '#CB8CAD',
                `VI Focus` = '#E41A1C',
                `VPI Focus` = '#A43E55')


# H&E with pathology as spots

##########
# FIG 2C_1
##########

fig2C_1_gg <- plotSurface(
  object = visium_spata_objects_path[[4]], 
  color_by = "pathology", 
  pt_clrp = NULL,
  display_image = TRUE,
  pt_alpha = 1,
  pt_size = 2.4) + 
  scale_color_manual(values = use_colors) + 
  theme(legend.position = "none",
        legend.title = element_text(size=20),
        legend.text = element_text(size=20),
        plot.title = element_text(hjust = 0.5, size = 38),
        text = element_text(family = 'Calibri')) +
  labs(color = "Pathology Features", title = 'Pathology annotations\n')

fig2C_1_gg$labels$y <- c("VI+ LUAD\n(Sample 4)\n")

fig2C_1_gg$theme$axis.title.y <- element_text(size=40, angle = 90)

fig2C_1_gg

# Get the current plot dimensions in inches

plot_dimensions <- dev.size("in")

# Calculate aspect ratio

aspect_ratio <- plot_dimensions[1] / plot_dimensions[2]

aspect_ratio <- 1.798058

# Define a reasonable height (in inches)

desired_height <- 6

# Calculate corresponding width to maintain aspect ratio

desired_width <- desired_height * aspect_ratio

# Define the resolution

TIFF_dpi <- 600

# Save the plot with adjusted dimensions

ggsave(here('figures','fig2C_1.tiff'),
       plot = fig2C_1_gg,
       device = 'tiff',
       width = desired_width,
       height = desired_height,
       units = 'in', dpi = TIFF_dpi)



##########
# FIG 2C_2
##########


fig2C_2_gg <- plotSurface(
  object = visium_spata_objects_path[[11]], 
  color_by = "pathology", 
  pt_clrp = NULL,
  display_image = TRUE,
  pt_alpha = 1,
  pt_size = 2.4) + 
  scale_color_manual(values = use_colors) + 
  theme(legend.position = "none",
        legend.title = element_text(size=20),
        legend.text = element_text(size=20),
        plot.title = element_text(hjust = 0.5, size = 38),
        text = element_text(family = 'Calibri')) +
  labs(color = "Pathology Features", title = '   \n')

fig2C_2_gg$labels$y <- c("VI+ LUAD\n(Sample 2B)\n")

fig2C_2_gg$theme$axis.title.y <- element_text(size=40, angle = 90)

fig2C_2_gg


# Get the current plot dimensions in inches

plot_dimensions <- dev.size("in")

# Calculate aspect ratio

aspect_ratio <- plot_dimensions[1] / plot_dimensions[2]

aspect_ratio <- 1.257796

# Define a reasonable height (in inches)

desired_height <- 6

# Calculate corresponding width to maintain aspect ratio

desired_width <- desired_height * aspect_ratio

# Save the plot with adjusted dimensions

ggsave(here('figures','fig2C_2.tiff'),
       plot = fig2C_2_gg,
       device = 'tiff',
       width = desired_width,
       height = desired_height,
       units = 'in', dpi = TIFF_dpi)



##########
# FIG 2C_3
##########


fig2C_3_gg <- plotSurface(
  object = visium_spata_objects_path[[15]], 
  color_by = "pathology", 
  pt_clrp = NULL,
  display_image = TRUE,
  pt_alpha = 1,
  pt_size = 2.6) + 
  scale_color_manual(values = use_colors) + 
  theme(legend.position = "none",
        legend.title = element_text(size=20),
        legend.text = element_text(size=20),
        plot.title = element_text(hjust = 0.5, size = 38),
        text = element_text(family = 'Calibri')) +
  labs(color = "Pathology Features", title = '   \n')

fig2C_3_gg$labels$y <- c("VI- LUAD\n(Sample 14)\n")

fig2C_3_gg$theme$axis.title.y <- element_text(size=40, angle = 90)

fig2C_3_gg


# Get the current plot dimensions in inches

plot_dimensions <- dev.size("in")

# Calculate aspect ratio

aspect_ratio <- plot_dimensions[1] / plot_dimensions[2]

aspect_ratio <- 1.798058

# Define a reasonable height (in inches)

desired_height <- 6

# Calculate corresponding width to maintain aspect ratio

desired_width <- desired_height * aspect_ratio

# Save the plot with adjusted dimensions
ggsave(here('figures','fig2C_3.tiff'),
       plot = fig2C_3_gg,
       device = 'tiff',
       width = desired_width,
       height = desired_height,
       units = 'in', dpi = TIFF_dpi)


# convert seurat objects to spata objects for plotting clusters or features

visium_spata_objects <- list()

visium_spata_objects <- lapply(X = visium_seurat_objects, FUN = function(x) {
  
  x <- SPATA2::asSPATA2(
    object = x,
    sample_name = x$orig.ident[1],
    img_name = "slice1", 
    platform = "VisiumSmall",
    assay_name = "Spatial",
    assay_modality = 'gene'
  )
  
})

sample_outline_VI <-
  ggpLayerGroupOutline(
    object = visium_spata_objects_path[[4]],
    grouping = "pathology",
    group_subset = c("VI Focus"),
    line_color = 'white',
    line_size = 1)



# plot VI cluster gene enrichment scores variables on spatial plot

##########
# FIG 2C_4
##########

custom_surface_plot <- function(spata_object, cluster, pt_size) {
  
  lowest_min <- min(sapply(lapply(X = visium_seurat_objects, 
                                  FUN = function(x) {range(x@meta.data[,cluster])}), function(x) x[1]))
  
  highest_max <- max(sapply(lapply(X = visium_seurat_objects, 
                                   FUN = function(x) {range(x@meta.data[,cluster])}), function(x) x[2]))
  
  cluster_sub <- gsub("[^0-9]", "", cluster)
  
  p <- SPATA2::plotSurface(
    object = spata_object, 
    color_by = cluster, 
    pt_clrsp = "inferno",
    smooth = FALSE, 
    display_image = FALSE,
    pt_size = 2.4,
    limits = c(lowest_min,highest_max)
  )  +
    sample_outline_VI +
    theme(legend.position = "bottom",
          legend.title=element_text(size=24),
          legend.text = element_text(size = 16),
          legend.key.size = unit(1, 'cm'),
          plot.title = element_text(hjust = 0.5, size = 38),
          text = element_text(family = 'Calibri')) +
    labs(color = "Enrichment score", title = paste('Cluster ',cluster_sub,'\n',sep=''))
  
  return(p)
  
}

fig2C_4_gg <- custom_surface_plot(spata_object = visium_spata_objects[[4]], 
                                  cluster = 'VI_gene_cluster_1',
                                  pt_size = 2.4)

fig2C_4_gg

# Save the plot with adjusted dimensions
ggsave(here('figures','fig2C_4.tiff'),
       plot = fig2C_4_gg,
       device = 'tiff',
       width = desired_width,
       height = desired_height,
       units = 'in', dpi = TIFF_dpi)

fig2C_5_gg <- custom_surface_plot(spata_object = visium_spata_objects[[4]], 
                                  cluster = 'VI_gene_cluster_2',
                                  pt_size = 2.4)

fig2C_5_gg

# Save the plot with adjusted dimensions
ggsave(here('figures','fig2C_5.tiff'),
       plot = fig2C_5_gg,
       device = 'tiff',
       width = desired_width,
       height = desired_height,
       units = 'in', dpi = TIFF_dpi)

fig2C_6_gg <- custom_surface_plot(spata_object = visium_spata_objects[[4]], 
                                  cluster = 'VI_gene_cluster_3',
                                  pt_size = 2.4)

fig2C_6_gg

# Save the plot with adjusted dimensions
ggsave(here('figures','fig2C_6.tiff'),
       plot = fig2C_6_gg,
       device = 'tiff',
       width = desired_width,
       height = desired_height,
       units = 'in', dpi = TIFF_dpi)

fig2C_7_gg <- custom_surface_plot(spata_object = visium_spata_objects[[4]], 
                                  cluster = 'VI_gene_cluster_4',
                                  pt_size = 2.4)

fig2C_7_gg

# Save the plot with adjusted dimensions
ggsave(here('figures','fig2C_7.tiff'),
       plot = fig2C_7_gg,
       device = 'tiff',
       width = desired_width,
       height = desired_height,
       units = 'in', dpi = TIFF_dpi)


custom_surface_plot <- function(spata_object, cluster, pt_size) {
  
  lowest_min <- min(sapply(lapply(X = visium_seurat_objects, FUN = function(x) {range(x@meta.data[,cluster])}), function(x) x[1]))
  
  highest_max <- max(sapply(lapply(X = visium_seurat_objects, FUN = function(x) {range(x@meta.data[,cluster])}), function(x) x[2]))
  
  cluster_sub <- gsub("[^0-9]", "", cluster)
  
  p <- SPATA2::plotSurface(
    object = spata_object, 
    color_by = cluster, 
    pt_clrsp = "inferno",
    smooth = FALSE, 
    display_image = FALSE,
    pt_size = 2.4,
    limits = c(lowest_min,highest_max)
  )  +
    theme(legend.position = "bottom",
          legend.title=element_text(size=24),
          legend.text = element_text(size = 16),
          legend.key.size = unit(1, 'cm'),
          plot.title = element_text(hjust = 0.5, size = 38),
          text = element_text(family = 'Calibri')) +
    labs(color = "Enrichment score", title = ' \n')
  
  return(p)
  
}


fig2C_8_gg <- custom_surface_plot(spata_object = visium_spata_objects[[11]], 
                                  cluster = 'VI_gene_cluster_1',
                                  pt_size = 2.4)

fig2C_8_gg

# Save the plot with adjusted dimensions

ggsave(here('figures','fig2C_8.tiff'),
       plot = fig2C_8_gg,
       device = 'tiff',
       width = desired_width,
       height = desired_height,
       units = 'in', dpi = TIFF_dpi)

fig2C_9_gg <- custom_surface_plot(spata_object = visium_spata_objects[[11]], 
                                  cluster = 'VI_gene_cluster_2',
                                  pt_size = 2.4)

fig2C_9_gg

# Save the plot with adjusted dimensions

ggsave(here('figures','fig2C_9.tiff'),
       plot = fig2C_9_gg,
       device = 'tiff',
       width = desired_width,
       height = desired_height,
       units = 'in', dpi = TIFF_dpi)

fig2C_10_gg <- custom_surface_plot(spata_object = visium_spata_objects[[11]], 
                                   cluster = 'VI_gene_cluster_3',
                                   pt_size = 2.4)

fig2C_10_gg

# Save the plot with adjusted dimensions

ggsave(here('figures','fig2C_10.tiff'),
       plot = fig2C_10_gg,
       device = 'tiff',
       width = desired_width,
       height = desired_height,
       units = 'in', dpi = TIFF_dpi)

fig2C_11_gg <- custom_surface_plot(spata_object = visium_spata_objects[[11]], 
                                   cluster = 'VI_gene_cluster_4',
                                   pt_size = 2.4)

fig2C_11_gg

# Save the plot with adjusted dimensions

ggsave(here('figures','fig2C_11.tiff'),
       plot = fig2C_11_gg,
       device = 'tiff',
       width = desired_width,
       height = desired_height,
       units = 'in', dpi = TIFF_dpi)


fig2C_12_gg <- custom_surface_plot(spata_object = visium_spata_objects[[15]], 
                                   cluster = 'VI_gene_cluster_1',
                                   pt_size = 2.6)

fig2C_12_gg

# Save the plot with adjusted dimensions

ggsave(here('figures','fig2C_12.tiff'),
       plot = fig2C_12_gg,
       device = 'tiff',
       width = desired_width,
       height = desired_height,
       units = 'in', dpi = TIFF_dpi)

fig2C_13_gg <- custom_surface_plot(spata_object = visium_spata_objects[[15]], 
                                   cluster = 'VI_gene_cluster_2',
                                   pt_size = 2.4)

fig2C_13_gg

# Save the plot with adjusted dimensions

ggsave(here('figures','fig2C_13.tiff'),
       plot = fig2C_13_gg,
       device = 'tiff',
       width = desired_width,
       height = desired_height,
       units = 'in', dpi = TIFF_dpi)

fig2C_14_gg <- custom_surface_plot(spata_object = visium_spata_objects[[15]], 
                                   cluster = 'VI_gene_cluster_3',
                                   pt_size = 2.4)

fig2C_14_gg

# Save the plot with adjusted dimensions

ggsave(here('figures','fig2C_14.tiff'),
       plot = fig2C_14_gg,
       device = 'tiff',
       width = desired_width,
       height = desired_height,
       units = 'in', dpi = TIFF_dpi)

fig2C_15_gg <- custom_surface_plot(spata_object = visium_spata_objects[[15]], 
                                   cluster = 'VI_gene_cluster_4',
                                   pt_size = 2.4)

fig2C_15_gg

# Save the plot with adjusted dimensions

ggsave(here('figures','fig2C_15.tiff'),
       plot = fig2C_15_gg,
       device = 'tiff',
       width = desired_width,
       height = desired_height,
       units = 'in', dpi = TIFF_dpi)



# divide spots into distal VI+ and proximal VI+

focus_circles <- list()

for (i in c(1,4,9)){
  
  VI_barcodes <- getMetaDf(visium_spata_objects_path[[i]]) %>%
    filter(pathology == "VI Focus") %>%
    pull(barcodes)
  
  visium_spata_objects_path[[i]] <- 
    identifyTissueOutline(
      visium_spata_objects_path[[i]],
      method = "obs",
      eps    = recDbscanEps(visium_spata_objects_path[[i]]),
      minPts = recDbscanMinPts(visium_spata_objects_path[[i]])
    )
  
  visium_spata_objects_path[[i]] <- createGroupAnnotations(
    object = visium_spata_objects_path[[i]],
    grouping = "pathology",   
    group = "VI Focus",  
    id = "VI_Focus",     
    method_outline= "concaveman",
    min_size = 1,
    overwrite = TRUE
  ) 
  
  getSpatAnnIds(visium_spata_objects_path[[i]])

  bin_df <- getCoordsDfSA(
    object       = visium_spata_objects_path[[i]],
    ids          = "VI_Focus_1",   
    distance     = "1.0mm",      
    resolution   = "1.0mm",     
    angle_span   = c(0,360),    
    n_bins_angle = 4,            
    format       = "wide"       
  )
  
  focus_circles[[i]] <- bin_df %>%
    select(barcodes, 
           bins_circle = bins_dist)

  
}


visium_seurat_objects <- readRDS(here('data','visium_seurat_objects_lognormal_clustered.rds'))

for (i in c(1,4,9)){
  
  seurat_obj <- visium_seurat_objects[[i]]
  
  seurat_obj@meta.data$barcodes <- rownames(seurat_obj@meta.data)
  
  merged <- seurat_obj@meta.data %>%
    as_tibble() %>%
    left_join(focus_circles[[i]], by = "barcodes") %>%
    mutate(
      `VI focus` = if_else(
        bins_circle == "Outside", 
        "Distal (>1 mm)", 
        "Proximal (<1 mm)"
      ) %>% 
        factor(levels = c("Proximal (<1 mm)", "Distal (>1 mm)"))
    ) %>%
    column_to_rownames("barcodes")
  
  seurat_obj@meta.data <- merged
  
  visium_seurat_objects[[i]] <- seurat_obj
  
}

# saveRDS(visium_seurat_objects, file = here('data','visium_seurat_objects_lognormal_clustered_annotated_distal.rds'))

# load for reproducibilty (issues with latest version of SPATA3)

visium_seurat_objects <- readRDS(here('data','visium_seurat_objects_lognormal_clustered_annotated_distal.rds'))

visium_spata_objects <- list()

visium_spata_objects <- lapply(X = visium_seurat_objects, FUN = function(x) {
  
  x <- SPATA2::asSPATA2(
    object = x,
    sample_name = x$orig.ident[1],
    img_name = "slice1", 
    platform = "VisiumSmall",
    assay_name = "Spatial",
    assay_modality = 'gene'
  )
  
})

for (i in 1:length(visium_spata_objects)){
  
  x <- names(visium_spata_objects)[i]
  
  visium_spata_objects[[i]]@meta_obs$pathology <- as.factor(visium_spata_objects[[i]]@meta_obs$pathology)
  
  # unannotated spots
  
  levels(visium_spata_objects[[i]]@meta_obs$pathology)[which(levels(visium_spata_objects[[i]]@meta_obs$pathology) == "")] <- 'Unannotated'
  
  # collapse levels (could also create a 'Other' level)
  
  visium_spata_objects[[i]]@meta_obs$pathology <- fct_collapse(visium_spata_objects[[i]]@meta_obs$pathology, Vessel = c("Hyal Vessel","Hyper Vessel","Normal Vessel"),
                                                                    `Normal Lung` = c("Normal Lung","Normal Pleura","Pulmonary Macrophage"),
                                                                    `Lymphoid Aggregate` = c("Lymphoid Aggregate","Lymphoid Aggregate (VI)"))
  
  levels(visium_spata_objects[[i]]@meta_obs$pathology)[which(levels(visium_spata_objects[[i]]@meta_obs$pathology) == "VI")] <- "VI Focus"
  
  levels(visium_spata_objects[[i]]@meta_obs$pathology)[which(levels(visium_spata_objects[[i]]@meta_obs$pathology) == "VPI")] <- "VPI Focus"
  
}


# sample 1A

# sample_outline_VI <- 
#   ggpLayerGroupOutline(
#     object = visium_spata_objects[[1]],
#     grouping = "pathology",
#     group_subset = c("VI Focus"),
#     line_color = 'white',
#     line_size = 1)

visium_spata_objects[[1]]@meta_obs$VI.focus <- factor(visium_spata_objects[[1]]@meta_obs$VI.focus, 
                                                                levels = c('Proximal (<1mm)','Distal (>1mm)'))

##########
# FIG 5A_1
##########

fig5A_1_gg <- plotSurface(
  object = visium_spata_objects[[1]], 
  color_by = "VI.focus", 
  pt_clrp = NULL,
  display_image = FALSE,
  pt_alpha = 1,
  pt_size = 2.4,
) + 
  scale_color_manual(values = c("#FF0054","#9E0059")) + 
  #sample_outline_VI +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size=20),
        plot.title = element_text(hjust = 0.5, size = 26),
        text = element_text(family = 'Calibri')) +
  ggtitle('Sample 1A VI focus') +
  guides(color=guide_legend(title="VI focus",
                            override.aes = list(size=8)))

fig5A_1_gg

# Save the plot with adjusted dimensions

ggsave(here('figures','fig5A_1.tiff'),
       plot = fig5A_1_gg,
       device = 'tiff',
       width = desired_width,
       height = desired_height,
       units = 'in', dpi = TIFF_dpi)




# sample 4 

sample_outline_VI <- 
  ggpLayerGroupOutline(
    object = visium_spata_objects[[4]],
    grouping = "pathology",
    group_subset = c("VI Focus"),
    line_color = 'white',
    line_size = 1)

visium_spata_objects[[4]]@meta_obs$VI.focus <- factor(visium_spata_objects[[4]]@meta_obs$VI.focus, 
                                                           levels = c('Proximal (<1mm)','Distal (>1mm)'))

##########
# FIG 5A_2
##########

fig5A_2_gg <- plotSurface(
  object = visium_spata_objects[[4]], 
  color_by = "VI.focus", 
  pt_clrp = NULL,
  display_image = FALSE,
  pt_alpha = 1,
  pt_size = 2.4,
) + 
  scale_color_manual(values = c("#FF0054","#9E0059")) + 
  sample_outline_VI +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size=20),
        plot.title = element_text(hjust = 0.5, size = 26),
        text = element_text(family = 'Calibri')) +
  ggtitle('Sample 4 VI focus\n') +
  guides(color=guide_legend(title="VI focus",
                            override.aes = list(size=8)))

fig5A_2_gg

# Save the plot with adjusted dimensions

ggsave(here('figures','fig5A_2.tiff'),
       plot = fig5A_2_gg,
       device = 'tiff',
       width = desired_width,
       height = desired_height,
       units = 'in', dpi = TIFF_dpi)



# sample 9

sample_outline_VI <- 
  ggpLayerGroupOutline(
    object = visium_spata_objects[[9]],
    grouping = "pathology",
    group_subset = c("VI Focus"),
    line_color = 'white',
    line_size = 1)

visium_spata_objects[[9]]@meta_obs$VI.focus <- factor(visium_spata_objects[[9]]@meta_obs$VI.focus, 
                                                           levels = c('Proximal (<1mm)','Distal (>1mm)'))

##########
# FIG 5A_3
##########

fig5A_3_gg <- plotSurface(
  object = visium_spata_objects[[9]], 
  color_by = "VI.focus", 
  pt_clrp = NULL,
  display_image = FALSE,
  pt_alpha = 1,
  pt_size = 2.4,
) + 
  scale_color_manual(values = c("#FF0054","#9E0059")) + 
  sample_outline_VI +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size=20),
        plot.title = element_text(hjust = 0.5, size = 26),
        text = element_text(family = 'Calibri')) +
  ggtitle('Sample 9 VI focus\n') +
  guides(color=guide_legend(title="VI focus",
                            override.aes = list(size=8)))

fig5A_3_gg

# Save the plot with adjusted dimensions

ggsave(here('figures','fig5A_3.tiff'),
       plot = fig5A_3_gg,
       device = 'tiff',
       width = desired_width,
       height = desired_height,
       units = 'in', dpi = TIFF_dpi)



# load scRNA-seq reference atlas

# load Salcher extended NSCLC scRNAseq atlas as sce object

Salcher_lung_cancer_atlas <- readRDS(here('data','extended_atlas.rds'))

names(assays(Salcher_lung_cancer_atlas)) <- 'counts'

# subset to stage I LUAD only

Salcher_lung_cancer_atlas_LUAD <- Salcher_lung_cancer_atlas[,which(Salcher_lung_cancer_atlas$condition == 'LUAD' &
                                                                     Salcher_lung_cancer_atlas$tissue == 'lung' &
                                                                     Salcher_lung_cancer_atlas$uicc_stage == 'I' &
                                                                     Salcher_lung_cancer_atlas$cell_type_tumor != 'Tumor cells LUSC' &
                                                                     Salcher_lung_cancer_atlas$cell_type_tumor != 'Tumor cells LUSC mitotic' &
                                                                     Salcher_lung_cancer_atlas$origin != 'tumor_metastasis' &
                                                                     Salcher_lung_cancer_atlas$origin != 'nan')]

dim(Salcher_lung_cancer_atlas_LUAD) # 295,813 cells

length(unique(Salcher_lung_cancer_atlas_LUAD$sample)) # 124 samples


assayNames(Salcher_lung_cancer_atlas_LUAD) <- 'counts'


# convert from sce to seurat object

Salcher_lung_cancer_atlas_LUAD_seurat <- CreateSeuratObject(counts = counts(Salcher_lung_cancer_atlas_LUAD), meta.data = as.data.frame(colData(Salcher_lung_cancer_atlas_LUAD)))

Salcher_lung_cancer_atlas_LUAD_seurat <- SetAssayData(object = Salcher_lung_cancer_atlas_LUAD_seurat, slot = "data", new.data = counts(Salcher_lung_cancer_atlas_LUAD))

Salcher_lung_cancer_atlas_LUAD_seurat <- RenameAssays(object = Salcher_lung_cancer_atlas_LUAD_seurat, RNA = 'RNA')



saveRDS(Salcher_lung_cancer_atlas_LUAD_seurat, file = here('data','Salcher_lung_cancer_atlas_LUAD.rds'))

Salcher_lung_cancer_atlas_LUAD_seurat <- readRDS(here('data','Salcher_lung_cancer_atlas_LUAD.rds'))


#####################################
# CELL TYPE ESTIMATION WITH CYTOSPACE
#####################################

# note: cytospace recommends estimating with seurat V3 with this script:

here('cytospace', 'get_cellfracs_seuratv3.R')

source(here('cytospace','cytospace','Prepare_input_files','generate_cytospace_from_seurat_object.R'))

# use non-normalized count matrix

DefaultAssay(Salcher_lung_cancer_atlas_LUAD_seurat) <- "RNA"

# randomly downsample seurat single cell object 

set.seed(123)

Salcher_lung_cancer_atlas_LUAD_seurat_subsampled <- Salcher_lung_cancer_atlas_LUAD_seurat[, sample(colnames(Salcher_lung_cancer_atlas_LUAD_seurat), size=25000, replace=F)]

dim(Salcher_lung_cancer_atlas_LUAD_seurat_subsampled)

Salcher_lung_cancer_atlas_LUAD_seurat_subsampled$cell_type <- gsub("[^A-Za-z0-9 ]", " ", Salcher_lung_cancer_atlas_LUAD_seurat_subsampled$cell_type)

Idents(Salcher_lung_cancer_atlas_LUAD_seurat_subsampled) <- Salcher_lung_cancer_atlas_LUAD_seurat_subsampled$cell_type

generate_cytospace_from_scRNA_seurat_object(Salcher_lung_cancer_atlas_LUAD_seurat_subsampled,
                                            dir_out=here('data'),
                                            fout_prefix='')

# use non-normalized count matrix to generate ST objects for cytospace

for (i in 1:length(visium_seurat_objects)) {
  
  DefaultAssay(visium_seurat_objects[[i]]) <- "Spatial"
  
  generate_cytospace_from_ST_seurat_object(visium_seurat_objects[[i]],
                                           dir_out=here('data'),
                                           fout_prefix=paste(names(visium_seurat_objects)[i],".",sep=""),slice='slice1')
  
}


# generate cytospace outputs with this script: 

here('cytospace.sh')

# read in cytospace absolute predicted cell counts

for (i in 1:length(visium_seurat_objects)) {
  
  if(i == 6) {
    next
  }
  
  cytospace_results <- read.csv(paste(here('data/'),
                                      names(visium_seurat_objects)[i],"_cell_type_assignments_by_spot.csv",sep = ""))
  
  visium_seurat_objects[[i]]@meta.data$SpotID <- rownames(visium_seurat_objects[[i]]@meta.data)
  
  # add spots back with 0s in all cell types (these had 0 nuclei predicted by stardist)
  
  merged_df <- merge(cytospace_results, visium_seurat_objects[[i]]@meta.data, by = "SpotID", all.y = TRUE)
  
  merged_df[is.na(merged_df)] <- 0
  
  merged_df <- merged_df[match(colnames(visium_seurat_objects[[i]]), merged_df$SpotID),]
  
  visium_seurat_objects[[i]]@meta.data <- merged_df
  
}

# salcher

visium_seurat_objects_pred <- visium_seurat_objects

# saveRDS(visium_seurat_objects_pred, here('data','visium_seurat_objects_clustered_named_predictions_cyto_salcher_stardist.rds'))

########
# FIG 3B
########

# cytospace salcher results

visium_seurat_objects_pred <- readRDS(here('data','visium_seurat_objects_clustered_named_predictions_cyto_salcher_stardist.rds'))


# remove sample 6

visium_seurat_objects_pred <- visium_seurat_objects_pred[-6]


# merge together

visium_seurat_merge_pred <- merge(x = visium_seurat_objects_pred[[1]], 
                                  y = visium_seurat_objects_pred[-1], 
                                  add.cell.ids = names(visium_seurat_objects_pred))


cell_type_data <- visium_seurat_merge_pred

# add additional metadata to Seurat object

cell_type_data@meta.data$column_name <- rownames(cell_type_data@meta.data)

tmp <- merge(cell_type_data@meta.data, geo_stRNAseq_metadata, by = "orig.ident")

tmp <- tmp[match(colnames(cell_type_data), tmp$column_name),]

cell_type_data@meta.data <- tmp

cell_type_data$predominant_VI <- ifelse(cell_type_data$`novel grade` == 'VI',1,0)

# change idents to variable of interest

Idents(cell_type_data) <- cell_type_data$predominant_VI

cell_type_data <- as.data.frame(cell_type_data@meta.data)

cell_type_all <- cell_type_data # save for downstream comparison of all cell types


cell_type_all$predominant_VI <- as.factor(cell_type_all$predominant_VI)

levels(cell_type_all$predominant_VI) <- c("VI-","VI+")

tmp <- cell_type_all %>% dplyr::select('Pericyte':'myeloid.dividing','Total.cells',
                                       'Ciliated':'Club','T.cell.CD4':'Neutrophils',
                                       'predominant_VI','orig.ident')

fig3B <- gather(tmp, cell_type, abundance, c(Pericyte:myeloid.dividing,
                                             Ciliated:Club,
                                             T.cell.CD4:Neutrophils), 
                factor_key=TRUE)

fig3B$abundance[is.na(fig3B$abundance)] <- 0

fig3B$Total.cells[is.na(fig3B$Total.cells)] <- 0

fig3B$cell_type <- as.character(fig3B$cell_type)

fig3B <- fig3B %>% group_by(predominant_VI, orig.ident, cell_type) %>% 
  dplyr::summarise(abundance = sum(abundance),Total.cells = sum(Total.cells)) %>%
  mutate(freq = abundance / Total.cells) %>% 
  mutate(`Cell type` = cell_type) %>%
  select(-cell_type) %>%
  mutate(`Cell type` = factor(`Cell type`))

levels(fig3B$`Cell type`) <-  str_replace_all(levels(fig3B$`Cell type`), "\\.", " ")

set.seed(100)

colorss <- sample(colorRampPalette(brewer.pal(name="Paired", n = 12))(length(unique(fig3B$`Cell type`))))

use_colors_VI <- c(`no VI` = "#390099",
                   `VI+` = "#FF0054")

# remove cell types with negligible proportions from visualization

fig3B <- fig3B %>% filter(freq > 0.01)

colorss <- sample(colorRampPalette(brewer.pal(name="Paired", n = 12))(length(unique(fig3B$`Cell type`))))


addWorksheet(source_data, "Figure 3B")

writeData(source_data, sheet = "Figure 3B", 
          x = "Figure 3B", 
          startCol =  1, startRow = 1)

writeData(source_data, sheet = 'Figure 3B',
          x = fig3B, startCol = 1, startRow = 3)

fig3B_gg <- ggplot(fig3B, aes(x = predominant_VI, y = freq, fill = `Cell type`))+
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = colorss) +
  labs(title = 'Cell proportions in LUAD stRNA-seq (n=15)\n',
       x = NULL, y = "CytoSPACE predicted proportion\n") +
  theme_classic() +
  theme(
    plot.title = element_text(size=22, hjust = 0.5),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size=20),
    axis.text.x = element_text(size=20,colour = 'white'),
    axis.text.y = element_text(size=20),
    legend.position="right", 
    legend.title = element_text(size=20),
    legend.text = element_text(size=20),
    text = element_text(family = 'Calibri')
  ) +
  scale_x_discrete(position = "top") +
  guides(fill=guide_legend(ncol=1)) +
  coord_cartesian(ylim = c(0, 1), 
                  clip = 'off') +
  annotate(
    xmin = c(.55,1.55), 
    xmax = c(1.45,2.45),
    ymin = 1.05,
    ymax = 1.1,
    geom = "rect",
    fill = c(unname(use_colors_VI))
  )

fig3B_gg

# Get the current plot dimensions in inches

plot_dimensions <- dev.size("in")

# Calculate aspect ratio

aspect_ratio <- plot_dimensions[1] / plot_dimensions[2]

aspect_ratio <- 0.9776952

# Define a reasonable height (in inches)

desired_height <- 8.5

# Calculate corresponding width to maintain aspect ratio

desired_width <- desired_height * aspect_ratio

# Define the resolution

TIFF_dpi <- 600

# Save the plot with adjusted dimensions

ggsave(here('figures','fig3B.tiff'),
       plot = fig3B_gg,
       device = 'tiff',
       width = desired_width,
       height = desired_height,
       units = 'in', dpi = TIFF_dpi)



plasma <- fig3B %>% filter(`Cell type` == 'Plasma cell')

wilcox.test(freq ~ predominant_VI, data = plasma)


fibroblast_peri <- fig3B %>% filter(`Cell type` == 'Fibroblast peribronchial')

wilcox.test(freq ~ predominant_VI, data = fibroblast_peri)



# cell type enrichment by pathology region

###########
# FIGURE 3A
###########

visium_seurat_merge_pred_bub <- visium_seurat_merge_pred


# add additional metadata to seurat object

visium_seurat_merge_pred_bub@meta.data$column_name <- rownames(visium_seurat_merge_pred_bub@meta.data)

tmp <- merge(visium_seurat_merge_pred_bub@meta.data, geo_stRNAseq_metadata, by = "orig.ident")

tmp <- tmp[match(colnames(visium_seurat_merge_pred_bub), tmp$column_name),]

visium_seurat_merge_pred_bub@meta.data <- cbind(visium_seurat_merge_pred_bub@meta.data, tmp)


# convert pathology annotations to factor, set normal lung as reference level

visium_seurat_merge_pred_bub$pathology <- factor(visium_seurat_merge_pred_bub$pathology, levels = c("Normal Lung",unique(visium_seurat_merge_pred_bub$pathology)[-which(unique(visium_seurat_merge_pred_bub$pathology) == "Normal Lung")]))


# set unannotated areas to NA to drop from factor levels

levels(visium_seurat_merge_pred_bub$pathology)[which(levels(visium_seurat_merge_pred_bub$pathology) == "")] <- "Unannotated"

# remove path annotations with < 200 spots

visium_seurat_merge_pred_bub$pathology <- visium_seurat_merge_pred_bub$pathology[visium_seurat_merge_pred_bub$pathology %!in% plyr::count(visium_seurat_merge_pred_bub$pathology)$x[which(plyr::count(visium_seurat_merge_pred_bub$pathology)$freq < 200)]]

visium_seurat_merge_pred_bub$pathology <- droplevels(visium_seurat_merge_pred_bub$pathology)

# downsample 

for(i in 1:length(levels(visium_seurat_merge_pred_bub$pathology))){
  
  all_path <- rownames(visium_seurat_merge_pred_bub@meta.data)[which(visium_seurat_merge_pred_bub@meta.data$pathology == levels(visium_seurat_merge_pred_bub$pathology)[i])]
  
  set.seed(123)
  
  downsampled_path <- sample(rownames(visium_seurat_merge_pred_bub@meta.data)[which(visium_seurat_merge_pred_bub@meta.data$pathology == levels(visium_seurat_merge_pred_bub$pathology)[i])], size=200, replace=F)
  
  remove <- all_path[all_path %!in% downsampled_path]
  
  visium_seurat_merge_pred_bub <- visium_seurat_merge_pred_bub[,-which(colnames(visium_seurat_merge_pred_bub) %in% remove)]
  
}



# plot cell type predictions by pathology region

# using mixed linear effects model with sample ID as random effect

mixed_effects_stats <- list()

levels(visium_seurat_merge_pred_bub@meta.data$pathology)[which(levels(visium_seurat_merge_pred_bub@meta.data$pathology) == 'Unannotated')] <- NA

visium_seurat_merge_pred_bub@meta.data$pathology <- droplevels(visium_seurat_merge_pred_bub@meta.data$pathology)


# change NAs to 0s

visium_seurat_merge_pred_bub@meta.data <- visium_seurat_merge_pred_bub@meta.data %>% 
  mutate_at(vars(Pericyte:myeloid.dividing,Ciliated:Club,
                 T.cell.CD4:Neutrophils), ~ ifelse(is.na(.), 0, .))



# get cell type names

tmp <- visium_seurat_merge_pred_bub@meta.data %>% dplyr::select('Pericyte':'myeloid.dividing',
                                                                'Ciliated':'Club','T.cell.CD4':'Neutrophils')

cell_types <- colnames(tmp)


for(i in 1:length(cell_types)){
  
  print(i)
  
  data <- visium_seurat_merge_pred_bub@meta.data
  
  scores <- aggregate(data[,cell_types[i]] ~ pathology, data = data, mean)
  
  colnames(scores)[2] <- cell_types[i]
  
  if (all(scores[,2] == 0)) {
    
    print(paste('No predicted cells in any region, creating NA values for',cell_types[i]))
    
    mixed_effects_stats[[i]] <- data.frame(cell_type = rep(0,length(levels(data$pathology))), 
                                           p_value = rep(NA, length(levels(data$pathology))))
    
    colnames(mixed_effects_stats[[i]]) <- c(cell_types[i], paste('p_value_',cell_types[i],sep=''))
    
    next
    
  }
  
  data <- data.frame(y = data[,cell_types[i]], x = data$pathology, sample = data$`sample id`)
  
  m <- lmer(y ~ 1 + x + (1 | sample), data = data)
  
  mixed_effects_stats[[i]] <- data.frame(summary(m)$coefficients)
  
  mixed_effects_stats[[i]]$pathology <- sapply(strsplit(rownames(mixed_effects_stats[[i]]),'x'),"[",2)
  
  mixed_effects_stats[[i]]$pathology[1] <- 'Normal Lung' # rename intercept
  
  mixed_effects_stats[[i]] <- merge(mixed_effects_stats[[i]], scores, by = 'pathology')
  
  mixed_effects_stats[[i]]$scores <- scores$`data[, cell_types[i]]`
  
  rownames(mixed_effects_stats[[i]]) <- mixed_effects_stats[[i]]$pathology
  
  mixed_effects_stats[[i]] <- mixed_effects_stats[[i]] %>% dplyr::select(cell_types[i], Pr...t..)
  
  colnames(mixed_effects_stats[[i]])[which(colnames(mixed_effects_stats[[i]]) == 'Pr...t..')] <- paste('p_value',cell_types[i],sep = '_')
  
}



mixed_effects_stats_all <- bind_cols(mixed_effects_stats)

# bubble plot heatmap to show scores and p values together

mixed_effects_stats_all$region <- rownames(mixed_effects_stats_all)

tmp <- mixed_effects_stats_all %>% dplyr::select(one_of('region',cell_types))

fig3A <- gather(tmp, cell_type, score, Pericyte:Neutrophils, factor_key=TRUE)

tmp <- mixed_effects_stats_all %>% dplyr::select(paste('p_value',cell_types,sep='_'))

data_long2 <- gather(tmp, cluster_p, p_value, p_value_Pericyte:p_value_Neutrophils, factor_key=TRUE)

fig3A <- cbind(fig3A, data_long2)


# need to adjust p values for multiple hypothesis testing

fig3A$FDR <- p.adjust(fig3A$p_value, n = length(fig3A$p_value))

fig3A$FDR <- -log10(fig3A$FDR)


# specify column order

fig3A$cell_type <- factor(fig3A$cell_type, levels = sort(levels(fig3A$cell_type)))


levels(fig3A$cell_type) <-  str_replace_all(levels(fig3A$cell_type), "\\.", " ")

# specify row order

fig3A$region <- factor(fig3A$region, levels = c('Normal Lung','Normal Vessel','Lymphoid Aggregate','Stroma','Desmoplastic Stroma','Normal Bronchus',
                                                'Collapsed Lung','Abnormal Pleura','Lepidic','Acinar','Papillary','Micropapillary','Cribriform','Solid',
                                                'VPI','VI'))

levels(fig3A$region)[which(levels(fig3A$region) == 'VI')] <- 'VI Focus'

levels(fig3A$region)[which(levels(fig3A$region) == 'VPI')] <- 'VPI Focus'

# plot cell type proportions

# remove from plot enrichment scores > 0.01 FDR

fig3A <- fig3A %>%
  mutate(FDR_filtered = ifelse(FDR > 2, FDR, NA))

addWorksheet(source_data, "Figure 3A")

writeData(source_data, sheet = "Figure 3A", 
          x = "Figure 3A", 
          startCol =  1, startRow = 1)

writeData(source_data, sheet = 'Figure 3A',
          x = fig3A, startCol = 1, startRow = 3)

fig3A_gg <- ggplot(data = fig3A, aes(x = cell_type, y = region)) +
  geom_point(aes(fill = score, size = FDR_filtered), shape = 21) +
  scale_size(breaks=c(3,20,50,100),labels=c(3,20,50,100), range = c(2,9)) +
  scale_fill_viridis() +
  theme_bw() +
  theme(
    legend.position = 'right',
    panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
    plot.title = element_text(size=26, hjust = 0.5),
    axis.title.x = element_text(size=20, ),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size=20, angle=45, vjust=1, hjust=1),
    axis.text.y = element_text(size=20),
    axis.text=element_text(size=20),
    axis.ticks=element_blank(),
    legend.title = element_text(size=20),
    legend.key.size = unit(0.6, "cm"),
    legend.text = element_text(size=20),
    text = element_text(family = 'Calibri')) +
  labs(title = 'Cell type association with LUAD pathology annotation (n=15)\n',
       x = element_blank(), y = 'Visium Region', fill = 'Average cells/spot') +
  scale_y_discrete(limits = rev) + 
  scale_x_discrete(position = "bottom") +
  guides(size=guide_legend(title="-log₁₀(p.adj)"))

fig3A_gg

plot_dimensions <- dev.size("in")

aspect_ratio <- plot_dimensions[1] / plot_dimensions[2]

aspect_ratio <- 1.759603

desired_height <- 8

desired_width <- desired_height * aspect_ratio

TIFF_dpi <- 600

ggsave(here('figures','fig3A.tiff'),
       plot = fig3A_gg,
       device = 'tiff',
       width = desired_width,
       height = desired_height,
       units = 'in', dpi = TIFF_dpi)





######################
# Extended Data Fig 5A
######################

visium_seurat_merge_pred_freq <- visium_seurat_merge_pred


# add additional metadata to seurat object

visium_seurat_merge_pred_freq@meta.data$column_name <- rownames(visium_seurat_merge_pred_freq@meta.data)

tmp <- merge(visium_seurat_merge_pred_freq@meta.data, geo_stRNAseq_metadata, by = "orig.ident")

tmp <- tmp[match(colnames(visium_seurat_merge_pred_freq), tmp$column_name),]

visium_seurat_merge_pred_freq@meta.data <- cbind(visium_seurat_merge_pred_freq@meta.data, tmp)


# convert pathology annotations to factor, set normal lung as reference level

visium_seurat_merge_pred_freq$pathology <- factor(visium_seurat_merge_pred_freq$pathology, levels = c("Normal Lung",unique(visium_seurat_merge_pred_freq$pathology)[-which(unique(visium_seurat_merge_pred_freq$pathology) == "Normal Lung")]))


# set unannotated areas to NA to drop from factor levels

levels(visium_seurat_merge_pred_freq$pathology)[which(levels(visium_seurat_merge_pred_freq$pathology) == "")] <- "Unannotated"

# remove path annotations with < 200 spots

visium_seurat_merge_pred_freq$pathology <- visium_seurat_merge_pred_freq$pathology[visium_seurat_merge_pred_freq$pathology %!in% plyr::count(visium_seurat_merge_pred_freq$pathology)$x[which(plyr::count(visium_seurat_merge_pred_freq$pathology)$freq < 200)]]

visium_seurat_merge_pred_freq$pathology <- droplevels(visium_seurat_merge_pred_freq$pathology)


# plot cell type predictions by pathology region

# using mixed linear effects model with sample ID as random effect

mixed_effects_stats <- list()

levels(visium_seurat_merge_pred_freq@meta.data$pathology)[which(levels(visium_seurat_merge_pred_freq@meta.data$pathology) == 'Unannotated')] <- NA

visium_seurat_merge_pred_freq@meta.data$pathology <- droplevels(visium_seurat_merge_pred_freq@meta.data$pathology)


# change NAs to 0s

visium_seurat_merge_pred_freq@meta.data <- visium_seurat_merge_pred_freq@meta.data %>% 
  mutate_at(vars(Pericyte:myeloid.dividing,Ciliated:Club,
                 T.cell.CD4:Neutrophils), ~ ifelse(is.na(.), 0, .))


df <- visium_seurat_merge_pred_freq@meta.data

df_long <- df %>%
  select(pathology, 'Pericyte':'myeloid.dividing',
         'Ciliated':'Club','T.cell.CD4':'Neutrophils') %>%
  pivot_longer(
    cols = Pericyte:Neutrophils,
    names_to = "cell_type",
    values_to = "count"
  ) %>%
  group_by(pathology, cell_type) %>%
  summarise(total_count = sum(count, na.rm = TRUE), .groups = "drop") %>%
  complete(pathology, cell_type, fill = list(total_count = 0)) %>%
  group_by(pathology) %>%
  mutate(pct = total_count / sum(total_count) * 100) %>%
  ungroup()

spot_counts <- df %>%
  group_by(pathology) %>%
  summarise(
    n_spots = n(),
    total_cell_count = sum(
      rowSums(across(c(Pericyte, Fibroblast.alveolar, Ciliated, Endothelial.cell.venous,
                       B.cell, transitional.club.AT2, Fibroblast.peribronchial,
                       Smooth.muscle.cell, Plasma.cell, Endothelial.cell.lymphatic,
                       Fibroblast.adventitial, Tumor.cells, Alveolar.cell.type.1,
                       Alveolar.cell.type.2, Macrophage, Macrophage.alveolar,
                       myeloid.dividing, Total.cells, Ciliated, T.cell.regulatory,
                       Club, T.cell.CD4, Endothelial.cell.capillary, 
                       Endothelial.cell.arterial, cDC2, ROS1..healthy.epithelial,
                       Mesothelial, B.cell.dividing, Plasma.cell.dividing,
                       Neutrophils)), na.rm = TRUE)
    )
  ) %>%
  mutate(facet_label = paste0(
    pathology, "\n(n = ", n_spots, " spots,\n",
    total_cell_count, " cells)"
  ))

df_long <- left_join(df_long, spot_counts, by = "pathology")

df_long$cell_type <-  str_replace_all(df_long$cell_type, "\\.", " ")

efig5A <- df_long %>% 
  filter(pct > 1) %>%
  filter(!is.na(pathology))

addWorksheet(source_data, "E. Figure 5A")

writeData(source_data, sheet = "E. Figure 5A", 
          x = "E. Figure 5A", 
          startCol =  1, startRow = 1)

writeData(source_data, sheet = 'E. Figure 5A',
          x = efig5A, startCol = 1, startRow = 3)

set.seed(111)

manual_colors <- c(
  "Alveolar cell type 1" = "#A281BD",
  "Alveolar cell type 2" = "#724899",
  "B cell" = "#E7372A",
  "Ciliated" = "#FDA03A",
  "Club" = "#4A96A7",
  "Endothelial cell arterial" = "#66a61e",
  "Endothelial cell lymphatic" = "#F88519",
  "Endothelial cell venous" = "#EE5656",
  "Fibroblast adventitial" = "#F8A160",
  "Fibroblast alveolar" = "#A6CEE3",
  "Fibroblast peribronchial" = "#B15928",
  "Macrophage" = "#E3C471",
  "Pericyte" = "#4E96C4",
  "Plasma cell" = "#A9D88C",
  "Smooth muscle cell" = "#D6A5A3",
  "T cell regulatory" = "#629E45",
  "transitional club AT2" = "#E39A8C",
  "Tumor cells" = "#D3C599"
)


efig5A_gg <- ggplot(efig5A, aes(x = factor(1), y = pct, fill = cell_type)) +
  geom_bar(width = 1, stat = "identity", color = "white") +
  coord_polar(theta = "y", start = 0, clip = 'off') +
  facet_wrap(~ facet_label) +
  theme_void() +
  theme(
    plot.title = element_text(size=18, hjust = 0.5),
    legend.position = "right",
    strip.text = element_text(size = 12),
    text = element_text(family = 'Calibri'),
    legend.title = element_text(size=12),
    legend.text = element_text(size=12),
    panel.spacing.x = unit(1, "cm"),
    strip.clip = "off"
  ) +
  geom_text(
    aes(label = ifelse(pct > 10, paste0(round(pct, 1), "%"), "")),
    position = position_stack(vjust = 0.5),
    size = 3
  ) +
  labs(
    title = "Cell type proportions across LUAD pathology annotations\n",
    fill = "Cell types (>1%)"
  ) +
  scale_fill_manual(values = manual_colors)

efig5A_gg

plot_dimensions <- dev.size("in")

aspect_ratio <- plot_dimensions[1] / plot_dimensions[2]

aspect_ratio <- 1.145455

desired_height <- 8

desired_width <- desired_height * aspect_ratio

TIFF_dpi <- 600

ggsave(here('figures','efig5A.tiff'),
       plot = efig5A_gg,
       device = 'tiff',
       width = desired_width,
       height = desired_height,
       units = 'in', dpi = TIFF_dpi)




# spatially weighted regression analysis of VI clusters and cell types

######################
# EXTENDED DATA FIG 5C
######################

# only include cell types that had cells in at least 20% of sections

cell_type_data <- visium_seurat_merge_pred

cell_type_data <- as.data.frame(cell_type_data@meta.data)

cell_type_all <- cell_type_data 

tmp <- cell_type_all %>% dplyr::select('Pericyte':'myeloid.dividing','Total.cells',
                                       'Ciliated':'Club','T.cell.CD4':'Neutrophils',
                                       'orig.ident')

data_long <- gather(tmp, cell_type, abundance, c(Pericyte:myeloid.dividing,
                                                 Ciliated:Club,
                                                 T.cell.CD4:Neutrophils), 
                    factor_key=TRUE)


data_long$abundance[is.na(data_long$abundance)] <- 0

data_long$Total.cells[is.na(data_long$Total.cells)] <- 0

data_long <- data_long %>%
  group_by(orig.ident, cell_type) %>%
  summarise(total = sum(abundance)) %>%
  ungroup()

# Filter out cell types that don't meet the criteria (total > 0 in at least 20% of orig.ident groups)

data_long <- data_long %>%
  group_by(cell_type) %>%
  filter(sum(total > 0) >= 0.2 * length(unique(data_long$orig.ident)))

remaining_cell_types <- unique(data_long$cell_type)

remaining_cell_types <- gsub("\\.", " ", remaining_cell_types)


# cell type markers from Salcher

options(future.globals.maxSize = 8000 * 1024^2)

Salcher_lung_cancer_atlas_LUAD_seurat <- readRDS(here('data','Salcher_lung_cancer_atlas_LUAD.rds'))

Idents(Salcher_lung_cancer_atlas_LUAD_seurat) <- Salcher_lung_cancer_atlas_LUAD_seurat$cell_type

# these steps take a very long time

# Salcher_lung_cancer_atlas_LUAD_seurat <- ScaleData(Salcher_lung_cancer_atlas_LUAD_seurat)

# saveRDS(Salcher_lung_cancer_atlas_LUAD_seurat, file = here('data','Salcher_lung_cancer_atlas_LUAD_scaled.rds'))

Salcher_lung_cancer_atlas_LUAD_seurat <- readRDS(here('data','Salcher_lung_cancer_atlas_LUAD_scaled.rds'))


# Salcher_lung_cancer_atlas_LUAD_seurat <- AddModuleScore(Salcher_lung_cancer_atlas_LUAD_seurat, features = VI_genesets, name = '_score')

# names(Salcher_lung_cancer_atlas_LUAD_seurat@meta.data)[grep("_score", names(Salcher_lung_cancer_atlas_LUAD_seurat@meta.data))] <- names(VI_genesets)


# salcher_markers <- FindAllMarkers(Salcher_lung_cancer_atlas_LUAD_seurat, assay = 'RNA')

# saveRDS(salcher_markers, file = here('data','salcher_markers.rds'))

salcher_markers <- readRDS(here('data','salcher_markers.rds'))

salcher_markers_top <- salcher_markers %>%
  group_by(cluster) %>%
  slice_head(n = 50) %>%
  ungroup()

subset.matrix <- Salcher_lung_cancer_atlas_LUAD_seurat@assays$RNA@scale.data[which(rownames(Salcher_lung_cancer_atlas_LUAD_seurat@assays$RNA@scale.data) %in% salcher_markers_top$gene), ] # Pull the raw expression matrix from the original Seurat object containing only the genes of interest

# Create a new Seurat object with just the genes of interest

object2 <- CreateSeuratObject(subset.matrix) 

object2$cell_type <- Salcher_lung_cancer_atlas_LUAD_seurat$cell_type   

Idents(object2) <- object2$cell_type

av.exp <- AverageExpression(object2)$RNA

cor.exp <- as.data.frame(cor(av.exp))

cor.exp$x <- rownames(cor.exp)


colnames(cor.exp) <- gsub("[^A-Za-z0-9]+", ".", colnames(cor.exp))

cor.exp$x <-  gsub("[^A-Za-z0-9]+", ".", cor.exp$x)


efig5C <- tidyr::gather(data = cor.exp, y, correlation, c(colnames(cor.exp)[colnames(cor.exp) %in% data_long$cell_type])) %>%
  filter(x %in% data_long$cell_type)

efig5C$x <- gsub("\\.", " ", efig5C$x)

efig5C$y <- gsub("\\.", " ", efig5C$y)

addWorksheet(source_data, "E. Figure 5C")

writeData(source_data, sheet = "E. Figure 5C", 
          x = "E. Figure 5C", 
          startCol =  1, startRow = 1)

writeData(source_data, sheet = 'E. Figure 5C',
          x = efig5C, startCol = 1, startRow = 3)


efig5C_gg <- ggplot(efig5C, aes(x, y, fill = correlation)) +
  geom_tile(color = "black") +
  scale_fill_viridis() +
  theme(axis.text.x = element_text(angle = 45, hjust=1))  +
  theme(
    panel.background = element_blank(),
    axis.text.y=element_text(size=20),
    axis.text.x=element_text(size=20),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.position = 'right',
    legend.title = element_text(size=20),
    legend.text = element_text(size=20),
    text = element_text(family = 'Calibri'),
    plot.title = element_text(size=24, hjust = 0.5)) +
  labs(title = 'Top 50 genes/cell type in stage I LUAD\nfrom Salcher et al atlas\n')

efig5C_gg

plot_dimensions <- dev.size("in")

aspect_ratio <- plot_dimensions[1] / plot_dimensions[2]

aspect_ratio <- 1.1

desired_height <- 10

desired_width <- desired_height * aspect_ratio

TIFF_dpi <- 600

ggsave(here('figures','efig5C.tiff'),
       plot = efig5C_gg,
       device = 'tiff',
       width = desired_width,
       height = desired_height,
       units = 'in', dpi = TIFF_dpi)



########
# FIG 3F
########

salcher_genesets <- list()

for (i in unique(salcher_markers_top$cluster)) {
  
  cluster_subset <- salcher_markers_top[salcher_markers_top$cluster == i, ]
  
  salcher_genesets[[as.character(i)]] <- cluster_subset$gene
  
}

names(salcher_genesets) <- gsub(' ','_',names(salcher_genesets))

names(salcher_genesets) <- gsub("[+/]", " ", names(salcher_genesets))

names(salcher_genesets) <- gsub('_',' ',names(salcher_genesets))


# subset to cell types recovered in deconvolution

salcher_genesets <- salcher_genesets[names(salcher_genesets) %in% remaining_cell_types]

names(salcher_genesets) <- gsub(' ','_',names(salcher_genesets))



visium_seurat_objects <- lapply(X = visium_seurat_objects, FUN = function(x) {
  
  x <- AddModuleScore(x, features = salcher_genesets, name = '_celltype_score')
  
})


for (i in 1:length(visium_seurat_objects)){
  
  names(visium_seurat_objects[[i]]@meta.data)[grep("_score", names(visium_seurat_objects[[i]]@meta.data))] <- names(salcher_genesets)
  
}

# this step takes a very long time

for (i in 1:length(visium_seurat_objects)){
  
  # get image coordinates
  
  cat('Getting coordinates', i, 'of', length(visium_seurat_objects),'\n')
  
  visium_coordinates <- visium_seurat_objects[[i]]@images$slice1@coordinates
  
  # add geneset scores for analysis
  
  visium_coordinates <- visium_coordinates %>%
    bind_cols(FetchData(visium_seurat_objects[[i]], c('VI_gene_cluster_1', 
                                                      'VI_gene_cluster_2', 
                                                      'VI_gene_cluster_3', 
                                                      'VI_gene_cluster_4')))
  
  
  visium_coordinates2 <- bind_cols(FetchData(visium_seurat_objects[[i]], c(names(salcher_genesets))))
  
  visium_coordinates <- cbind(visium_coordinates, visium_coordinates2)
  
  # convert to spatial points dataframe
  
  coordinates(visium_coordinates) <- c('row','col')
  
  # perform spatially weighted regression analysis
  
  gws_results <- gwss(visium_coordinates, vars = names(visium_coordinates)[-c(1:3)], bw = 5)
  
  gws_results <- data.frame(gws_results$SDF)
  
  # extract spearman correlation summary stats
  
  gws_results <- gws_results %>% 
    select(starts_with('Spearman'))

  gws_results <- gws_results %>% mutate_all(~ replace(., is.nan(.),0))
  
  visium_seurat_objects[[i]]@meta.data <- cbind(visium_seurat_objects[[i]]@meta.data, gws_results)
  
  
}

# saveRDS(visium_seurat_objects, file = here('data','visium_seurat_objects_gws_subset.rds'))

visium_seurat_objects <- readRDS(here('data','visium_seurat_objects_gws_subset.rds'))

visium_gw_cor_means_all <- list()

for (i in 1:length(visium_seurat_objects)){
  
  # take mean of each spatially weighted correlation pair
  
  visium_gw_cor_means <- visium_seurat_objects[[i]]@meta.data %>%
    select(starts_with("Spearman")) %>%
    dplyr::summarise(across(everything(), ~mean(., na.rm = TRUE)))
  
  names(visium_gw_cor_means) <- paste(unique(visium_seurat_objects[[i]]$orig.ident), names(visium_gw_cor_means), sep = '_')
  
  visium_gw_cor_means_all[[i]] <- data.frame(t(visium_gw_cor_means))
  
}

# plot global correlation

visium_gw_cor_means_all_df <- bind_rows(visium_gw_cor_means_all, .id = "orig.ident")

# remove sample 6

visium_gw_cor_means_all_df <- visium_gw_cor_means_all_df %>%
  filter(orig.ident != 6)

geo_stRNAseq_metadata$Library <- as.integer(sub(".*_", "", geo_stRNAseq_metadata$`*library name`))

visium_gw_cor_means_all_df <- visium_gw_cor_means_all_df %>% 
  rownames_to_column(var = 'corr') %>%
  mutate(clusters = sapply(strsplit(corr,'rho_'),"[",2)) %>%
  dplyr::rename(Library = orig.ident) %>%
  mutate(Library = as.integer(Library)) %>%
  inner_join(geo_stRNAseq_metadata, by = 'Library')

# calculate mean over all samples

visium_gw_cor_means_all_df_mean <- visium_gw_cor_means_all_df %>%
  group_by(clusters) %>%
  dplyr::summarise(mean_clusters = mean(t.visium_gw_cor_means.))

# Extract the gene cluster names

unique_clusters <- unique(unlist(strsplit(visium_gw_cor_means_all_df_mean$clusters,"\\.")))

# Create a square matrix to store the correlations

num_clusters <- length(unique_clusters)

cor_matrix <- matrix(0, nrow = num_clusters, ncol = num_clusters,
                     dimnames = list(unique_clusters, unique_clusters))

# Fill in the correlation matrix

for (i in seq_len(length(visium_gw_cor_means_all_df_mean$mean_clusters))) {
  
  row_cluster <- unique(unlist(strsplit(visium_gw_cor_means_all_df_mean$clusters[i],"\\.")))[1]
  
  col_cluster <- unique(unlist(strsplit(visium_gw_cor_means_all_df_mean$clusters[i],"\\.")))[2]
  
  cor_value <- visium_gw_cor_means_all_df_mean$mean_clusters[i]
  
  cor_matrix[row_cluster, col_cluster] <- cor_value
  
  cor_matrix[col_cluster, row_cluster] <- cor_value
  
  print(i)
  
}

cor_matrix[cor_matrix == 0] <- 1

# plot mean spatially weighted correlation across samples

colnames(cor_matrix) <- gsub('_',' ',colnames(cor_matrix))

rownames(cor_matrix) <- gsub('_',' ',rownames(cor_matrix))


cor_matrix <- cor_matrix[which(rownames(cor_matrix) %in% c(remaining_cell_types, 
                                                           'VI gene cluster 1',
                                                           'VI gene cluster 2',
                                                           'VI gene cluster 3',
                                                           'VI gene cluster 4')),]

cor_matrix <- cor_matrix[,which(colnames(cor_matrix) %in% c(remaining_cell_types, 
                                                            'VI gene cluster 1',
                                                            'VI gene cluster 2',
                                                            'VI gene cluster 3',
                                                            'VI gene cluster 4'))]


tmp <- colnames(cor_matrix)

tmp[-grep('VI',tmp, invert = FALSE)] <- ' '

annotation_col <- data.frame(tmp)

rownames(annotation_col) <- colnames(cor_matrix)

colnames(annotation_col) <- 'VI gene cluster'

annotation_row <- data.frame(tmp)

rownames(annotation_row) <- colnames(cor_matrix)

colnames(annotation_row) <- 'VI gene cluster'

annotation_colors = list(
  `VI gene cluster` = c("VI gene cluster 1" = "#cab2d6",
                        "VI gene cluster 2" = "#b2df8a",
                        "VI gene cluster 3" = "#E31A1C",
                        "VI gene cluster 4" = "#A6CEE3",
                        ' ' = 'white'))

bluered <- colorRampPalette(c("#007FFF","white","red"))(length(seq(-1, 1,
                                                                   by = 0.01)))


fig3F <- cor_matrix

addWorksheet(source_data, "Figure 3F")

writeData(source_data, sheet = "Figure 3F", 
          x = "Figure 3F", 
          startCol =  1, startRow = 1)

writeData(source_data, sheet = 'Figure 3F',
          x = fig3F, startCol = 1, startRow = 3)


fig3F_gg <- pheatmap(fig3F,
                     color = diverging_hcl(length(seq(-1, 1, by = 0.01)), 
                                           palette = 'Blue-Red 3'), 
                     border_color = 'white', scale = "none", show_rownames = TRUE,
                     show_colnames = TRUE, cluster_rows = T, cluster_cols = T,
                     clustering_distance_rows = 'euclidean', 
                     clustering_distance_cols = 'euclidean',
                     legend = TRUE, annotation_legend = FALSE,
                     annotation_col = annotation_col,
                     annotation_row = annotation_row,
                     angle_col = 45,
                     annotation_colors = annotation_colors,
                     fontsize = 12,
                     clustering_method = "ward.D2", 
                     main = "Mean spatially weighted correlation (n=15)\n",
                     breaks = seq(-1,1, by = 0.01),
                     fontfamily = 'Calibri')


dev.off()

fig3F_gg$gtable$grobs[[1]]$gp = gpar(fontface = 'plain', fontsize = 20)

fig3F_gg$gtable$grobs[[5]]$gp = gpar(fontface = 'plain')

fig3F_gg


plot_dimensions <- dev.size("in")

aspect_ratio <- plot_dimensions[1] / plot_dimensions[2]

aspect_ratio <- 1.008879

desired_height <- 8.5  

desired_width <- desired_height * aspect_ratio

TIFF_dpi <- 600

ggsave(here('figures','fig3F.tiff'),
       plot = fig3F_gg,
       device = 'tiff',
       width = desired_width,
       height = desired_height,
       units = 'in', dpi = TIFF_dpi)



##########
# FIG 3C_1
##########

# peribronchial fibroblasts are myofibroblast CAFs

caf_genesets <- list()

# CAF signatures from Hanley et al 2023

hanley_caf_sigs <- read_excel(here('data','Hanley_suppData2_CAF_signatures.xlsx'),
                              sheet = 2, col_names = TRUE, skip = 1)

myo_caf <- hanley_caf_sigs %>% 
  filter(Cluster == 'Myo' & p_val_adj.Sample < 0.01 & avg_log2FC > 1) %>%
  select(gene)

caf_genesets$myo_caf <- as.vector(myo_caf)$gene

adv_caf <- hanley_caf_sigs %>% 
  filter(Cluster == 'Adventitial' & p_val_adj.Sample < 0.01 & avg_log2FC > 1) %>%
  select(gene)

caf_genesets$adv_caf <- as.vector(adv_caf)$gene

alv_caf <- hanley_caf_sigs %>% 
  filter(Cluster == 'Alveolar' & p_val_adj.Sample < 0.01 & avg_log2FC > 1) %>%
  select(gene)

caf_genesets$alv_caf <- as.vector(alv_caf)$gene


# score salcher LUAD atlas for caf signatures

Salcher_lung_cancer_atlas_LUAD_seurat <- readRDS(here('data','Salcher_lung_cancer_atlas_LUAD.rds'))

Idents(Salcher_lung_cancer_atlas_LUAD_seurat) <- Salcher_lung_cancer_atlas_LUAD_seurat$cell_type


Salcher_lung_cancer_atlas_LUAD_seurat <- AddModuleScore(Salcher_lung_cancer_atlas_LUAD_seurat, features = caf_genesets, name = '_score')


names(Salcher_lung_cancer_atlas_LUAD_seurat@meta.data)[grep("_score", names(Salcher_lung_cancer_atlas_LUAD_seurat@meta.data))] <- names(caf_genesets)

tmp <- subset(Salcher_lung_cancer_atlas_LUAD_seurat, idents = c('Fibroblast alveolar','Fibroblast peribronchial','Fibroblast adventitial'))


data_long <- gather(tmp@meta.data, fibroblast, score, myo_caf:alv_caf, factor_key=TRUE)


fig3C_1 <- data_long %>%
  filter(fibroblast == 'myo_caf') %>%
  group_by(sample, cell_type) %>%
  dplyr::summarise(mean_score = mean(score)) %>%
  group_by(sample) %>%
  mutate(mean_score = (mean_score - mean(mean_score)) / sd(mean_score))

fig3C_1 <- tibble(fig3C_1)

fig3C_1$cell_type <- droplevels(fig3C_1$cell_type)

addWorksheet(source_data, "Figure 3B-1")

writeData(source_data, sheet = "Figure 3B-1", 
          x = "Figure 3B-1", 
          startCol =  1, startRow = 1)

writeData(source_data, sheet = 'Figure 3B-1',
          x = fig3C_1, startCol = 1, startRow = 3)


pval <- fig3C_1 %>% 
  wilcox_test(
    mean_score ~ cell_type
  ) %>% 
  add_xy_position() %>% 
  mutate(p = sub("e","%.% 10^",p) )


fig3C_1_gg <- ggplot(data = fig3C_1, aes(x = cell_type, y = mean_score)) +
  theme_classic()+
  scale_fill_brewer(palette = "Dark2") +
  geom_violin(lwd=0.5, width = 1, aes(fill=cell_type)) +
  geom_boxplot(color = "black", fill = 'white', width = .15, 
               outlier.alpha = 0, lwd=0.5) +
  geom_jitter(colour="black",fill='white', size=1, alpha=1, 
              width = 0.2, shape=21) +
  theme(plot.title = element_text(size=22, hjust = 0.5),
        axis.title.x = element_text(size=22, face="bold"),
        axis.title.y = element_text(size=22),
        axis.text.x = element_text(colour="black", size=22),
        axis.text.y = element_text(colour="black", size=22),
        legend.position="none", axis.text=element_text(size=12),
        legend.title = element_text(size=14),
        legend.key.size = unit(0.6, "cm"),
        legend.text = element_text(size=14),
        text = element_text(family = 'Calibri')) +
  scale_x_discrete(labels = c("Fibroblast adventitial" = "Adventitial", 
                              "Fibroblast alveolar" = "Alveolar", 
                              "Fibroblast peribronchial" = "Peribronchial")) +
  labs(title = "Salcher et al fibroblasts (n=65)\n",
       x = NULL, y = "Hanley et al myoCAF enrichment score\n") +
  add_pvalue(pval,label = 'p', step.increase = 0.1,
             family = 'Calibri', label.size = 6, parse = TRUE) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))

fig3C_1_gg$layers[[4]]$aes_params$family <- "Calibri"

fig3C_1_gg

plot_dimensions <- dev.size("in")

aspect_ratio <- plot_dimensions[1] / plot_dimensions[2]

aspect_ratio <- 1.148601

desired_height <- 6 

desired_width <- desired_height * aspect_ratio

TIFF_dpi <- 600

ggsave(here('figures','fig3C_1.tiff'),
       plot = fig3C_1_gg,
       device = 'tiff',
       width = desired_width,
       height = desired_height,
       units = 'in', dpi = TIFF_dpi)


##########
# FIG 3C_2
##########

fig3C_2 <- data_long %>%
  filter(fibroblast == 'adv_caf') %>%
  group_by(sample, cell_type) %>%
  dplyr::summarise(mean_score = mean(score)) %>%
  group_by(sample) %>%
  mutate(mean_score = (mean_score - mean(mean_score)) / sd(mean_score))

fig3C_2 <- tibble(fig3C_2)

fig3C_2$cell_type <- droplevels(fig3C_2$cell_type)

addWorksheet(source_data, "Figure 3B-2")

writeData(source_data, sheet = "Figure 3B-2", 
          x = "Figure 3B-2", 
          startCol =  1, startRow = 1)

writeData(source_data, sheet = 'Figure 3B-2',
          x = fig3C_2, startCol = 1, startRow = 3)

pval <- fig3C_2 %>% 
  wilcox_test(
    mean_score ~ cell_type
  ) %>% 
  add_xy_position() %>% 
  mutate(p = sub("e","%.% 10^",p) )

fig3C_2_gg <- ggplot(data = fig3C_2, aes(x = cell_type, y = mean_score)) +
  theme_classic()+
  scale_fill_brewer(palette = "Dark2") +
  geom_violin(lwd=0.5, width = 1, aes(fill=cell_type)) +
  geom_boxplot(color = "black", fill = 'white', width = .15, 
               outlier.alpha = 0, lwd=0.5) +
  geom_jitter(colour="black",fill='white', size=1, alpha=1, 
              width = 0.2, shape=21) +
  theme(plot.title = element_text(size=22, hjust = 0.5),
        axis.title.x = element_text(size=22, face="bold"),
        axis.title.y = element_text(size=22),
        axis.text.x = element_text(colour="black", size=22),
        axis.text.y = element_text(colour="black", size=22),
        legend.position="none", axis.text=element_text(size=12),
        legend.title = element_text(size=14),
        legend.key.size = unit(0.6, "cm"),
        legend.text = element_text(size=14),
        text = element_text(family = 'Calibri')) +
  scale_x_discrete(labels = c("Fibroblast adventitial" = "Adventitial", 
                              "Fibroblast alveolar" = "Alveolar", 
                              "Fibroblast peribronchial" = "Peribronchial")) +
  labs(title = "Salcher et al fibroblasts (n=65)\n",
       x = NULL, y = "Hanley et al adv fib enrichment score\n") +
  add_pvalue(pval,label = 'p', step.increase = 0.1,
             family = 'Calibri', label.size = 6, parse = TRUE,
             bracket.nudge.y = 0.1) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))

fig3C_2_gg$layers[[4]]$aes_params$family <- "Calibri"

fig3C_2_gg

ggsave(here('figures','fig3C_2.tiff'),
       plot = fig3C_2_gg,
       device = 'tiff',
       width = desired_width,
       height = desired_height,
       units = 'in', dpi = TIFF_dpi)

##########
# FIG 3C_3
##########

fig3C_3 <- data_long %>%
  filter(fibroblast == 'alv_caf') %>%
  group_by(sample, cell_type) %>%
  dplyr::summarise(mean_score = mean(score)) %>%
  group_by(sample) %>%
  mutate(mean_score = (mean_score - mean(mean_score)) / sd(mean_score))

fig3C_3 <- tibble(fig3C_3)

fig3C_3$cell_type <- droplevels(fig3C_3$cell_type)

addWorksheet(source_data, "Figure 3B-3")

writeData(source_data, sheet = "Figure 3B-3", 
          x = "Figure 3B-3", 
          startCol =  1, startRow = 1)

writeData(source_data, sheet = 'Figure 3B-3',
          x = fig3C_3, startCol = 1, startRow = 3)

pval <- fig3C_3 %>% 
  wilcox_test(
    mean_score ~ cell_type
  ) %>% 
  add_xy_position() %>% 
  mutate(p = sub("e","%.% 10^",p) )


fig3C_3_gg <- ggplot(data = fig3C_3, aes(x = cell_type, y = mean_score)) +
  theme_classic()+
  scale_fill_brewer(palette = "Dark2") +
  geom_violin(lwd=0.5, width = 1, aes(fill=cell_type)) +
  geom_boxplot(color = "black", fill = 'white', width = .15, 
               outlier.alpha = 0, lwd=0.5) +
  geom_jitter(colour="black",fill='white', size=1, alpha=1, 
              width = 0.2, shape=21) +
  theme(plot.title = element_text(size=22, hjust = 0.5),
        axis.title.x = element_text(size=22, face="bold"),
        axis.title.y = element_text(size=22),
        axis.text.x = element_text(colour="black", size=22),
        axis.text.y = element_text(colour="black", size=22),
        legend.position="none", axis.text=element_text(size=12),
        legend.title = element_text(size=14),
        legend.key.size = unit(0.6, "cm"),
        legend.text = element_text(size=14),
        text = element_text(family = 'Calibri')) +
  scale_x_discrete(labels = c("Fibroblast adventitial" = "Adventitial", 
                              "Fibroblast alveolar" = "Alveolar", 
                              "Fibroblast peribronchial" = "Peribronchial")) +
  labs(title = "Salcher et al fibroblasts (n=65)\n",
       x = NULL, y = "Hanley et al alv fib enrichment score\n") +
  add_pvalue(pval,label = 'p', step.increase = 0.1,
             family = 'Calibri', label.size = 6, parse = TRUE,
             bracket.nudge.y = 0.1) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))

fig3C_3_gg$layers[[4]]$aes_params$family <- "Calibri"

fig3C_3_gg

ggsave(here('figures','fig3C_3.tiff'),
       plot = fig3C_3_gg,
       device = 'tiff',
       width = desired_width,
       height = desired_height,
       units = 'in', dpi = TIFF_dpi)





########
# FIG 3D
########

# load Hanley et al data (downloaded from https://zenodo.org/records/7400873)

load(here('data', 'IntegratedFibs_Zenodo.Rdata'))

class(Fibs.integrated)

Fibs.integrated <- subset(x = Fibs.integrated, subset = Sample.Subtype == 'LUAD')

Fibs.integrated <- subset(x = Fibs.integrated, subset = Sample.type != 'Control')

DefaultAssay(Fibs.integrated) <- 'RNA'

Fibs.integrated@meta.data$COL1A1_log_counts <- as.numeric(Fibs.integrated@assays$RNA['COL1A1',])


# VI cluster 2 expression 

Fibs.integrated <- AddModuleScore(Fibs.integrated, features = VI_genesets,
                                  name = '_score')

names(Fibs.integrated@meta.data)[grep("_score", names(Fibs.integrated@meta.data))] <- names(VI_genesets)


Idents(Fibs.integrated) <- Fibs.integrated$Fibs_MajorClusters

tmp <- subset(Fibs.integrated, idents = c('Myo','Alveolar','Adventitial'))


fig3D <- tmp@meta.data %>% select(SampleID,            
                                      cell_type = Fibs_MajorClusters,
                                      score = VI_gene_cluster_2)


fig3D <- fig3D %>% 
  group_by(SampleID, cell_type) %>%                     
  summarise(mean_score = mean(score, na.rm = TRUE), .groups = "drop") %>% 
  group_by(SampleID) %>%                               
  mutate(mean_score = as.numeric(scale(mean_score))) %>% 
  ungroup()

fig3D <- tibble(fig3D)

pval <- fig3D %>% 
  wilcox_test(mean_score ~ cell_type) %>% 
  add_xy_position() %>% 
  mutate(p = sub("e","%.% 10^",p) )

fig3D$cell_type <- factor(fig3D$cell_type, levels = c('Adventitial',
                                                              'Alveolar',
                                                              'Myo'))

addWorksheet(source_data, "Figure 3D")

writeData(source_data, sheet = "Figure 3D", 
          x = "Figure 3D", 
          startCol =  1, startRow = 1)

writeData(source_data, sheet = 'Figure 3D',
          x = fig3D, startCol = 1, startRow = 3)

fig3D_gg <- ggplot(fig3D,
       aes(x = cell_type, y = mean_score, fill = cell_type)) +
  geom_violin(width = 1, lwd = 0.5) +
  geom_boxplot(width = .15, fill = "white", colour = "black",
               lwd = 0.5, outlier.alpha = 0) +
  geom_jitter(shape = 21, size = 1, width = 0.2,
              colour = "black", fill = "white") +
  scale_x_discrete(labels = c(
    "Myo"="Myofibroblast",
    "Adventitial"="Adventitial",
    "Alveolar"="Alveolar")) +
  scale_fill_manual(
    values = c("Myo" = "#7570B3",
               "Alveolar"      = "#D95F02",
               "Adventitial"   = "#1B9E77")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  labs(title = "Hanley et al LUAD scRNA-seq dataset (n=46)\n",
       x = NULL,
       y = "Scaled VI gene cluster 2 enrichment score\n") +
  theme_classic(base_family = "Calibri") +
  theme(plot.title   = element_text(size = 22, hjust = 0.5),
        axis.title.x = element_text(size = 22, face = "bold"),
        axis.title.y = element_text(size = 22),
        axis.text.x  = element_text(size = 22, colour = "black"),
        axis.text.y  = element_text(size = 22, colour = "black"),
        legend.position = "none") +
  add_pvalue(pval,label = 'p', step.increase = 0.1,
             family = 'Calibri', label.size = 6, parse = TRUE, inherit.aes=FALSE)

fig3D_gg

ggsave(here('figures','fig3D.tiff'),
       plot = fig3D_gg,
       device = 'tiff',
       width = desired_width,
       height = desired_height,
       units = 'in', dpi = TIFF_dpi)


######################
# EXTENDED DATA FIG 6A
######################

efig6A <- tmp@meta.data %>% select(SampleID,            
                                      cell_type = Fibs_MajorClusters,
                                      score = COL1A1_log_counts)


efig6A <- efig6A %>% 
  group_by(SampleID, cell_type) %>%                     
  summarise(mean_score = mean(score, na.rm = TRUE), .groups = "drop") %>% 
  group_by(SampleID) %>%                               
  mutate(mean_score = as.numeric(scale(mean_score))) %>% 
  ungroup()

efig6A <- tibble(efig6A)

pval <- efig6A %>% 
  wilcox_test(mean_score ~ cell_type) %>% 
  add_xy_position() %>% 
  mutate(p = sub("e","%.% 10^",p) )

efig6A$cell_type <- factor(efig6A$cell_type, levels = c('Adventitial',
                                                              'Alveolar',
                                                              'Myo'))


addWorksheet(source_data, "E. Figure 6A")

writeData(source_data, sheet = "E. Figure 6A", 
          x = "E. Figure 6A", 
          startCol =  1, startRow = 1)

writeData(source_data, sheet = 'E. Figure 6A',
          x = efig6A, startCol = 1, startRow = 3)


efig6A_gg <- ggplot(efig6A,
       aes(x = cell_type, y = mean_score, fill = cell_type)) +
  geom_violin(width = 1, lwd = 0.5) +
  geom_boxplot(width = .15, fill = "white", colour = "black",
               lwd = 0.5, outlier.alpha = 0) +
  geom_jitter(shape = 21, size = 1, width = 0.2,
              colour = "black", fill = "white") +
  scale_x_discrete(labels = c(
    "Myo"="Myofibroblast",
    "Adventitial"="Adventitial",
    "Alveolar"="Alveolar")) +
  scale_fill_manual(
    values = c("Myo" = "#7570B3",
               "Alveolar"      = "#D95F02",
               "Adventitial"   = "#1B9E77")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  labs(title = "Hanley et al LUAD scRNA-seq dataset (n=46)\n",
       x = NULL,
       y = "Scaled COL1A1 log counts\n") +
  theme_classic(base_family = "Calibri") +
  theme(plot.title   = element_text(size = 22, hjust = 0.5),
        axis.title.x = element_text(size = 22, face = "bold"),
        axis.title.y = element_text(size = 22),
        axis.text.x  = element_text(size = 22, colour = "black"),
        axis.text.y  = element_text(size = 22, colour = "black"),
        legend.position = "none") +
  add_pvalue(pval,label = 'p', step.increase = 0.1,
             family = 'Calibri', label.size = 6, parse = TRUE, inherit.aes=FALSE)


efig6A_gg

ggsave(here('figures','efig6A.tiff'),
       plot = efig6A_gg,
       device = 'tiff',
       width = desired_width,
       height = desired_height,
       units = 'in', dpi = TIFF_dpi)




######################
# EXTENDED DATA FIG 5D
######################

# heatmap enrichment of the VI cluster signatures across all the cell types 

efig5D <- Salcher_lung_cancer_atlas_LUAD_seurat

efig5D <- AddModuleScore(efig5D, features = VI_genesets, name = '_score')

names(efig5D@meta.data)[grep("_score", names(efig5D@meta.data))] <- names(VI_genesets)


data_long <- gather(efig5D@meta.data, cluster, score, VI_gene_cluster_1, VI_gene_cluster_2,
                    VI_gene_cluster_3, VI_gene_cluster_4, factor_key=TRUE)

# Compute the mean score for each cell type and cluster
data_summary <- data_long %>%
  group_by(cluster, cell_type) %>%
  summarize(mean_score = mean(score, na.rm = TRUE), .groups = "drop")

efig5D <- data_summary %>%
  group_by(cluster) %>%
  ungroup() %>%
  filter(cell_type %in% remaining_cell_types)

addWorksheet(source_data, "E. Figure 5E")

writeData(source_data, sheet = "E. Figure 5E", 
          x = "E. Figure 5E", 
          startCol =  1, startRow = 1)

writeData(source_data, sheet = 'E. Figure 5E',
          x = efig5D, startCol = 1, startRow = 3)


efig5D_gg <- ggplot(efig5D, aes(x = cluster, y = cell_type, fill = mean_score)) +
    geom_tile(color = 'white', lwd = 1, linetype = 1) +
    scale_fill_gradient2(
      low     = "steelblue",  # color at your min (‑10)
      mid     = "white",      # color at 0
      high    = "darkred",    # color at your max (40)
      midpoint = 0,
      limits  = c(-10, 55),
      space   = "Lab"
    ) +
    theme(plot.title = element_text(size=22, hjust = 0.5),
          axis.title.x = element_text(size=22, face="bold"),
          axis.title.y = element_text(size=22),
          axis.text.x = element_text(colour='black',angle = 45, hjust = 1, size=14),
          axis.text.y = element_text(colour="black", size=14),
          legend.position="right", axis.text=element_text(size=12),
          legend.title = element_text(size=14),
          legend.key.size = unit(0.6, "cm"),
          legend.text = element_text(size=14),
          text = element_text(family = 'Calibri')) +
    labs(
      title = "Enrichment of VI gene clusters in scRNA-seq data\nfrom stage I LUAD (Salcher et al)",
      x = NULL,
      y = NULL,
      fill = "Mean module score"
    ) +
    coord_flip()

efig5D_gg

plot_dimensions <- dev.size("in")

aspect_ratio <- plot_dimensions[1] / plot_dimensions[2]

aspect_ratio <- 2.945

desired_height <- 4  

desired_width <- desired_height * aspect_ratio

TIFF_dpi <- 600

ggsave(here('figures','efig5D.tiff'),
       plot = efig5D_gg,
       device = 'tiff',
       width = desired_width,
       height = desired_height,
       units = 'in', dpi = TIFF_dpi)





#####################
# CIBERSORTX ANALYSIS
#####################

# use cibersortx on bulk data to validate cell deconvolution findings of stRNA-seq data

# create signature genes file

set.seed(123)

Salcher_lung_cancer_atlas_LUAD_seurat_subsampled <- Salcher_lung_cancer_atlas_LUAD_seurat[, sample(colnames(Salcher_lung_cancer_atlas_LUAD_seurat), size=25000, replace=F)]

Salcher_lung_cancer_atlas_LUAD_seurat_subsampled$cell_type_sample <- paste(Salcher_lung_cancer_atlas_LUAD_seurat_subsampled$sample,Salcher_lung_cancer_atlas_LUAD_seurat_subsampled$cell_type,sep='$')

Salcher_lung_cancer_atlas_LUAD_seurat_sce <- as.SingleCellExperiment(Salcher_lung_cancer_atlas_LUAD_seurat_subsampled)

Salcher_lung_cancer_atlas_LUAD_seurat_sce <- aggregateAcrossCells(Salcher_lung_cancer_atlas_LUAD_seurat_sce, 
                                                                  id=colData(Salcher_lung_cancer_atlas_LUAD_seurat_sce)[,c("cell_type_sample")],
                                                                  statistics = 'mean')

cibersortx_ref <- as.matrix(Salcher_lung_cancer_atlas_LUAD_seurat_sce@assays@data$counts)

colnames(cibersortx_ref) <- sub(".*\\$", "", colnames(cibersortx_ref))

cibersortx_ref <- as.data.frame(cibersortx_ref, colnames = TRUE)

cibersortx_ref <- cbind(GeneSymbol = rownames(cibersortx_ref), cibersortx_ref)

write.table(cibersortx_ref,
            here('data','cibersortx_salcher_ref.txt'), 
            sep = "\t", row.names = FALSE)


# create mixture file from discovery cohort

discovery_dge <- readRDS(here('data','discovery_dge.rds'))

cibersortx_mixture <- discovery_dge$counts %>%
  as.data.frame(.) %>%
  mutate(GeneSymbol = rownames(.)) %>%
  select(GeneSymbol, everything())

write.table(cibersortx_mixture,
            here('data','cibersortx_mixture_discovery.txt'), 
            sep = "\t", row.names = FALSE)


# read in cibersortx results for discovery cohort

######################
# EXTENDED DATA FIG 5B
######################

cibersortx_results <- read.csv(here('data','CIBERSORTx_Job11_Results.csv'))

rownames(cibersortx_results) <- cibersortx_results$Mixture

# put in same order as discovery dge

cibersortx_results <- cibersortx_results[rownames(discovery_dge$samples),]

cibersortx_results <- cibersortx_results %>%
  select(-c(Mixture, P.value, Correlation, RMSE))

colnames(cibersortx_results) <- paste(colnames(cibersortx_results),'cibersort',sep='_')

discovery_dge$samples <- cbind(discovery_dge$samples, cibersortx_results)

discovery_dge$samples$Plasma <- as.factor(discovery_dge$samples$Plasma)

efig5B <- discovery_dge$samples %>%
  filter(!is.na(Plasma)) %>%
  select(Plasma, Plasma.cell_cibersort)

levels(efig5B$Plasma) <- c('None','Scattered','Clusters','Sheets')

addWorksheet(source_data, "E. Figure 5B")

writeData(source_data, sheet = "E. Figure 5B", 
          x = "E. Figure 5B", 
          startCol =  1, startRow = 1)

writeData(source_data, sheet = 'E. Figure 5B',
          x = efig5B, startCol = 1, startRow = 3)

p_val <- kruskal.test(efig5B$Plasma.cell_cibersort ~ efig5B$Plasma)$p.value

pval <- grobTree(textGrob('Kruskal-Wallis, p =',
                          x=0.4,  y=0.9, 
                          gp=gpar(fontsize=17, 
                                  fontfamily = 'Calibri')),
                 textGrob(eval(bquote(scientific_10(p_val))),
                          x=0.7, y=0.9,
                          gp=gpar(fontsize=17, 
                                  fontfamily = 'Calibri')))

efig5B_gg <- ggplot(data = efig5B, aes(x = Plasma, y = Plasma.cell_cibersort, fill = Plasma)) +
  theme_classic() + 
  scale_fill_manual(values = c( "#B294C7",'#B294C7',"#B294C7","#B294C7")) +
  geom_violin(lwd=0.5, width = 1.) +
  geom_boxplot(color = "black", fill = 'white', width = .1, outlier.alpha = 0, lwd=0.5) +
  geom_jitter(colour="black",fill='white', size=1, alpha=1, width = 0.2, shape=21) +
  labs(title="Plasma cells in discovery cohort (n=102)\n", 
       x = '\nPathologist plasma cell grade', 
       y = "CIBERSORTx plasma cell proportion\n") +
  theme(
    plot.title = element_text(size=20, hjust = 0.5),
    axis.title.x = element_text(size=20),
    axis.title.y = element_text(size=20),
    axis.text.x = element_text(size=20,color = 'black'),
    legend.position="none", 
    axis.text=element_text(size=20),
    legend.title = element_text(),
    text = element_text(family = 'Calibri')) +  
  annotation_custom(pval)

efig5B_gg

# Get the current plot dimensions in inches

plot_dimensions <- dev.size("in")

# Calculate aspect ratio

aspect_ratio <- plot_dimensions[1] / plot_dimensions[2]

aspect_ratio <- 1.057692

# Define a reasonable height (in inches)

desired_height <- 6 

# Calculate corresponding width to maintain aspect ratio

desired_width <- desired_height * aspect_ratio

# Define the resolution

TIFF_dpi <- 600

# Save the plot with adjusted dimensions

ggsave(here('figures','efig5B.tiff'),
       plot = efig5B_gg,
       device = 'tiff',
       width = desired_width,
       height = desired_height,
       units = 'in', dpi = TIFF_dpi)



###########
# FIGURE 3E
###########

discovery_dge$samples$VI <- as.factor(discovery_dge$samples$VI)

levels(discovery_dge$samples$VI) <- c('VI-','VI+')

fig3E <- discovery_dge$samples %>%
  select(VI, Fibroblast.peribronchial_cibersort)

addWorksheet(source_data, "Figure 3E")

writeData(source_data, sheet = "Figure 3E", 
          x = "Figure 3E", 
          startCol =  1, startRow = 1)

writeData(source_data, sheet = 'Figure 3E',
          x = fig3E, startCol = 1, startRow = 3)

my_comparisons = list(c(1,2),c(2,3),c(1,3))

fig3E_gg <- ggplot(data = fig3E, 
                     aes(x = VI, y = Fibroblast.peribronchial_cibersort, fill = VI)) +
  theme_classic() + 
  scale_fill_manual(values = c('#390099','#FF0054')) +
  geom_violin(lwd=0.5, width = 0.5) +
  geom_boxplot(color = "black", fill = 'white', width = .15, 
               outlier.alpha = 0, lwd=0.5) +
  geom_jitter(colour="black",fill='white', size=1, 
              alpha=1, width = 0.2, shape=21) +
  labs(title="Peribronchial fibroblasts in bulk RNA-seq\ndiscovery cohort (n=103)\n", x = NULL, 
       y = "CIBERSORTx predicted proportion\n") +
  theme(
    plot.title = element_text(size=22, hjust = 0.5),
    axis.title.x = element_text(size=20),
    axis.title.y = element_text(size=20),
    axis.text.x = element_text(size=20),
    axis.text.y = element_text(size=20),
    legend.position="none", 
    legend.title = element_text(),
    text = element_text(family = 'Calibri')) +
  ylim(min(discovery_dge$samples$Fibroblast.peribronchial_cibersort),
       max(discovery_dge$samples$Fibroblast.peribronchial_cibersort)+0.02) +
  geom_bracket(
    xmin = 'VI-', xmax = "VI+", 
    y.position = max(discovery_dge$samples$Fibroblast.peribronchial_cibersort)+0.01,
    label = paste('p.adj =',
                  round(p.adjust(wilcox.test(Fibroblast.peribronchial_cibersort ~ VI, 
                                             data = discovery_dge$samples)$p.value, 
                                 n = ncol(cibersortx_results)),3)), 
    tip.length = c(0.02, 0.02), vjust = -0.25,
    size = 0.5,
    label.size = 6,
    family = 'Calibri',
    inherit.aes = FALSE)

fig3E_gg

# Save the plot with adjusted dimensions

ggsave(here('figures','fig3E.tiff'),
       plot = fig3E_gg,
       device = 'tiff',
       width = desired_width,
       height = desired_height,
       units = 'in', dpi = TIFF_dpi)


################################
# IHC / ISH validation of COL1A1
################################

# read in COL1A1+ ISH cells/spot from adjacent section

LM_SD_4_ISH_spot_measurements <- read.csv(here('data','LM_SD_4_ISH_spot_measurements_per_cell.csv'),
                                          sep = '\t')

# subset to spot IDs that passed visium qc

LM_SD_4_ISH_spot_measurements <- LM_SD_4_ISH_spot_measurements[which(LM_SD_4_ISH_spot_measurements$Name %in% rownames(visium_seurat_objects[[4]]@meta.data)),]

LM_SD_4_ISH_spot_measurements <- LM_SD_4_ISH_spot_measurements %>%
  rename(spot_id = Name) %>%
  rename(COL1A1_ISH_positive_cells = Positive.COL1A1.Cells) %>%
  rename(ISH_total_cells = Num.Detections) %>%
  select(spot_id, COL1A1_ISH_positive_cells, ISH_total_cells, ISH.staining..Positive.area.µm.2)

# read in COL1A1+ IHC cells/spot from adjacent section

LM_SD_4_IHC_spot_measurements <- read.csv(here('data','LM_SD_4_IHC_spot_measurements_per_cell.csv'),
                                          sep = '\t')

LM_SD_4_IHC_spot_measurements <- LM_SD_4_IHC_spot_measurements[which(LM_SD_4_IHC_spot_measurements$Name %in% rownames(visium_seurat_objects[[4]]@meta.data)),]

LM_SD_4_IHC_spot_measurements <- LM_SD_4_IHC_spot_measurements %>%
  rename(spot_id = Name) %>%
  rename(COL1A1_IHC_positive_cells = Positive.COL1A1.Cells) %>%
  rename(IHC_total_cells = Num.Detections) %>%
  select(spot_id, COL1A1_IHC_positive_cells, IHC_total_cells, IHC.staining..Positive.area.µm.2)


# read in THY+ IHC cells/spot from adjacent section

LM_SD_4_IHC_THY_spot_measurements <- read.csv(here('data','LM_SD_4_IHC_THY1_spot_measurements_per_cell.csv'),
                                              sep = '\t')

LM_SD_4_IHC_THY_spot_measurements <- LM_SD_4_IHC_THY_spot_measurements[which(LM_SD_4_IHC_THY_spot_measurements$Name %in% rownames(visium_seurat_objects[[4]]@meta.data)),]

LM_SD_4_IHC_THY_spot_measurements <- LM_SD_4_IHC_THY_spot_measurements %>%
  rename(spot_id = Name) %>%
  rename(THY1_IHC_positive_cells = Num.Positive) %>%
  rename(IHC_THY_total_cells = Num.Detections) %>%
  select(spot_id, THY1_IHC_positive_cells, IHC_THY_total_cells)


# read in H&E total cells/spot from original visium section 

LM_SD_4_HE_spot_measurements <- read.csv(here('data','LM_SD_4_HE_spot_measurements_per_cell.csv'),
                                         sep = '\t')

LM_SD_4_HE_spot_measurements <- LM_SD_4_HE_spot_measurements[which(LM_SD_4_HE_spot_measurements$Name %in% rownames(visium_seurat_objects[[4]]@meta.data)),]

LM_SD_4_HE_spot_measurements <- LM_SD_4_HE_spot_measurements %>%
  rename(spot_id = Name) %>%
  rename(HE_total_cells = Total.cells) %>%
  select(spot_id, HE_total_cells)


visium_seurat_objects[[4]]@meta.data$spot_id <- rownames(visium_seurat_objects[[4]]@meta.data)


LM_SD_4_all_stain_measurements <- merge(visium_seurat_objects[[4]]@meta.data,
                                        LM_SD_4_ISH_spot_measurements, 
                                        by = 'spot_id',
                                        all.x = TRUE)

LM_SD_4_all_stain_measurements <- merge(LM_SD_4_all_stain_measurements,
                                        LM_SD_4_IHC_spot_measurements, 
                                        by = 'spot_id',
                                        all.x = TRUE)

LM_SD_4_all_stain_measurements <- merge(LM_SD_4_all_stain_measurements,
                                        LM_SD_4_IHC_THY_spot_measurements, 
                                        by = 'spot_id',
                                        all.x = TRUE)

LM_SD_4_all_stain_measurements <- merge(LM_SD_4_all_stain_measurements,
                                        LM_SD_4_HE_spot_measurements, 
                                        by = 'spot_id',
                                        all.x = TRUE)




# read in cytospace cell type predictions

visium_seurat_objects_pred <- readRDS(here('data','visium_seurat_objects_clustered_named_predictions_cyto_salcher_stardist.rds'))


LM_SD_4_all_stain_measurements$Fibroblast.peribronchial <- visium_seurat_objects_pred[[4]]$Fibroblast.peribronchial

LM_SD_4_all_stain_measurements$Fibroblast.adventitial <- visium_seurat_objects_pred[[4]]$Fibroblast.adventitial

LM_SD_4_all_stain_measurements$Fibroblast.alveolar <- visium_seurat_objects_pred[[4]]$Fibroblast.alveolar

LM_SD_4_all_stain_measurements$cytospace_total_cells <- visium_seurat_objects_pred[[4]]$Total.cells



LM_SD_4_all_stain_measurements$COL1A1_ISH_positive_cells_frac <- LM_SD_4_all_stain_measurements$COL1A1_ISH_positive_cells / LM_SD_4_all_stain_measurements$ISH_total_cells

LM_SD_4_all_stain_measurements$COL1A1_IHC_positive_cells_frac <- LM_SD_4_all_stain_measurements$COL1A1_IHC_positive_cells / LM_SD_4_all_stain_measurements$IHC_total_cells

LM_SD_4_all_stain_measurements$THY1_IHC_positive_cells_frac <- LM_SD_4_all_stain_measurements$THY1_IHC_positive_cells / LM_SD_4_all_stain_measurements$IHC_THY_total_cells

LM_SD_4_all_stain_measurements$Fibroblast.peribronchial_frac <- LM_SD_4_all_stain_measurements$Fibroblast.peribronchial / LM_SD_4_all_stain_measurements$cytospace_total_cells

LM_SD_4_all_stain_measurements$Fibroblast.adventitial_frac <- LM_SD_4_all_stain_measurements$Fibroblast.adventitial / LM_SD_4_all_stain_measurements$cytospace_total_cells

LM_SD_4_all_stain_measurements$Fibroblast.alveolar_frac <- LM_SD_4_all_stain_measurements$Fibroblast.alveolar / LM_SD_4_all_stain_measurements$cytospace_total_cells



# 7 spot sliding window
coords <- visium_seurat_objects[[4]]@images[["slice1"]]@coordinates
head(coords)

coords$spot_id <- rownames(coords)

LM_SD_4_all_stain_measurements_fnn <- merge(LM_SD_4_all_stain_measurements, 
                                            coords, by = "spot_id")


knn_result <- get.knnx(
  data  = LM_SD_4_all_stain_measurements_fnn[, c("imagerow", "imagecol")],
  query = LM_SD_4_all_stain_measurements_fnn[, c("imagerow", "imagecol")],
  k = 7  # the spot itself + 6 neighbors
)

# Initialize
LM_SD_4_all_stain_measurements_fnn$COL1A1_ISH_positive_cells_frac_mean <- NA
LM_SD_4_all_stain_measurements_fnn$VI_gene_cluster_2_mean   <- NA

for (i in seq_len(nrow(LM_SD_4_all_stain_measurements_fnn))) {
  # Indices of the spot itself + 6 neighbors
  nbr_idx <- knn_result$nn.index[i, ]
  
  # If you want the mean:
  LM_SD_4_all_stain_measurements_fnn$COL1A1_ISH_positive_cells_frac_mean[i] <- mean(LM_SD_4_all_stain_measurements_fnn$COL1A1_ISH_positive_cells_frac[nbr_idx])
  LM_SD_4_all_stain_measurements_fnn$VI_gene_cluster_2_mean[i]   <- mean(LM_SD_4_all_stain_measurements_fnn$VI_gene_cluster_2[nbr_idx])
  
}


######################
# EXTENDED DATA FIG 6B
######################

efig6B <- LM_SD_4_all_stain_measurements_fnn %>%
  select(COL1A1_ISH_positive_cells_frac_mean, VI_gene_cluster_2_mean)

addWorksheet(source_data, "E. Figure 6B")

writeData(source_data, sheet = "E. Figure 6B", 
          x = "E. Figure 6B", 
          startCol =  1, startRow = 1)

writeData(source_data, sheet = 'E. Figure 6B',
          x = efig6B, startCol = 1, startRow = 3)


efig6B_gg <- ggplot(data = efig6B, 
       aes(x = COL1A1_ISH_positive_cells_frac_mean, y = VI_gene_cluster_2_mean)) +
  geom_point() +
  theme_classic() +
  labs(title="Sample 4 (VI+)\n", 
       x = "\nMean COL1A1 ISH+ cell fraction\n(7-spot smoothed)", 
       y = "VI gene cluster 2 Enrichment score\n(7-spot smoothed\n") +
  theme(
    plot.title = element_text(size=22, hjust = 0.5),
    axis.title.x = element_text(size=20),
    axis.title.y = element_text(size=20),
    legend.position="bottom", 
    axis.text=element_text(size=20),
    legend.title = element_text(size=20),
    text = element_text(family = 'Calibri')) +
  geom_smooth(method = "lm", se = TRUE, formula = y ~ x, color = "#FF0054") + 
  stat_cor(method="spearman",
           aes(label = paste(after_stat(r.label),
                             sub("e","%.% 10^",
                                 after_stat(p.label)), 
                             sep = "~`,`~")), 
           size = 6, family = 'Calibri', label.y = 0.9) 

efig6B_gg

plot_dimensions <- dev.size("in")

aspect_ratio <- plot_dimensions[1] / plot_dimensions[2]

aspect_ratio <- 0.9561934

desired_height <- 7  

desired_width <- desired_height * aspect_ratio

TIFF_dpi <- 600

ggsave(here('figures','efig6B.tiff'),
       plot = efig6B_gg,
       device = 'tiff',
       width = desired_width,
       height = desired_height,
       units = 'in', dpi = TIFF_dpi)



LM_SD_4_all_stain_measurements_fnn$COL1A1_IHC_positive_cells_frac_mean <- NA
LM_SD_4_all_stain_measurements_fnn$VI_gene_cluster_2_mean   <- NA

for (i in seq_len(nrow(LM_SD_4_all_stain_measurements_fnn))) {
  nbr_idx <- knn_result$nn.index[i, ]
  
  LM_SD_4_all_stain_measurements_fnn$COL1A1_IHC_positive_cells_frac_mean[i] <- mean(LM_SD_4_all_stain_measurements_fnn$COL1A1_IHC_positive_cells_frac[nbr_idx])
  LM_SD_4_all_stain_measurements_fnn$VI_gene_cluster_2_mean[i]   <- mean(LM_SD_4_all_stain_measurements_fnn$VI_gene_cluster_2[nbr_idx])
  
}

efig6C <- LM_SD_4_all_stain_measurements_fnn %>%
  select(COL1A1_IHC_positive_cells_frac_mean, VI_gene_cluster_2_mean)

addWorksheet(source_data, "E. Figure 6C")

writeData(source_data, sheet = "E. Figure 6C", 
          x = "E. Figure 6C", 
          startCol =  1, startRow = 1)

writeData(source_data, sheet = 'E. Figure 6C',
          x = efig6C, startCol = 1, startRow = 3)

######################
# EXTENDED DATA FIG 6C
######################

efig6C_gg <- ggplot(data = efig6C, aes(x = COL1A1_IHC_positive_cells_frac_mean, y = VI_gene_cluster_2_mean)) +
  geom_point() +
  theme_classic() +
  labs(title="Sample 4 (VI+)\n", x = "\nMean COL1A1 IHC+ cell fraction\n(7-spot smoothed)", y = "VI gene cluster 2 Enrichment score\n(7-spot smoothed\n") +
  theme(
    plot.title = element_text(size=22, hjust = 0.5),
    axis.title.x = element_text(size=20),
    axis.title.y = element_text(size=20),
    legend.position="bottom", 
    axis.text=element_text(size=20),
    legend.title = element_text(size=20),
    text = element_text(family = 'Calibri')) +
  geom_smooth(method = "lm", se = TRUE, formula = y ~ x, color = "#FF0054") + 
  stat_cor(method="spearman",
           aes(label = paste(after_stat(r.label),
                             sub("e","%.% 10^",
                                 after_stat(p.label)), 
                             sep = "~`,`~")), 
           size = 6, family = 'Calibri', label.y = 0.9) 


efig6C_gg

plot_dimensions <- dev.size("in")

aspect_ratio <- plot_dimensions[1] / plot_dimensions[2]

aspect_ratio <- 0.9561934

desired_height <- 7  

desired_width <- desired_height * aspect_ratio

TIFF_dpi <- 600

ggsave(here('figures','efig6C.tiff'),
       plot = efig6C_gg,
       device = 'tiff',
       width = desired_width,
       height = desired_height,
       units = 'in', dpi = TIFF_dpi)




LM_SD_4_all_stain_measurements_fnn$COL1A1_IHC_positive_cells_frac_mean <- NA
LM_SD_4_all_stain_measurements_fnn$THY1_IHC_positive_cells_frac_mean   <- NA

for (i in seq_len(nrow(LM_SD_4_all_stain_measurements_fnn))) {
  nbr_idx <- knn_result$nn.index[i, ]
  
  LM_SD_4_all_stain_measurements_fnn$COL1A1_IHC_positive_cells_frac_mean[i] <- mean(LM_SD_4_all_stain_measurements_fnn$COL1A1_ISH_positive_cells_frac[nbr_idx])
  LM_SD_4_all_stain_measurements_fnn$THY1_IHC_positive_cells_frac_mean[i]   <- mean(LM_SD_4_all_stain_measurements_fnn$THY1_IHC_positive_cells_frac[nbr_idx])
  
}

efig6D <- LM_SD_4_all_stain_measurements_fnn %>%
  select(COL1A1_IHC_positive_cells_frac_mean, THY1_IHC_positive_cells_frac_mean)

addWorksheet(source_data, "E. Figure 6D")

writeData(source_data, sheet = "E. Figure 6D", 
          x = "E. Figure 6D", 
          startCol =  1, startRow = 1)

writeData(source_data, sheet = 'E. Figure 6D',
          x = efig6D, startCol = 1, startRow = 3)

######################
# EXTENDED DATA FIG 6D
######################

efig6D_gg <- ggplot(data = LM_SD_4_all_stain_measurements_fnn, aes(x = COL1A1_IHC_positive_cells_frac_mean, y = THY1_IHC_positive_cells_frac_mean)) +
  geom_point() +
  theme_classic() +
  labs(title="Sample 4 (VI+)\n", x = "\nMean COL1A1 IHC+ cell fraction\n(7-spot smoothed)", y = "Mean THY1 IHC+ cell fraction\n(7-spot smoothed)\n") +
  theme(
    plot.title = element_text(size=22, hjust = 0.5),
    axis.title.x = element_text(size=20),
    axis.title.y = element_text(size=20),
    legend.position="bottom", 
    axis.text=element_text(size=20),
    legend.title = element_text(size=20),
    text = element_text(family = 'Calibri')) +
  geom_smooth(method = "lm", se = TRUE, formula = y ~ x, color = "#FF0054") + 
  stat_cor(method="spearman",
           aes(label = paste(after_stat(r.label),
                             sub("e","%.% 10^",
                                 after_stat(p.label)), 
                             sep = "~`,`~")), 
           size = 6, family = 'Calibri', label.y = 1.1) 

efig6D_gg

plot_dimensions <- dev.size("in")

aspect_ratio <- plot_dimensions[1] / plot_dimensions[2]

aspect_ratio <- 0.9561934

desired_height <- 7  

desired_width <- desired_height * aspect_ratio

TIFF_dpi <- 600

ggsave(here('figures','efig6D.tiff'),
       plot = efig6D_gg,
       device = 'tiff',
       width = desired_width,
       height = desired_height,
       units = 'in', dpi = TIFF_dpi)




######################
# EXTENDED DATA FIG 6F
######################

LM_SD_4_all_stain_measurements$desmoplastic_stroma <- as.factor(ifelse(LM_SD_4_all_stain_measurements$pathology == 'Desmoplastic Stroma',
                                                                       'Yes','No'))

efig6F <- LM_SD_4_all_stain_measurements %>%
  select(COL1A1_IHC_positive_cells_frac, desmoplastic_stroma)

pval <- efig6F %>% 
  wilcox_test(
    COL1A1_IHC_positive_cells_frac ~ desmoplastic_stroma
  ) %>% 
  add_xy_position() %>% 
  mutate(p = sub("e","%.% 10^",p) )

addWorksheet(source_data, "E. Figure 6F")

writeData(source_data, sheet = "E. Figure 6F", 
          x = "E. Figure 6F", 
          startCol =  1, startRow = 1)

writeData(source_data, sheet = 'E. Figure 6F',
          x = efig6F, startCol = 1, startRow = 3)


efig6F_gg <- ggplot(data = efig6F, aes(x = desmoplastic_stroma, y = COL1A1_IHC_positive_cells_frac)) +
  theme_classic() + 
  scale_fill_manual(values = c("#999999","#B16C29")) +
  geom_violin(lwd=0.5, width = 1, aes(fill = desmoplastic_stroma)) +
  geom_boxplot(color = "black", fill = 'white', width = .15, 
               outlier.alpha = 0, lwd=0.5) +
  labs(title="Sample 4\n", 
       x = "\nDesmoplastic Stroma spots", y = "COL1A1 IHC+ cell fraction\n") +
  theme(
    plot.title = element_text(size=22, hjust = 0.5),
    axis.title.x = element_text(size=20),
    axis.title.y = element_text(size=20),
    axis.text.x = element_text(size=20),
    legend.position="none", 
    axis.text=element_text(size=20),
    legend.title = element_text(),
    text = element_text(family = 'Calibri')) +
  scale_y_continuous(
    limits = c(0,1.2),
    breaks = seq(0,1,0.2),
    expand = expansion(mult = c(0,0))
  ) +
  add_pvalue(pval,label = 'p', step.increase = 0.1,
             family = 'Calibri', label.size = 6, parse = TRUE, y.position = 1.1)


efig6F_gg

plot_dimensions <- dev.size("in")

aspect_ratio <- plot_dimensions[1] / plot_dimensions[2]

aspect_ratio <- 0.9566004

desired_height <- 7

desired_width <- desired_height * aspect_ratio

TIFF_dpi <- 600

ggsave(here('figures','efig6F.tiff'),
       plot = efig6F_gg,
       device = 'tiff',
       width = desired_width,
       height = desired_height,
       units = 'in', dpi = TIFF_dpi)





# correlation of COL1A1 IHC/ISH with fibroblast types

LM_SD_4_all_stain_measurements_fnn$Fibroblast.peribronchial_frac_mean <- NA
LM_SD_4_all_stain_measurements_fnn$Fibroblast.adventitial_frac_mean   <- NA
LM_SD_4_all_stain_measurements_fnn$Fibroblast.alveolar_frac_mean      <- NA
LM_SD_4_all_stain_measurements_fnn$COL1A1_IHC_positive_cells_frac_mean <- NA
LM_SD_4_all_stain_measurements_fnn$COL1A1_ISH_positive_cells_frac_mean <- NA


LM_SD_4_all_stain_measurements_fnn <- 
  LM_SD_4_all_stain_measurements_fnn[complete.cases(LM_SD_4_all_stain_measurements_fnn[ ,c("imagerow","imagecol")]), ]

knn_result <- get.knnx(
  data  = LM_SD_4_all_stain_measurements_fnn[, c("imagerow", "imagecol")],
  query = LM_SD_4_all_stain_measurements_fnn[, c("imagerow", "imagecol")],
  k = 7
)


for (i in seq_len(nrow(LM_SD_4_all_stain_measurements_fnn))) {
  nbr_idx <- knn_result$nn.index[i, ]
  
  LM_SD_4_all_stain_measurements_fnn$Fibroblast.peribronchial_frac_mean[i] <- 
    mean(LM_SD_4_all_stain_measurements_fnn$Fibroblast.peribronchial_frac[nbr_idx], na.rm = TRUE)
  
  LM_SD_4_all_stain_measurements_fnn$Fibroblast.adventitial_frac_mean[i] <- 
    mean(LM_SD_4_all_stain_measurements_fnn$Fibroblast.adventitial_frac[nbr_idx], na.rm = TRUE)
  
  LM_SD_4_all_stain_measurements_fnn$Fibroblast.alveolar_frac_mean[i] <- 
    mean(LM_SD_4_all_stain_measurements_fnn$Fibroblast.alveolar_frac[nbr_idx], na.rm = TRUE)
  
  # Smoothed IHC fraction
  LM_SD_4_all_stain_measurements_fnn$COL1A1_IHC_positive_cells_frac_mean[i] <-
    mean(LM_SD_4_all_stain_measurements_fnn$COL1A1_IHC_positive_cells_frac[nbr_idx], na.rm = TRUE)
  
  # Smoothed ISH fraction
  LM_SD_4_all_stain_measurements_fnn$COL1A1_ISH_positive_cells_frac_mean[i] <-
    mean(LM_SD_4_all_stain_measurements_fnn$COL1A1_ISH_positive_cells_frac[nbr_idx], na.rm = TRUE)
  
}

efig6E <- LM_SD_4_all_stain_measurements_fnn %>%
  select(COL1A1_IHC_positive_cells_frac_mean, COL1A1_ISH_positive_cells_frac_mean,
         Fibroblast.peribronchial_frac_mean, Fibroblast.alveolar_frac_mean,
         Fibroblast.adventitial_frac_mean)

addWorksheet(source_data, "E. Figure 6E")

writeData(source_data, sheet = "E. Figure 6E", 
          x = "E. Figure 6E", 
          startCol =  1, startRow = 1)

writeData(source_data, sheet = 'E. Figure 6E',
          x = efig6E, startCol = 1, startRow = 3)

df_cor <- data.frame(
  comparison = c("IHC_COL1A1 vs. Peribronchial",
                 "IHC_COL1A1 vs. Adventitial",
                 "IHC_COL1A1 vs. Alveolar",
                 "ISH_COL1A1 vs. Peribronchial",
                 "ISH_COL1A1 vs. Adventitial",
                 "ISH_COL1A1 vs. Alveolar"),
  spearman_rho = NA_real_,
  p_value       = NA_real_
)

df <- LM_SD_4_all_stain_measurements_fnn

get_spearman_stats <- function(x, y) {
  test_res <- cor.test(x, y, method = "spearman", use = "complete.obs")
  list(rho = test_res$estimate, p = test_res$p.value)
}

res1 <- get_spearman_stats(df$COL1A1_IHC_positive_cells_frac_mean,
                           df$Fibroblast.peribronchial_frac_mean)
df_cor$spearman_rho[1] <- res1$rho
df_cor$p_value[1]      <- res1$p

res2 <- get_spearman_stats(df$COL1A1_IHC_positive_cells_frac_mean,
                           df$Fibroblast.adventitial_frac_mean)
df_cor$spearman_rho[2] <- res2$rho
df_cor$p_value[2]      <- res2$p

res3 <- get_spearman_stats(df$COL1A1_IHC_positive_cells_frac_mean,
                           df$Fibroblast.alveolar_frac_mean)
df_cor$spearman_rho[3] <- res3$rho
df_cor$p_value[3]      <- res3$p

res4 <- get_spearman_stats(df$COL1A1_ISH_positive_cells_frac_mean,
                           df$Fibroblast.peribronchial_frac_mean)
df_cor$spearman_rho[4] <- res4$rho
df_cor$p_value[4]      <- res4$p

res5 <- get_spearman_stats(df$COL1A1_ISH_positive_cells_frac_mean,
                           df$Fibroblast.adventitial_frac_mean)
df_cor$spearman_rho[5] <- res5$rho
df_cor$p_value[5]      <- res5$p

res6 <- get_spearman_stats(df$COL1A1_ISH_positive_cells_frac_mean,
                           df$Fibroblast.alveolar_frac_mean)
df_cor$spearman_rho[6] <- res6$rho
df_cor$p_value[6]      <- res6$p

df_cor

df_cor <- df_cor %>%
  mutate(
    assay = ifelse(grepl("^IHC", comparison), "IHC", "ISH")
  ) %>%
  arrange(assay, desc(spearman_rho)) %>%
  mutate(
    comparison = factor(comparison, levels = unique(comparison)),
    sig_label  = ifelse(p_value < 0.01, "*", ""),
    fibro_type = case_when(
      grepl("Peribronchial", comparison) ~ "Peribronchial",
      grepl("Alveolar",      comparison) ~ "Alveolar",
      grepl("Adventitial",   comparison) ~ "Adventitial"
    ),
    offset_y   = ifelse(spearman_rho >= 0,
                        spearman_rho + 0.03, 
                        spearman_rho - 0.03)
  )

df_cor <- df_cor %>%
  mutate(
    assay = ifelse(grepl("^IHC", comparison), "COL1A1 IHC", "COL1A1 ISH"),
    comparison = factor(comparison, levels = comparison)  
  )


my_labels <- c(
  "IHC_COL1A1 vs. Peribronchial" = "COL1A1 IHC+ ~ Peri",
  "IHC_COL1A1 vs. Adventitial"   = "COL1A1 IHC+ ~ Adv",
  "IHC_COL1A1 vs. Alveolar"      = "COL1A1 IHC+ ~ Alv",
  "ISH_COL1A1 vs. Peribronchial" = "COL1A1 ISH+ ~ Peri",
  "ISH_COL1A1 vs. Adventitial"   = "COL1A1 ISH+ ~ Adv",
  "ISH_COL1A1 vs. Alveolar"      = "COL1A1 ISH+ ~ Alv"
)

df_cor <- df_cor %>%
  mutate(
    comparison = factor(
      comparison,
      levels = names(my_labels),  # full, ordered set
      labels = my_labels           # the “pretty” names
    )
  )


n <- nrow(df)

# three IHC correlations
rIHC <- list(
  peri = cor(df$COL1A1_IHC_positive_cells_frac_mean,
             df$Fibroblast.peribronchial_frac_mean,
             method="spearman", use="complete.obs"),
  adv  = cor(df$COL1A1_IHC_positive_cells_frac_mean,
             df$Fibroblast.adventitial_frac_mean,
             method="spearman", use="complete.obs"),
  alv  = cor(df$COL1A1_IHC_positive_cells_frac_mean,
             df$Fibroblast.alveolar_frac_mean,
             method="spearman", use="complete.obs")
)

# three ISH correlations
rISH <- list(
  peri = cor(df$COL1A1_ISH_positive_cells_frac_mean,
             df$Fibroblast.peribronchial_frac_mean,
             method="spearman", use="complete.obs"),
  adv  = cor(df$COL1A1_ISH_positive_cells_frac_mean,
             df$Fibroblast.adventitial_frac_mean,
             method="spearman", use="complete.obs"),
  alv  = cor(df$COL1A1_ISH_positive_cells_frac_mean,
             df$Fibroblast.alveolar_frac_mean,
             method="spearman", use="complete.obs")
)

# the three overlaps (same for both assays)
r23 <- list(
  peri_adv = cor(df$Fibroblast.peribronchial_frac_mean,
                 df$Fibroblast.adventitial_frac_mean,
                 method="spearman", use="complete.obs"),
  peri_alv = cor(df$Fibroblast.peribronchial_frac_mean,
                 df$Fibroblast.alveolar_frac_mean,
                 method="spearman", use="complete.obs"),
  adv_alv  = cor(df$Fibroblast.adventitial_frac_mean,
                 df$Fibroblast.alveolar_frac_mean,
                 method="spearman", use="complete.obs")
)

pv <- function(r1, r2, r12) {
  r.test(n, r1, r2, r12)$p
}

p_IHC_peri_adv <- pv(rIHC$peri, rIHC$adv,  r23$peri_adv)
p_IHC_peri_alv <- pv(rIHC$peri, rIHC$alv,  r23$peri_alv)
p_IHC_adv_alv  <- pv(rIHC$adv,  rIHC$alv,  r23$adv_alv)

p_ISH_peri_adv <- pv(rISH$peri, rISH$adv,  r23$peri_adv)
p_ISH_peri_alv <- pv(rISH$peri, rISH$alv,  r23$peri_alv)
p_ISH_adv_alv  <- pv(rISH$adv,  rISH$alv,  r23$adv_alv)

my_labels <- c(
  "IHC_COL1A1 vs. Peribronchial" = "COL1A1 IHC+ ~ Peri",
  "IHC_COL1A1 vs. Adventitial"   = "COL1A1 IHC+ ~ Alv",
  "IHC_COL1A1 vs. Alveolar"      = "COL1A1 IHC+ ~ Adv",
  "ISH_COL1A1 vs. Peribronchial" = "COL1A1 ISH+ ~ Peri",
  "ISH_COL1A1 vs. Adventitial"   = "COL1A1 ISH+ ~ Alv",
  "ISH_COL1A1 vs. Alveolar"      = "COL1A1 ISH+ ~ Adv"
)

plot_lvls <- unname(my_labels)

max_rho <- df_cor %>% group_by(assay) %>% summarise(max_rho = max(spearman_rho))

# assemble
stat.test <- tribble(
  ~assay, ~comp1, ~comp2, ~p_val,
  "COL1A1 IHC", "COL1A1 IHC+ ~ Peri", "COL1A1 IHC+ ~ Alv", p_IHC_peri_adv,
  "COL1A1 IHC", "COL1A1 IHC+ ~ Peri", "COL1A1 IHC+ ~ Adv", p_IHC_peri_alv,
  "COL1A1 IHC", "COL1A1 IHC+ ~ Alv", "COL1A1 IHC+ ~ Adv", p_IHC_adv_alv,
  "COL1A1 ISH", "COL1A1 ISH+ ~ Peri", "COL1A1 ISH+ ~ Alv", p_ISH_peri_adv,
  "COL1A1 ISH", "COL1A1 ISH+ ~ Peri", "COL1A1 ISH+ ~ Adv", p_ISH_peri_alv,
  "COL1A1 ISH", "COL1A1 ISH+ ~ Alv", "COL1A1 ISH+ ~ Adv", p_ISH_adv_alv
) %>%
  left_join(max_rho, by = "assay") %>%
  group_by(assay) %>%
  mutate(
    y.position = max_rho + seq(0.02, by = 0.03, length.out = 3)
  ) %>%
  ungroup() %>%
  mutate(
    comp1 = factor(comp1, levels = plot_lvls),
    comp2 = factor(comp2, levels = plot_lvls),
    assay = factor(assay, levels = c("COL1A1 IHC", "COL1A1 ISH"))
  )

sup_digits <- c(
  "0"="⁰","1"="¹","2"="²","3"="³","4"="⁴",
  "5"="⁵","6"="⁶","7"="⁷","8"="⁸","9"="⁹"
)

make_superscript <- function(e) {
  digs <- strsplit(as.character(abs(e)), "")[[1]]
  s    <- paste0(sup_digits[digs], collapse = "")
  if (e < 0) s <- paste0("⁻", s)
  s
}

stat.test <- stat.test %>%
  rename(
    group1 = comp1,
    group2 = comp2
  ) %>%
  mutate(
    sci = format(signif(p_val, 3), scientific = TRUE),
    parts = strsplit(sci, "e"),
    mant  = vapply(parts, `[`, 1, FUN.VALUE = ""),
    expo  = as.integer(vapply(parts, `[`, 2, FUN.VALUE = "")),
    p_label = ifelse(
      p_val < 0.001,
      paste0(mant, " × 10", vapply(expo, make_superscript, FUN.VALUE="")),
      as.character(signif(p_val, 3))
    )
  )

new_labels <- c(
  "COL1A1 IHC+ ~ Peri" = "~ Peri",
  "COL1A1 IHC+ ~ Alv"  = "~ Alv",
  "COL1A1 IHC+ ~ Adv"  = "~ Adv",
  "COL1A1 ISH+ ~ Peri" = "~ Peri",
  "COL1A1 ISH+ ~ Alv"  = "~ Alv",
  "COL1A1 ISH+ ~ Adv"  = "~ Adv"
)

efig6E_gg <- ggplot(df_cor, aes(x = comparison, y = spearman_rho, fill = fibro_type)) +
  geom_col() +
  scale_x_discrete(labels = my_labels) +
  facet_grid(. ~ assay,
             switch = "x",
             scales = "free_x",
             space  = "free_x") +
  coord_cartesian(ylim = c(-0.15, max(df_cor$spearman_rho) + 0.15),
                  clip = "off") +
  scale_y_continuous(breaks = seq(-0.1, 0.4, 0.1)) +
  scale_fill_manual(values = c(
    "Peribronchial" = "#7570B3",
    "Alveolar"      = "#D95F02",
    "Adventitial"   = "#1B9E77"
  )) +
  theme_classic() +
  theme(
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text = element_text(size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
    axis.text.y  = element_text(size = 16),
    axis.title.x = element_text(size = 16, margin = margin(t = 10)),
    axis.title.y = element_text(size = 16, margin = margin(r = 10)),
    plot.margin = margin(t=10, r=10, b=80, l=10),
    legend.title = element_text(size = 14),
    legend.text  = element_text(size = 14),
    legend.position = 'right',
    legend.box = 'horizontal'
  ) +
  scale_x_discrete(
    breaks = names(new_labels),  
    labels = new_labels          
  ) +
  labs(
    x = "\nMarker-positive cell fraction per spot comparison",
    y = "Spearman correlation coefficient\n(7‑spot smoothed)\n",
    fill = "Deconvoluted\nfibroblast"
  ) +
  stat_pvalue_manual(
    data    = stat.test,
    label   = "p_label",        # now contains e.g. "7.9 × 10⁻¹⁷"
    mapping = aes(
      xmin       = group1,
      xmax       = group2,
      y.position = y.position
    ),
    tip.length   = 0.01,
    bracket.size = 0.5,
    size         = 5,
    inherit.aes  = FALSE,
    step.increase = 0.02
  )

efig6E_gg

plot_dimensions <- dev.size("in")

aspect_ratio <- plot_dimensions[1] / plot_dimensions[2]

aspect_ratio <- 0.9178082

desired_height <- 8

desired_width <- desired_height * aspect_ratio

TIFF_dpi <- 600

ggsave(here('figures','efig6E.tiff'),
       plot = efig6E_gg,
       device = 'tiff',
       width = desired_width,
       height = desired_height,
       units = 'in', dpi = TIFF_dpi)








# comparison of VI predictor expression between VI- and distal VI 

###########
# FIGURE 5B
###########

visium_seurat_objects <- readRDS(file = here('data','visium_seurat_objects_lognormal_clustered_annotated_distal.rds'))


visium_seurat_merge_tmp <- merge(x = visium_seurat_objects[[1]], 
                                 y = visium_seurat_objects[-1], 
                                 add.cell.ids = names(visium_seurat_objects))


# add additional metadata to seurat object

visium_seurat_merge_tmp@meta.data$column_name <- rownames(visium_seurat_merge_tmp@meta.data)

geo_stRNAseq_metadata$predominant_VI <- ifelse(geo_stRNAseq_metadata$`novel grade` == 'VI',1,0)

geo_stRNAseq_metadata_sub <- geo_stRNAseq_metadata[,which(colnames(geo_stRNAseq_metadata) %!in% colnames(visium_seurat_merge_tmp@meta.data)[-which(colnames(visium_seurat_merge_tmp@meta.data) == 'orig.ident')])]

tmp <- merge(visium_seurat_merge_tmp@meta.data, geo_stRNAseq_metadata_sub, by = "orig.ident")

tmp <- tmp[match(colnames(visium_seurat_merge_tmp), tmp$column_name),]

visium_seurat_merge_tmp@meta.data <- cbind(visium_seurat_merge_tmp@meta.data, tmp)


# remove sample LM_SD_6 for low mean counts due to tissue loss during workflow

Idents(visium_seurat_merge_tmp) <- visium_seurat_merge_tmp$orig.ident

visium_seurat_merge_tmp <- subset(visium_seurat_merge_tmp, idents = 'LM_SD_6', invert = TRUE)

visium_seurat_merge_tmp <- AddModuleScore(visium_seurat_merge_tmp, features = VI_genesets,
                                          name = '_score')

names(visium_seurat_merge_tmp@meta.data)[grep("_score", names(visium_seurat_merge_tmp@meta.data))] <- names(VI_genesets)


tmp <- visium_seurat_merge_tmp

tmp@meta.data$VI.focus <- as.factor(tmp@meta.data$VI.focus)

tmp <- tmp@meta.data %>% dplyr::filter(VI.focus %!in% c('Proximal (<1mm)'))


tmp$VI_predictor_genes <- tmp$VI_predictor_genes_up - tmp$VI_predictor_genes_dn

fig5B <- gather(tmp, cluster, score, VI_predictor_genes, factor_key=TRUE)

names(fig5B)[names(fig5B) == "predominant_VI"] <- "VI Status"



fig5B <- fig5B %>% filter(cluster == 'VI_predictor_genes')

data <- data.frame(y = fig5B$score, x = fig5B$`VI Status`, 
                   sample = fig5B$`sample id`)

m <- lmer(y ~ 1 + x + (1 | sample), data = data)

# type II anova 
Anova(m, type = 2)

mixed_effects_stats <- data.frame(Anova(m, type = 2, test.statistic = 'Chisq'))

mixed_effects_stats <- mixed_effects_stats %>% select(Pr..Chisq.)

fig5B$mixed_effects_p <- mixed_effects_stats$Pr..Chisq.


# scale by cluster (for visualization purposes)

fig5B <- fig5B %>%
  group_by(cluster) %>%
  mutate(score = (score - mean(score)) / sd(score))

fig5B$`VI Status` <- as.factor(fig5B$`VI Status`)

levels(fig5B$`VI Status`) <- c('VI-','VI+')


addWorksheet(source_data, "Figure 5B")

writeData(source_data, sheet = "Figure 5B", 
          x = "Figure 5B", 
          startCol =  1, startRow = 1)

writeData(source_data, sheet = 'Figure 5B',
          x = fig5B, startCol = 1, startRow = 3)


pval <- grobTree(textGrob(paste('Anova, p =',
                                round(mixed_effects_stats$Pr..Chisq., 
                                      digits = 3)),
                          x=0.3,  y=0.9, hjust=0,
                          gp=gpar(fontsize=20, 
                                  fontfamily = 'Calibri')))

fig5B_gg <- ggplot(data = fig5B, aes(x = `VI Status`, y = score, fill = `VI Status`)) +
  theme_classic() + 
  coord_cartesian(clip = "off") +
  scale_fill_manual(values = c('#390099','#9E0059')) +
  geom_violin(lwd=0.5, width = 0.5) +
  geom_boxplot(color = "black", fill = 'white', width = .15, outlier.alpha = 0, lwd=0.5) +
  theme(plot.title = element_text(size=20, hjust = 0.5),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        legend.position="none", 
        text = element_text(family = 'Calibri')) +
  labs(title = 'VI predictor in stRNA-seq spots\n',x = NULL, y = "VI predictor enrichment score\n") +
  ylim(min(fig5B$score),max(fig5B$score)+1) +
  annotation_custom(pval) +
  scale_x_discrete(labels = c("VI-\n(n=8 samples)","Distal VI+\n(n=7 samples)"))

fig5B_gg

ggsave(here('figures','fig5B.tiff'),
       plot = fig5B_gg,
       device = 'tiff',
       width = desired_width,
       height = desired_height,
       units = 'in', dpi = TIFF_dpi)



# save source data file

saveWorkbook(source_data, 
             here('data','steiner_source_data_stRNAseq.xlsx'), 
             overwrite = TRUE)


sessionInfo()
