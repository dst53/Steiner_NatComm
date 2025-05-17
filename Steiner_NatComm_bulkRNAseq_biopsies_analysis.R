#' ---
#' title: Bulk RNA-seq analysis of biopsies for "Identification of a gene expression 
#' signature of vascular invasion and recurrence in stage I lung adenocarcinoma 
#' via bulk and spatial transcriptomics"
#' output: github_document
#' ---

# expected run time < 1 min

#######
# SETUP
#######

# load R packages (all installed under R 4.2.1)

library(here)
library(openxlsx)
library(h2o)
library(readxl)
library(dplyr)
library(edgeR)
library(tableone)
library(tibble)
library(sva)
library(pROC)
library(grid)
library(ggplot2)
library(stringr)
library(tidyr)
library(ggpubr)

source(here("utils.R"))

# define source data file for figures and tables

source_data <- createWorkbook()

# Start the H2O cluster (locally) for machine learning

h2o.init()

###########################################
# VALIDATION IN BIOPSIES/MATCHED RESECTIONS
###########################################

# RNA-seq of 24 retrospective stage I LUAD biopsies sequenced from Inova
# 12 of which had RNA_seq on matched resections
# 36 samples total

# read processed count matrix from GEO (placeholder)

biopsies_and_resections_se <- read.table(here('data','biopsies_exp_count.txt'),
                                 header = TRUE, row.names = 1,
                                 check.names = FALSE)

# read GEO metadata for biopsies - placeholder until GEO series is public

geo_bulkRNAseq_metadata <- read_excel(here('data','Steiner_NatComm_bulkRNAseq_biopsies_GEO_submission.xlsx'), 
                                      sheet = "Metadata", range = c('A33:P69'))


# load supp clinical data for biopsy cohort

steiner_patient_metadata_biopsies <- read.csv(here('data','steiner_patient_metadata_biopsies.csv'), 
                                     row.names = 1, check.names = FALSE)


geo_bulkRNAseq_metadata <- merge(geo_bulkRNAseq_metadata, 
                                 steiner_patient_metadata_biopsies, 
                                 by.x = '*library name',
                                 by.y = 'sample_id')

# put samples in same order

rownames(geo_bulkRNAseq_metadata) <- geo_bulkRNAseq_metadata$`*library name`

geo_bulkRNAseq_metadata <- geo_bulkRNAseq_metadata[colnames(biopsies_and_resections_se),]

biopsies_and_resections_dge <- DGEList(counts = biopsies_and_resections_se,
                                       samples = geo_bulkRNAseq_metadata) 


head(biopsies_and_resections_dge$samples$lib.size)


###############
# NORMALIZATION
###############

# Trimmed mean of m values (TMM) normalization

biopsies_and_resections_dge <- calcNormFactors(biopsies_and_resections_dge) 

head(biopsies_and_resections_dge$samples$lib.size)

biopsies_and_resections_dge$samples <- biopsies_and_resections_dge$samples %>%
  mutate(`Smoking status` = Smoking.status,
         `TNM stage (8th edition)` = TNM.8th.edition.stage,
         LMPVI = novel.grade
  )


##########
# TABLE S4
##########

# table 1 for biopsy cohort (summarized by LMPVI grade)

biopsies_and_resections_cohort_tableone <- CreateTableOne(data = biopsies_and_resections_dge$samples, 
                                                          vars = c("TNM stage (8th edition)","Age","Gender",
                                                                   "Race","Smoking status",
                                                                   "Procedure"),
                                                          strata = 'LMPVI')

tableS4 <- print(biopsies_and_resections_cohort_tableone)

write.csv(tableS4, file = here('data','tableS4.csv'))


################
# GENE FILTERING
################

gene_annotations <- readRDS(here('data','gene_annotations.rds'))

# convert gene names to hgnc

biopsies_and_resections_dge <- ensemble_to_hgnc(biopsies_and_resections_dge)

# Filter uppsala counts to filtered discovery set counts

discovery_dge <- readRDS(here('data','discovery_dge.rds'))

biopsies_and_resections_dge$counts <- biopsies_and_resections_dge$counts[which(rownames(biopsies_and_resections_dge$counts) %in% 
                                                                                 rownames(discovery_dge$counts)),]

dim(biopsies_and_resections_dge$counts)

discovery_dge$counts <- discovery_dge$counts[which(rownames(discovery_dge$counts) %in% 
                                                     rownames(biopsies_and_resections_dge$counts)),]

dim(discovery_dge$counts)


# put genes in same order as train

biopsies_and_resections_dge$counts <- biopsies_and_resections_dge$counts[rownames(discovery_dge$counts),]


biopsies_and_resections_dge <- DGEList(counts = biopsies_and_resections_dge$counts, 
                                       samples = biopsies_and_resections_dge$samples)

head(biopsies_and_resections_dge$samples$lib.size)

biopsies_and_resections_dge <- calcNormFactors(biopsies_and_resections_dge)


# Reference combat uppsala count matrix with discovery count matrix

biopsies_and_resections_dge$counts <- ComBat(cbind(cpm(discovery_dge$counts,log=T), 
                                                   cpm(biopsies_and_resections_dge$counts,log=T)),
                                             c(rep(1,ncol(discovery_dge$counts)),
                                               rep(2,ncol(biopsies_and_resections_dge$counts))),
                                             ref.batch = 1)

biopsies_and_resections_dge$counts <- biopsies_and_resections_dge$counts[,(ncol(discovery_dge$counts)+1):ncol(biopsies_and_resections_dge$counts)]

VI_predictor_genes <- readRDS(here('data','VI_predictor_genes.rds'))

# 48/48 genes present

length(which(VI_predictor_genes %in% rownames(biopsies_and_resections_dge$counts))) 


############################
# GENERATE MODEL PREDICTIONS
############################

# use final VI predictor model to predict

test_biopsies_and_resection <- data.frame(t(biopsies_and_resections_dge$counts[which(rownames(biopsies_and_resections_dge$counts) %in% 
                                                                                       VI_predictor_genes),]))

# convert to h2o format and add dummy class variable

test_biopsies_and_resection$response <- rep(c(1,0),length.out = ncol(biopsies_and_resections_dge$counts))

test_biopsies_and_resection <- as.h2o(test_biopsies_and_resection)

# Identify predictors and response

y <- "response"

x <- setdiff(names(test_biopsies_and_resection), y)

# For binary classification, response should be a factor

# test_biopsies_and_resection[, y] <- as.factor(test_biopsies_and_resection[, y])

# to ensure absolute reproducibility, load exact model trained

aml_final_leader <- h2o.loadModel(here('data','GLM_1_AutoML_1_20230802_104952'))

pred_test <- h2o.predict(aml_final_leader, test_biopsies_and_resection)  

pred_test <- data.frame(as.matrix(pred_test))

biopsies_and_resections_dge$samples$VI_predictor <- as.numeric(pred_test$p1)


# AUROC of predicting VI status at resection from biopsies using VI predictor

########
# FIG 5F
########

fig5F <- biopsies_and_resections_dge$samples %>%
  mutate(VI = ifelse(LMPVI == "VI", 1, 0)) %>%
  select(VI, VI_predictor)

addWorksheet(source_data, "Figure 5F")

writeData(source_data, sheet = "Figure 5F", 
          x = "Figure 5F", 
          startCol =  1, startRow = 1)

writeData(source_data, sheet = 'Figure 5F',
          x = fig5F, startCol = 1, startRow = 3,
          rowNames = TRUE)

gene_roc <- roc(as.factor(fig5F$VI), fig5F$VI_predictor,
                levels=c(0, 1), ci = TRUE)

gene_roc


ciobj <- ci.se(gene_roc, specificities=seq(0, 1, l=25))

dat.ci <- data.frame(x = as.numeric(rownames(ciobj)),
                     lower = ciobj[, 1],
                     upper = ciobj[, 3])

rocobj <- pROC::plot.roc(gene_roc, print.thres = TRUE, print.auc = TRUE)

auc <- round(pROC::auc(gene_roc),2)


p_val <- round(wilcox.test(fig5F$VI_predictor ~ fig5F$VI)$p.value,9)

pval <- grobTree(textGrob(paste0('VI predictor ',
                                 '(AUROC = ', auc, ')',sep=''),
                          x=0.7,  y=0.25, 
                          gp=gpar(fontsize=17, 
                                  fontfamily = 'Calibri')),
                 textGrob("p = ", 
                          x = 0.65, y = 0.15, 
                          gp = gpar(fontsize = 17, fontfamily = 'Calibri')),
                 textGrob(eval(bquote(scientific_10(p_val))),
                          x=0.8, y=0.15,
                          gp=gpar(fontsize=17, 
                                  fontfamily = 'Calibri')))

fig5F_gg <- pROC::ggroc(rocobj, colour = '#FF0054', size = 1, alpha = 1) +
  theme_classic() +
  annotation_custom(pval) +
  labs(title='Stage I LUAD presurgical biopsy\ncohort (n=24)\n', 
       x = "\nSpecificity", y = "Sensitivity\n") +
  theme(plot.title = element_text(size=22, hjust = 0.5),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        text = element_text(family = 'Calibri')) +
  geom_abline(slope=1, intercept = 1, linetype = "dashed", 
              alpha=0.7, color = "grey") + 
  geom_ribbon(data = dat.ci, aes(x = x, ymin = lower, ymax = upper), 
              fill = '#FF0054', alpha= 0.1) 


fig5F_gg

plot_dimensions <- dev.size("in")

aspect_ratio <- plot_dimensions[1] / plot_dimensions[2]

aspect_ratio <- 0.9122807

desired_height <- 6.5 

desired_width <- desired_height * aspect_ratio

TIFF_dpi <- 600

ggsave(here('figures','fig5F.tiff'),
       plot = fig5F_gg,
       device = 'tiff',
       width = desired_width,
       height = desired_height,
       units = 'in', dpi = TIFF_dpi)



# correlation of VI predictor score between matched biopsies and resections

########
# FIG 5G
########

biopsies_and_resections_dge$samples$sample_id <- rownames(biopsies_and_resections_dge$samples)

biopsies_and_resections_dge$samples <- biopsies_and_resections_dge$samples %>% 
  mutate(base = str_remove(sample_id, "[BR]$")) %>% 
  group_by(base) %>% 
  mutate(
    hasB = any(str_ends(sample_id, "B")),
    hasR = any(str_ends(sample_id, "R")),
    matched = if_else(hasB & hasR, "yes", "no")
  ) %>% 
  ungroup() %>% 
  select(-base, -hasB, -hasR)

fig5G <- biopsies_and_resections_dge$samples %>% 
  mutate(
    base = str_remove(sample_id, "[BR]$"),
    suffix = str_sub(sample_id, -1)                   
  ) %>% 
  group_by(base) %>% 
  filter(any(suffix == "B") & any(suffix == "R")) %>% 
  ungroup() %>% 
  select(base, suffix, VI_predictor) %>% 
  pivot_wider(
    names_from = suffix,
    values_from = VI_predictor,
    names_glue = "{suffix}_VI_predictor"                   
  )

addWorksheet(source_data, "Figure 5G")

writeData(source_data, sheet = "Figure 5G", 
          x = "Figure 5G", 
          startCol =  1, startRow = 1)

writeData(source_data, sheet = 'Figure 5G',
          x = fig5G, startCol = 1, startRow = 3)


fig5G_gg <- ggplot(data = fig5G, aes(x = B_VI_predictor, y = R_VI_predictor)) +
  geom_point() +
  theme_classic() +
  labs(title="Matched stage I LUAD biopsies &\nresections (n=12)\n", 
       x = "\nBiospy VI predictor score", y = "Resection VI predictor score\n") +
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
           size = 6, family = 'Calibri') 

fig5G_gg


plot_dimensions <- dev.size("in")

aspect_ratio <- plot_dimensions[1] / plot_dimensions[2]

aspect_ratio <- 0.9815436

desired_height <- 6.5

desired_width <- desired_height * aspect_ratio

TIFF_dpi <- 600

ggsave(here('figures','fig5G.tiff'),
       plot = fig5G_gg,
       device = 'tiff',
       width = desired_width,
       height = desired_height,
       units = 'in', dpi = TIFF_dpi)

h2o.shutdown(prompt = FALSE)

# save source data file

saveWorkbook(source_data, 
             here('data','steiner_source_data_bulkRNAseq_biopsies.xlsx'), 
             overwrite = TRUE)


sessionInfo()



