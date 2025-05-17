#' ---
#' title: Bulk RNA-seq analysis for "Identification of a gene expression 
#' signature of vascular invasion and recurrence in stage I lung adenocarcinoma 
#' via bulk and spatial transcriptomics"
#' output: github_document
#' ---

# run with Rstudio Source

# expected run time is 10 min
# or 1-2 hours for re-running cross validation (recommend min of 12 CPUs and 16gb RAM per core)

#######
# SETUP
#######

# load R packages (all installed under R 4.2.1)

library(readxl)
library(dplyr)
library(survival)
library(survminer)
library(sva)
library(pheatmap)
library(edgeR)
library(GSVA)
library(pROC)
library(plotROC)
library(ggsci)
library(rms)
library(SummarizedExperiment)
library(ggpubr)
library(caret)
library(h2o)
library(broom)
library(msigdbr)
library(ggcorrplot)
library(ggfortify)
library(data.table)
library(robustbase)
library(plyr)
library(tidyr)
library(enrichR)
library(phenoTest)
library(tibble)
library(estimate)
library(CePa)
library(ggrepel)
library(tableone)
library(RColorBrewer)
library(corrplot)
library(rcompanion)
library(ComplexUpset)
library(clusterProfiler)
library(cmprsk)
library(modelsummary)
library(forcats)
library(enrichplot)
library(grid)
library(scales)
library(rstatix)
library(ggprism)
library(TCGAbiolinks)
library(fst)
library(stringr)
library(GEOquery)
library(openxlsx)
library(here)
library(tidycmprsk)
library(ggsurvfit)
library(ConsensusClusterPlus)

source(here("utils.R"))

# define source data file for figures and tables

source_data <- createWorkbook()


##############
# DATA LOADING
##############

# load supp clinical data for discovery and validation cohorts

steiner_patient_metadata <- read.csv(here('data','steiner_patient_metadata.csv'), 
                                     row.names = 1, check.names = FALSE)

# split supp clinical annotations by discovery and validation cohorts

discovery_patient_metadata <- steiner_patient_metadata %>% 
  filter(cohort == 'discovery')

# 192 discovery tumors

dim(discovery_patient_metadata)

validation_patient_metadata <- steiner_patient_metadata %>% 
  filter(cohort == 'validation')

# 61 validation tumors from 60 patients

dim(validation_patient_metadata)




# outcome analysis of pathology features across stage I LUAD discovery cohort

######################
# EXTENDED DATA FIG 1A
######################

efig1A <- discovery_patient_metadata %>%
  select(Sample, time_7_year_RFS, RFS_7_year, LMPVI)

addWorksheet(source_data, "E. Figure 1A")

writeData(source_data, sheet = "E. Figure 1A", 
          x = "E. Figure 1A", 
          startCol =  1, startRow = 1)

writeData(source_data, sheet = 'E. Figure 1A',
          x = efig1A, startCol = 1, startRow = 3)

mod.cox <- coxph(formula = Surv(time_7_year_RFS, RFS_7_year) ~ LMPVI,
                 data=efig1A)

summary(mod.cox)

surv_object <- Surv(time = efig1A$time_7_year_RFS,event = efig1A$RFS_7_year)

fit1 <- survfit(surv_object ~ LMPVI, data = efig1A, conf.type = "log-log")

diff.1 <- survdiff(surv_object ~ LMPVI, data = efig1A)

pval <- grobTree(textGrob('Log-rank  p =',
                          x=0.2,  y=0.15, 
                          gp=gpar(fontsize=17, 
                                  fontfamily = 'Calibri')),
                 textGrob(eval(bquote(scientific_10(round(diff.1$pvalue,11)))),
                          x=0.5, y=0.16,
                          gp=gpar(fontsize=17, 
                                  fontfamily = 'Calibri')))

efig1A_gg <- ggsurvplot(fit1, data = efig1A, pval = F, censored=T, 
                        ggtheme = theme_classic(base_family = 'Calibri',
                                                base_size = 20),
                        legend = "none", title = "Novel grade\n",
                        legend.title = "",
                        legend.labs = c("LMP","NST","VI"),
                        xlab = "\nMonths",
                        ylab = "RFS (%)\n",
                        risk.table = TRUE,
                        tables.theme = theme_cleantable(base_family = 'Calibri',
                                                        base_size = 20,
                                                        plot.title = element_text(size = 2)),
                        pval.coord = c(20,0.1),
                        pval.size = 6,
                        pval.method = TRUE,
                        pval.method.coord = c(0.5,0.1),
                        font.family = 'Calibri',
                        fontsize = 6,
                        palette = c("#1f78b4","#33a02c","#FF0054"))

efig1A_gg$table$theme$text$size <- 20

efig1A_gg$table <- efig1A_gg$table + theme(plot.title = element_text(size=20))

efig1A_gg$plot <- efig1A_gg$plot + theme(plot.title = element_text(size=22,
                                                                   hjust=0.5)) +
  annotation_custom(pval)

efig1A_gg 

# add method to grid.draw for ggsurvplot

grid.draw.ggsurvplot <- function(x){
  survminer:::print.ggsurvplot(x, newpage = FALSE)
}

# Get the current plot dimensions in inches
plot_dimensions <- dev.size("in")

# Calculate aspect ratio
aspect_ratio <- plot_dimensions[1] / plot_dimensions[2]

aspect_ratio <- 0.7628399

# Define a reasonable height (in inches)
desired_height <- 7 

# Calculate corresponding width to maintain aspect ratio
desired_width <- desired_height * aspect_ratio

# Define the resolution
TIFF_dpi <- 600

# Save the plot with adjusted dimensions
ggsave(here('figures','efig1A.tiff'),
       plot = efig1A_gg,
       device = 'tiff',
       width = desired_width,
       height = desired_height,
       units = 'in', dpi = TIFF_dpi)


######################
# EXTENDED DATA FIG 1B
######################

efig1B <- discovery_patient_metadata %>% 
  mutate(WHO = factor(WHO, levels = c('M','AISMIA','1','2','3'))) %>%
  filter(WHO != 'M') %>%
  select(Sample, time_7_year_RFS, RFS_7_year, WHO)

addWorksheet(source_data, "E. Figure 1B")

writeData(source_data, sheet = "E. Figure 1B", 
          x = "E. Figure 1B", 
          startCol =  1, startRow = 1)

writeData(source_data, sheet = 'E. Figure 1B',
          x = efig1B, startCol = 1, startRow = 3)


mod.cox <- coxph(formula = Surv(time_7_year_RFS, RFS_7_year) ~ WHO,
                 data=efig1B)

summary(mod.cox)

surv_object <- Surv(time = efig1B$time_7_year_RFS,event = efig1B$RFS_7_year)

fit1 <- survfit(surv_object ~ WHO, data = efig1B, conf.type = "log-log")

diff.1 <- survdiff(surv_object ~ WHO, data = efig1B)

efig1B_gg <- ggsurvplot(fit1, data = efig1B, pval = T, censored=T, 
                        ggtheme = theme_classic(base_family = 'Calibri',
                                                base_size = 20),
                        legend = "none", title = "WHO 2021 grade\n",
                        legend.title = "",
                        legend.labs = c("AIS/MIA","G1","G2","G3"),
                        xlab = "\nMonths",
                        ylab = "RFS (%)\n",
                        risk.table = TRUE,
                        tables.theme = theme_cleantable(base_family = 'Calibri',
                                                        base_size = 20,
                                                        plot.title = element_text(size = 2)),
                        pval.coord = c(20,0.1),
                        pval.size = 6,
                        pval.method = TRUE,
                        pval.method.coord = c(0.5,0.1),
                        font.family = 'Calibri',
                        fontsize = 6,
                        palette = c("grey","#1f78b4","#33a02c","#FF0054"))

efig1B_gg$table$theme$text$size <- 20

efig1B_gg$table <- efig1B_gg$table + theme(plot.title = element_text(size=20))

efig1B_gg$plot <- efig1B_gg$plot + theme(plot.title = element_text(size=22,
                                                                   hjust=0.5))

efig1B_gg

plot_dimensions <- dev.size("in")

aspect_ratio <- plot_dimensions[1] / plot_dimensions[2]

aspect_ratio <- 0.7921348

desired_height <- 7.5

desired_width <- desired_height * aspect_ratio

TIFF_dpi <- 600

ggsave(here('figures','efig1B.tiff'),
       plot = efig1B_gg,
       device = 'tiff',
       width = desired_width,
       height = desired_height,
       units = 'in', dpi = TIFF_dpi)


######################
# EXTENDED DATA FIG 1C
######################

efig1C <- discovery_patient_metadata %>%
  select(Sample, time_7_year_RFS, RFS_7_year, VI)

addWorksheet(source_data, "E. Figure 1C")

writeData(source_data, sheet = "E. Figure 1C", 
          x = "E. Figure 1C", 
          startCol =  1, startRow = 1)

writeData(source_data, sheet = 'E. Figure 1C',
          x = efig1C, startCol = 1, startRow = 3)

surv_object <- Surv(time = efig1C$time_7_year_RFS,event = efig1C$RFS_7_year)

fit1 <- survfit(surv_object ~ VI, data = efig1C, conf.type = "log-log")

diff.1 <- survdiff(surv_object ~ VI, data = efig1C)

pval <- grobTree(textGrob('Log-rank  p =',
                          x=0.2,  y=0.15, 
                          gp=gpar(fontsize=17, 
                                  fontfamily = 'Calibri')),
                 textGrob(eval(bquote(scientific_10(round(diff.1$pvalue,11)))),
                          x=0.5, y=0.16,
                          gp=gpar(fontsize=17, 
                                  fontfamily = 'Calibri')))

efig1C_gg <- ggsurvplot(fit1, data = efig1C, pval = F, censored=T, 
                        ggtheme = theme_classic(base_family = 'Calibri',
                                                base_size = 20),
                        legend = "none", title = "Vascular invasion\n",
                        legend.title = "",
                        legend.labs = c("VI-","VI+"),
                        xlab = "\nMonths",
                        ylab = "RFS (%)\n",
                        risk.table = TRUE,
                        tables.theme = theme_cleantable(base_family = 'Calibri',
                                                        base_size = 20,
                                                        plot.title = element_text(size = 2)),
                        pval.coord = c(20,0.1),
                        pval.size = 6,
                        pval.method = TRUE,
                        pval.method.coord = c(0.5,0.1),
                        font.family = 'Calibri',
                        fontsize = 6,
                        palette = c("#999999","#FF0054"))

efig1C_gg$table$theme$text$size <- 20

efig1C_gg$table <- efig1C_gg$table + theme(plot.title = element_text(size=20))

efig1C_gg$plot <- efig1C_gg$plot + theme(plot.title = element_text(size=22,
                                                                   hjust=0.5)) +
  annotation_custom(pval)

efig1C_gg


plot_dimensions <- dev.size("in")

aspect_ratio <- plot_dimensions[1] / plot_dimensions[2]

aspect_ratio <- 0.7628399

desired_height <- 7 

desired_width <- desired_height * aspect_ratio

TIFF_dpi <- 600

ggsave(here('figures','efig1C.tiff'),
       plot = efig1C_gg,
       device = 'tiff',
       width = desired_width,
       height = desired_height,
       units = 'in', dpi = TIFF_dpi)



# Cox regression result reported in text

mod.cox <- coxph(formula = Surv(time_7_year_RFS, RFS_7_year) ~ VI,
                 data=efig1C)

summary(mod.cox)


######################
# EXTENDED DATA FIG 1D
######################

efig1D <- discovery_patient_metadata %>%
  select(Sample, time_7_year_RFS, RFS_7_year, STAS)

addWorksheet(source_data, "E. Figure 1D")

writeData(source_data, sheet = "E. Figure 1D", 
          x = "E. Figure 1D", 
          startCol =  1, startRow = 1)

writeData(source_data, sheet = 'E. Figure 1D',
          x = efig1D, startCol = 1, startRow = 3)

surv_object <- Surv(time = efig1D$time_7_year_RFS,event = efig1D$RFS_7_year)

fit1 <- survfit(surv_object ~ STAS, data = efig1D, conf.type = "log-log")

diff.1 <- survdiff(surv_object ~ STAS, data = efig1D)

efig1D_gg <- ggsurvplot(fit1, data = efig1D, pval = T, censored=T, 
                        ggtheme = theme_classic(base_family = 'Calibri',
                                                base_size = 20),
                        legend = "none", title = "Spread through air spaces\n",
                        legend.title = "",
                        legend.labs = c("STAS-","STAS+"),
                        xlab = "\nMonths",
                        ylab = "RFS (%)\n",
                        risk.table = TRUE,
                        tables.theme = theme_cleantable(base_family = 'Calibri',
                                                        base_size = 20,
                                                        plot.title = element_text(size = 2)),
                        pval.coord = c(20,0.1),
                        pval.size = 6,
                        pval.method = TRUE,
                        pval.method.coord = c(0.5,0.1),
                        font.family = 'Calibri',
                        fontsize = 6,
                        palette = c("#999999","#CB8CAD"))

efig1D_gg$table$theme$text$size <- 20

efig1D_gg$table <- efig1D_gg$table + theme(plot.title = element_text(size=20))

efig1D_gg$plot <- efig1D_gg$plot + theme(plot.title = element_text(size=22,
                                                                   hjust=0.5))
efig1D_gg

plot_dimensions <- dev.size("in")

aspect_ratio <- plot_dimensions[1] / plot_dimensions[2]

aspect_ratio <- 0.7921348

desired_height <- 7.5 

desired_width <- desired_height * aspect_ratio

TIFF_dpi <- 600

ggsave(here('figures','efig1D.tiff'),
       plot = efig1D_gg,
       device = 'tiff',
       width = desired_width,
       height = desired_height,
       units = 'in', dpi = TIFF_dpi)


######################
# EXTENDED DATA FIG 1E
######################

efig1E <- discovery_patient_metadata %>%
  select(Sample, time_7_year_RFS, RFS_7_year, VPI)

addWorksheet(source_data, "E. Figure 1E")

writeData(source_data, sheet = "E. Figure 1E", 
          x = "E. Figure 1E", 
          startCol =  1, startRow = 1)

writeData(source_data, sheet = 'E. Figure 1E',
          x = efig1E, startCol = 1, startRow = 3)

surv_object <- Surv(time = efig1E$time_7_year_RFS,event = efig1E$RFS_7_year)

fit1 <- survfit(surv_object ~ VPI, data = efig1E, conf.type = "log-log")

diff.1 <- survdiff(surv_object ~ VPI, data = efig1E)

efig1E_gg <- ggsurvplot(fit1, data = efig1E, pval = T, censored=T, 
                        ggtheme = theme_classic(base_family = 'Calibri',
                                                base_size = 20),
                        legend = "none", title = "Visceral pleural invasion\n",
                        legend.title = "",
                        legend.labs = c("VPI-","VPI+"),
                        xlab = "\nMonths",
                        ylab = "RFS (%)\n",
                        risk.table = TRUE,
                        tables.theme = theme_cleantable(base_family = 'Calibri',
                                                        base_size = 20,
                                                        plot.title = element_text(size = 2)),
                        pval.coord = c(20,0.1),
                        pval.size = 6,
                        pval.method = TRUE,
                        pval.method.coord = c(0.5,0.1),
                        font.family = 'Calibri',
                        fontsize = 6,
                        palette = c("#999999","#874F6F"))

efig1E_gg$table$theme$text$size <- 20

efig1E_gg$table <- efig1E_gg$table + theme(plot.title = element_text(size=20))

efig1E_gg$plot <- efig1E_gg$plot + theme(plot.title = element_text(size=22,
                                                                   hjust=0.5))
efig1E_gg

plot_dimensions <- dev.size("in")

aspect_ratio <- plot_dimensions[1] / plot_dimensions[2]

aspect_ratio <- 0.7921348

desired_height <- 7.5  

desired_width <- desired_height * aspect_ratio

TIFF_dpi <- 600

ggsave(here('figures','efig1E.tiff'),
       plot = efig1E_gg,
       device = 'tiff',
       width = desired_width,
       height = desired_height,
       units = 'in', dpi = TIFF_dpi)


######################
# EXTENDED DATA FIG 1F
######################

efig1F <- discovery_patient_metadata %>%
  select(Sample, time_7_year_RFS, RFS_7_year, LI)

addWorksheet(source_data, "E. Figure 1F")

writeData(source_data, sheet = "E. Figure 1F", 
          x = "E. Figure 1F", 
          startCol =  1, startRow = 1)

writeData(source_data, sheet = 'E. Figure 1F',
          x = efig1F, startCol = 1, startRow = 3)

surv_object <- Surv(time = efig1F$time_7_year_RFS,event = efig1F$RFS_7_year)

fit1 <- survfit(surv_object ~ LI, data = efig1F, conf.type = "log-log")

diff.1 <- survdiff(surv_object ~ LI, data = efig1F)

efig1F_gg <- ggsurvplot(fit1, data = efig1F, pval = T, censored=T, 
                        ggtheme = theme_classic(base_family = 'Calibri',
                                                base_size = 20),
                        legend = "none", title = "Lymphatic invasion\n",
                        legend.title = "",
                        legend.labs = c("LI-","LI+"),
                        xlab = "\nMonths",
                        ylab = "RFS (%)\n",
                        risk.table = TRUE,
                        tables.theme = theme_cleantable(base_family = 'Calibri',
                                                        base_size = 20,
                                                        plot.title = element_text(size = 2)),
                        pval.coord = c(20,0.1),
                        pval.size = 6,
                        pval.method = TRUE,
                        pval.method.coord = c(0.5,0.1),
                        font.family = 'Calibri',
                        fontsize = 6,
                        palette = c("#999999","#1A5276"))

efig1F_gg$table$theme$text$size <- 20

efig1F_gg$table <- efig1F_gg$table + theme(plot.title = element_text(size=20))

efig1F_gg$plot <- efig1F_gg$plot + theme(plot.title = element_text(size=22,
                                                                   hjust=0.5))
efig1F_gg

plot_dimensions <- dev.size("in")

aspect_ratio <- plot_dimensions[1] / plot_dimensions[2]

aspect_ratio <- 0.7921348

desired_height <- 7.5

desired_width <- desired_height * aspect_ratio

TIFF_dpi <- 600

ggsave(here('figures','efig1F.tiff'),
       plot = efig1F_gg,
       device = 'tiff',
       width = desired_width,
       height = desired_height,
       units = 'in', dpi = TIFF_dpi)


######################
# EXTENDED DATA FIG 1G
######################

efig1G <- discovery_patient_metadata %>%
  mutate(LVI = case_when(VI == 1 | LI == 1 ~ 1, TRUE ~ 0)) %>%
  select(Sample, time_7_year_RFS, RFS_7_year, LVI)

addWorksheet(source_data, "E. Figure 1G")

writeData(source_data, sheet = "E. Figure 1G", 
          x = "E. Figure 1G", 
          startCol =  1, startRow = 1)

writeData(source_data, sheet = 'E. Figure 1G',
          x = efig1G, startCol = 1, startRow = 3)

surv_object <- Surv(time = efig1G$time_7_year_RFS,event = efig1G$RFS_7_year)

fit1 <- survfit(surv_object ~ LVI, data = efig1G, conf.type = "log-log")

diff.1 <- survdiff(surv_object ~ LVI, data = efig1G)

pval <- grobTree(textGrob('Log-rank  p =',
                          x=0.2,  y=0.15, 
                          gp=gpar(fontsize=17, 
                                  fontfamily = 'Calibri')),
                 textGrob(eval(bquote(scientific_10(round(diff.1$pvalue,11)))),
                          x=0.5, y=0.16,
                          gp=gpar(fontsize=17, 
                                  fontfamily = 'Calibri')))

efig1G_gg <- ggsurvplot(fit1, data = efig1G, pval = F, censored=T, 
                        ggtheme = theme_classic(base_family = 'Calibri',
                                                base_size = 20),
                        legend = "none", title = "Lymphovascular invasion\n",
                        legend.title = "",
                        legend.labs = c("LVI-","LVI+"),
                        xlab = "\nMonths",
                        ylab = "RFS (%)\n",
                        risk.table = TRUE,
                        tables.theme = theme_cleantable(base_family = 'Calibri',
                                                        base_size = 20,
                                                        plot.title = element_text(size = 2)),
                        pval.coord = c(20,0.1),
                        pval.size = 6,
                        pval.method = TRUE,
                        pval.method.coord = c(0.5,0.1),
                        font.family = 'Calibri',
                        fontsize = 6,
                        palette = c("#999999","#8c2965"))

efig1G_gg$table$theme$text$size <- 20

efig1G_gg$table <- efig1G_gg$table + 
  theme(plot.title = element_text(size=20))

efig1G_gg$plot <- efig1G_gg$plot + 
  theme(plot.title = element_text(size=22, hjust=0.5)) +
  annotation_custom(pval)

efig1G_gg


plot_dimensions <- dev.size("in")

aspect_ratio <- plot_dimensions[1] / plot_dimensions[2]

aspect_ratio <- 0.7921348

desired_height <- 7.5

desired_width <- desired_height * aspect_ratio

TIFF_dpi <- 600

ggsave(here('figures','efig1G.tiff'),
       plot = efig1G_gg,
       device = 'tiff',
       width = desired_width,
       height = desired_height,
       units = 'in', dpi = TIFF_dpi)




# Visualize co-occurrence of invasion types across all samples

######################
# EXTENDED DATA FIG 1H
######################

efig1H <- discovery_patient_metadata %>% 
  dplyr::select(c('VI','STAS','VPI','LI'))

# VI prevalence result reported in text

table(efig1H$VI)

prevalence <- efig1H %>%
  summarise(
    total = n(),  
    VI_count = sum(VI == 1),  
    VI_percent = VI_count / total * 100 
  )

print(prevalence)


efig1H <- as.data.frame(lapply(efig1H, as.factor))

efig1H <- lapply(efig1H, function(x) if (is.factor(x)) as.numeric(x) - 1 else x)

efig1H <- as.data.frame(lapply(efig1H, as.integer))

efig1H$`No Invasion` <- ifelse(efig1H$VI + efig1H$STAS + efig1H$VPI + efig1H$LI == 0,1,0)

efig1H$Sample <- discovery_patient_metadata$Sample

efig1H <- efig1H[, c("Sample", names(efig1H)[names(efig1H) != "Sample"])]

addWorksheet(source_data, "E. Figure 1H")

writeData(source_data, sheet = "E. Figure 1H", 
          x = "E. Figure 1H", 
          startCol =  1, startRow = 1)

writeData(source_data, sheet = 'E. Figure 1H',
          x = efig1H, startCol = 1, startRow = 3)

efig1H_gg <- upset(
  efig1H, c('VI','STAS','VPI','LI','No Invasion'), name = NULL, width_ratio=0.1,
  themes=ComplexUpset::upset_default_themes(text=element_text(family='Calibri',
                                                              size = 30)),
  sort_sets = FALSE,
  base_annotations = list(
    'Intersection size'=(
      intersection_size(text=list(size=7)) + 
        theme_classic() +
        theme(text=element_text(size=25),
              axis.text.x = element_blank()) +
        ylab('Invasion Co-occurence') +
        xlab(NULL)
    )),
  set_sizes=FALSE,
  matrix=(
    intersection_matrix(
      geom=geom_point(
        shape='circle',
        size=3.5
      )
    )
  )
) + patchwork::plot_layout(heights=c(4, 1),widths = c(1,3))


efig1H_gg

plot_dimensions <- dev.size("in")

aspect_ratio <- plot_dimensions[1] / plot_dimensions[2]

aspect_ratio <- 1.161765

desired_height <- 7 

desired_width <- desired_height * aspect_ratio

TIFF_dpi <- 600

ggsave(here('figures','efig1H.tiff'),
       plot = efig1H_gg,
       device = 'tiff',
       width = desired_width,
       height = desired_height,
       units = 'in', dpi = TIFF_dpi)


# VI multivariate coxph model

######################
# EXTENDED DATA FIG 1I
######################

efig1I <- discovery_patient_metadata %>%
  mutate(STAS = factor(STAS)) %>%
  mutate(LI = factor(LI)) %>%
  mutate(VPI = factor(VPI)) %>%
  mutate(Gender = factor(Gender)) %>%
  mutate(X8thStg = factor(X8thStg)) %>%
  filter(WHO != 'M') %>%
  mutate(WHO = factor(WHO, levels = c('AISMIA','1','2','3'))) %>%
  mutate(WHO = fct_collapse(WHO,
                            "AIS/MIA/1/2" = c("AISMIA","1","2"),
                            "3" = c("3"))) %>%
  mutate(X8thStg = fct_collapse(X8thStg,
                                "0/IA" = c("0","IA1","IA2","IA3"),
                                "IB" = c("IB"))) %>%
  rename("Surgical Procedure" = Surgical.Procedure) %>%
  mutate(`Surgical Procedure` = factor(`Surgical Procedure`, levels = c('Wedge','Segment','Lobe'))) %>%
  mutate(`Surgical Procedure` = fct_collapse(`Surgical Procedure`,
                                             'Sublobar' = c('Wedge','Segment'),
                                             'Lobectomy' = c('Lobe'))) %>%
  mutate(Race = factor(Race, levels = c('Asian','Black / African American',
                                        'Hispanic / Latino','Unknown','White'))) %>%
  mutate(Race = fct_collapse(Race,
                             'White' = 'White',
                             'Non-white' = c('Asian','Black / African American',
                                             'Hispanic / Latino'),
                             'Unknown' = c('Unknown'))) %>%
  filter(Race != 'Unknown') %>%
  mutate(Race = fct_drop(Race)) %>%
  select(Sample, time_7_year_RFS, RFS_7_year, VI, Gender, Site, Age, Pack_yrs,
         X8thStg, WHO, Invasive_size, Total_size, `Surgical Procedure`,
         LI, STAS, VPI, Race)

addWorksheet(source_data, "E. Figure 1I")

writeData(source_data, sheet = "E. Figure 1I", 
          x = "E. Figure 1I", 
          startCol =  1, startRow = 1)

writeData(source_data, sheet = 'E. Figure 1I',
          x = efig1I, startCol = 1, startRow = 3)


multivar_cox.res <- coxph(formula = Surv(time_7_year_RFS, RFS_7_year) ~ VI + 
                            Gender + Race + Site + Age + Pack_yrs + WHO + 
                            Invasive_size + Total_size + `Surgical Procedure` + 
                            LI + STAS + VPI + Race + X8thStg,
                          data=efig1I)

summary(multivar_cox.res)


# Extract HR and confidence interval

efig1I <- data.frame(
  HR = summary(multivar_cox.res)$coefficients[, 2], 
  lower_ci = summary(multivar_cox.res)$conf.int[, 3], 
  upper_ci = summary(multivar_cox.res)$conf.int[, 4], 
  P = as.character(summary(multivar_cox.res)$coefficients[, 5])
) %>% 
  rownames_to_column("Predictor") %>% 
  mutate(Predictor = case_when(
    Predictor == "GenderM" ~ "Male",
    Predictor == "RaceWhite" ~ "White",
    Predictor == "SiteLHMC" ~ "LHMC",
    Predictor == "Pack_yrs" ~ "Pack years",
    Predictor == "X8thStgIB" ~ "Stage IB",
    Predictor == "WHO3" ~ "WHO G3",
    Predictor == "Invasive_size" ~ "Invasive size",
    Predictor == "Total_size" ~ "Total size",
    Predictor == "`Surgical Procedure`Lobectomy" ~ "Lobectomy",
    Predictor == "LI1" ~ "LI",
    Predictor == "STAS1" ~ "STAS",
    Predictor == "VPI1" ~ "VPI",
    TRUE ~ Predictor
  )) %>% 
  bind_rows(
    tibble(Predictor = c("Female", "Non-white", "BMC","Stage 0/IA", "WHO G1/2", "Sublobar"),
           HR = 1, lower_ci = NA, upper_ci = NA, P = "Reference")
  ) %>% 
  mutate(Predictor = factor(Predictor, levels = c('VI', 'Female', 'Male', 
                                                  'Non-white', 'White', 
                                                  'BMC', 'LHMC', 'Age', 
                                                  'Pack years', 'Stage 0/IA', 
                                                  'Stage IB', 'WHO G1/2',
                                                  'WHO G3', 'Invasive size', 
                                                  'Total size', 'Sublobar', 
                                                  'Lobectomy', 'LI', 'STAS', 
                                                  'VPI'))) %>% 
  arrange(Predictor) %>% 
  mutate(across(c(HR, lower_ci, upper_ci), as.numeric)) %>% 
  mutate(across(c(HR, lower_ci, upper_ci), ~signif(., digits = 3))) %>%
  mutate(across(-1, ~ signif(as.numeric(.), digits = 3))) %>%
  mutate(P = as.character(signif(as.numeric(P), digits = 3))) %>%
  mutate(P = sub("e", " %*% 10^", P)) %>%
  mutate(P = case_when(
    Predictor %in% c("Female", "Non-white", "BMC", "Stage 0/IA", "WHO G1/2", "Sublobar") ~ "Reference",
    TRUE ~ P
  ))


# Generate forest plot range

hr_max <- ceiling(max(efig1I$upper_ci, na.rm = TRUE) / 10) * 10

efig1I_gg <- ggplot(efig1I, aes(x=Predictor, y=HR)) + coord_flip() +
  geom_hline(yintercept = 1, linetype="dotted") + 
  geom_errorbar(aes(x=Predictor, ymin=lower_ci, ymax=upper_ci), 
                width=0.5, color="#FF0054") +
  geom_point(pch=15, size=3) +
  theme_minimal() + 
  theme(legend.position = "none", 
        plot.title = element_text(size=14, face="bold"),
        aspect.ratio = 0.8,
        panel.grid = element_blank(),
        axis.line.x = element_line(),
        axis.ticks.x = element_line(),
        axis.text.x=element_text(size=15),
        axis.text.y=element_text(size=15),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        text = element_text(family = 'Calibri')) + 
  scale_y_continuous(trans = "log", expand = c(0,0), breaks = c(0.05,2,8,32), 
                     limits = c(0.05, 1200)) + 
  scale_x_discrete(limits=c(rev(as.character(efig1I$Predictor)),' ')) + 
  xlab("Predictor\n") + ylab("\nHazard Ratio (95% CI)") + 
  geom_text(aes(x=Predictor, y=hr_max + 300, label=P), 
            size=5, family = 'Calibri', parse = TRUE) + 
  annotate("text", x = 21, y = 330, label = "p value", 
           size = 5, fontface = 'bold', family = 'Calibri')

efig1I_gg

plot_dimensions <- dev.size("in")

aspect_ratio <- plot_dimensions[1] / plot_dimensions[2]

aspect_ratio <- 0.8632597

desired_height <- 8  

desired_width <- desired_height * aspect_ratio

TIFF_dpi <- 600

ggsave(here('figures','efig1I.tiff'),
       plot = efig1I_gg,
       device = 'tiff',
       width = desired_width,
       height = desired_height,
       units = 'in', dpi = TIFF_dpi)




# association of invasion types with histologic patterns in discovery cohort

######################
# EXTENDED DATA FIG 1J
######################

efig1J <- discovery_patient_metadata %>% 
  dplyr::select(c('Sample','VI','VPI','LI','STAS',
                  'Solid','Micropapillary',
                  'Cribriform','Acinar',
                  'Lepidic','Papillary'))

addWorksheet(source_data, "E. Figure 1J")

writeData(source_data, sheet = "E. Figure 1J", 
          x = "E. Figure 1J", 
          startCol =  1, startRow = 1)

writeData(source_data, sheet = 'E. Figure 1J',
          x = efig1J, startCol = 1, startRow = 3)

LUAD_invasion_types <- discovery_patient_metadata %>% 
  dplyr::select(c('VI','VPI','LI','STAS')) %>%
  mutate(across(everything(), as.numeric)) %>%
  mutate(VI = ifelse(VI == 1, 2, 1)) %>%
  mutate(VPI = ifelse(VPI == 1, 2, 1)) %>%
  mutate(LI = ifelse(LI == 1, 2, 1)) %>%
  mutate(STAS = ifelse(STAS == 1, 2, 1))


LUAD_patterns <- discovery_patient_metadata %>% 
  dplyr::select(c('Solid','Micropapillary',
                  'Cribriform','Acinar',
                  'Lepidic','Papillary')) %>%
  mutate(across(everything(), as.numeric))

# calculate associations

invasion_p_values <- function(LUAD_pattern, LUAD_invasion_type){
  
  tmpy <- data.frame('pattern' = LUAD_pattern, 'invasion' = LUAD_invasion_type)
  
  p_value <- wilcox.test(x = tmpy[which(tmpy$invasion == 1),][,1], 
                         y = tmpy[which(tmpy$invasion == 2),][,1])$p.value
  
  return(p_value)
  
}

# get signed z statistic

invasion_zvalues <- function(LUAD_pattern, LUAD_invasion_type){
  
  tmpy <- data.frame('pattern' = LUAD_pattern, 'invasion' = LUAD_invasion_type)
  
  z_value <- wilcoxonZ(x = tmpy[which(tmpy$invasion == 1),][,1], 
                       y = tmpy[which(tmpy$invasion == 2),][,1], exact = FALSE)*-1
  
  return(z_value)
  
}

# convert to long format 

p_values <- data.frame(sapply(as.list(LUAD_invasion_types), 
                              function(x) sapply(LUAD_patterns, invasion_p_values, 
                                                 LUAD_invasion_type = x)))

p_values <- p_values %>% rownames_to_column('pattern')

p_values <- tidyr::gather(p_values, invasion, p_value, VI:STAS, factor_key=TRUE)

z_values <- data.frame(sapply(as.list(LUAD_invasion_types), 
                              function(x) sapply(LUAD_patterns, invasion_zvalues, 
                                                 LUAD_invasion_type = x)))

z_values <- tidyr::gather(z_values, invasion, z_value, VI:STAS, factor_key=TRUE)

data_long <- cbind(p_values, z_values)[,-2]

data_long$pattern <- factor(data_long$pattern, levels = c('Lepidic','Papillary',
                                                          'Acinar','Cribriform',
                                                          'Micropapillary','Solid'))

# convert to -log10 p values

data_long$FDR <- p.adjust(data_long$p_value, n = length(data_long$p_value))

data_long$FDR <- -log10(data_long$FDR)

# dot plot heatmap

data_long$pattern <- paste('%',data_long$pattern,sep=' ')

data_long$pattern <- factor(data_long$pattern, 
                            levels = c('% Lepidic','% Papillary','% Acinar',
                                       '% Cribriform','% Micropapillary','% Solid'))

efig1J <- data_long %>% select(-p_value)

efig1J_gg <- ggplot(data = efig1J, aes(x = invasion, y = pattern)) +
  geom_point(aes(fill = z_value, size = FDR), shape = 21) +
  scale_size(breaks=c(1,2,5),labels=c(1,2,5), range = c(0,13)) +
  scale_fill_gradient2(low = "#007FFF", mid = "white", high = "red", 
                       limits = c(range(efig1J$z_value))) +
  theme_classic() +
  theme(
    legend.position = 'right',
    text = element_text(family = 'Calibri'),
    panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
    plot.title = element_text(size=22, hjust = 0.5),
    axis.text=element_text(size=18),
    axis.ticks=element_blank(),
    legend.title = element_text(size=20),
    legend.key.size = unit(0.6, "cm"),
    legend.text = element_text(size=16)) +
  labs(title = 'Association of invasion types\n with LUAD growth patterns\n',
       x = element_blank(), y = NULL, fill = 'Z Statistic') +
  scale_y_discrete(limits = rev) + 
  scale_x_discrete(position = "top") + 
  guides(size=guide_legend(title="-log₁₀(p.adj)"))


efig1J_gg

plot_dimensions <- dev.size("in")

aspect_ratio <- plot_dimensions[1] / plot_dimensions[2]

aspect_ratio <- 0.8622995

desired_height <- 6 

desired_width <- desired_height * aspect_ratio

TIFF_dpi <- 600

ggsave(here('figures','efig1J.tiff'),
       plot = efig1J_gg,
       device = 'tiff',
       width = desired_width,
       height = desired_height,
       units = 'in', dpi = TIFF_dpi)


# site specific recurrence by invasion status

######################
# EXTENDED DATA FIG 1K
######################

efig1K <- discovery_patient_metadata %>%
  select(Sample, VI, Recurrence, RFS) %>%
  mutate(VI = factor(VI)) %>%
  mutate(Recurrence = addNA(Recurrence))

levels(efig1K$Recurrence) <- c('distant (+/- loco-regional)',
                               'loco-regional only','unknown','none')

efig1K$Recurrence <- as.character(efig1K$Recurrence)

efig1K <- efig1K %>% 
  mutate(Recurrence = ifelse(Recurrence == 'none' & RFS == 1, 
                             'unknown', Recurrence)) %>%
  mutate(Recurrence = ifelse(Recurrence == 'none' & RFS == 0, 
                             'none', Recurrence)) %>%
  filter(Recurrence != 'unknown')

efig1K$Recurrence <- as.factor(efig1K$Recurrence)

efig1K$Recurrence <- droplevels(efig1K$Recurrence)

levels(efig1K$VI) <- c('VI-','VI+')

addWorksheet(source_data, "E. Figure 1K")

writeData(source_data, sheet = "E. Figure 1K", 
          x = "E. Figure 1K", 
          startCol =  1, startRow = 1)

writeData(source_data, sheet = 'E. Figure 1K',
          x = efig1K, startCol = 1, startRow = 3)

efig1K_tbl <- table(efig1K$VI, efig1K$Recurrence) / rowSums(table(efig1K$VI, efig1K$Recurrence))

df_counts <- as.data.frame(table(efig1K$VI, efig1K$Recurrence))

# Convert proportions to a data frame

efig1K_tbl <- as.data.frame(efig1K_tbl)

efig1K_tbl$count <- df_counts$Freq

colnames(efig1K_tbl) <- c('VI status','Site of recurrence','Freq','count')

efig1K_gg <- ggplot(efig1K_tbl, aes(fill=`Site of recurrence`, x=`VI status`, y = Freq)) + 
  theme_classic() + 
  coord_flip() +
  scale_fill_manual(values = c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF")) +
  geom_bar(position="fill", stat="identity", colour='black') +
  theme(plot.title = element_text(size=18, hjust = 0.5),
        axis.title.x = element_text(size=22),
        axis.title.y = element_text(size=22),
        axis.text.x = element_text(size=22),
        axis.text.y = element_text(size=22),
        legend.position="bottom",
        legend.title = element_text(size=22),
        legend.key.size = unit(0.6, "cm"),
        legend.text = element_text(size=20),
        text = element_text(family = 'Calibri')) +
  labs(x = '\nVI Status',
       y = 'Proportion\n') +
  ylim(0,1) + 
  guides(fill = guide_legend(nrow = 2, byrow = TRUE, title.position = 'top')) +
  geom_text(aes(label = count, group = `Site of recurrence`),
            position = position_stack(vjust = 0.5), 
            color = "white", size = 8, family = 'Calibri')

efig1K_gg

plot_dimensions <- dev.size("in")

aspect_ratio <- plot_dimensions[1] / plot_dimensions[2]

aspect_ratio <- 1.430195

desired_height <- 5

desired_width <- desired_height * aspect_ratio

TIFF_dpi <- 600

ggsave(here('figures','efig1K.tiff'),
       plot = efig1K_gg,
       device = 'tiff',
       width = desired_width,
       height = desired_height,
       units = 'in', dpi = TIFF_dpi)


# competing risks regression

######################
# EXTENDED DATA FIG 1L
######################

# include non-lung cancer deaths as competing events (OS 1 but not RFS 1)

efig1L <- discovery_patient_metadata %>%
  mutate(VI = factor(VI)) %>%
  mutate(Recurrence = addNA(Recurrence))

levels(efig1L$Recurrence) <- c('distant (+/- loco-regional)',
                               'loco-regional only','unknown','none')

efig1L$Recurrence <- as.character(efig1L$Recurrence)

efig1L <- efig1L %>% 
  mutate(Recurrence = ifelse(Recurrence == 'none' & RFS == 1, 
                             'unknown', Recurrence)) %>%
  mutate(Recurrence = ifelse(Recurrence == 'none' & RFS == 0, 
                             'none', Recurrence)) %>%
  filter(Recurrence != 'unknown')

efig1L$Recurrence <- as.factor(efig1L$Recurrence)

efig1L$Recurrence <- droplevels(efig1L$Recurrence)

levels(efig1L$VI) <- c('VI-','VI+')

efig1L$Recurrence <- as.character(efig1L$Recurrence)

efig1L$Recurrence[which(efig1L$OS_7_year == 1 & efig1L$RFS_7_year == 0 & efig1L$Recurrence == 'none')] <- 'non-cancer death'

efig1L <- efig1L %>% 
  rename("Surgical Procedure" = Surgical.Procedure) %>%
  dplyr::select('Sample','time_7_year_RFS','Recurrence','VI','Surgical Procedure',
                'Age','Pack_yrs','Gender','X8thStg','Site') %>%
  filter(Recurrence != 'unknown') %>%
  mutate(Recurrence = as.factor(Recurrence)) %>%
  filter(!is.na(X8thStg)) %>%
  filter(X8thStg != 'IIA') %>%
  filter(!is.na(time_7_year_RFS)) %>%
  filter(!is.na(`Surgical Procedure`)) %>%
  mutate(`Surgical Procedure` = as.factor(`Surgical Procedure`)) %>%
  mutate(`Surgical Procedure` = ifelse(`Surgical Procedure` == 'Lobe','Lobectomy','Sublobar')) %>%
  filter(!is.na(Age)) %>%
  filter(!is.na(Pack_yrs)) %>%
  filter(!is.na(Gender)) %>%
  mutate(Gender = as.factor(Gender))

addWorksheet(source_data, "E. Figure 1L")

writeData(source_data, sheet = "E. Figure 1L", 
          x = "E. Figure 1L", 
          startCol =  1, startRow = 1)

writeData(source_data, sheet = 'E. Figure 1L',
          x = efig1L, startCol = 1, startRow = 3)


recurrence_data <- lapply(efig1L, function(x) if(is.factor(x)) droplevels(x) else x)

# create model matrix

covs <- model.matrix(~ VI + Gender + Age + Pack_yrs + `Surgical Procedure` + Site,
                     data = recurrence_data)[, -1]

ftime <- recurrence_data$time_7_year_RFS

recurrence_data$Recurrence <- droplevels(recurrence_data$Recurrence)

# recode factor levels, 3 is distant, 2 loco, 1 non-cancer death, 0 none

levels(recurrence_data$Recurrence) <- c(3,2,1,0)

fstatus <- recurrence_data$Recurrence

# Fit the cumulative incidence model of Fine and Gray
# return the subdistribution hazard 
# (effect of VI and covariates on recurrence that is distant, with recurrence in other sites as competing event)

# subdistribution hazard ratio of VI for distant recurrence

fit_crr_distant <- cmprsk::crr(ftime=ftime, fstatus=fstatus, cov1 = covs, failcode = 3, 
                       cencode = 0, maxiter=20)

summary(fit_crr_distant)

# subdistribution hazard ratio of VI for locoregional recurrence

fit_crr_loco <- cmprsk::crr(ftime=ftime, fstatus=fstatus, cov1 = covs, failcode = 2, 
                    cencode = 0, maxiter=20)

summary(fit_crr_loco)

# plot sub distribution hazard ratios

efig1L <- list('Distant (+/- loco-regional)' = fit_crr_distant, 
               'Loco-regional only' = fit_crr_loco)

efig1L_gg <- modelplot(efig1L, coef_omit = 'Interc', exponentiate = TRUE,
                       coef_map = c('VIVI+' = 'VI'), 
                       size = 1.5, linewidth = 1.5,
                       conf_level = 0.95) + 
  theme_classic() +
  coord_flip() +
  ggtitle("VI\n") +
  scale_color_manual(values = c("#E64B35FF","#4DBBD5FF")) +
  aes(shape = ifelse(p.value < 0.01, "p < 0.01", "p > 0.01")) +
  scale_shape_manual(values = c(16,1)) +
  labs(x = '\nSub-distribution HR (95% CI)\n',
       y = NULL,
       color = 'Site of recurrence',
       shape = 'P value') +
  theme(plot.title = element_text(size=22, hjust = 0.5),
        aspect.ratio = 1,
        axis.title.x = element_text(),
        axis.title.y = element_text(size=20),
        axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        legend.title = element_text(size=20),
        legend.key.size = unit(0.6, "cm"),
        legend.text = element_text(size=20),
        legend.position = 'bottom',
        legend.box = "horizontal",
        text = element_text(family = 'Calibri')) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE, title.position = 'top'),
         shape = guide_legend(nrow = 2, byrow = TRUE, title.position = 'top')) +
  geom_vline(xintercept = 1, color = 'grey', linetype = 'dashed') +
  scale_y_discrete(breaks=c("0.5","1","2"))


efig1L_gg

plot_dimensions <- dev.size("in")

aspect_ratio <- plot_dimensions[1] / plot_dimensions[2]

aspect_ratio <- 1.868881

desired_height <- 6 

desired_width <- desired_height * aspect_ratio

TIFF_dpi <- 600

ggsave(here('figures','efig1L.tiff'),
       plot = efig1L_gg,
       device = 'tiff',
       width = desired_width,
       height = desired_height,
       units = 'in', dpi = TIFF_dpi)


# site of recurrence is not related to time to recurrence

######################
# EXTENDED DATA FIG 1M
######################

efig1M <- discovery_patient_metadata %>% 
  dplyr::select('time_7_year_RFS','Recurrence','VI','Surgical.Procedure',
                'Age','Pack_yrs','Gender','X8thStg','Site') %>%
  mutate(Recurrence = addNA(Recurrence)) %>%
  filter(!is.na(X8thStg)) %>%
  filter(X8thStg != 'IIA') %>%
  filter(!is.na(time_7_year_RFS)) %>%
  filter(!is.na(Surgical.Procedure)) %>%
  mutate(Surgical.Procedure = as.factor(Surgical.Procedure)) %>%
  mutate(Surgical.Procedure = ifelse(Surgical.Procedure == 'Lobe',
                                     'Lobectomy','Sublobar')) %>%
  filter(!is.na(Age)) %>%
  filter(!is.na(Pack_yrs)) %>%
  filter(!is.na(Gender)) %>%
  mutate(Gender = as.factor(Gender))

levels(efig1M$Recurrence) <- c('distant (+/- loco-regional)',
                               'loco-regional only','unknown','none')


efig1M <- efig1M %>% filter(Recurrence != 'unknown')

efig1M$Recurrence <- droplevels(efig1M$Recurrence)


efig1M$Recurrence <- factor(efig1M$Recurrence,
                            levels = c('none','loco-regional only',
                                       'distant (+/- loco-regional)'))

ftime <- efig1M$time_7_year_RFS


table(efig1M$Recurrence)

fstatus <- efig1M$Recurrence

addWorksheet(source_data, "E. Figure 1M")

writeData(source_data, sheet = "E. Figure 1M", 
          x = "E. Figure 1M", 
          startCol =  1, startRow = 1)

writeData(source_data, sheet = 'E. Figure 1M',
          x = efig1M, startCol = 1, startRow = 3)


# Create the cumulative incidence object with cencode = 0
ci_obj <- tidycmprsk::cuminc(
  formula = Surv(time_7_year_RFS, Recurrence, type = "mstate") ~ 1,
  data = efig1M
)

names(ci_obj)

efig1M_gg <- ci_obj %>%
  ggcuminc(outcome = c("loco-regional only", "distant (+/- loco-regional)")) +
  add_confidence_interval(alpha=0.1, type = 'ribbon') +
  scale_ggsurvfit() +
  labs(title = NULL,
       x = "\nTime (Months)",
       y = "Cumulative Incidence\n") +
  theme_classic(base_family = "Calibri", base_size = 20) +
  theme(plot.title    = element_text(size = 22, hjust = 0.5),
        axis.title.x  = element_text(size = 20),
        axis.title.y  = element_text(size = 20),
        axis.text.x   = element_text(size = 20),
        axis.text.y   = element_text(size = 20),
        legend.position = "bottom",
        legend.direction = "vertical",
        legend.title    = element_text(size = 22),
        legend.key.size = unit(0.6, "cm"),
        legend.text     = element_text(size = 20))

efig1M_gg <- ggplot_build(efig1M_gg)

efig1M_gg$data[[1]] <- efig1M_gg$data[[1]] %>%
  mutate(colour = ifelse(group == 1, "#E64B35FF", "#4DBBD5FF"))
efig1M_gg$data[[2]] <- efig1M_gg$data[[2]] %>%
  mutate(colour = ifelse(group == 1, "#E64B35FF", "#4DBBD5FF"),
         fill = ifelse(group == 1, "#E64B35FF", "#4DBBD5FF"))

efig1M_gg$data[[1]]$linewidth <- 1
efig1M_gg$data[[2]]$linewidth <- 0

efig1M_gg <- ggplot_gtable(efig1M_gg)

plot(efig1M_gg)


plot_dimensions <- dev.size("in")

aspect_ratio <- plot_dimensions[1] / plot_dimensions[2]

aspect_ratio <- 0.9545455

desired_height <- 6 

desired_width <- desired_height * aspect_ratio

TIFF_dpi <- 600

ggsave(here('figures','efig1M.tiff'),
       plot = efig1M_gg,
       device = 'tiff',
       width = desired_width,
       height = desired_height,
       units = 'in', dpi = TIFF_dpi)

###################
# LOAD RNA-SEQ DATA
###################

# read processed count matrix from GEO
# first download GSE273377_exp_count.txt.gz	file from subseries GSE273378 within superseries GSE273528
# unzip and place in Steiner_NatComm/data

all_stageI_LUAD_se <- read.table(here('data','GSE273377_exp_count.txt'),
                                 header = TRUE, row.names = 1,
                                 check.names = FALSE)


# read GEO metadata - placeholder until GEO series is public

# gse <- getGEO("GSE273528", GSEMatrix = TRUE)

geo_bulkRNAseq_metadata <- read_excel(here('data','Steiner_NatComm_bulkRNAseq_GEO_submission.xlsx'), 
                                      sheet = "Metadata", range = c('A33:S202'))

geo_bulkRNAseq_metadata_discovery <- geo_bulkRNAseq_metadata %>% 
  filter(cohort == 'discovery cohort')

geo_bulkRNAseq_metadata_validation <- geo_bulkRNAseq_metadata %>% 
  filter(cohort == 'validation cohort')


# create discovery cohort DGE list object

# put samples in same order

discovery_se <- all_stageI_LUAD_se[,which(colnames(all_stageI_LUAD_se) %in% geo_bulkRNAseq_metadata_discovery$`*library name`)]

discovery_se <- discovery_se[,geo_bulkRNAseq_metadata_discovery$`*library name`]

discovery_dge <- DGEList(counts = discovery_se, 
                         samples = geo_bulkRNAseq_metadata_discovery) 

colnames(discovery_dge$samples)[4:length(colnames(discovery_dge$samples))] <- colnames(geo_bulkRNAseq_metadata_discovery)

# create validation cohort DGE list object

validation_se <- all_stageI_LUAD_se[,which(colnames(all_stageI_LUAD_se) %in% geo_bulkRNAseq_metadata_validation$`*library name`)]

validation_se <- validation_se[,geo_bulkRNAseq_metadata_validation$`*library name`]

validation_dge <- DGEList(counts = validation_se, 
                          samples = geo_bulkRNAseq_metadata_validation) 

colnames(validation_dge$samples)[4:length(colnames(validation_dge$samples))] <- colnames(geo_bulkRNAseq_metadata_validation)

###############
# NORMALIZATION
###############

# Trimmed mean of m values (TMM) normalization

discovery_dge <- calcNormFactors(discovery_dge) 

head(discovery_dge$samples$lib.size)

validation_dge <- calcNormFactors(validation_dge) 

head(validation_dge$samples$lib.size)

###############
# QC EVALUATION
###############

# remove outlier samples (> 2 sd from mean TIN)

criteria <- mean(discovery_dge$samples$`mean TIN`) - sd(discovery_dge$samples$`mean TIN`)*2

remove <- discovery_dge$samples[which(discovery_dge$samples$`mean TIN` < criteria),]$`*library name`

discovery_dge$samples <- discovery_dge$samples[!(rownames(discovery_dge$samples) %in% remove), ]  

discovery_dge$counts <- discovery_dge$counts[,!(colnames(discovery_dge$counts) %in% remove)]

dim(discovery_dge)


criteria <- mean(validation_dge$samples$`mean TIN`) - sd(validation_dge$samples$`mean TIN`)*2

remove <- validation_dge$samples[which(validation_dge$samples$`mean TIN` < criteria),]$`*library name`

validation_dge$samples <- validation_dge$samples[!(rownames(validation_dge$samples) %in% remove), ]  

validation_dge$counts <- validation_dge$counts[,!(colnames(validation_dge$counts) %in% remove)]

dim(validation_dge)

# save validation cohort dge object for comparison with matched Visium data (includes 1 stage II)

# saveRDS(validation_dge, file = here('data','validation_dge.rds'))


# keep only stage I (8th TNM edition) bulk RNA-seq validation cohort
# note that 13034 has matching stRNA-seq data and will be used in that analysis

remove <- validation_dge$samples[which(validation_dge$samples$`TNM 8th edition stage` == 'IIA'),]$`*library name`

validation_dge$samples <- validation_dge$samples[!(rownames(validation_dge$samples) %in% remove), ]  

validation_dge$counts <- validation_dge$counts[,!(colnames(validation_dge$counts) %in% remove)]

dim(validation_dge)



# combine geo metadata with additional supplementary clinical data

# for now remove any duplicated columns from geo metadata

common_columns <- intersect(names(discovery_dge$samples), names(discovery_patient_metadata))

discovery_dge$samples <- discovery_dge$samples %>% select(-all_of(common_columns))

discovery_dge$samples <- merge(discovery_dge$samples, 
                               discovery_patient_metadata, 
                               by.x = "*library name",
                               by.y = "Sample")

rownames(discovery_dge$samples) <- discovery_dge$samples$`*library name`

discovery_dge$samples <- discovery_dge$samples[colnames(discovery_dge$counts),]


common_columns <- intersect(names(validation_dge$samples), names(validation_patient_metadata))

validation_dge$samples <- validation_dge$samples %>% select(-all_of(common_columns))

# need to update validation_patient_metadata to have same ids as geo metadata

validation_dge$samples <- merge(validation_dge$samples, 
                                validation_patient_metadata, 
                                by.x = "*library name",
                                by.y = "BU_ID")

rownames(validation_dge$samples) <- validation_dge$samples$`*library name`

validation_dge$samples <- validation_dge$samples[colnames(validation_dge$counts),]

##########
# TABLE S1
##########

# discovery cohort (summarized by LMPVI grade)

discovery_cohort_tableone <- CreateTableOne(data = discovery_dge$samples, 
                                            vars = c("X8thStg","Age","Gender","Race",
                                                     "Pack_yrs","Smoker",
                                                     "Total_size","Invasive_size",
                                                     "Lepidic","Acinar",
                                                     "Papillary", "Solid",
                                                     "Micropapillary","Cribriform",
                                                     "LI","VPI","STAS","Site",
                                                     "Surgical.Procedure"),
                                            strata = 'LMPVI')

tableS1 <- print(discovery_cohort_tableone)

write.csv(tableS1, file = here('data','tableS1.csv'))


##########
# TABLE S3
##########

# validation cohort (summarized by LMPVI grade)

validation_cohort_tableone <- CreateTableOne(data = validation_dge$samples, 
                                             vars = c("X8thStg","Age","Gender",
                                                      "Pack_yrs","Smoker",
                                                      "Total_size","Invasive_size",
                                                      "Lepidic","Acinar","Papillary",
                                                      "Solid","Micropapillary","Cribriform",
                                                      "LI","VPI","STAS","Site",
                                                      "Surgical.Procedure"),
                                             strata = 'LMPVI')

summary(validation_cohort_tableone)

tableS3 <- print(validation_cohort_tableone)

write.csv(tableS3, file = here('data','tableS3.csv'))




############################
# GENE FILTERING (DISCOVERY)
############################

gene_annotations <- readRDS(here('data','gene_annotations.rds'))

# convert ensemble gene names to hgnc

discovery_dge <- ensemble_to_hgnc(discovery_dge)

rownames(discovery_dge$samples) <- discovery_dge$samples$`*library name`

colnames(discovery_dge$counts) <- discovery_dge$samples$`*library name`


cpm.tmp <- edgeR::cpm(discovery_dge$counts)

# at least 10% of samples should express a particular gene

keep <- rowSums(cpm.tmp>1) > 0.1*ncol(cpm.tmp) 

counts.tmp <- discovery_dge$counts[keep,]

dim(counts.tmp)

discovery_dge <- DGEList(counts = counts.tmp, samples = discovery_dge$samples) 

head(discovery_dge$samples$lib.size)

# Trimmed mean of m values (TMM) normalization

discovery_dge <- calcNormFactors(discovery_dge)

######################
# ASSESS BATCH EFFECTS
######################

# PCA dimensionality reduction

pca_batch <- prcomp(t(edgeR::cpm(discovery_dge,log=T)),scale=T, center = T)

# Visualize possible batch effects due to collection site in discovery cohort

autoplot(pca_batch, data = discovery_dge$samples, 
         colour = 'Site', label = TRUE) 

# perform batch correction with combat-seq

discovery_dge$counts <- ComBat_seq(discovery_dge$counts, discovery_dge$samples$Site)

# visualize batch-corrected samples

pca_batch <- prcomp(t(edgeR::cpm(discovery_dge,log=T)),scale=T, center = T)

autoplot(pca_batch, data = discovery_dge$samples, 
         colour = 'Site', label = TRUE) 


# save discovery cohort dge object

# saveRDS(discovery_dge, here('data','discovery_dge.rds'))


###################################################
# DIFFERENTIALLY EXPRESSED GENES ASSOCIATED WITH VI
###################################################

design <- model.matrix(~0+LMPVI, discovery_dge$samples)

discovery_dge <- estimateDisp(discovery_dge, design)

colnames(design) = c('LMP','NST','VI') # rename column names of design matrix

VI_contrasts = makeContrasts(VI-LMP,levels = design)

fit <- glmFit(discovery_dge, design) 

VI_lrt <- glmLRT(fit, contrast = VI_contrasts)

VI_toptags <- topTags(VI_lrt, n = dim(discovery_dge)[1])

print(summary(decideTests(VI_lrt, adjust.method = "fdr", p.value = 0.01)))

VI_genes_up <- rownames(VI_toptags$table[which((VI_toptags$table$FDR)<0.01 & (VI_toptags$table$logFC)>0),])

VI_genes_dn <- rownames(VI_toptags$table[which((VI_toptags$table$FDR)<0.01 & (VI_toptags$table$logFC)<0),])

VI_genes <- c(VI_genes_up, VI_genes_dn)


VI_toptags$table$rank <- -log10(VI_toptags$table$PValue) * sign(VI_toptags$table$logFC)

VI_ranked_list <- sort(setNames(VI_toptags$table$rank, rownames(VI_toptags$table)), decreasing = TRUE)

# look at genes in heatmap

discovery_dge$samples$VI <- as.factor(discovery_dge$samples$VI)

discovery_dge$samples$NST <- as.factor(discovery_dge$samples$NST)

discovery_dge$samples$LMP <- as.factor(discovery_dge$samples$LMP)

discovery_dge$samples$LI <- as.factor(discovery_dge$samples$LI)

levels(discovery_dge$samples$VI) <- c('-','+')

levels(discovery_dge$samples$NST) <- c('-','+')

levels(discovery_dge$samples$LMP) <- c('-','+')

levels(discovery_dge$samples$LI) <- c('-','+')

annotation_col = subset(discovery_dge$samples, select=c("Site","LI","LMPVI"))

names(annotation_col)[3] <- 'Novel grade'

annotation_col <- annotation_col %>%
  mutate_all(as.factor)

rownames(annotation_col) <- discovery_dge$samples$SAMPLE_ID

colnames(discovery_dge) <- rownames(annotation_col)

annotation_colors = list(
  `Novel grade` = c( "LMP" = "#1f78b4",
                     "NST" = "#33a02c",
                     "VI" = "#FF0054"),
  Cluster = c("1" = "#cab2d6","2" = "#b2df8a",
              "3" = "#E31A1C","4" = "#A6CEE3"),
  Site = c("LHMC" = "#ff7f00","BMC" = "#E5E4E2"),
  VPI = c("+" = "#6a3d9a","-" = "#E5E4E2"),
  LI = c("+" = "#1A5276","-" = "#E5E4E2"))

# plot heatmap at FDR < 0.01

VIheatmap <- pheatmap::pheatmap((edgeR::cpm(discovery_dge,log=T)[VI_genes,]),
                                color = bluered, border_color = NA, scale = "row", show_rownames = FALSE,
                                show_colnames = FALSE, cluster_rows = T, cluster_cols = T,
                                clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean",
                                legend = TRUE, annotation_legend = TRUE,
                                annotation_col = annotation_col,
                                annotation_colors = annotation_colors,
                                clustering_method = "ward.D2", main = "VI vs. LMP, FDR < 0.01")


# Perform consensus clustering to identify best k

efig2A <- edgeR::cpm(discovery_dge, log = TRUE)[VI_genes,]

addWorksheet(source_data, "E. Figure 2A")

writeData(source_data, sheet = "E. Figure 2A", 
          x = "E. Figure 2A", 
          startCol =  1, startRow = 1)

writeData(source_data, sheet = 'E. Figure 2A',
          x = efig2A, startCol = 1, startRow = 3)

results <- ConsensusClusterPlus(efig2A,
                                maxK = 10,          
                                reps = 10,         
                                pItem = 0.8,       
                                pFeature = 1,    
                                clusterAlg = "hc",
                                distance = "euclidean",
                                seed = 222)

# navigate 1 plot pane to the left before saving

dev.list()
dev.cur()
dev.set(2)  

aspect_ratio <- plot_dimensions[1] / plot_dimensions[2]
aspect_ratio <- 0.9627193
desired_height <- 4  
desired_width <- desired_height * aspect_ratio
TIFF_dpi <- 600

dev.copy(
  tiff,
  filename = here("figures","efig2A.tiff"),
  width = desired_width,    
  height = desired_height,  
  units = "in",
  res = TIFF_dpi          
)

dev.off()

dev.off()


# Derive gene clusters from heatmap:

VI_gene_cluster <- as.character(cutree(as.hclust(VIheatmap$tree_row), k = 4)) # derive gene clusters from row clustering

VI_gene_cluster <- data.frame(Cluster = VI_gene_cluster, row.names = NULL)

rownames(VI_gene_cluster) <- make.names(VIheatmap$tree_row$labels, unique = TRUE)

# Replot heatmap with gene clusters annotated:


###########
# FIGURE 1B
###########

fig1B <- edgeR::cpm(discovery_dge,log=T)[VI_genes,]

addWorksheet(source_data, "Figure 1B")

writeData(source_data, sheet = "Figure 1B", 
          x = "Figure 1B", 
          startCol =  1, startRow = 1)

writeData(source_data, sheet = 'Figure 1B',
          x = fig1B, startCol = 1, startRow = 3,
          rowNames = TRUE)


fig1B_gg <- pheatmap::pheatmap(fig1B,
                               color = bluered, border_color = NA, scale = "row", show_rownames = FALSE,
                               show_colnames = FALSE, cluster_rows = T, cluster_cols = T,
                               clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean",
                               legend = TRUE, annotation_legend = TRUE,
                               annotation_row = VI_gene_cluster,
                               annotation_col = annotation_col,
                               annotation_colors = annotation_colors,
                               treeheight_row = 40,
                               angle_col = 0,
                               treeheight_col = 0,
                               fontsize = 16,
                               fontfamily = 'Calibri',
                               clustering_method = "ward.D2", main = "")

dev.off()

fig1B_gg$gtable$grobs[[5]]$gp = gpar(fontface = 'plain')

fig1B_gg$gtable$grobs[[7]]$gp = gpar(fontsize = 16, fontface = 'plain')

fig1B_gg



plot_dimensions <- dev.size("in")

aspect_ratio <- plot_dimensions[1] / plot_dimensions[2]

aspect_ratio <- 0.8481013

desired_height <- 8.5  

desired_width <- desired_height * aspect_ratio

TIFF_dpi <- 600

ggsave(here('figures','fig1B.tiff'),
       plot = fig1B_gg,
       device = 'tiff',
       width = desired_width,
       height = desired_height,
       units = 'in', dpi = TIFF_dpi)


# save VI gene clusters

VI_gene_cluster_1 <- VI_gene_cluster %>% filter(Cluster == 1)

VI_gene_cluster_2 <- VI_gene_cluster %>% filter(Cluster == 2)

VI_gene_cluster_3 <- VI_gene_cluster %>% filter(Cluster == 3)

VI_gene_cluster_4 <- VI_gene_cluster %>% filter(Cluster == 4)

addWorksheet(source_data, "Figure 1B-2")

writeData(source_data, sheet = "Figure 1B-2", 
          x = "Figure 1B-2", 
          startCol =  1, startRow = 1)

writeData(source_data, sheet = 'Figure 1B-2',
          x = VI_gene_cluster, startCol = 1, startRow = 3,
          rowNames = TRUE)

# saveRDS(VI_gene_cluster, file = here('data','VI_gene_clusters.rds'))


levels(discovery_dge$samples$VI) <- c(0,1)

levels(discovery_dge$samples$NST) <- c(0,1)

levels(discovery_dge$samples$LMP) <- c(0,1)

levels(discovery_dge$samples$LI) <- c(0,1)


# examine individual cluster predictions of VI vs no VI

discovery_dge$samples$cluster1_scores <- colMeans(zscore.rows(cpm(discovery_dge$counts,log=T)[rownames(VI_gene_cluster_1),]))

discovery_dge$samples$cluster2_scores <- colMeans(zscore.rows(cpm(discovery_dge$counts,log=T)[rownames(VI_gene_cluster_2),]))

discovery_dge$samples$cluster3_scores <- colMeans(zscore.rows(cpm(discovery_dge$counts,log=T)[rownames(VI_gene_cluster_3),]))

discovery_dge$samples$cluster4_scores <- colMeans(zscore.rows(cpm(discovery_dge$counts,log=T)[rownames(VI_gene_cluster_4),]))

# cluster 1

efig2B <- discovery_dge$samples %>% 
  select(VI, cluster1_scores, cluster2_scores, cluster3_scores, cluster4_scores)

addWorksheet(source_data, "E. Figure 2B")

writeData(source_data, sheet = "E. Figure 2B", 
          x = "E. Figure 2B", 
          startCol =  1, startRow = 1)

writeData(source_data, sheet = 'E. Figure 2B',
          x = efig2B, startCol = 1, startRow = 3,
          rowNames = TRUE)

gene_roc <- roc(as.factor(efig2B$VI), efig2B$cluster1_scores,
                levels=c(0, 1), ci = TRUE)

gene_roc

ciobj <- ci.se(gene_roc, specificities=seq(0, 1, l=25))

dat.ci <- data.frame(x = as.numeric(rownames(ciobj)),
                     lower = ciobj[, 1],
                     upper = ciobj[, 3])

rocobj <- pROC::plot.roc(gene_roc, print.thres = TRUE, print.auc = TRUE)

p_val <- round(wilcox.test(efig2B$cluster1_scores ~ efig2B$VI)$p.value,7)

auc <- round(pROC::auc(gene_roc),2)

######################
# EXTENDED DATA FIG 2B
######################


pval <- grobTree(textGrob(paste0('Cluster 1 Z-score\n ',
                                 '(AUROC = ', auc, ')',sep=''),
                          x=0.75,  y=0.25, 
                          gp=gpar(fontsize=17, 
                                  fontfamily = 'Calibri')),
                 textGrob("p = ", 
                          x = 0.65, y = 0.15, 
                          gp = gpar(fontsize = 17, fontfamily = 'Calibri')),
                 textGrob(eval(bquote(scientific_10(p_val))),
                          x=0.8, y=0.15,
                          gp=gpar(fontsize=17, 
                                  fontfamily = 'Calibri')))

efig2B_gg <- pROC::ggroc(rocobj, colour = '#cab2d6', size = 1, alpha = 1) +
  theme_classic() +
  annotation_custom(pval) +
  labs(title='Discovery cohort (n=103)\n VI vs. no VI\n', 
       x = "\nSpecificity", y = "Sensitivity\n") +
  theme(plot.title = element_text(size=22, hjust = 0.5),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        legend.title = element_text(size=20),
        text = element_text(family = 'Calibri')) +
  geom_abline(slope=1, intercept = 1, linetype = "dashed", 
              alpha=0.7, color = "grey") + 
  geom_ribbon(data = dat.ci, 
              aes(x = x, ymin = lower, ymax = upper), 
              fill = '#cab2d6', alpha= 0.1) 


efig2B_gg

plot_dimensions <- dev.size("in")

aspect_ratio <- plot_dimensions[1] / plot_dimensions[2]

aspect_ratio <- 0.9122807

desired_height <- 6.5  

desired_width <- desired_height * aspect_ratio

TIFF_dpi <- 600

ggsave(here('figures','efig2B.tiff'),
       plot = efig2B_gg,
       device = 'tiff',
       width = desired_width,
       height = desired_height,
       units = 'in', dpi = TIFF_dpi)





# cluster 2

gene_roc <- roc(as.factor(efig2B$VI), efig2B$cluster2_scores,
                levels=c(0, 1), ci = TRUE)

gene_roc

ciobj <- ci.se(gene_roc, specificities=seq(0, 1, l=25))

dat.ci <- data.frame(x = as.numeric(rownames(ciobj)),
                     lower = ciobj[, 1],
                     upper = ciobj[, 3])

rocobj <- pROC::plot.roc(gene_roc, print.thres = TRUE, print.auc = TRUE)

p_val <- round(wilcox.test(efig2B$cluster2_scores ~ efig2B$VI)$p.value,7)

auc <- round(pROC::auc(gene_roc),2)

pval <- grobTree(textGrob(paste0('Cluster 2 Z-score\n ',
                                 '(AUROC = ', auc, ')',sep=''),
                          x=0.75,  y=0.25, 
                          gp=gpar(fontsize=17, 
                                  fontfamily = 'Calibri')),
                 textGrob("p = ", 
                          x = 0.65, y = 0.15, 
                          gp = gpar(fontsize = 17, fontfamily = 'Calibri')),
                 textGrob(eval(bquote(scientific_10(p_val))),
                          x=0.8, y=0.15,
                          gp=gpar(fontsize=17, 
                                  fontfamily = 'Calibri')))

efig2B_2_gg <- pROC::ggroc(rocobj, colour = '#b2df8a', size = 1, alpha = 1) +
  theme_classic() +
  annotation_custom(pval) +
  labs(title='Discovery cohort (n=103)\n VI vs. no VI\n', 
       x = "\nSpecificity", y = "Sensitivity\n") +
  theme(plot.title = element_text(size=22, hjust = 0.5),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        legend.title = element_text(size=20),
        text = element_text(family = 'Calibri')) +
  geom_abline(slope=1, intercept = 1, linetype = "dashed", 
              alpha=0.7, color = "grey") + 
  geom_ribbon(data = dat.ci, 
              aes(x = x, ymin = lower, ymax = upper), 
              fill = '#b2df8a', alpha= 0.1) 


plot_dimensions <- dev.size("in")

aspect_ratio <- plot_dimensions[1] / plot_dimensions[2]

aspect_ratio <- 0.9122807

desired_height <- 6.5

desired_width <- desired_height * aspect_ratio

TIFF_dpi <- 600

ggsave(here('figures','efig2B_2.tiff'),
       plot = efig2B_2_gg,
       device = 'tiff',
       width = desired_width,
       height = desired_height,
       units = 'in', dpi = TIFF_dpi)





# cluster 3

gene_roc <- roc(as.factor(efig2B$VI), efig2B$cluster3_scores,
                levels=c(0, 1), ci = TRUE)

gene_roc

ciobj <- ci.se(gene_roc, specificities=seq(0, 1, l=25))

dat.ci <- data.frame(x = as.numeric(rownames(ciobj)),
                     lower = ciobj[, 1],
                     upper = ciobj[, 3])

rocobj <- pROC::plot.roc(gene_roc, print.thres = TRUE, print.auc = TRUE)

p_val <- round(wilcox.test(efig2B$cluster3_scores ~ efig2B$VI)$p.value,9)

wilcox.test(efig2B$cluster3_scores ~ efig2B$VI)

auc <- round(pROC::auc(gene_roc),2)

pval <- grobTree(textGrob(paste0('Cluster 3 Z-score\n ',
                                 '(AUROC = ', auc, ')',sep=''),
                          x=0.75,  y=0.25, 
                          gp=gpar(fontsize=17, 
                                  fontfamily = 'Calibri')),
                 textGrob("p = ", 
                          x = 0.65, y = 0.15, 
                          gp = gpar(fontsize = 17, fontfamily = 'Calibri')),
                 textGrob(eval(bquote(scientific_10(p_val))),
                          x=0.8, y=0.15,
                          gp=gpar(fontsize=17, 
                                  fontfamily = 'Calibri')))

efig2B_3_gg <- pROC::ggroc(rocobj, colour = '#E31A1C', size = 1, alpha = 1) +
  theme_classic() +
  annotation_custom(pval) +
  labs(title='Discovery cohort (n=103)\n VI vs. no VI\n', 
       x = "\nSpecificity", y = "Sensitivity\n") +
  theme(plot.title = element_text(size=22, hjust = 0.5),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        legend.title = element_text(size=20),
        text = element_text(family = 'Calibri')) +
  geom_abline(slope=1, intercept = 1, linetype = "dashed", 
              alpha=0.7, color = "grey") + 
  geom_ribbon(data = dat.ci, 
              aes(x = x, ymin = lower, ymax = upper), 
              fill = '#E31A1C', alpha= 0.1) 


plot_dimensions <- dev.size("in")

aspect_ratio <- plot_dimensions[1] / plot_dimensions[2]

aspect_ratio <- 0.9122807

desired_height <- 6.5

desired_width <- desired_height * aspect_ratio

TIFF_dpi <- 600

ggsave(here('figures','efig2B_3.tiff'),
       plot = efig2B_3_gg,
       device = 'tiff',
       width = desired_width,
       height = desired_height,
       units = 'in', dpi = TIFF_dpi)




# cluster 4

gene_roc <- roc(as.factor(efig2B$VI), efig2B$cluster4_scores,
                levels=c(0, 1), ci = TRUE)

gene_roc

ciobj <- ci.se(gene_roc, specificities=seq(0, 1, l=25))

dat.ci <- data.frame(x = as.numeric(rownames(ciobj)),
                     lower = ciobj[, 1],
                     upper = ciobj[, 3])

rocobj <- pROC::plot.roc(gene_roc, print.thres = TRUE, print.auc = TRUE)

# H0 = 'AUC is equal to 0.5' is equivalent to H0 = 'the distribution of ranks in the two groups are equal'
# therefore can use the mann-whitney-wilcoxon test

p_val <- round(wilcox.test(efig2B$cluster4_scores ~ efig2B$VI)$p.value,6)

auc <- round(pROC::auc(gene_roc),2)

pval <- grobTree(textGrob(paste0('Cluster 4 Z-score\n ',
                                 '(AUROC = ', auc, ')',sep=''),
                          x=0.75,  y=0.25, 
                          gp=gpar(fontsize=17, 
                                  fontfamily = 'Calibri')),
                 textGrob("p = ", 
                          x = 0.65, y = 0.15, 
                          gp = gpar(fontsize = 17, fontfamily = 'Calibri')),
                 textGrob(eval(bquote(scientific_10(p_val))),
                          x=0.8, y=0.15,
                          gp=gpar(fontsize=17, 
                                  fontfamily = 'Calibri')))


efig2B_4_gg <- pROC::ggroc(rocobj, colour = '#A6CEE3', size = 1, alpha = 1) +
  theme_classic() +
  annotation_custom(pval) +
  labs(title='Discovery cohort (n=103)\n VI vs. no VI\n', 
       x = "\nSpecificity", y = "Sensitivity\n") +
  theme(plot.title = element_text(size=22, hjust = 0.5),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        legend.title = element_text(size=20),
        text = element_text(family = 'Calibri')) +
  geom_abline(slope=1, intercept = 1, linetype = "dashed", 
              alpha=0.7, color = "grey") + 
  geom_ribbon(data = dat.ci, 
              aes(x = x, ymin = lower, ymax = upper), 
              fill = '#A6CEE3', alpha= 0.1) 




plot_dimensions <- dev.size("in")

aspect_ratio <- plot_dimensions[1] / plot_dimensions[2]

aspect_ratio <- 0.9122807

desired_height <- 6.5  

desired_width <- desired_height * aspect_ratio

TIFF_dpi <- 600

ggsave(here('figures','efig2B_4.tiff'),
       plot = efig2B_4_gg,
       device = 'tiff',
       width = desired_width,
       height = desired_height,
       units = 'in', dpi = TIFF_dpi)






# examine individual cluster predictions of VI vs NST

efig2C <- discovery_dge$samples %>% 
  filter(LMPVI != 'LMP') %>%
  select(VI, cluster1_scores, cluster2_scores, cluster3_scores, cluster4_scores)

addWorksheet(source_data, "E. Figure 2C")

writeData(source_data, sheet = "E. Figure 2C", 
          x = "E. Figure 2C", 
          startCol =  1, startRow = 1)

writeData(source_data, sheet = 'E. Figure 2C',
          x = efig2C, startCol = 1, startRow = 3,
          rowNames = TRUE)



# cluster 1

gene_roc <- roc(as.factor(efig2C$VI), efig2C$cluster1_scores,
                levels=c(0, 1), ci = TRUE)

gene_roc

ciobj <- ci.se(gene_roc, specificities=seq(0, 1, l=25))

dat.ci <- data.frame(x = as.numeric(rownames(ciobj)),
                     lower = ciobj[, 1],
                     upper = ciobj[, 3])

rocobj <- pROC::plot.roc(gene_roc, print.thres = TRUE, print.auc = TRUE)

p_val <- round(wilcox.test(efig2C$cluster1_scores ~ efig2C$VI)$p.value,3)

auc <- round(pROC::auc(gene_roc),2)

######################
# EXTENDED DATA FIG 2C
######################

pval <- grobTree(textGrob(paste('Cluster 1 Z-score\n ',
                                '(AUROC = ', auc, ')','\np = ',p_val,sep=''),
                          x=0.75,  y=0.25, 
                          gp=gpar(fontsize=17, 
                                  fontfamily = 'Calibri')))

efig2C_1_gg <- pROC::ggroc(rocobj, colour = '#cab2d6', size = 1, alpha = 1) +
  theme_classic() +
  annotation_custom(pval) +
  labs(title='Discovery cohort (n=80)\n VI vs. NST only\n', 
       x = "\nSpecificity", y = "Sensitivity\n") +
  theme(plot.title = element_text(size=22, hjust = 0.5),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        legend.title = element_text(size=20),
        text = element_text(family = 'Calibri')) +
  geom_abline(slope=1, intercept = 1, linetype = "dashed", 
              alpha=0.7, color = "grey") + 
  geom_ribbon(data = dat.ci, 
              aes(x = x, ymin = lower, ymax = upper), 
              fill = '#cab2d6', alpha= 0.1) 



plot_dimensions <- dev.size("in")

aspect_ratio <- plot_dimensions[1] / plot_dimensions[2]

aspect_ratio <- 0.9122807

desired_height <- 6.5

desired_width <- desired_height * aspect_ratio

TIFF_dpi <- 600

ggsave(here('figures','efig2C_1.tiff'),
       plot = efig2C_1_gg,
       device = 'tiff',
       width = desired_width,
       height = desired_height,
       units = 'in', dpi = TIFF_dpi)


# cluster 2

gene_roc <- roc(as.factor(efig2C$VI), efig2C$cluster2_scores,
                levels=c(0, 1), ci = TRUE)

gene_roc

ciobj <- ci.se(gene_roc, specificities=seq(0, 1, l=25))

dat.ci <- data.frame(x = as.numeric(rownames(ciobj)),
                     lower = ciobj[, 1],
                     upper = ciobj[, 3])

rocobj <- pROC::plot.roc(gene_roc, print.thres = TRUE, print.auc = TRUE)

p_val <- round(wilcox.test(efig2C$cluster2_scores ~ efig2C$VI)$p.value,3)

wilcox.test(efig2C$cluster2_scores ~ efig2C$VI)

auc <- round(pROC::auc(gene_roc),2)

pval <- grobTree(textGrob(paste('Cluster 2 Z-score\n ',
                                '(AUROC = ', auc, ')','\np = ',p_val,sep=''),
                          x=0.75,  y=0.25, 
                          gp=gpar(fontsize=17, 
                                  fontfamily = 'Calibri')))

efig2C_2_gg <- pROC::ggroc(rocobj, colour = '#b2df8a', size = 1, alpha = 1) +
  theme_classic() +
  annotation_custom(pval) +
  labs(title='Discovery cohort (n=80)\n VI vs. NST only\n', 
       x = "\nSpecificity", y = "Sensitivity\n") +
  theme(plot.title = element_text(size=22, hjust = 0.5),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        legend.title = element_text(size=20),
        text = element_text(family = 'Calibri')) +
  geom_abline(slope=1, intercept = 1, linetype = "dashed", 
              alpha=0.7, color = "grey") + 
  geom_ribbon(data = dat.ci, 
              aes(x = x, ymin = lower, ymax = upper), 
              fill = '#b2df8a', alpha= 0.1) 



plot_dimensions <- dev.size("in")

aspect_ratio <- plot_dimensions[1] / plot_dimensions[2]

aspect_ratio <- 0.9122807

desired_height <- 6.5 

desired_width <- desired_height * aspect_ratio

TIFF_dpi <- 600

ggsave(here('figures','efig2C_2.tiff'),
       plot = efig2C_2_gg,
       device = 'tiff',
       width = desired_width,
       height = desired_height,
       units = 'in', dpi = TIFF_dpi)



# cluster 3

gene_roc <- roc(as.factor(efig2C$VI), efig2C$cluster3_scores,
                levels=c(0, 1), ci = TRUE)

gene_roc

ciobj <- ci.se(gene_roc, specificities=seq(0, 1, l=25))

dat.ci <- data.frame(x = as.numeric(rownames(ciobj)),
                     lower = ciobj[, 1],
                     upper = ciobj[, 3])

rocobj <- pROC::plot.roc(gene_roc, print.thres = TRUE, print.auc = TRUE)

p_val <- round(wilcox.test(efig2C$cluster3_scores ~ efig2C$VI)$p.value,7)

auc <- round(pROC::auc(gene_roc),2)

pval <- grobTree(textGrob(paste0('Cluster 3 Z-score\n ',
                                 '(AUROC = ', auc, ')',sep=''),
                          x=0.75,  y=0.25, 
                          gp=gpar(fontsize=17, 
                                  fontfamily = 'Calibri')),
                 textGrob("p = ", 
                          x = 0.65, y = 0.15, 
                          gp = gpar(fontsize = 17, fontfamily = 'Calibri')),
                 textGrob(eval(bquote(scientific_10(p_val))),
                          x=0.8, y=0.15,
                          gp=gpar(fontsize=17, 
                                  fontfamily = 'Calibri')))

efig2C_3_gg <- pROC::ggroc(rocobj, colour = '#E31A1C', size = 1, alpha = 1) +
  theme_classic() +
  annotation_custom(pval) +
  labs(title='Discovery cohort (n=80)\n VI vs. NST only\n', 
       x = "\nSpecificity", y = "Sensitivity\n") +
  theme(plot.title = element_text(size=22, hjust = 0.5),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        legend.title = element_text(size=20),
        text = element_text(family = 'Calibri')) +
  geom_abline(slope=1, intercept = 1, linetype = "dashed", 
              alpha=0.7, color = "grey") + 
  geom_ribbon(data = dat.ci, 
              aes(x = x, ymin = lower, ymax = upper), 
              fill = '#E31A1C', alpha= 0.1) 



plot_dimensions <- dev.size("in")

aspect_ratio <- plot_dimensions[1] / plot_dimensions[2]

aspect_ratio <- 0.9122807

desired_height <- 6.5

desired_width <- desired_height * aspect_ratio

TIFF_dpi <- 600

ggsave(here('figures','efig2C_3.tiff'),
       plot = efig2C_3_gg,
       device = 'tiff',
       width = desired_width,
       height = desired_height,
       units = 'in', dpi = TIFF_dpi)



# cluster 4

gene_roc <- roc(as.factor(efig2C$VI), efig2C$cluster4_scores,
                levels=c(0, 1), ci = TRUE)

gene_roc

ciobj <- ci.se(gene_roc, specificities=seq(0, 1, l=25))

dat.ci <- data.frame(x = as.numeric(rownames(ciobj)),
                     lower = ciobj[, 1],
                     upper = ciobj[, 3])

rocobj <- pROC::plot.roc(gene_roc, print.thres = TRUE, print.auc = TRUE)

p_val <- round(wilcox.test(efig2C$cluster4_scores ~ efig2C$VI)$p.value,3)

auc <- round(pROC::auc(gene_roc),2)


pval <- grobTree(textGrob(paste('Cluster 4 Z-score\n ',
                                '(AUROC = ', auc, ')','\np = ',p_val,sep=''),
                          x=0.75,  y=0.25, 
                          gp=gpar(fontsize=17, 
                                  fontfamily = 'Calibri')))

efig2C_4_gg <- pROC::ggroc(rocobj, colour = '#A6CEE3', size = 1, alpha = 1) +
  theme_classic() +
  annotation_custom(pval) +
  labs(title='Discovery cohort (n=80)\n VI vs. NST only\n', 
       x = "\nSpecificity", y = "Sensitivity\n") +
  theme(plot.title = element_text(size=22, hjust = 0.5),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        legend.title = element_text(size=20),
        text = element_text(family = 'Calibri')) +
  geom_abline(slope=1, intercept = 1, linetype = "dashed", 
              alpha=0.7, color = "grey") + 
  geom_ribbon(data = dat.ci, 
              aes(x = x, ymin = lower, ymax = upper), 
              fill = '#A6CEE3', alpha= 0.1) 



plot_dimensions <- dev.size("in")

aspect_ratio <- plot_dimensions[1] / plot_dimensions[2]

aspect_ratio <- 0.9122807

desired_height <- 6.5

desired_width <- desired_height * aspect_ratio

TIFF_dpi <- 600

ggsave(here('figures','efig2C_4.tiff'),
       plot = efig2C_4_gg,
       device = 'tiff',
       width = desired_width,
       height = desired_height,
       units = 'in', dpi = TIFF_dpi)



#######################################
# BIOLOGICAL ENRICHMENT OF VI SIGNATURE
#######################################

# Identify enriched biological processes in VI signature clusters

websiteLive <- getOption("enrichR.live")
if (websiteLive) {
  listEnrichrSites()
  setEnrichrSite("Enrichr")
}

if (websiteLive) dbs <- listEnrichrDbs()

if (websiteLive) head(dbs)



# Cluster 1 enrichment

#############
# FIGURE 1C_1
#############

signature_enriched_1 <- enrichr(c(rownames(VI_gene_cluster_1)), 
                                c("MSigDB_Hallmark_2020",
                                  "GO_Biological_Process_2023"))

signature_enriched_1 <- bind_rows(signature_enriched_1, .id = "column_label")

signature_enriched_1 <- signature_enriched_1[order(-signature_enriched_1$Combined.Score),]

signature_enriched_1 <- subset(signature_enriched_1, 
                               select=c("column_label","Term","Combined.Score",
                                        "Adjusted.P.value","Odds.Ratio",
                                        "Overlap"))

signature_enriched_1$`VI Cluster` <- '1'

# Plot significantly enriched biological processes

signature_enriched_1 <- signature_enriched_1[order(-signature_enriched_1$Combined.Score),]

# first select all pathways that pass FDR < 0.05

signature_enriched_1 <- signature_enriched_1[which(signature_enriched_1$Adjusted.P.value < 0.05),]

# then select top 10 pathways ranked by combined score

signature_enriched_1_top10 <- signature_enriched_1[1:10,]

signature_enriched_1_top10$Adjusted.P.value <- -log(signature_enriched_1_top10$Adjusted.P.value, base = 10)

signature_enriched_1_top10$Term <- sapply(strsplit(signature_enriched_1_top10$Term, 
                                                   split='(', fixed=TRUE), function(x) (x[1]))

signature_enriched_1_top10$column_label <- as.factor(signature_enriched_1_top10$column_label)

levels(signature_enriched_1_top10$column_label) <- c('BP','H')

fig1C_1 <- signature_enriched_1_top10

# shorten some names for plotting purposes

fig1C_1$Term[3] <- 'Attachment of Spindle Microtubules to Kinetochore'

fig1C_1$Term[10] <- 'Negative Regulation of Mitotic Metaphase/Anaphase'


fig1C_1_gg <- ggplot(fig1C_1, aes(x=reorder(Term,Adjusted.P.value),
                                  y=Adjusted.P.value, group = 1, 
                                  fill=`VI Cluster`)) + 
  geom_bar(stat = "identity", width=0.8, colour = 'black', show.legend = FALSE) +
  scale_fill_manual(values = VI_cluster_colors, 
                    guide = guide_legend(title = NULL, label = FALSE)) +
  geom_hline(yintercept = 1.3,linetype = "dotted",linewidth=1.5) +
  theme_minimal() +
  geom_text(aes(0,1.5,label = "FDR 0.05", vjust = -0.5, hjust= -0.02),
            family = "Calibri", size = 6) +
  labs(title=NULL, x = NULL, y = "\n-log₁₀(FDR)") +
  theme(plot.title = element_text(size=28, hjust = -1, vjust = 2),
        axis.title.x = element_text(size=24),
        axis.title.y = element_text(size=24, 
                                    colour = 'black'),
        axis.text.x.top = element_text(color = 'black'), 
        axis.title.x.bottom = element_text(color = 'black'),
        axis.text.x.bottom = element_text(color = "black"), 
        axis.text.x = element_text(colour="black", size=24),
        axis.text.y = element_text(colour="black", size=24),
        axis.ticks = element_line(),
        legend.position="none", 
        axis.text=element_text(size=12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        text = element_text(family = 'Calibri')) +
  scale_y_continuous(breaks =  seq(0,56, 10),labels =  seq(0,56, 10), 
                     limits = c(-2.5,56)) +
  coord_flip() +
  geom_hline(yintercept=0) +
  geom_segment(aes(x=0, xend=0,y=0,yend=56), color = 'black') +
  annotate(
    "text",
    x = 1:10, 
    y = -2.5,
    label = rev(fig1C_1[order(-fig1C_1$Adjusted.P.value),]$column_label),
    fontface = 'plain',
    family = 'Calibri',
    size = 7.5) 

fig1C_1_gg

plot_dimensions <- dev.size("in")

aspect_ratio <- plot_dimensions[1] / plot_dimensions[2]

aspect_ratio <- 2.458101

desired_height <- 5 

desired_width <- desired_height * aspect_ratio

TIFF_dpi <- 600

ggsave(here('figures','fig1C_1.tiff'),
       plot = fig1C_1_gg,
       device = 'tiff',
       width = desired_width,
       height = desired_height,
       units = 'in', dpi = TIFF_dpi)



# Cluster 2 enrichment

#############
# FIGURE 1C_2
#############

signature_enriched_2 <- enrichr(c(rownames(VI_gene_cluster_2)), 
                                c("MSigDB_Hallmark_2020",
                                  "GO_Biological_Process_2023")) 

signature_enriched_2 <- bind_rows(signature_enriched_2, .id = "column_label")

signature_enriched_2 <- signature_enriched_2[order(-signature_enriched_2$Combined.Score),]

signature_enriched_2 <- subset(signature_enriched_2, 
                               select=c("column_label","Term","Combined.Score",
                                        "Adjusted.P.value","Odds.Ratio",
                                        "Overlap"))

signature_enriched_2$`VI Cluster` <- '2'

# Plot significantly enriched biological processes

signature_enriched_2 <- signature_enriched_2[order(-signature_enriched_2$Combined.Score),]

signature_enriched_2 <- signature_enriched_2[which(signature_enriched_2$Adjusted.P.value < 0.11),]

signature_enriched_2_top10 <- signature_enriched_2[1:10,]

signature_enriched_2_top10$Adjusted.P.value <- -log(signature_enriched_2_top10$Adjusted.P.value, base = 10)

signature_enriched_2_top10$Term <- sapply(strsplit(signature_enriched_2_top10$Term, 
                                                   split='(', fixed=TRUE), function(x) (x[1]))


signature_enriched_2_top10$column_label <- as.factor(signature_enriched_2_top10$column_label)

levels(signature_enriched_2_top10$column_label) <- c('BP','H')

fig1C_2 <- signature_enriched_2_top10

fig1C_2$Term[10] <- 'Negative Regulation of IGF-R Signaling Pathway'

fig1C_2_gg <- ggplot(fig1C_2, aes(x=reorder(Term,Adjusted.P.value),
                                  y=Adjusted.P.value, group = 1, 
                                  fill=`VI Cluster`)) + 
  geom_bar(stat = "identity", width=0.8, colour = 'black', show.legend = FALSE) +
  scale_fill_manual(values = VI_cluster_colors, 
                    guide = guide_legend(title = NULL, label = FALSE)) +
  geom_hline(yintercept = 1.3,linetype = "dotted",linewidth=1.5) +
  theme_minimal() +
  geom_text(aes(0,1.5,label = "FDR 0.05", vjust = -0.5, hjust= -0.02),
            family = "Calibri", size = 6) +
  labs(title=NULL, x = NULL, y = "\n-log₁₀(FDR)") +
  theme(plot.title = element_text(size=28, hjust = -1, vjust = 2),
        axis.title.x = element_text(size=24),
        axis.title.y = element_text(size=24, 
                                    colour = 'black'),
        axis.text.x.top = element_text(color = 'black'), 
        axis.title.x.bottom = element_text(color = 'black'),
        axis.text.x.bottom = element_text(color = "black"), 
        axis.text.x = element_text(colour="black", size=24),
        axis.text.y = element_text(colour="black", size=24),
        axis.ticks = element_line(),
        legend.position="none", 
        axis.text=element_text(size=12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        text = element_text(family = 'Calibri')) +
  scale_y_continuous(breaks =  seq(0,32, 5),labels =  seq(0,32, 5), 
                     limits = c(-1.5,32)) +
  coord_flip() +
  geom_hline(yintercept=0) +
  geom_segment(aes(x=0, xend=0,y=0,yend=32), color = 'black') +
  annotate(
    "text",
    x = 1:10, 
    y = -1.5,
    label = rev(fig1C_2[order(-fig1C_2$Adjusted.P.value),]$column_label),
    fontface = 'plain',
    family = 'Calibri',
    size = 7.5) 

fig1C_2_gg

plot_dimensions <- dev.size("in")

aspect_ratio <- plot_dimensions[1] / plot_dimensions[2]

aspect_ratio <- 2.361266

desired_height <- 5

desired_width <- desired_height * aspect_ratio

TIFF_dpi <- 600

ggsave(here('figures','fig1C_2.tiff'),
       plot = fig1C_2_gg,
       device = 'tiff',
       width = desired_width,
       height = desired_height,
       units = 'in', dpi = TIFF_dpi)



# Cluster 3 enrichment

#############
# FIGURE 1C_3
#############

signature_enriched_3 <- enrichr(c(rownames(VI_gene_cluster_3)), 
                                c("MSigDB_Hallmark_2020",
                                  "GO_Biological_Process_2023")) 

signature_enriched_3 <- bind_rows(signature_enriched_3, .id = "column_label")

signature_enriched_3 <- signature_enriched_3[order(-signature_enriched_3$Combined.Score),]

signature_enriched_3 <- subset(signature_enriched_3, 
                               select=c("column_label","Term","Combined.Score",
                                        "Adjusted.P.value","Odds.Ratio",
                                        "Overlap"))

signature_enriched_3$`VI Cluster` <- '3'

# Plot significantly enriched biological processes

signature_enriched_3 <- signature_enriched_3[order(-signature_enriched_3$Combined.Score),]

signature_enriched_3 <- signature_enriched_3[which(signature_enriched_3$Adjusted.P.value < 0.1),]

signature_enriched_3_top10 <- signature_enriched_3[1:10,]

signature_enriched_3_top10$Adjusted.P.value <- -log(signature_enriched_3_top10$Adjusted.P.value, base = 10)

signature_enriched_3_top10$Term <- sapply(strsplit(signature_enriched_3_top10$Term, split='(', fixed=TRUE), function(x) (x[1]))

signature_enriched_3_top10$column_label <- as.factor(signature_enriched_3_top10$column_label)

levels(signature_enriched_3_top10$column_label) <- c('BP','H')


fig1C_3 <- signature_enriched_3_top10

fig1C_3_gg <- ggplot(fig1C_3, aes(x=reorder(Term,Adjusted.P.value),
                                  y=Adjusted.P.value, group = 1, 
                                  fill=`VI Cluster`)) + 
  geom_bar(stat = "identity", width=0.8, colour = 'black', show.legend = FALSE) +
  scale_fill_manual(values = VI_cluster_colors, 
                    guide = guide_legend(title = NULL, label = FALSE)) +
  geom_hline(yintercept = 1.3,linetype = "dotted",linewidth=1.5) +
  theme_minimal() +
  geom_text(aes(0,1.5,label = "FDR 0.05", vjust = -0.5, hjust= -0.02),
            family = "Calibri", size = 6) +
  labs(title=NULL, x = NULL, y = "\n-log₁₀(FDR)") +
  theme(plot.title = element_text(size=28, hjust = -1, vjust = 2),
        axis.title.x = element_text(size=24),
        axis.title.y = element_text(size=24, 
                                    colour = 'black'),
        axis.title.x.bottom = element_text(color = 'black'),
        axis.text.x.bottom = element_text(color = "black"), 
        axis.text.x = element_text(colour="black", size=24),
        axis.text.y = element_text(colour="black", size=24),
        axis.ticks = element_line(),
        legend.position="none", 
        axis.text=element_text(size=12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        text = element_text(family = 'Calibri')) +
  scale_y_continuous(breaks =  seq(0,10, 2),labels =  seq(0,10, 2), 
                     limits = c(-0.5,10)) +
  coord_flip() +
  geom_hline(yintercept=0) +
  geom_segment(aes(x=0, xend=0,y=0,yend=10), color = 'black') +
  annotate(
    "text",
    x = 1:10, 
    y = -0.5,
    label = rev(fig1C_3[order(-fig1C_3$Adjusted.P.value),]$column_label),
    fontface = 'plain',
    family = 'Calibri',
    size = 7.5) 

fig1C_3_gg

plot_dimensions <- dev.size("in")

aspect_ratio <- plot_dimensions[1] / plot_dimensions[2]

aspect_ratio <- 2.117318

desired_height <- 5

desired_width <- desired_height * aspect_ratio

TIFF_dpi <- 600

ggsave(here('figures','fig1C_3.tiff'),
       plot = fig1C_3_gg,
       device = 'tiff',
       width = desired_width,
       height = desired_height,
       units = 'in', dpi = TIFF_dpi)



# Cluster 4 enrichment

#############
# FIGURE 1C_4
#############

signature_enriched_4 <- enrichr(c(rownames(VI_gene_cluster_4)), 
                                c("MSigDB_Hallmark_2020",
                                  "GO_Biological_Process_2023")) 

signature_enriched_4 <- bind_rows(signature_enriched_4, .id = "column_label")

signature_enriched_4 <- signature_enriched_4[order(-signature_enriched_4$Combined.Score),]

signature_enriched_4 <- subset(signature_enriched_4, 
                               select=c("column_label","Term","Combined.Score",
                                        "Adjusted.P.value","Odds.Ratio",
                                        "Overlap"))

signature_enriched_4$`VI Cluster` <- '4'

# Plot significantly enriched biological processes

signature_enriched_4 <- signature_enriched_4[order(signature_enriched_4$Adjusted.P.value),]

signature_enriched_4_top10 <- signature_enriched_4[1:10,]

signature_enriched_4 <- signature_enriched_4[order(-signature_enriched_4$Combined.Score),]

signature_enriched_4_top10$Adjusted.P.value <- -log(signature_enriched_4_top10$Adjusted.P.value, base = 10)

signature_enriched_4_top10$Term <- sapply(strsplit(signature_enriched_4_top10$Term, split='(', fixed=TRUE), function(x) (x[1]))

signature_enriched_4_top10$column_label <- as.factor(signature_enriched_4_top10$column_label)

levels(signature_enriched_4_top10$column_label) <- c('BP','H')

fig1C_4 <- signature_enriched_4_top10

fig1C_4$Term[5] <- 'Pathway-Restricted SMAD Protein Phosphorylation'

fig1C_4$Term[2] <- 'Small GTPase Mediated Signal Transduction'


fig1C <- rbind(fig1C_1, fig1C_2, fig1C_3, fig1C_4)

addWorksheet(source_data, "Figure 1C")

writeData(source_data, sheet = "Figure 1C", 
          x = "Figure 1C", 
          startCol =  1, startRow = 1)

writeData(source_data, sheet = 'Figure 1C',
          x = fig1C, startCol = 1, startRow = 3)


fig1C_4_gg <- ggplot(fig1C_4, aes(x=reorder(Term,Adjusted.P.value),
                                  y=Adjusted.P.value, group = 1, 
                                  fill=`VI Cluster`)) + 
  geom_bar(stat = "identity", width=0.8, colour = 'black', show.legend = FALSE) +
  scale_fill_manual(values = VI_cluster_colors, 
                    guide = guide_legend(title = NULL, label = FALSE)) +
  geom_hline(yintercept = 1.3,linetype = "dotted",linewidth=1.5) +
  theme_minimal() +
  geom_text(aes(0,1.5,label = "FDR 0.05", vjust = -0.5, hjust= -0.02),
            family = "Calibri", size = 6) +
  labs(title=NULL, x = NULL, y = "\n-log₁₀(FDR)") +
  theme(plot.title = element_text(size=28, hjust = -1, vjust = 2),
        axis.title.x = element_text(size=24),
        axis.title.y = element_text(size=24, 
                                    colour ='black'),
        axis.text.x.top = element_text(color = 'black'), 
        axis.title.x.bottom = element_text(color = 'black'),
        axis.text.x.bottom = element_text(color = "black"), 
        axis.text.x = element_text(colour="black", size=24),
        axis.text.y = element_text(colour="black", size=24),
        axis.ticks = element_line(),
        legend.position="none", 
        axis.text=element_text(size=12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        text = element_text(family = 'Calibri')) +
  scale_y_continuous(breaks =  seq(0,2, 0.5),labels =  seq(0,2, 0.5), 
                     limits = c(-0.1,2)) +
  coord_flip() +
  geom_hline(yintercept=0) +
  geom_segment(aes(x=0, xend=0,y=0,yend=2), color = 'black') +
  annotate(
    "text",
    x = 1:10, 
    y = -0.1,
    label = rev(fig1C_4[order(-fig1C_4$Adjusted.P.value),]$column_label),
    fontface = 'plain',
    family = 'Calibri',
    size = 7.5) 

fig1C_4_gg

plot_dimensions <- dev.size("in")

# Calculate aspect ratio

aspect_ratio <- plot_dimensions[1] / plot_dimensions[2]

aspect_ratio <- 2.394786

# Define a reasonable height (in inches)

desired_height <- 5

# Calculate corresponding width to maintain aspect ratio

desired_width <- desired_height * aspect_ratio

TIFF_dpi <- 600

ggsave(here('figures','fig1C_4.tiff'),
       plot = fig1C_4_gg,
       device = 'tiff',
       width = desired_width,
       height = desired_height,
       units = 'in', dpi = TIFF_dpi)


# differentially expressed genes associated with LI

# first re-derive VI signature with LI as a covariate

######################
# EXTENDED DATA FIG 3A
######################

# Feature filtering on all data

design_withLI <- model.matrix(~1+LMPVI+LI, discovery_dge$samples)

discovery_dge_withLI <- estimateDisp(discovery_dge, design)

colnames(design_withLI) = c('intercept','NST','VI','LI') # rename column names of design matrix

VI_contrasts_withLI = makeContrasts(VI,levels = design_withLI)

fit_withLI <- glmFit(discovery_dge_withLI, design_withLI) 

VI_lrt_withLI <- glmLRT(fit_withLI, contrast = VI_contrasts_withLI)

VI_toptags_withLI <- topTags(VI_lrt_withLI, n = dim(discovery_dge_withLI)[1])

print(summary(decideTests(VI_lrt_withLI, adjust.method = "fdr", p.value = 0.01)))

VI_genes_up_withLI <- rownames(VI_toptags_withLI$table[which((VI_toptags_withLI$table$FDR)<0.01 & (VI_toptags_withLI$table$logFC)>0),])

VI_genes_dn_withLI <- rownames(VI_toptags_withLI$table[which((VI_toptags_withLI$table$FDR)<0.01 & (VI_toptags_withLI$table$logFC)<0),])

VI_genes_withLI <- c(VI_genes_up_withLI, VI_genes_dn_withLI)

efig3A <- edgeR::cpm(discovery_dge,log=T)[VI_genes_withLI,]

addWorksheet(source_data, "E. Figure 3A")

writeData(source_data, sheet = "E. Figure 3A", 
          x = "E. Figure 3A", 
          startCol =  1, startRow = 1)

writeData(source_data, sheet = 'E. Figure 3A',
          x = efig3A, startCol = 1, startRow = 3,
          rowNames = TRUE)

efig3A_gg <- pheatmap::pheatmap(efig3A,
                                color = bluered, border_color = NA, scale = "row", show_rownames = FALSE,
                                show_colnames = FALSE, cluster_rows = T, cluster_cols = T,
                                clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean",
                                legend = TRUE, annotation_legend = TRUE,
                                annotation_row = VI_gene_cluster,
                                annotation_col = annotation_col,
                                annotation_colors = annotation_colors,
                                treeheight_row = 40,
                                angle_col = 0,
                                treeheight_col = 0,
                                fontsize = 16,
                                fontfamily = 'Calibri',
                                clustering_method = "ward.D2", main = "")


dev.off()

efig3A_gg$gtable$grobs[[5]]$gp = gpar(fontface = 'plain')

efig3A_gg$gtable$grobs[[7]]$gp = gpar(fontsize = 16, fontface = 'plain')

efig3A_gg$gtable$grobs[[8]]$gp = gpar(fontface = 'plain')

efig3A_gg


plot_dimensions <- dev.size("in")

aspect_ratio <- plot_dimensions[1] / plot_dimensions[2]

aspect_ratio <- 1.022727

desired_height <- 8.5

desired_width <- desired_height * aspect_ratio

TIFF_dpi <- 600

ggsave(here('figures','efig3A.tiff'),
       plot = efig3A_gg,
       device = 'tiff',
       width = desired_width,
       height = desired_height,
       units = 'in', dpi = TIFF_dpi)





# genes associated with LI contrasts

######################
# EXTENDED DATA FIG 3B
######################

# Feature filtering on all data

design_withLI <- model.matrix(~1+LMPVI+LI, discovery_dge$samples)

discovery_dge_withLI <- estimateDisp(discovery_dge, design)

colnames(design_withLI) = c('intercept','NST','VI','LI') 

LI_contrasts = makeContrasts(LI,levels = design_withLI)

LI_fit <- glmFit(discovery_dge_withLI, design_withLI) 

LI_lrt <- glmLRT(LI_fit, contrast = LI_contrasts)

LI_toptags <- topTags(LI_lrt, n = dim(discovery_dge_withLI)[1])

print(summary(decideTests(LI_lrt, adjust.method = "fdr", p.value = 0.01)))

LI_genes_up <- rownames(LI_toptags$table[which((LI_toptags$table$FDR)<0.01 & (LI_toptags$table$logFC)>0),])

LI_genes_dn <- rownames(LI_toptags$table[which((LI_toptags$table$FDR)<0.01 & (LI_toptags$table$logFC)<0),])

LI_genes <- c(LI_genes_up, LI_genes_dn)


LI_toptags$table$rank <- -log10(LI_toptags$table$PValue) * sign(LI_toptags$table$logFC)

LI_ranked_list <- sort(setNames(LI_toptags$table$rank, rownames(LI_toptags$table)), decreasing = TRUE)

efig3B <- edgeR::cpm(discovery_dge_withLI,log=T)[LI_genes,]

addWorksheet(source_data, "E. Figure 3B")

writeData(source_data, sheet = "E. Figure 3B", 
          x = "E. Figure 3B", 
          startCol =  1, startRow = 1)

writeData(source_data, sheet = 'E. Figure 3B',
          x = efig3B, startCol = 1, startRow = 3,
          rowNames = TRUE)

efig3B_gg <- pheatmap::pheatmap(edgeR::cpm(discovery_dge_withLI,log=T)[LI_genes,],
                                color = bluered, border_color = NA, scale = "row", show_rownames = FALSE,
                                show_colnames = FALSE, cluster_rows = T, cluster_cols = T,
                                clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean",
                                legend = TRUE, annotation_legend = TRUE,
                                annotation_row = VI_gene_cluster,
                                annotation_col = annotation_col,
                                annotation_colors = annotation_colors,
                                treeheight_row = 40,
                                angle_col = 0,
                                treeheight_col = 0,
                                fontsize = 16,
                                fontfamily = 'Calibri',
                                clustering_method = "ward.D2", main = "")


dev.off()

efig3B_gg$gtable$grobs[[5]]$gp = gpar(fontface = 'plain')

efig3B_gg$gtable$grobs[[7]]$gp = gpar(fontsize = 16, fontface = 'plain')

efig3B_gg$gtable$grobs[[8]]$gp = gpar(fontface = 'plain')

efig3B_gg


plot_dimensions <- dev.size("in")

aspect_ratio <- plot_dimensions[1] / plot_dimensions[2]

aspect_ratio <- 1.022727

desired_height <- 8.5

desired_width <- desired_height * aspect_ratio

TIFF_dpi <- 600

ggsave(here('figures','efig3B.tiff'),
       plot = efig3B_gg,
       device = 'tiff',
       width = desired_width,
       height = desired_height,
       units = 'in', dpi = TIFF_dpi)



# look at VI cluster enrichment on LI ranked list


VI_biomarker_genes_df <- rbind(data.frame(term = 'VI_gene_cluster_1',
                                          gene = rownames(VI_gene_cluster_1)),
                               data.frame(term = 'VI_gene_cluster_2',
                                          gene = rownames(VI_gene_cluster_2)),
                               data.frame(term = 'VI_gene_cluster_3',
                                          gene = rownames(VI_gene_cluster_3)),
                               data.frame(term = 'VI_gene_cluster_4',
                                          gene = rownames(VI_gene_cluster_4)))

set.seed(111)

LI_gsea <- GSEA(LI_ranked_list, TERM2GENE = VI_biomarker_genes_df, 
                pvalueCutoff = 1, seed = TRUE)

# build dataframe for plot

LI_ranked_df <- data.frame('gene' = names(LI_ranked_list),
                           'rank' = unname(LI_ranked_list))

LI_ranked_df <- LI_ranked_df %>%
  mutate('VI gene cluster 1' = ifelse(gene %in% rownames(VI_gene_cluster_1), 'Yes','No')) %>%
  mutate('VI gene cluster 2' = ifelse(gene %in% rownames(VI_gene_cluster_2), 'Yes','No')) %>%
  mutate('VI gene cluster 3' = ifelse(gene %in% rownames(VI_gene_cluster_3), 'Yes','No')) %>%
  mutate('VI gene cluster 4' = ifelse(gene %in% rownames(VI_gene_cluster_4), 'Yes','No')) %>%
  arrange(-rank) %>%
  mutate('Rank' = row_number())

# correct p values

gsea_pvals <- data.frame('Cluster' = LI_gsea@result$ID,
                         'p_value' = LI_gsea@result$pvalue)

gsea_pvals$p_value <- p.adjust(gsea_pvals$p_value, method = 'bonferroni',
                               n=length(gsea_pvals$p_value))


# function to get running enrichment score for custom plot

gsInfo <- function(object, geneSetID) {
  
  geneList <- object@geneList
  
  if (is.numeric(geneSetID))
    
    geneSetID <- object@result[geneSetID, "ID"]
  
  geneSet <- object@geneSets[[geneSetID]]
  
  exponent <- object@params[["exponent"]]
  
  df <- gseaScores(geneList, geneSet, exponent, fortify=TRUE)
  
  df$ymin <- 0
  
  df$ymax <- 0
  
  pos <- df$position == 1
  
  h <- diff(range(df$runningScore))/20
  
  df$ymin[pos] <- -h
  
  df$ymax[pos] <- h
  
  df$geneList <- geneList
  
  df$Description <- object@result[geneSetID, "Description"]
  
  return(df)
  
}


gseaScores <- getFromNamespace("gseaScores", "DOSE")

# re-scale ranking metric for plotting purposes

LI_ranked_df$rank_scaled <- scale(LI_ranked_df$rank, 
                                  center = mean(LI_ranked_df$rank), 
                                  scale = max(LI_ranked_df$rank) - min(LI_ranked_df$rank))

LI_ranked_df$rank_scaled <- LI_ranked_df$rank_scaled[,1]


efig3C_F <- LI_ranked_df

addWorksheet(source_data, "E. Figure 3C-F")

writeData(source_data, sheet = "E. Figure 3C-F", 
          x = "E. Figure 3C-F", 
          startCol =  1, startRow = 1)

writeData(source_data, sheet = 'E. Figure 3C-F',
          x = efig3C_F, startCol = 1, startRow = 3)


# function for custom gsea plot

custom_gsea_plot <- function(gsea_cluster, pval_rnd) {
  
  VI_cluster_color <- unname(VI_cluster_colors[as.numeric(gsub("\\D", "", gsea_cluster))])
  
  LI_ranked_df$runningScore <- gsInfo(LI_gsea, geneSetID = c(gsea_cluster))$runningScore
  
  gsea_cluster_sub <- gsub('_',' ',gsea_cluster)
  
  gsea_cluster_sub <- ensym(gsea_cluster_sub)
  
  direction <- ifelse(abs(min(gsInfo(LI_gsea, geneSetID = c(gsea_cluster))$runningScore)) < abs(max(gsInfo(LI_gsea, geneSetID = c(gsea_cluster))$runningScore)),0.4,1.4)
  
  if (gsea_cluster_sub == 'VI gene cluster 2') {
    scientific_10 <- function(input_value){
      return(input_value)
    }
  }
  
  pval <- grobTree(textGrob("p =", 
                            x = 0.53, y = 0.95, hjust=0,
                            gp = gpar(fontsize = 17, fontfamily = 'Calibri')),
                   textGrob(eval(bquote(scientific_10(round(subset(gsea_pvals, 
                                                                   Cluster == gsea_cluster)$p_value,pval_rnd)))),
                            x=0.6, y=0.95, hjust=0,
                            gp=gpar(fontsize=17, 
                                    fontfamily = 'Calibri')))
  
  
  p <- ggplot(LI_ranked_df, aes(x=Rank, y=0.25, group=1)) + 
    geom_line(aes(y = rank_scaled-1), linetype = 'solid', linewidth = 1) +
    geom_line(aes(y = (runningScore+direction)), 
              color=VI_cluster_color, linewidth=1) +
    geom_bar(stat = "identity", alpha = 1, width=0.5, 
             aes(colour = !!gsea_cluster_sub, 
                 fill = !!gsea_cluster_sub), 
             show.legend = TRUE) +
    scale_fill_manual(values = c('No' = 'transparent',
                                 'Yes' = VI_cluster_color)) +
    scale_colour_manual(values = c('No' = 'transparent',
                                   'Yes' = VI_cluster_color)) +
    theme_classic() +
    labs(title=NULL, x = "\nGene rank", 
         y = "Association with LI           Running ES\n") +
    theme(plot.title = element_text(size=22, hjust = 0.5),
          axis.title.x = element_text(size=20),
          axis.title.y = element_text(size=20, hjust = 0.8),
          axis.text.x = element_text(colour="black", size=20),
          axis.text.y = element_text(colour="black", size=20),
          legend.position="right", 
          axis.text=element_text(size=20),
          legend.title = element_text(size=20),
          legend.key.size = unit(0.6, "cm"),
          legend.text = element_text(size=20),
          text = element_text(family = 'Calibri')) +
    scale_x_continuous(breaks = c(4000,8000,12000), 
                       labels = c(4000,8000,12000)) +
    scale_y_continuous(breaks = c(min(LI_ranked_df$rank_scaled)-1,
                                  max(LI_ranked_df$rank_scaled)-1,
                                  0.4,max(LI_ranked_df$runningScore)+direction),
                       labels = c(round(min(LI_ranked_df$rank),1),
                                  round(max(LI_ranked_df$rank),1),
                                  0,1)) +
    annotation_custom(pval) +
    geom_hline(yintercept = 0.4, linetype = "dashed") +
    geom_hline(yintercept = max(LI_ranked_df$rank_scaled)-1, 
               linetype = "dashed") +
    guides(size = 'none') 
  
  return(p)
  
}

######################
# EXTENDED DATA FIG 3C
######################

efig3C_gg <- custom_gsea_plot(gsea_cluster = 'VI_gene_cluster_1', pval_rnd=10)

efig3C_gg

plot_dimensions <- dev.size("in")

aspect_ratio <- plot_dimensions[1] / plot_dimensions[2]

aspect_ratio <- 1.616034

desired_height <- 5

desired_width <- desired_height * aspect_ratio

TIFF_dpi <- 600

# Save the plot with adjusted dimensions
ggsave(here('figures','efig3C.tiff'),
       plot = efig3C_gg,
       device = 'tiff',
       width = desired_width,
       height = desired_height,
       units = 'in', dpi = TIFF_dpi)

######################
# EXTENDED DATA FIG 3D
######################

efig3D_gg <- custom_gsea_plot(gsea_cluster = 'VI_gene_cluster_2', pval_rnd=3)

efig3D_gg

# Save the plot with adjusted dimensions
ggsave(here('figures','efig3D.tiff'),
       plot = efig3D_gg,
       device = 'tiff',
       width = desired_width,
       height = desired_height,
       units = 'in', dpi = TIFF_dpi)

######################
# EXTENDED DATA FIG 3E
######################

efig3E_gg <- custom_gsea_plot(gsea_cluster = 'VI_gene_cluster_3', pval_rnd=10)

efig3E_gg

# Save the plot with adjusted dimensions
ggsave(here('figures','efig3E.tiff'),
       plot = efig3E_gg,
       device = 'tiff',
       width = desired_width,
       height = desired_height,
       units = 'in', dpi = TIFF_dpi)

######################
# EXTENDED DATA FIG 3F
######################

efig3F_gg <- custom_gsea_plot(gsea_cluster = 'VI_gene_cluster_4', pval_rnd=10)

efig3F_gg

# Save the plot with adjusted dimensions
ggsave(here('figures','efig3F.tiff'),
       plot = efig3F_gg,
       device = 'tiff',
       width = desired_width,
       height = desired_height,
       units = 'in', dpi = TIFF_dpi)


######################
# EXTENDED DATA FIG 3G
######################

# compare VI predictor with LVI breast cancer signature from Kurozumi 2019
# and with VI endometrial cancer signature from Mannelqvist 2011


# 26 up genes from 35 gene signature from Mannelqvist et al 2011

`Mannelqvist 2011` <- c('SERPINE1','MMP1','IL8','GPR109B','I_963021','TNFAIP6',
                        'HSD11B1','IL6','I_962844','FN1','BCL2A1','IL1B',
                        'I_961179','PHLDA2','GBP5','PLEX','ANGPTL4','VCAM1',
                        'INHBA','CA2','FPR2','C3orf54','CTHRC1',
                        'MT2A','THBS2','G0S2')


# 42 up genes from 99 gene signature from Kurozumi et al 2019


`Kurozumi 2019` <- c("APOC1","APOE","CALML5","CCNB2","CDCA5","COX6C","DNAJA4",
                     "EEF1A2","ELF3","ERBB2","GNAS","HMGA1","HMGB3","HSPB1","IDH2",
                     "IFI27","ISG15","KRT18","KRT18P55","KRT19","KRT7","KRT8",
                     "LAPTM4B","LRRC26","LY6E","MMP11","MX1","NME1","NOP56","PGAP3",
                     "PITX1","PTTG1","S100P","SCD","SLC52A2","SLC9A3R1","SPDEF",
                     "TM7SF2","UBE2C","UBE2S","UCP2","YWHAZ")

VI_pancancer_genes_df <- rbind(data.frame(term = "Mannelqvist 2011",
                                          gene = `Mannelqvist 2011`),
                               data.frame(term = "Kurozumi 2019",
                                          gene = `Kurozumi 2019`))

set.seed(111)

VI_pancancer_gsea <- GSEA(VI_ranked_list, TERM2GENE = VI_pancancer_genes_df, 
                          pvalueCutoff = 1, seed = TRUE)

VI_ranked_df <- data.frame('gene' = names(VI_ranked_list),
                           'rank' = unname(VI_ranked_list))

VI_ranked_df <- VI_ranked_df %>%
  mutate(`Mannelqvist 2011` = ifelse(gene %in% `Mannelqvist 2011`, 'Yes','No')) %>%
  mutate(`Kurozumi 2019` = ifelse(gene %in% `Kurozumi 2019`, 'Yes','No')) %>%
  arrange(-rank) %>%
  mutate('Rank' = row_number())


efigG_H <- VI_ranked_df

addWorksheet(source_data, "E. Figure 3G-H")

writeData(source_data, sheet = "E. Figure 3G-H", 
          x = "E. Figure 3G-H", 
          startCol =  1, startRow = 1)

writeData(source_data, sheet = 'E. Figure 3G-H',
          x = efigG_H, startCol = 1, startRow = 3)


# adjust p values

gsea_pvals <- data.frame('Cluster' = VI_pancancer_gsea@result$ID,
                         'p_value' = VI_pancancer_gsea@result$pvalue)

gsea_pvals$p_value <- p.adjust(gsea_pvals$p_value, method = 'bonferroni',
                               n=length(gsea_pvals$p_value))


# function to get running enrichment score for custom plot

gsInfo <- function(object, geneSetID) {
  
  geneList <- object@geneList
  
  if (is.numeric(geneSetID))
    
    geneSetID <- object@result[geneSetID, "ID"]
  
  geneSet <- object@geneSets[[geneSetID]]
  
  exponent <- object@params[["exponent"]]
  
  df <- gseaScores(geneList, geneSet, exponent, fortify=TRUE)
  
  df$ymin <- 0
  
  df$ymax <- 0
  
  pos <- df$position == 1
  
  h <- diff(range(df$runningScore))/20
  
  df$ymin[pos] <- -h
  
  df$ymax[pos] <- h
  
  df$geneList <- geneList
  
  df$Description <- object@result[geneSetID, "Description"]
  
  return(df)
  
}


gseaScores <- getFromNamespace("gseaScores", "DOSE")

# re-scale ranking metric for plotting purposes

efigG_H$rank_scaled <- scale(efigG_H$rank, 
                                  center = mean(efigG_H$rank), 
                                  scale = max(efigG_H$rank) - min(efigG_H$rank))

efigG_H$rank_scaled <- efigG_H$rank_scaled[,1]

# function for custom gsea plot

custom_gsea_plot <- function(gsea_cluster, pval_rnd, plot_title = "plot title\n\n") {
  
  VI_cluster_color <- c('black')
  
  efigG_H$runningScore <- gsInfo(VI_pancancer_gsea, geneSetID = c(gsea_cluster))$runningScore
  
  gsea_cluster_sub <- ensym(gsea_cluster)
  
  direction <- ifelse(
    abs(min(gsInfo(VI_pancancer_gsea, geneSetID = gsea_cluster)$runningScore)) < 
      abs(max(gsInfo(VI_pancancer_gsea, geneSetID = gsea_cluster)$runningScore)),
    0.4, 1.4
  )
  
  pval <- grobTree(
    textGrob("p =", x = 0.53, y = 0.95, hjust = 0,
             gp = gpar(fontsize = 17, fontfamily = 'Calibri')),
    textGrob(
      eval(bquote(scientific_10(
        round(subset(gsea_pvals, Cluster == gsea_cluster)$p_value, pval_rnd)
      ))),
      x = 0.6, y = 0.95, hjust = 0,
      gp = gpar(fontsize = 17, fontfamily = 'Calibri')
    )
  )
  
  p <- ggplot(efigG_H, aes(x = Rank, y = 0.25, group = 1)) + 
    geom_line(aes(y = rank_scaled - 1), linetype = 'solid', linewidth = 1) +
    geom_line(aes(y = runningScore + direction),
              color = VI_cluster_color, linewidth = 1) +
    geom_bar(stat = "identity", alpha = 1, width = 0.5,
             aes(colour = !!gsea_cluster_sub, fill = !!gsea_cluster_sub),
             show.legend = TRUE) +
    scale_fill_manual(values = c('No' = 'transparent', 'Yes' = VI_cluster_color)) +
    scale_colour_manual(values = c('No' = 'transparent', 'Yes' = VI_cluster_color)) +
    theme_classic() +
    labs(
      title = plot_title,
      x     = "\nGene rank",
      y     = "Association with VI           Running ES\n"
    ) +
    theme(
      plot.title   = element_text(size = 22, hjust = 0.5),
      axis.title.x = element_text(size = 20),
      axis.title.y = element_text(size = 20, hjust = 0.8),
      axis.text.x  = element_text(colour = "black", size = 20),
      axis.text.y  = element_text(colour = "black", size = 20),
      legend.position  = "right",
      axis.text        = element_text(size = 20),
      legend.title     = element_text(size = 20),
      legend.key.size  = unit(0.6, "cm"),
      legend.text      = element_text(size = 20),
      text             = element_text(family = 'Calibri')
    ) +
    scale_x_continuous(
      breaks = c(4000, 8000, 12000),
      labels = c(4000, 8000, 12000)
    ) +
    scale_y_continuous(
      breaks = c(
        min(efigG_H$rank_scaled) - 1,
        max(efigG_H$rank_scaled) - 1,
        0.4,
        max(efigG_H$runningScore) + direction
      ),
      labels = c(
        round(min(efigG_H$rank), 1),
        round(max(efigG_H$rank), 1),
        0, 1
      )
    ) +
    annotation_custom(pval) +
    geom_hline(yintercept = 0.4, linetype = "dashed") +
    geom_hline(yintercept = max(efigG_H$rank_scaled) - 1, linetype = "dashed") +
    guides(size = 'none')
  
  return(p)
}


efig3G_gg <- custom_gsea_plot(
  gsea_cluster = "Mannelqvist 2011",
  pval_rnd = 10,
  plot_title = "26 genes upregulated with VI+ endometrial cancer\n"
)

efig3G_gg

aspect_ratio <- 1.616034

desired_height <- 5

desired_width <- desired_height * aspect_ratio

TIFF_dpi <- 600

ggsave(here('figures','efig3G.tiff'),
       plot = efig3G_gg,
       device = 'tiff',
       width = desired_width,
       height = desired_height,
       units = 'in', dpi = TIFF_dpi)


efig3H_gg <- custom_gsea_plot(
    gsea_cluster = "Kurozumi 2019",
    pval_rnd    = 10,
    plot_title  = "42 genes upregulated with LVI+ breast cancer\n"
  )


efig3H_gg

aspect_ratio <- 1.616034

desired_height <- 5

desired_width <- desired_height * aspect_ratio

TIFF_dpi <- 600

ggsave(here('figures','efig3H.tiff'),
       plot = efig3H_gg,
       device = 'tiff',
       width = desired_width,
       height = desired_height,
       units = 'in', dpi = TIFF_dpi)






#########################
# VI PREDICTOR DERIVATION
#########################

# nested cross validation in discovery cohort to identify best ML model

set.seed(111)

# 100x 70/30 train/test split, balanced for VI prevalence

trainIndex <- createDataPartition(discovery_dge$samples$VI, p = .7, 
                                  list = FALSE, 
                                  times = 100)


# Start the H2O cluster (locally) for machine learning

h2o.init()

if (interactive()) {
  
  yn <- utils::askYesNo("Skip the nested CV loop and load saved results?")
  skip_loop <- isTRUE(yn)
  
} else {
  
  skip_loop <- FALSE
  
}


if (skip_loop) {
  
  message("→ loading precomputed results…")

  } else {
  
  message("→ Running the full nested CV loop (this will take a while)…")
  
  aml_list <- list()
  
  aml_leader_list <- list()
  
  test_performance <- list()
  
  selected_features <- list()
  
  sizes <- list()
  
  # skip loop and load saved model for faster figure generation
  
  # outer loop
  
  for (i in 1:ncol(trainIndex)) {
    
    print(paste('Starting outer loop',i))
    
    # create outer train
    
    discovery_dge_train <- discovery_dge
    
    discovery_dge_train$samples <- discovery_dge$samples[trainIndex[,i],]
    
    discovery_dge_train$counts <- discovery_dge$counts[,trainIndex[,i]]
    
    # create outer test
    
    discovery_dge_test <- discovery_dge
    
    discovery_dge_test$samples  <- discovery_dge$samples[-trainIndex[,i],]
    
    discovery_dge_test$counts  <- discovery_dge$counts[,-trainIndex[,i]]
    
    # Feature filtering - derive differentially expressed genes for VI within train split
    
    design <- model.matrix(~0+LMPVI, discovery_dge_train$samples)
    
    discovery_dge_train <- estimateDisp(discovery_dge_train, design)
    
    colnames(design) = c('LMP','NST','VI')
    
    VI_contrasts = makeContrasts(VI-LMP,levels = design)
    
    fit <- glmFit(discovery_dge_train, design) 
    
    VI_lrt <- glmLRT(fit, contrast = VI_contrasts)
    
    VI_toptags <- topTags(VI_lrt, n = dim(discovery_dge_train)[1])
    
    VI_genes_cv <- rownames(VI_toptags$table[which((VI_toptags$table$FDR)<0.01),])
    
    bluered <- colorRampPalette(c("#007FFF", "white", "red"))(256)
    
    discovery_dge_train$samples$VI <- as.factor(discovery_dge_train$samples$VI)
    
    annotation_col = subset(discovery_dge_train$samples, select=c("VI","NST","LMP"))
    
    rownames(annotation_col) <- discovery_dge_train$samples$SAMPLE_ID
    
    colnames(discovery_dge_train) <- rownames(annotation_col)
    
    annotation_colors = list(
      VI = c("0" = "#390099", "1" = "#FF0054"))
    
    # plot heatmap at FDR < 0.01
    
    VIheatmap <- pheatmap::pheatmap((edgeR::cpm(discovery_dge_train,log=T)[VI_genes_cv,]),
                                    color = bluered, border_color = NA, scale = "row", show_rownames = FALSE,
                                    show_colnames = FALSE, cluster_rows = T, cluster_cols = T,
                                    clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean",
                                    legend = TRUE, annotation_legend = TRUE,
                                    annotation_col = annotation_col,
                                    annotation_colors = annotation_colors,
                                    clustering_method = "ward.D2", main = "~LMP+NST+VI, VI-LMP FDR < 0.01")
    
    # Derive gene clusters from heatmap:
    
    VI.gene.cluster <- as.character(cutree(as.hclust(VIheatmap$tree_row), k = 4)) # derive gene clusters from row clustering
    
    VI.gene.cluster <- data.frame(Cluster = VI.gene.cluster, row.names = NULL)
    
    rownames(VI.gene.cluster) <- make.names(VIheatmap$tree_row$labels, unique = TRUE)
    
    # Replot heatmap with gene clusters annotated:
    
    VIheatmapclust <- pheatmap::pheatmap((edgeR::cpm(discovery_dge_train,log=T)[VI_genes_cv,]),
                                         color = bluered, border_color = NA, scale = "row", show_rownames = FALSE,
                                         show_colnames = FALSE, cluster_rows = T, cluster_cols = T,
                                         clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean",
                                         legend = TRUE, annotation_legend = TRUE,
                                         annotation_row = VI.gene.cluster,
                                         annotation_col = annotation_col,
                                         annotation_colors = annotation_colors,
                                         treeheight_row = 40,
                                         treeheight_col = 40,
                                         fontsize = 6,
                                         clustering_method = "ward.D2", main = "~LMP+NST+VI, VI-LMP FDR < 0.01")
    
    VI.gene.cluster1 <- VI.gene.cluster %>% filter(Cluster == 1)
    
    VI.gene.cluster2 <- VI.gene.cluster %>% filter(Cluster == 2)
    
    VI.gene.cluster3 <- VI.gene.cluster %>% filter(Cluster == 3)
    
    VI.gene.cluster4 <- VI.gene.cluster %>% filter(Cluster == 4)
    
    
    
    # cluster 1 feature selection
    
    # get top genes by mean expression
    
    cluster_1_genes <- rowMeans(rbind(cpm(discovery_dge_train$counts,log=T)[rownames(VI.gene.cluster1),], NA))
    
    names(cluster_1_genes) <- c(rownames(VI.gene.cluster1), NA)
    
    cluster_1_genes <- cluster_1_genes[1:nrow(VI.gene.cluster1)]
    
    # proportional to size of cluster
    
    cluster_1_genes <- names(sort(cluster_1_genes,decreasing=TRUE))[1:round(length(cluster_1_genes)/(length(VI_genes_cv)/48))]
    
    
    # cluster 2 feature selection
    
    # get top genes by mean expression
    
    cluster_2_genes <- rowMeans(rbind(cpm(discovery_dge_train$counts,log=T)[rownames(VI.gene.cluster2),], NA))
    
    names(cluster_2_genes) <- c(rownames(VI.gene.cluster2), NA)
    
    cluster_2_genes <- cluster_2_genes[1:nrow(VI.gene.cluster2)]
    
    # proportional to size of cluster
    
    cluster_2_genes <- names(sort(cluster_2_genes,decreasing=TRUE))[1:round(length(cluster_2_genes)/(length(VI_genes_cv)/48))]
    
    
    # cluster 3 feature selection
    
    # get top genes by mean expression
    
    cluster_3_genes <- rowMeans(rbind(cpm(discovery_dge_train$counts,log=T)[rownames(VI.gene.cluster3),], NA))
    
    names(cluster_3_genes) <- c(rownames(VI.gene.cluster3), NA)
    
    cluster_3_genes <- cluster_3_genes[1:nrow(VI.gene.cluster3)]
    
    # proportional to size of cluster
    
    cluster_3_genes <- names(sort(cluster_3_genes,decreasing=TRUE))[1:round(length(cluster_3_genes)/(length(VI_genes_cv)/48))]
    
    
    # cluster 4 feature selection
    
    # get top genes by mean expression
    
    cluster_4_genes <- rowMeans(rbind(cpm(discovery_dge_train$counts,log=T)[rownames(VI.gene.cluster4),], NA))
    
    names(cluster_4_genes) <- c(rownames(VI.gene.cluster4), NA)
    
    cluster_4_genes <- cluster_4_genes[1:nrow(VI.gene.cluster4)]
    
    # proportional to size of cluster
    
    cluster_4_genes <- names(sort(cluster_4_genes,decreasing=TRUE))[1:round(length(cluster_4_genes)/(length(VI_genes_cv)/48))]
    
    
    
    VI_predictor_genes_cv <- c(cluster_1_genes, cluster_2_genes, cluster_3_genes, cluster_4_genes)
    
    VI_predictor_genes_cv <- VI_predictor_genes_cv[nzchar(VI_predictor_genes_cv)]
    
    sizes[[i]] <- length(na.omit(VI_predictor_genes_cv))
    
    print(paste0('Gene set size is ', sizes[[i]]))
    
    print(paste('Starting inner loop'))
    
    # re-subset train and subset test to selected features
    
    train <- data.frame(t(edgeR::cpm(discovery_dge_train,log=T)[na.omit(VI_predictor_genes_cv),]))
    
    
    test <- data.frame(t(edgeR::cpm(discovery_dge_test,log=T)[na.omit(VI_predictor_genes_cv),]))
    
    # store top features from each loop
    
    selected_features[[i]] <- VI_predictor_genes_cv
    
    # convert train and test to h2o format
    
    train$response <- discovery_dge_train$samples$VI
    
    train <- as.h2o(train)
    
    test$response <- discovery_dge_test$samples$VI
    
    test <- as.h2o(test)
    
    # Identify predictors and response
    
    y <- "response"
    
    x <- setdiff(names(train), y)
    
    # Run AutoML to tune parameters for 20 base models in 5-fold cross validation 
    # inner loop
    
    aml <- h2o.automl(x = x, y = y,
                      training_frame = train,
                      exclude_algos = 'DeepLearning',
                      max_models = 20,
                      seed = 1,
                      nfolds = 5,
                      project_name = paste('loop',i,sep = '_'))
    
    # print the AutoML leaderboard
    
    lb <- aml@leaderboard
    
    print(lb, n = nrow(lb))
    
    # store leading model from each feature cutoff
    
    aml_leader_list[[i]] <- aml@leader
    
    # store all models from each train loop
    
    aml_list[[i]] <- aml@leaderboard
    
    # predict on outer test using leading model from 5-fold cross validation
    
    test_perf <- h2o.performance(aml@leader, test)
    
    # store test performance across all outer loops
    
    test_performance[[i]] <- test_perf
    
  }
}

# remove all objects on the h2o cluster

h2o.removeAll(timeout_secs = 0, retained_elements = c())


# plot average AUC for model category across all outer train and test

# to ensure absolute reproducibility, load exact training metrics

aml_leader_list <- readRDS(here('data','VI_biomarker_CV_train_leader_metrics.rds'))

test_performance <- readRDS(here('data','VI_biomarker_CV_test_metrics.rds'))


top_models_train <- data.frame('Model' = unlist(lapply(aml_leader_list,function(x) x@algorithm)),
                               'AUC' = unlist(lapply(aml_leader_list,function(x) x@model$cross_validation_metrics_summary[2,1])),
                               'Dataset' = 'Train')

top_models_test <- data.frame('Model' = unlist(lapply(test_performance,function(x) x@algorithm)),
                              'AUC' = unlist(lapply(test_performance,function(x) x@metrics$AUC)),
                              'Dataset' = 'Test')

top_models <- rbind(top_models_train, top_models_test)

top_models$Dataset <- ordered(top_models$Dataset, levels = c('Train','Test'))

efig7A <- top_models %>%
  dplyr::group_by(Dataset, Model) %>%
  dplyr::summarize(Mean = round(mean(AUC, na.rm=TRUE),digits=2),
                   Median = round(median(AUC, na.rm=TRUE),digits=2),
                   sd = sd(AUC, na.rm=TRUE),
                   AUC = AUC) %>%
  filter(Dataset != 'Train')


addWorksheet(source_data, "E. Figure 7A")

writeData(source_data, sheet = "E. Figure 7A", 
          x = "E. Figure 6A", 
          startCol =  1, startRow = 1)

writeData(source_data, sheet = 'E. Figure 7A',
          x = efig7A, startCol = 1, startRow = 3)


######################
# EXTENDED DATA FIG 7A
######################

efig7A_gg <- ggplot(efig7A, aes(fill = Dataset,y=Mean, x=Model)) + 
  theme_classic() +
  geom_bar(position="dodge", stat="summary",fun='mean') +
  geom_errorbar(position=position_dodge(width=0.9),
                aes(x=Model, ymin=Mean-sd, ymax=Mean+sd),
                colour = 'black',width=0.1,linewidth=0.5) +
  scale_fill_manual(values=c("grey")) +
  labs(title="Cross validation outer test performance\n", 
       x = "\nModel", y = "Mean AUROC (sd)\n") +
  geom_text(aes(label = Mean), vjust = -0.5, hjust=1.2, 
            position = position_dodge(.9),
            size = 6, family = 'binomial') +
  theme(plot.title = element_text(size=22, hjust = 0.5),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        legend.position="none", 
        text = element_text(family = 'Calibri')) +
  ylim(0,1)

efig7A_gg

plot_dimensions <- dev.size("in")

aspect_ratio <- plot_dimensions[1] / plot_dimensions[2]

aspect_ratio <- 1.11828

desired_height <- 5.5

desired_width <- desired_height * aspect_ratio

TIFF_dpi <- 600

ggsave(here('figures','efig7A.tiff'),
       plot = efig7A_gg,
       device = 'tiff',
       width = desired_width,
       height = desired_height,
       units = 'in', dpi = TIFF_dpi)


# comparing best model across different gene set sizes

######################
# EXTENDED DATA FIG 7B
######################

top_models_12 <- readRDS(here('data','top_models_12.rds'))

top_models_24 <- readRDS(here('data','top_models_24.rds'))

top_models_48 <- rbind(top_models_train, top_models_test)

top_models_96 <- readRDS(here('data','top_models_96.rds'))

top_models_192 <- readRDS(here('data','top_models_192.rds'))

model_list <- list(top_models_12, top_models_24, 
                   top_models_48, top_models_96, top_models_192)

target_sizes <- c('2.5%','5%','10%','20%','40%')

efig7B <- do.call(rbind, lapply(seq_along(model_list), function(i) {
  df <- model_list[[i]]
  df$target_size <- target_sizes[i]  
  subset(df, Model == "glm")         
}))

addWorksheet(source_data, "E. Figure 7B")

writeData(source_data, sheet = "E. Figure 7B", 
          x = "E. Figure 7B", 
          startCol =  1, startRow = 1)

writeData(source_data, sheet = 'E. Figure 7B',
          x = efig7B, startCol = 1, startRow = 3)


efig7B <- efig7B %>%
  dplyr::group_by(Dataset, target_size) %>%
  dplyr::summarize(Mean = round(mean(AUC, na.rm=TRUE),digits=2),
                   Median = round(median(AUC, na.rm=TRUE),digits=2),
                   sd = sd(AUC, na.rm=TRUE),
                   AUC = round(AUC,)) %>%
  filter(Dataset != 'Train')

efig7B$target_size <- factor(
  efig7B$target_size,
  levels = c("2.5%", "5%", "10%", "20%", "40%")  # desired order
)

efig7B_gg <- ggplot(efig7B, aes(x = factor(target_size), y = Mean, fill = Dataset)) +
    theme_classic() +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
    geom_errorbar(aes(ymin = Mean - sd, ymax = Mean + sd),
                  position = position_dodge(width = 0.9),
                  width = 0.1, color = "black", linewidth = 0.5) +
    geom_text(aes(label = Mean),
              vjust = -0.5, hjust = 1.2,
              position = position_dodge(0.9),
              size = 6, family = "binomial") +
    scale_fill_manual(values = c("grey", "lightblue")) +
    scale_x_discrete(labels = c('12\n(2.5%)','24\n(5%)','48\n(10%)','96\n(20%)','192\n(40%)')) +
    labs(title = "Cross validation GLM performance\nby gene set size",
         x = "\n Median gene set size (% of total signature)", 
         y = "Mean AUC (sd)\n") +
    theme(plot.title = element_text(size = 22, hjust = 0.5),
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20),
          axis.text.x = element_text(size = 20),
          axis.text.y = element_text(size = 20),
          legend.position = "none",
          text = element_text(family = "Calibri")) +
    ylim(0, 1)

efig7B_gg

plot_dimensions <- dev.size("in")

aspect_ratio <- plot_dimensions[1] / plot_dimensions[2]

aspect_ratio <- 1.25

desired_height <- 5.5

desired_width <- desired_height * aspect_ratio

TIFF_dpi <- 600

ggsave(here('figures','efig7B.tiff'),
       plot = efig7B_gg,
       device = 'tiff',
       width = desired_width,
       height = desired_height,
       units = 'in', dpi = TIFF_dpi)






###################################
# CONFIGURE FINAL MODEL ON ALL DATA
###################################

# Final feature selection

# cluster 1

cluster_1_genes <- rowMeans(cpm(discovery_dge$counts,log=T)[rownames(VI_gene_cluster_1),])

cluster_1_genes <- names(sort(cluster_1_genes,decreasing=TRUE))[1:round(length(cluster_1_genes)/(length(VI_genes)/48))]


# cluster 2

cluster_2_genes <- rowMeans(cpm(discovery_dge$counts,log=T)[rownames(VI_gene_cluster_2),])

cluster_2_genes <- names(sort(cluster_2_genes,decreasing=TRUE))[1:round(length(cluster_2_genes)/(length(VI_genes)/48))]


# cluster 3

cluster_3_genes <- rowMeans(cpm(discovery_dge$counts,log=T)[rownames(VI_gene_cluster_3),])

cluster_3_genes <- names(sort(cluster_3_genes,decreasing=TRUE))[1:round(length(cluster_3_genes)/(length(VI_genes)/48))]


# cluster 4

cluster_4_genes <- rowMeans(cpm(discovery_dge$counts,log=T)[rownames(VI_gene_cluster_4),])

cluster_4_genes <- names(sort(cluster_4_genes,decreasing=TRUE))[1:round(length(cluster_4_genes)/(length(VI_genes)/48))]


# final predictor genes

VI_predictor_genes <- c(cluster_1_genes, cluster_2_genes, cluster_3_genes, cluster_4_genes)

VI_predictor_genes_up <- c(cluster_1_genes, cluster_2_genes, cluster_3_genes)

VI_predictor_genes_dn <- c(cluster_4_genes)

# save predictor genes

# saveRDS(VI_predictor_genes, file = here('data','VI_predictor_genes.rds'))

# saveRDS(VI_predictor_genes_up, file = here('data','VI_predictor_genes_up.rds'))

# saveRDS(VI_predictor_genes_dn, file = here('data','VI_predictor_genes_dn.rds'))



# subset discovery dge to final predictor genes

train_final <- data.frame(t(edgeR::cpm(discovery_dge,log=T)[VI_predictor_genes,]))

# convert to h2o format

train_final$response <- discovery_dge$samples$VI

train_final <- as.h2o(train_final)

# Identify predictors and response

y <- "response"

x <- setdiff(names(train_final), y)

# re-train best performing model on entire discovery cohort

aml_final <- h2o.automl(x = x, y = y,
                        training_frame = train_final,
                        include_algos = 'GLM',
                        max_models = 20,
                        seed = 1,
                        nfolds = 5,
                        project_name = 'final_VI_model')

aml_final@leader@model$cross_validation_metrics_summary

# to ensure absolute reproducibility, load exact model trained

aml_final_leader <- h2o.loadModel(here('data','GLM_1_AutoML_1_20230802_104952'))

pred_train <- h2o.predict(aml_final_leader, train_final)  

pred_train <- data.frame(as.matrix(pred_train))

discovery_dge$samples$VI_predictor <- as.numeric(pred_train$p1)

########
# FIG 4B
########

fig4B <- discovery_dge$samples %>% 
  select(VI, VI_predictor)

addWorksheet(source_data, "Figure 4B")

writeData(source_data, sheet = "Figure 4B", 
          x = "Figure 4B", 
          startCol =  1, startRow = 1)

writeData(source_data, sheet = 'Figure 4B',
          x = fig4B, startCol = 1, startRow = 3,
          rowNames = TRUE)

gene_roc <- roc(as.factor(fig4B$VI), fig4B$VI_predictor,
                levels=c(0, 1), ci = TRUE)

gene_roc


ciobj <- ci.se(gene_roc, specificities=seq(0, 1, l=25))

dat.ci <- data.frame(x = as.numeric(rownames(ciobj)),
                     lower = ciobj[, 1],
                     upper = ciobj[, 3])

rocobj <- pROC::plot.roc(gene_roc, print.thres = TRUE, print.auc = TRUE)

auc <- round(pROC::auc(gene_roc),2)

p_val <- round(wilcox.test(fig4B$VI_predictor ~ fig4B$VI)$p.value,9)

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

fig4B_gg <- pROC::ggroc(rocobj, colour = '#FF0054', size = 1, alpha = 1) +
  theme_classic() +
  annotation_custom(pval) +
  labs(title='Discovery cohort (n=103)\n', 
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

fig4B_gg

plot_dimensions <- dev.size("in")

aspect_ratio <- plot_dimensions[1] / plot_dimensions[2]

aspect_ratio <- 0.9122807

desired_height <- 6.5 

desired_width <- desired_height * aspect_ratio

TIFF_dpi <- 600

ggsave(here('figures','fig4B.tiff'),
       plot = fig4B_gg,
       device = 'tiff',
       width = desired_width,
       height = desired_height,
       units = 'in', dpi = TIFF_dpi)


# plot feature importance in model by original cluster membership

########
# FIG 4C
########

cluster_importance <- VI_gene_cluster %>% dplyr::mutate(variable = rownames(.))

tmp <- aml_final_leader@model$variable_importances

names(tmp)[1] <- 'variable'

fig4C <- merge(tmp, cluster_importance, by = 'variable') %>%
  select(variable, scaled_importance, Cluster)

addWorksheet(source_data, "Figure 4C")

writeData(source_data, sheet = "Figure 4C", 
          x = "Figure 4C", 
          startCol =  1, startRow = 1)

writeData(source_data, sheet = 'Figure 4C',
          x = fig4C, startCol = 1, startRow = 3)


fig4C_gg <- ggplot(fig4C, aes(x = reorder(variable,scaled_importance), 
                              y = scaled_importance, fill = Cluster)) + 
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = VI_cluster_colors) +
  labs(title=NULL, x = "VI predictor features\n",
       y = "\nGLM standardized coefficients") +
  theme_classic() +
  theme(plot.title = element_text(size=20, hjust = 0.5),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=18),
        legend.title = element_text(size=20),
        legend.key.size = unit(0.6, "cm"),
        legend.text = element_text(size=20),
        text = element_text(family = 'Calibri'))

fig4C_gg

plot_dimensions <- dev.size("in")

aspect_ratio <- plot_dimensions[1] / plot_dimensions[2]

aspect_ratio <- 0.65

desired_height <- 10

desired_width <- desired_height * aspect_ratio

TIFF_dpi <- 600

ggsave(here('figures','fig4C.tiff'),
       plot = fig4C_gg,
       device = 'tiff',
       width = desired_width,
       height = desired_height,
       units = 'in', dpi = TIFF_dpi)



# association of predictor with recurrence in discovery cohort

discovery_cox <- coxph(formula = Surv(time_7_year_RFS, RFS_7_year) ~ scale(VI_predictor),
                       data=discovery_dge$samples)

summary(discovery_cox)

names(discovery_cox$coefficients) <- 'Discovery (n=103)'




# multivariate analysis of VI predictor

######################
# EXTENDED DATA FIG 7C
######################

efig7C <- discovery_dge$samples %>%
  filter(Invasive_size <= 4) %>%
  mutate(STAS = factor(STAS)) %>%
  mutate(LI = factor(LI)) %>%
  mutate(VPI = factor(VPI)) %>%
  mutate(Gender = factor(Gender)) %>%
  mutate(X8thStg = factor(X8thStg)) %>%
  filter(WHO != 'M') %>%
  mutate(WHO = factor(WHO, levels = c('AISMIA','1','2','3'))) %>%
  mutate(WHO = fct_collapse(WHO,
                            "AIS/MIA/1/2" = c("AISMIA","1","2"),
                            "3" = c("3"))) %>%
  mutate(X8thStg = fct_collapse(X8thStg,
                                "0/IA" = c("0","IA1","IA2","IA3"),
                                "IB" = c("IB"))) %>%
  mutate(`Surgical Procedure` = factor(`Surgical.Procedure`, levels = c('Wedge','Segment','Lobe'))) %>%
  mutate(`Surgical Procedure` = fct_collapse(`Surgical Procedure`,
                                             'Sublobar' = c('Wedge','Segment'),
                                             'Lobectomy' = c('Lobe'))) %>%
  mutate(Race = factor(Race, levels = c('Asian','Black / African American',
                                        'Hispanic / Latino','Unknown','White'))) %>%
  mutate(Race = fct_collapse(Race,
                             'White' = 'White',
                             'Non-white' = c('Asian','Black / African American',
                                             'Hispanic / Latino'),
                             'Unknown' = c('Unknown'))) %>%
  filter(Race != 'Unknown') %>%
  mutate(Race = fct_drop(Race)) %>%
  select(time_7_year_RFS, RFS_7_year, Gender, Site, Age, Pack_yrs,
         X8thStg, WHO, Invasive_size, Total_size, `Surgical Procedure`,
         Race, X8thStg, VI_predictor)


addWorksheet(source_data, "E. Figure 7C")

writeData(source_data, sheet = "E. Figure 7C", 
          x = "E. Figure 7C", 
          startCol =  1, startRow = 1)

writeData(source_data, sheet = 'E. Figure 7C',
          x = efig7C, startCol = 1, startRow = 3)


discovery_cox_multi <- coxph(formula = Surv(time_7_year_RFS, RFS_7_year) ~ 
                               scale(VI_predictor) + Gender + Race + Site + 
                               Age + Pack_yrs + WHO + Invasive_size + 
                               Total_size + `Surgical Procedure` + X8thStg,
                             data=efig7C)

summary(discovery_cox_multi)

efig7C <- data.frame(
  HR = summary(discovery_cox_multi)$coefficients[, 2], 
  lower_ci = summary(discovery_cox_multi)$conf.int[, 3], 
  upper_ci = summary(discovery_cox_multi)$conf.int[, 4], 
  P = as.character(summary(discovery_cox_multi)$coefficients[, 5])
) %>% 
  rownames_to_column("Predictor") %>% 
  mutate(Predictor = case_when(
    Predictor == "scale(VI_predictor)" ~ "VI predictor",
    Predictor == "GenderM" ~ "Male",
    Predictor == "RaceWhite" ~ "White",
    Predictor == "SiteLHMC" ~ "LHMC",
    Predictor == "Pack_yrs" ~ "Pack years",
    Predictor == "X8thStgIB" ~ "Stage IB",
    Predictor == "WHO3" ~ "WHO G3",
    Predictor == "Invasive_size" ~ "Invasive size",
    Predictor == "Total_size" ~ "Total size",
    Predictor == "`Surgical Procedure`Lobectomy" ~ "Lobectomy",
    TRUE ~ Predictor
  )) %>% 
  bind_rows(
    tibble(Predictor = c("Female", "Non-white", "BMC","Stage 0/IA", "WHO G1/2", "Sublobar"),
           HR = 1, lower_ci = NA, upper_ci = NA, P = "Reference")
  ) %>% 
  mutate(Predictor = factor(Predictor, levels = c('VI predictor', 'Female', 'Male', 
                                                  'Non-white', 'White', 
                                                  'BMC', 'LHMC', 'Age', 'Pack years',
                                                  'Stage 0/IA', 
                                                  'Stage IB', 'WHO G1/2',
                                                  'WHO G3', 'Invasive size', 
                                                  'Total size', 'Sublobar', 
                                                  'Lobectomy'))) %>% 
  arrange(Predictor) %>% 
  mutate(across(c(HR, lower_ci, upper_ci), as.numeric)) %>% 
  mutate(across(c(HR, lower_ci, upper_ci), ~signif(., digits = 3))) %>%
  mutate(across(-1, ~ signif(as.numeric(.), digits = 3))) %>%
  mutate(P = as.character(signif(as.numeric(P), digits = 3))) %>%
  mutate(P = sub("e", " %*% 10^", P)) %>%
  mutate(P = case_when(
    Predictor %in% c("Female", "Non-white", "BMC", "Stage 0/IA", "WHO G1/2", "Sublobar") ~ "Reference",
    TRUE ~ P
  ))

hr_max <- ceiling(max(efig7C$upper_ci, na.rm = TRUE) / 10) * 10

efig7C_gg <- ggplot(efig7C, aes(x=Predictor, y=HR)) + coord_flip() +
  geom_hline(yintercept = 1, linetype="dotted") + 
  geom_errorbar(aes(x=Predictor, ymin=lower_ci, ymax=upper_ci), 
                width=0.5, color="#FF0054") +
  geom_point(pch=15, size=3) +
  theme_minimal() + 
  theme(legend.position = "none", 
        plot.title = element_text(size=20),
        aspect.ratio = 0.8,
        panel.grid = element_blank(),
        axis.line.x = element_line(),
        axis.ticks.x = element_line(),
        axis.text.x=element_text(size=15),
        axis.text.y=element_text(size=15),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        text = element_text(family = 'Calibri')) + 
  scale_y_continuous(trans = "log", expand = c(0,0), breaks = c(0.02,2,8,32), 
                     limits = c(0.02, 1200)) + 
  scale_x_discrete(limits=c(rev(as.character(efig7C$Predictor)),' ')) + 
  xlab("Predictor\n") + ylab("\nHazard Ratio (95% CI)") + 
  geom_text(aes(x=Predictor, y=hr_max + 300, label=P), 
            size=5, family = 'Calibri', parse = TRUE) + 
  ggtitle("Discovery cohort (n=99) 7-Year RFS\n") + 
  annotate("text", x = 18.2, y = 330, label = "p value", 
           size = 5, fontface = 'bold', family = 'Calibri')

efig7C_gg

plot_dimensions <- dev.size("in")

aspect_ratio <- plot_dimensions[1] / plot_dimensions[2]

aspect_ratio <- 0.8632597

desired_height <- 8  

desired_width <- desired_height * aspect_ratio

TIFF_dpi <- 600

ggsave(here('figures','efig7C.tiff'),
       plot = efig7C_gg,
       device = 'tiff',
       width = desired_width,
       height = desired_height,
       units = 'in', dpi = TIFF_dpi)


h2o.removeAll(timeout_secs = 0, retained_elements = c())

#########################
# VI PREDICTOR VALIDATION
#########################

# convert gene names to hgnc

validation_dge <- ensemble_to_hgnc(validation_dge)


#######################
# GENE FILTERING (LHMC)
#######################

# subset to filtered genes in discovery cohort

validation_dge$counts <- validation_dge$counts[which(rownames(validation_dge$counts) %in% rownames(discovery_dge$counts)),]

dim(validation_dge$counts)

validation_dge <- DGEList(counts = validation_dge$counts, samples = validation_dge$samples)

head(validation_dge$samples$lib.size)

# Trimmed mean of m values (TMM) normalization

validation_dge <- calcNormFactors(validation_dge) 

##########################################
# ASSESS BATCH EFFECTS (VALIDATION COHORT)
##########################################

# PCA dimensionality reduction

pca_batch <- prcomp(t(cpm(validation_dge,log=T)),scale=T, center = T)

# Visualize possible batch effects

autoplot(pca_batch, data = validation_dge$samples, 
         colour = 'Sequence_batch', label = TRUE) 


# correct batch effect within validation cohort

validation_dge$counts <- ComBat_seq(validation_dge$counts, 
                                    validation_dge$samples$Sequence_batch)

pca_batch <- prcomp(t(cpm(validation_dge,log=T)),scale=T, center = T)

# Visualize possible batch effects

autoplot(pca_batch, data = validation_dge$samples, 
         colour = 'Sequence_batch', label = TRUE) 


# reference combat to discovery cohort

validation_dge$counts <- ComBat(cbind(cpm(discovery_dge$counts,log=T), 
                                      cpm(validation_dge$counts,log=T)),
                                c(rep(1,ncol(discovery_dge$counts)),
                                  rep(2,ncol(validation_dge$counts))), 
                                ref.batch = 1)

validation_dge$counts <- validation_dge$counts[,(ncol(discovery_dge$counts)+1):ncol(validation_dge$counts)]


# external validation of predictor

test_final <- data.frame(t(validation_dge$counts[VI_predictor_genes,]))

# convert to h2o format

validation_dge$samples$VI <- as.factor(validation_dge$samples$VI)

test_final$response <- validation_dge$samples$VI

test_final <- as.h2o(test_final)

# Identify predictors and response

y <- "response"

x <- setdiff(names(test_final), y)

aml_final_leader <- h2o.loadModel(here('data','GLM_1_AutoML_1_20230802_104952'))

pred_test <- h2o.predict(aml_final_leader, test_final)  

pred_test <- data.frame(as.matrix(pred_test))

validation_dge$samples$VI_predictor <- as.numeric(pred_test$p1)

########
# FIG 4D
########

fig4D <- validation_dge$samples %>%
  select(VI, VI_predictor)

addWorksheet(source_data, "Figure 4D")

writeData(source_data, sheet = "Figure 4D", 
          x = "Figure 4D", 
          startCol =  1, startRow = 1)

writeData(source_data, sheet = 'Figure 4D',
          x = fig4D, startCol = 1, startRow = 3,
          rowNames = TRUE)

gene_roc <- roc(as.factor(fig4D$VI), fig4D$VI_predictor,
                levels=c(0, 1), ci = TRUE)

gene_roc

ciobj <- ci.se(gene_roc, specificities=seq(0, 1, l=25))

dat.ci <- data.frame(x = as.numeric(rownames(ciobj)),
                     lower = ciobj[, 1],
                     upper = ciobj[, 3])

rocobj <- pROC::plot.roc(gene_roc, print.thres = TRUE, print.auc = TRUE)

auc <- round(pROC::auc(gene_roc),2)

p_val <- round(wilcox.test(fig4D$VI_predictor ~ fig4D$VI)$p.value,8)

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

fig4D_gg <- pROC::ggroc(rocobj, colour = '#FF0054', size = 1, alpha = 1) +
  theme_classic() +
  annotation_custom(pval) +
  labs(title='Validation cohort (n=59)\n', 
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

fig4D_gg

plot_dimensions <- dev.size("in")

aspect_ratio <- plot_dimensions[1] / plot_dimensions[2]

aspect_ratio <- 0.9122807

desired_height <- 6.5

desired_width <- desired_height * aspect_ratio

TIFF_dpi <- 600

ggsave(here('figures','fig4D.tiff'),
       plot = fig4D_gg,
       device = 'tiff',
       width = desired_width,
       height = desired_height,
       units = 'in', dpi = TIFF_dpi)




# performance classifying VI in NST + VI samples only (removing LMP)

########
# FIG 4E
########

fig4E <- validation_dge$samples %>% 
  filter(LMPVI != 'LMP') %>%
  select(VI, VI_predictor)

addWorksheet(source_data, "Figure 4E")

writeData(source_data, sheet = "Figure 4E", 
          x = "Figure 4E", 
          startCol =  1, startRow = 1)

writeData(source_data, sheet = 'Figure 4E',
          x = fig4E, startCol = 1, startRow = 3,
          rowNames = TRUE)

gene_roc <- roc(as.factor(fig4E$VI), fig4E$VI_predictor,
                levels=c(0, 1), ci = TRUE)

gene_roc

ciobj <- ci.se(gene_roc, specificities=seq(0, 1, l=25))

dat.ci <- data.frame(x = as.numeric(rownames(ciobj)),
                     lower = ciobj[, 1],
                     upper = ciobj[, 3])

rocobj <- pROC::plot.roc(gene_roc, print.thres = TRUE, print.auc = TRUE)

auc <- round(pROC::auc(gene_roc),2)

p_val <- round(wilcox.test(fig4E$VI_predictor ~ fig4E$VI)$p.value,8)

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

fig4E_gg <- pROC::ggroc(rocobj, colour = '#FF0054', size = 1, alpha = 1) +
  theme_classic() +
  annotation_custom(pval) +
  labs(title='Validation cohort (n=52, NST+VI only)\n', 
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

fig4E_gg

plot_dimensions <- dev.size("in")

aspect_ratio <- plot_dimensions[1] / plot_dimensions[2]

aspect_ratio <- 0.9122807

desired_height <- 6.5 

desired_width <- desired_height * aspect_ratio

TIFF_dpi <- 600

ggsave(here('figures','fig4E.tiff'),
       plot = fig4E_gg,
       device = 'tiff',
       width = desired_width,
       height = desired_height,
       units = 'in', dpi = TIFF_dpi)



# look at association of predictor scores with novel grading system

validation_dge$samples$LMPVI <- as.factor(validation_dge$samples$LMPVI)

levels(validation_dge$samples$LMPVI) <- c('LMP (G1)','NST (G2)','VI (G3)')

efig7D_1 <- validation_dge$samples %>%
  select(LMPVI, VI_predictor)

addWorksheet(source_data, "E. Figure 7D-1")

writeData(source_data, sheet = "E. Figure 7D-1", 
          x = "E. Figure 7D-1", 
          startCol =  1, startRow = 1)

writeData(source_data, sheet = 'E. Figure 7D-1',
          x = efig7D_1, startCol = 1, startRow = 3,
          rowNames = TRUE)

########################
# EXTENDED DATA FIG 7D_1
########################

pval <- efig7D_1 %>% 
  wilcox_test(
    VI_predictor ~ LMPVI
  ) %>% 
  add_xy_position() %>% 
  mutate(p = sub("e","%.% 10^",p) )

efig7D_1_gg <- ggplot(data = efig7D_1, aes(x = LMPVI, y = VI_predictor)) +
  theme_classic() + 
  scale_fill_manual(values = c("#1f78b4","#33a02c","#FF0054")) +
  geom_violin(lwd=0.5, width = 1, aes(fill = LMPVI)) +
  geom_boxplot(color = "black", fill = 'white', width = .15, 
               outlier.alpha = 0, lwd=0.5) +
  geom_jitter(colour="black",fill='white', size=1, alpha=1, 
              width = 0.2, shape=21) +
  labs(title="Validation cohort (n=59)\n", 
       x = "\nNovel Grading System", y = "VI predictor score\n") +
  theme(
    plot.title = element_text(size=22, hjust = 0.5),
    axis.title.x = element_text(size=20),
    axis.title.y = element_text(size=20),
    axis.text.x = element_text(size=20),
    legend.position="none", 
    axis.text=element_text(size=20),
    legend.title = element_text(),
    text = element_text(family = 'Calibri')) +
  ylim(0,max(efig7D_1$VI_predictor)+0.4) +
  ggprism::add_pvalue(pval,label = 'p', step.increase = 0.1,
             family = 'Calibri', label.size = 6, parse = TRUE)


efig7D_1_gg$layers[[4]]$aes_params$family <- "Calibri"

efig7D_1_gg

plot_dimensions <- dev.size("in")

aspect_ratio <- plot_dimensions[1] / plot_dimensions[2]

aspect_ratio <- 0.9644513

desired_height <- 6.5 

desired_width <- desired_height * aspect_ratio

TIFF_dpi <- 600

ggsave(here('figures','efig7D_1.tiff'),
       plot = efig7D_1_gg,
       device = 'tiff',
       width = desired_width,
       height = desired_height,
       units = 'in', dpi = TIFF_dpi)





# association with WHO 2015 grading system

validation_dge$samples$WHO <- factor(validation_dge$samples$WHO, 
                                     levels = c('AISMIA','M','1','2','3'))

levels(validation_dge$samples$WHO) <- c('AISMIA','M','G1','G2','G3')

efig7D_2 <- validation_dge$samples %>% 
  filter(WHO != "M") %>%
  select(WHO, VI_predictor)

efig7D_2$WHO <- droplevels(efig7D_2$WHO)

addWorksheet(source_data, "E. Figure 7D-2")

writeData(source_data, sheet = "E. Figure 7D-2", 
          x = "E. Figure 7D-2", 
          startCol =  1, startRow = 1)

writeData(source_data, sheet = 'E. Figure 7D-2',
          x = efig7D_2, startCol = 1, startRow = 3,
          rowNames = TRUE)

########################
# EXTENDED DATA FIG 7D_2
########################

pval <- efig7D_2 %>% 
  wilcox_test(
    VI_predictor ~ WHO
  ) %>% 
  add_xy_position() %>% 
  mutate(p = sub("e","%.% 10^",p) )

efig7D_2_gg <- ggplot(data = efig7D_2, aes(x = WHO, y = VI_predictor)) +
  theme_classic() + 
  scale_fill_manual(values = c("#1f78b4","#33a02c","#FF0054")) +
  geom_violin(lwd=0.5, width = 1, aes(fill=WHO)) +
  geom_boxplot(color = "black", fill = 'white', width = .15, 
               outlier.alpha = 0, lwd=0.5) +
  geom_jitter(colour="black",fill='white', size=1, alpha=1, 
              width = 0.2, shape=21) +
  labs(title="Validation cohort (n=57)\n", 
       x = "\nWHO 2015 Grading System", y = "VI predictor score\n") +
  theme(
    plot.title = element_text(size=22, hjust = 0.5),
    axis.title.x = element_text(size=20),
    axis.title.y = element_text(size=20),
    axis.text.x = element_text(size=20),
    legend.position="none", 
    axis.text=element_text(size=20),
    legend.title = element_text(),
    text = element_text(family = 'Calibri')) +  
  ylim(0,max(efig7D_2$VI_predictor)+0.4) +
  ggprism::add_pvalue(pval,label = 'p', step.increase = 0.05,
             family = 'Calibri', label.size = 6, parse = TRUE)

efig7D_2_gg$layers[[4]]$aes_params$family <- "Calibri"

efig7D_2_gg

ggsave(here('figures','efig7D_2.tiff'),
       plot = efig7D_2_gg,
       device = 'tiff',
       width = desired_width,
       height = desired_height,
       units = 'in', dpi = TIFF_dpi)



################################
# OUTCOME ANALYSIS (VALIDATION)
################################

# univariate coxph analyis of VI predictor association with outcome

fig4F_1 <- validation_dge$samples %>%
  select(time_7_year_RFS, RFS_7_year, VI_predictor)

addWorksheet(source_data, "Figure 4F-1")

writeData(source_data, sheet = "Figure 4F-1", 
          x = "Figure 4F-1", 
          startCol =  1, startRow = 1)

writeData(source_data, sheet = 'Figure 4F-1',
          x = fig4F_1, startCol = 1, startRow = 3,
          rowNames = TRUE)

validation_cox <- coxph(formula = Surv(time_7_year_RFS, RFS_7_year) ~ scale(VI_predictor),
                        data=fig4F_1)

summary(validation_cox)

names(validation_cox$coefficients) <- 'Validation (n=59)'


# predicting outcome in VI- samples only

fig4F_2 <- validation_dge$samples %>%
  dplyr::filter(LMPVI != 'VI (G3)') %>%
  select(time_7_year_RFS, RFS_7_year, VI_predictor)

addWorksheet(source_data, "Figure 4F-2")

writeData(source_data, sheet = "Figure 4F-2", 
          x = "Figure 4F-2", 
          startCol =  1, startRow = 1)

writeData(source_data, sheet = 'Figure 4F-2',
          x = fig4F_2, startCol = 1, startRow = 3,
          rowNames = TRUE)

# univariate coxph analyis of VI biomarker association with outcome in VI- samples

validation_cox_VI_neg <- coxph(formula = Surv(time_7_year_RFS, RFS_7_year) ~ scale(VI_predictor),
                               data=fig4F_2)

summary(validation_cox_VI_neg)

names(validation_cox_VI_neg$coefficients) <- 'Validation, VI- (n=42)'

h2o.removeAll(timeout_secs = 0, retained_elements = c())



##############
# LOAD TRACERX
##############

# TRACERx 421 cohort
# originally downloaded from https://doi.org/10.5281/zenodo.7683605 and https://doi.org/10.5281/zenodo.76033862

TRACERx_all_patient <- readRDS(here('data','20221109_TRACERx421_all_patient_df.rds'))

dim(TRACERx_all_patient) # 421 patients

TRACERx_all_tumor <- readRDS(here('data','20221109_TRACERx421_all_tumour_df.rds'))

dim(TRACERx_all_tumor) # 432 tumors

# TRACERx tumor evolutionary metrics

TRACERx_tumor_evo_metrics <- read.table(here('data','20221110_TRACERx421_evolutionary_metrics.tsv'),
                                        sep ='\t', header = TRUE)

dim(TRACERx_tumor_evo_metrics) # 432 x 172

TRACERx_tumor_evo_metrics <- TRACERx_tumor_evo_metrics %>%
  dplyr::rename(tumour_id_muttable_cruk = tumour_id)

TRACERx_all_tumor <- TRACERx_all_tumor %>%
  left_join(TRACERx_tumor_evo_metrics, 
            by = 'tumour_id_muttable_cruk', 
            multiple = 'all')

dim(TRACERx_all_tumor)

TRACERx_all_regions <- read_fst(here('data','2022-10-14_clinicohistopathological_data.fst'))

dim(TRACERx_all_regions) # 1515 x 10

# TRACERx I-TED intratumor heterogeneity metric

TRACERx_ITED <- read_fst(here('data','2022-10-21_ited_primaries.fst'))

dim(TRACERx_ITED) # 813 per region

TRACERx_ITED <- TRACERx_ITED %>%
  dplyr::rename(sample_name_cruk = region)

TRACERx_all_regions <- TRACERx_all_regions %>%
  left_join(TRACERx_ITED, by = 'sample_name_cruk', multiple = 'all')

dim(TRACERx_all_regions)


# TRACERx 421 cohort RNA-seq data

TRACERx_rsem_counts <- read_fst(here('data','2022-10-17_rsem_counts_mat.fst'))

dim(TRACERx_rsem_counts) # 28073 x 1052


# TRACERx additional data (e.g. ctDNA) about LUADs

TRACERx_LUAD_patient <- read_excel(here('data','TRACERx_LUAD_path_data_STAS_necrosis_ctDNA.xlsx'),
                                   sheet = 3)

dim(TRACERx_LUAD_patient) # 242 LUAD patients

TRACERx_LUAD_patient <- TRACERx_LUAD_patient %>%
  dplyr::rename(tumour_id_muttable_cruk = Patient_ID)

TRACERx_LUAD_tumor <- read_excel(here('data','TRACERx_LUAD_path_data_STAS_necrosis_ctDNA.xlsx'),
                                 sheet = 4)

dim(TRACERx_LUAD_tumor) # 248 LUAD tumors

TRACERx_LUAD_tumor <- TRACERx_LUAD_tumor %>%
  dplyr::rename(tumour_id_muttable_cruk = `Tumour ID`)

TRACERx_LUAD_regions <- read_excel(here('data','TRACERx_LUAD_path_data_STAS_necrosis_ctDNA.xlsx'),
                                   sheet = 5)

dim(TRACERx_LUAD_regions) # 805 LUAD regions

TRACERx_LUAD_regions <- TRACERx_LUAD_regions %>%
  dplyr::rename(tumour_id_muttable_cruk = Tumour_ID)

# combine patient, tumor, region metadata

TRACERx_all_pheno <- TRACERx_all_regions %>%
  left_join(TRACERx_all_tumor, by = 'tumour_id_muttable_cruk', multiple = 'all') %>%
  left_join(TRACERx_all_patient, by = 'tumour_id_muttable_cruk', multiple = 'all') %>%
  mutate(tumour_id_muttable_cruk = str_replace(str_replace(tumour_id_muttable_cruk, "Cluster", "Tumour"), "-", "_"))


dim(TRACERx_all_pheno) # 1515 x 78

# combine patient, tumor, region LUAD specific metadata

TRACERx_LUAD_pheno <- TRACERx_LUAD_regions %>%
  left_join(TRACERx_LUAD_tumor, by = 'tumour_id_muttable_cruk', multiple = 'all') %>%
  left_join(TRACERx_LUAD_patient, by = 'tumour_id_muttable_cruk', multiple = 'all') %>%
  dplyr::rename(sample_name_cruk = Sample_ID)

dim(TRACERx_LUAD_pheno) # 805 x 40

# add LUAD specific metadata to overall metadata

TRACERx_all_pheno <- TRACERx_all_pheno %>%
  left_join(TRACERx_LUAD_pheno, by = 'sample_name_cruk', multiple = 'all')

dim(TRACERx_all_pheno) # 1515 x 115

# TRACERx used 7th edition TNM staging

TRACERx_all_pheno <- TRACERx_all_pheno %>% 
  filter(histology_3 == 'LUAD') %>%
  filter(pathologyTNM == 'IA' | pathologyTNM == 'IB' | pathologyTNM == 'IIA' | pathologyTNM == 'IIB') #%>%
#filter(Size <= 40) # this is total size not invasive size

dim(TRACERx_all_pheno) # 297 stage I LUAD regions


# outcome analysis

TRACERx_all_pheno_patients <- TRACERx_all_pheno %>%
  distinct(tumour_id_muttable_cruk.x, .keep_all = TRUE)

dim(TRACERx_all_pheno_patients) # 118 stage I LUAD patients

# Converting to 5 year OS, DSS, RFS data below 

TRACERx_all_pheno_patients <- TRACERx_all_pheno_patients %>%
  mutate(os_5_year = ifelse(cens_os == 1 & os_time > 1825, 0 , cens_os)) %>%
  mutate(time_os_5_year = ((os_time / 365) * 12)) %>%
  mutate(time_os_5_year = ifelse(time_os_5_year >= 60, 60, time_os_5_year)) %>%
  mutate(dfs_5_year = ifelse(cens_dfs == 1 & dfs_time > 1825, 0 , cens_dfs)) %>%
  mutate(time_dfs_5_year = ((dfs_time / 365) * 12)) %>%
  mutate(time_dfs_5_year = ifelse(time_dfs_5_year >= 60, 60, time_dfs_5_year)) %>%
  mutate(rfs = ifelse(Relapse_cat_new == 'No rec',0,1)) %>%
  mutate(rfs_5_year = ifelse(rfs == 1 & lung_event_time > 1825, 0 , rfs)) %>%
  mutate(time_rfs_5_year = ((lung_event_time / 365) * 12)) %>%
  mutate(time_rfs_5_year = ifelse(time_rfs_5_year >= 60, 60, time_rfs_5_year)) %>%
  mutate(ctDNA = dplyr::coalesce(`Pre-op ctD (Tx100 assay)`, `Pre-op ctD (Tx421 assay)`))


# subset to regions with rna-seq data

rownames(TRACERx_rsem_counts) <- TRACERx_rsem_counts$gene_id

TRACERx_rsem_counts <- TRACERx_rsem_counts %>%
  dplyr::select(-gene_id) %>%
  dplyr::select(which(colnames(.) %in% TRACERx_all_pheno$sample_name_cruk))

dim(TRACERx_rsem_counts)

TRACERx_all_pheno <- TRACERx_all_pheno %>%
  filter(sample_name_cruk %in% colnames(TRACERx_rsem_counts))

dim(TRACERx_all_pheno)



# put counts in same order as clinical data

TRACERx_rsem_counts <- TRACERx_rsem_counts[,TRACERx_all_pheno$sample_name_cruk]


# create DGE list object

TRACERx_dge <- DGEList(counts = TRACERx_rsem_counts, samples = TRACERx_all_pheno) # create DGE list

TRACERx_dge_stageI <- TRACERx_dge

TRACERx_dge_stageI$samples <- TRACERx_dge_stageI$samples %>%
  filter(pathologyTNM == 'IA' | pathologyTNM == 'IB') %>%
  filter(Size <= 40)

TRACERx_dge_stageI$counts <- TRACERx_dge_stageI$counts[,which(colnames(TRACERx_dge_stageI$counts) %in% rownames(TRACERx_dge_stageI$samples))]


################
# GENE FILTERING
################

# Filter TRACERx counts to filtered discovery cohort counts

TRACERx_dge_stageI$counts <- TRACERx_dge_stageI$counts[which(rownames(TRACERx_dge_stageI$counts) %in% rownames(discovery_dge$counts)),]

dim(TRACERx_dge_stageI$counts)

discovery_dge$counts <- discovery_dge$counts[which(rownames(discovery_dge$counts) %in% rownames(TRACERx_dge_stageI$counts)),]

dim(discovery_dge$counts)


# put genes in same order as train

TRACERx_dge_stageI$counts <- TRACERx_dge_stageI$counts[rownames(discovery_dge$counts),]



TRACERx_dge_stageI <- DGEList(counts = TRACERx_dge_stageI$counts, samples = TRACERx_dge_stageI$samples)

head(TRACERx_dge_stageI$samples$lib.size)

TRACERx_dge_stageI <- calcNormFactors(TRACERx_dge_stageI) # Trimmed mean of m values (TMM) normalization


#######################
# CORRECT BATCH EFFECTS
#######################

# Reference combat TRACERx count matrix with discovery cohort count matrix


TRACERx_dge_stageI$counts <- ComBat(cbind(edgeR::cpm(discovery_dge$counts,log=T), edgeR::cpm(TRACERx_dge_stageI$counts,log=T)),
                                    c(rep(1,ncol(discovery_dge$counts)),rep(2,ncol(TRACERx_dge_stageI$counts))),
                                    ref.batch = 1)

TRACERx_dge_stageI$counts <- TRACERx_dge_stageI$counts[,(ncol(discovery_dge$counts)+1):ncol(TRACERx_dge_stageI$counts)]

length(which(VI_predictor_genes %in% rownames(TRACERx_dge_stageI$counts))) # 48/48 genes present



############################
# GENERATE MODEL PREDICTIONS
############################


# use final VI predictor model to predict

test_TRACERx <- data.frame(t(TRACERx_dge_stageI$counts[which(rownames(TRACERx_dge_stageI$counts) %in% VI_predictor_genes),]))

# convert to h2o format and add dummy class variable

test_TRACERx$response <- rep(c(0,1),length.out = ncol(TRACERx_dge_stageI$counts))

test_TRACERx <- as.h2o(test_TRACERx)

# Identify predictors and response

y <- "response"

x <- setdiff(names(test_TRACERx), y)

aml_final_leader <- h2o.loadModel(here('data','GLM_1_AutoML_1_20230802_104952'))

pred_test <- h2o.predict(aml_final_leader, test_TRACERx)  

pred_test <- data.frame(as.matrix(pred_test))

TRACERx_dge_stageI$samples$VI_predictor <- as.numeric(pred_test$p1)


##############
# ITH ANALYSIS
##############



TRACERx_ITH <- as.data.frame(table(TRACERx_dge_stageI$samples$tumour_id_per_patient)) %>%
  dplyr::rename(tumour_id_per_patient = Var1) %>%
  dplyr::rename(tumour_regions = Freq) %>%
  left_join(TRACERx_dge_stageI$samples, by = 'tumour_id_per_patient', multiple = 'all')

# remove tumors with only 1 region sampled

TRACERx_ITH <- TRACERx_ITH %>% filter(tumour_regions > 1) %>%
  filter(pathologyTNM == 'IA' | pathologyTNM == 'IB') %>%
  filter(Size <= 40)

# downsample all tumors to 2 regions each

set.seed(111)

TRACERx_ITH_2 <- TRACERx_ITH %>%
  group_by(tumour_id_per_patient) %>%
  mutate(region_number = row_number()) %>%
  slice_sample(n = 2) %>%
  mutate(region_number = row_number()) %>%
  mutate_at(vars(c('region_number')), factor) %>% 
  ungroup()

dim(TRACERx_ITH_2)

TRACERx_ITH_2 <- TRACERx_ITH_2 %>%
  mutate(tumour_regions = as.factor(tumour_regions))



# look at correlation in scores between region 1 and region 2 (randomly chosen if > 2 regions)
# create new 'region_number' as not all regions may have been sequenced
# if more than 2 regions are available, randomly sample to 2

###########
# FIGURE 5C
###########

set.seed(111)

data_long <- TRACERx_ITH %>%
  group_by(tumour_id_per_patient) %>%
  mutate(region_number = row_number()) %>%
  slice_sample(n = 2) %>%
  mutate(region_number = row_number()) %>%
  mutate_at(vars(c('region_number')), factor) %>% 
  ungroup() %>%
  dplyr::select(tumour_id_per_patient, region_number, VI_predictor)

fig5C <- tidyr::spread(data_long, region_number, VI_predictor)

addWorksheet(source_data, "Figure 5C")

writeData(source_data, sheet = "Figure 5C", 
          x = "Figure 5C", 
          startCol =  1, startRow = 1)

writeData(source_data, sheet = 'Figure 5C',
          x = fig5C, startCol = 1, startRow = 3)

fig5C_gg <- ggplot(data = fig5C, aes(x = `1`, y = `2`)) +
  geom_point() +
  theme_classic() +
  labs(title="TRACERx stage I LUAD (n=63 tumors)\n", x = "\nRegion 1 VI predictor score", y = "Region 2 VI predictor score\n") +
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

fig5C_gg

plot_dimensions <- dev.size("in")

aspect_ratio <- plot_dimensions[1] / plot_dimensions[2]

aspect_ratio <- 0.9815436

desired_height <- 6.5

desired_width <- desired_height * aspect_ratio

TIFF_dpi <- 600

ggsave(here('figures','fig5C.tiff'),
       plot = fig5C_gg,
       device = 'tiff',
       width = desired_width,
       height = desired_height,
       units = 'in', dpi = TIFF_dpi)



# rank genes by how correlated they are between region 1 and region 2

###########
# FIGURE 5D
###########

TRACERx_ITH_counts <- TRACERx_dge_stageI$counts[,which(colnames(TRACERx_dge_stageI$counts) %in% TRACERx_ITH_2$sample_name_cruk)]

# put counts in same order as clinical data

TRACERx_ITH_counts <- TRACERx_ITH_counts[,TRACERx_ITH_2$sample_name_cruk]

colnames(TRACERx_ITH_counts) <- TRACERx_ITH_2$sample_name_cruk

data_long <- data_long[order(match(data_long$tumour_id_per_patient, TRACERx_ITH_2$tumour_id_per_patient)),]

data_long$sample_name_cruk <- TRACERx_ITH_2$sample_name_cruk

TRACERx_ITH_counts <- TRACERx_ITH_counts[,data_long$sample_name_cruk]

data_long <- cbind(data_long, t(TRACERx_ITH_counts))

data_long <- data_long %>%
  dplyr::select(-sample_name_cruk)

gene_data_long <- data_long[, -c(1:3)]

gene_correlations <- diag(cor(gene_data_long[data_long$region_number == 1, ], 
                              gene_data_long[data_long$region_number == 2, ],
                              method = 'spearman'))

gene_correlations <- na.omit(gene_correlations)

gene_correlations_ranked <- gene_correlations[order(gene_correlations, decreasing = TRUE)]


VI_predictor_genes_df <- rbind(data.frame(term = 'VI_predictor_genes',
                                          gene = VI_predictor_genes))

set.seed(123)

gsea_region_correlation <- GSEA(gene_correlations_ranked, 
                                TERM2GENE = VI_predictor_genes_df,
                                pvalueCutoff = 1, seed = TRUE)

fig5D <- data.frame('gene' = names(gene_correlations_ranked),
                    'cor' = unname(gene_correlations_ranked))

fig5D <- fig5D %>%
  mutate('Predictor gene' = ifelse(gene %in% VI_predictor_genes, 'Yes','No')) %>%
  arrange(-cor) %>%
  mutate('Rank' = row_number())

# get running enrichment score for custom plot

gsInfo <- function(object, geneSetID) {
  geneList <- object@geneList
  
  if (is.numeric(geneSetID))
    geneSetID <- object@result[geneSetID, "ID"]
  
  geneSet <- object@geneSets[[geneSetID]]
  exponent <- object@params[["exponent"]]
  df <- gseaScores(geneList, geneSet, exponent, fortify=TRUE)
  df$ymin <- 0
  df$ymax <- 0
  pos <- df$position == 1
  h <- diff(range(df$runningScore))/20
  df$ymin[pos] <- -h
  df$ymax[pos] <- h
  df$geneList <- geneList
  
  df$Description <- object@result[geneSetID, "Description"]
  return(df)
}


gseaScores <- getFromNamespace("gseaScores", "DOSE")

fig5D$runningScore <- gsInfo(gsea_region_correlation, geneSetID = c(1))$runningScore

addWorksheet(source_data, "Figure 5D")

writeData(source_data, sheet = "Figure 5D", 
          x = "Figure 5D", 
          startCol =  1, startRow = 1)

writeData(source_data, sheet = 'Figure 5D',
          x = fig5D, startCol = 1, startRow = 3)

pval <- grobTree(textGrob("p =", 
                          x = 0.53, y = 0.95, hjust=0,
                          gp = gpar(fontsize = 17, fontfamily = 'Calibri')),
                 textGrob(eval(bquote(scientific_10(round(gsea_region_correlation$p.adjust, 
                                                          digits = 8)))),
                          x=0.6, y=0.95, hjust=0,
                          gp=gpar(fontsize=17, 
                                  fontfamily = 'Calibri')))

fig5D_gg <- ggplot(fig5D, aes(x=Rank, y=cor, group=1)) + 
  geom_line(linetype = 'solid', linewidth = 2) +
  annotation_custom(pval) +
  geom_line(aes(y = runningScore+1.2), color='black',linewidth=1) +
  geom_bar(stat = "identity", alpha = 1, width=0.5, 
           aes(colour = `Predictor gene`, fill = `Predictor gene`), 
           show.legend = TRUE) +
  scale_fill_manual(values = c('No' = 'transparent', 'Yes' = '#FF0054')) +
  scale_colour_manual(values = c('No' = 'transparent', 'Yes' = '#FF0054')) +
  theme_classic() +
  labs(title="TRACERx stage I LUAD (n=63 tumors)\n", x = "\nGene rank", 
       y = "Region 1 and 2 Correlation     Running ES\n") +
  theme(plot.title = element_text(size=22, hjust = 0.5),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20, hjust = 0.8),
        axis.text.x = element_text(colour="black", size=20),
        axis.text.y = element_text(colour="black", size=20),
        legend.position="right", 
        axis.text=element_text(size=20),
        legend.title = element_text(size=20),
        legend.key.size = unit(0.6, "cm"),
        legend.text = element_text(size=20),
        text = element_text(family = 'Calibri')) +
  scale_x_continuous(breaks = c(4000,8000,12000), labels = c(4000,8000,12000)) +
  scale_y_continuous(breaks = c(0,0.5,1,1.2,1.62), labels = c(0,0.5,1,0,1)) +
  geom_hline(yintercept = 1.2, linetype = "dashed") +
  guides(size = 'none') 

fig5D_gg

plot_dimensions <- dev.size("in")

aspect_ratio <- plot_dimensions[1] / plot_dimensions[2]

aspect_ratio <- 1.213087

desired_height <- 6.5  

desired_width <- desired_height * aspect_ratio

TIFF_dpi <- 600

ggsave(here('figures','fig5D.tiff'),
       plot = fig5D_gg,
       device = 'tiff',
       width = desired_width,
       height = desired_height,
       units = 'in', dpi = TIFF_dpi)


# intra-tumor variation in scores (absolute difference between two regions of same tumor)

###########
# FIGURE 5E
###########

TRACERx_ITH_2 <- TRACERx_ITH_2 %>%
  group_by(tumour_id_per_patient) %>%
  dplyr::summarise(VI_predictor_diff = abs(last(VI_predictor) - dplyr::first(VI_predictor))) %>%
  right_join(TRACERx_ITH_2, by = 'tumour_id_per_patient', multiple = 'all')

# inter-tumor variation in scores (randomly repeated drawing from regions)

set.seed(111)

# Function to randomly sample 1 region from two two different tumors and calculate absolute difference

calculate_inter_diff <- function(df) {
  
  sampled_df <- df %>%
    group_by(tumour_id_per_patient) %>%
    sample_n(1) %>%
    ungroup() %>%
    sample_n(2)
  
  diff_score <- abs(last(sampled_df$VI_predictor) - dplyr::first(sampled_df$VI_predictor))
  
  return(diff_score)
  
}

# Perform the sampling and calculation n times and store in a vector

VI_predictor_inter_diff <- replicate(length(unique(TRACERx_ITH_2$tumour_id_per_patient)), 
                                     calculate_inter_diff(TRACERx_ITH_2))


tmp <- TRACERx_ITH_2 %>%
  distinct(tumour_id_per_patient, .keep_all = TRUE) %>%
  select(VI_predictor_diff)

# plot intra and inter TH together

fig5E <- data.frame('Diff' = c(VI_predictor_inter_diff, 
                               unique(tmp$VI_predictor_diff)),
                    'Tumor_heterogeneity' = c(rep('InterTH',length(VI_predictor_inter_diff)),
                                              rep('IntraTH',length(unique(TRACERx_ITH_2$tumour_id_per_patient)))))

fig5E$Tumor_heterogeneity <- factor(fig5E$Tumor_heterogeneity, levels = c('IntraTH','InterTH'))

addWorksheet(source_data, "Figure 5E")

writeData(source_data, sheet = "Figure 5E", 
          x = "Figure 5E", 
          startCol =  1, startRow = 1)

writeData(source_data, sheet = 'Figure 5E',
          x = fig5E, startCol = 1, startRow = 3)

pval <- fig5E %>% 
  wilcox_test(
    Diff ~ Tumor_heterogeneity
  ) %>% 
  add_xy_position() %>% 
  mutate(p = sub("e","%.% 10^",p) )

fig5E_gg <- ggplot(data = fig5E, aes(x = Tumor_heterogeneity, y = Diff)) +
  theme_classic() +
  geom_violin(lwd=0.5, width = 1, aes(fill = Tumor_heterogeneity)) +
  geom_boxplot(color = "black", fill = 'white', width = .15, outlier.alpha = 0, lwd=0.5) +
  geom_jitter(colour="black",fill='white', size=1, alpha=1, width = 0.2, shape=21) +
  scale_fill_manual(values = c("#0B3954","#BFD7EA")) +
  labs(title="TRACERx stage I LUAD (n=63 tumors)\n", x = NULL, y = "Absolute difference in \n predictor scores, region 1 vs region 2\n") +
  theme(plot.title = element_text(size=22, hjust = 0.5),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        legend.position="none", 
        axis.text=element_text(size=12),
        text = element_text(family = 'Calibri')) +
  ylim(0,max(fig5E$Diff)+0.05) +
  ggprism::add_pvalue(pval, step.increase = 0.1,
                      family = 'Calibri', label.size = 6, parse = TRUE)

fig5E_gg


plot_dimensions <- dev.size("in")

aspect_ratio <- plot_dimensions[1] / plot_dimensions[2]

aspect_ratio <- 1.14094

desired_height <- 6.5

desired_width <- desired_height * aspect_ratio

TIFF_dpi <- 600

ggsave(here('figures','fig5E.tiff'),
       plot = fig5E_gg,
       device = 'tiff',
       width = desired_width,
       height = desired_height,
       units = 'in', dpi = TIFF_dpi)


##################
# OUTCOME ANALYSIS
##################

# select mean regional VI biomarker score for each tumor

TRACERx_dge_stageI$samples <- TRACERx_dge_stageI$samples %>%
  group_by(tumour_id_muttable_cruk.x) %>%
  dplyr::mutate(mean_VI_predictor = mean(VI_predictor, na.rm=TRUE)) %>%
  distinct(tumour_id_muttable_cruk.x, .keep_all = TRUE)


dim(TRACERx_dge_stageI$samples) # 82 stage I LUAD patients

# Converting to 5 year OS, DSS, RFS data below 

TRACERx_dge_stageI$samples <- TRACERx_dge_stageI$samples %>%
  mutate(os_5_year = ifelse(cens_os == 1 & os_time > 1825, 0 , cens_os)) %>%
  mutate(time_os_5_year = ((os_time / 365) * 12)) %>%
  mutate(time_os_5_year = ifelse(time_os_5_year >= 60, 60, time_os_5_year)) %>%
  mutate(dfs_5_year = ifelse(cens_dfs == 1 & dfs_time > 1825, 0 , cens_dfs)) %>%
  mutate(time_dfs_5_year = ((dfs_time / 365) * 12)) %>%
  mutate(time_dfs_5_year = ifelse(time_dfs_5_year >= 60, 60, time_dfs_5_year)) %>%
  mutate(rfs = ifelse(Relapse_cat_new == 'No rec',0,1)) %>%
  mutate(rfs_5_year = ifelse(rfs == 1 & lung_event_time > 1825, 0 , rfs)) %>%
  mutate(time_rfs_5_year = ((lung_event_time / 365) * 12)) %>%
  mutate(time_rfs_5_year = ifelse(time_rfs_5_year >= 60, 60, time_rfs_5_year))


# convert variables of interest to factors

TRACERx_dge_stageI$samples <- TRACERx_dge_stageI$samples %>%
  mutate(ctDNA = dplyr::coalesce(Pre.op.ctD..Tx100.assay., Pre.op.ctD..Tx421.assay.)) %>%
  mutate_at(vars(c('ctDNA','STAS','vascular_invasion_per_lesion','pleural_invasion_per_lesion','Necrosis')), factor) %>%
  mutate(ctDNA = fct_recode(ctDNA,'0' = 'neg','1' = 'pos')) %>%
  mutate(STAS = fct_recode(STAS,'0' = 'absent','1' = 'present')) %>%
  mutate(Necrosis = fct_recode(Necrosis,'0' = 'absent','1' = 'present')) %>%
  mutate(vascular_invasion_per_lesion = fct_recode(vascular_invasion_per_lesion,'0' = 'No','1' = 'Yes')) %>%
  mutate(pleural_invasion_per_lesion = fct_recode(pleural_invasion_per_lesion,'0' = 'No','1' = 'Yes'))



levels(TRACERx_dge_stageI$samples$Relapse_cat_new) <- c('No rec','Intrathoracic','Extrathoracic(+/- Intra)','Extrathoracic(+/- Intra)','Unknown Site')

TRACERx_dge_stageI$samples <- TRACERx_dge_stageI$samples %>%
  mutate(LVI_per_patient = factor(LVI_per_patient, levels = c('No','Yes')))

levels(TRACERx_dge_stageI$samples$LVI_per_patient) <- c('LVI-','LVI+')


# association with LVI status

######################
# EXTENDED DATA FIG 7E
######################

efig7E <- TRACERx_dge_stageI$samples %>%
  select(LVI_per_patient, VI_predictor)

addWorksheet(source_data, "E. Figure 7E")

writeData(source_data, sheet = "E. Figure 7E", 
          x = "E. Figure 7E", 
          startCol =  1, startRow = 1)

writeData(source_data, sheet = 'E. Figure 7E',
          x = efig7E, startCol = 1, startRow = 3)

efig7E_gg <- ggplot(data = efig7E, aes(x = LVI_per_patient, y = VI_predictor, 
                                       fill = LVI_per_patient)) +
  theme_classic() +
  geom_violin(lwd=0.5, width = 0.5) +
  geom_boxplot(color = "black", fill = 'white', width = .15, outlier.alpha = 0, lwd=0.5) +
  geom_jitter(colour="black",fill='white', size=1, alpha=1, width = 0.2, shape=21) +
  scale_fill_manual(values = c('#999999','#8D2965')) +
  labs(title="TRACERx stage I LUAD (n=82 tumors)\n", 
       x = NULL, y = "Mean VI predictor score across regions\n") +
  theme(plot.title = element_text(size=22, hjust = 0.5),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        legend.position="none", 
        axis.text=element_text(size=20),
        text = element_text(family = 'Calibri')) +
  stat_compare_means(aes(group=LVI_per_patient, label = 'p.format'),
                     family = 'Calibri',
                     comparisons = NULL,
                     method = "wilcox.test",
                     label.x = 1.3,
                     label.y = max(efig7E$VI_predictor)+0.05,
                     size = 6,
                     step.increase = 0.05, tip.length = 0.01)

efig7E_gg

plot_dimensions <- dev.size("in")

aspect_ratio <- plot_dimensions[1] / plot_dimensions[2]

aspect_ratio <- 1.07047

desired_height <- 6.5  

desired_width <- desired_height * aspect_ratio

TIFF_dpi <- 600

# Save the plot with adjusted dimensions
ggsave(here('figures','efig7E.tiff'),
       plot = efig7E_gg,
       device = 'tiff',
       width = desired_width,
       height = desired_height,
       units = 'in', dpi = TIFF_dpi)


# association with ctDNA status

######################
# EXTENDED DATA FIG 7F
######################

efig7F <- TRACERx_dge_stageI$samples %>% filter(!is.na(ctDNA)) %>%
  select(ctDNA, VI_predictor)

levels(efig7F$ctDNA) <- c('ctDNA-','ctDNA+')

addWorksheet(source_data, "E. Figure 7F")

writeData(source_data, sheet = "E. Figure 7F", 
          x = "E. Figure 7F", 
          startCol =  1, startRow = 1)

writeData(source_data, sheet = 'E. Figure 7F',
          x = efig7F, startCol = 1, startRow = 3)

efig7F_gg <- ggplot(data = efig7F, aes(x = ctDNA, y = VI_predictor, fill = ctDNA)) +
  theme_classic() +
  geom_violin(lwd=0.5, width = 0.5) +
  geom_boxplot(color = "black", fill = 'white', width = .15, outlier.alpha = 0, lwd=0.5) +
  geom_jitter(colour="black",fill='white', size=1, alpha=1, width = 0.2, shape=21) +
  scale_fill_manual(values = c('#999999','#880808')) +
  labs(title="TRACERx stage I LUAD (n=51 tumors)\n", x = NULL, y = "Mean VI predictor score across regions\n") +
  theme(plot.title = element_text(size=22, hjust = 0.5),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=16),
        legend.position="none", 
        axis.text=element_text(size=12),
        text = element_text(family = 'Calibri')) +
  stat_compare_means(aes(group=ctDNA, label = 'p.format'),
                     family = 'Calibri',
                     comparisons = NULL,
                     method = "wilcox.test",
                     label.x = 1.3,
                     label.y = max(efig7F$VI_predictor)+0.05,
                     size = 6,
                     step.increase = 0.05, tip.length = 0.01)

efig7F_gg

ggsave(here('figures','efig7F.tiff'),
       plot = efig7F_gg,
       device = 'tiff',
       width = desired_width,
       height = desired_height,
       units = 'in', dpi = TIFF_dpi)



# VI predictor is independent of TNM stage I and II in TRACERx

# here we repeat the VI predictor score generation on the stage I/IIs for reproducbility of previously reported results above

######################
# EXTENDED DATA FIG 7G
######################

# Filter TRACERx stage I/II counts to filtered discovery cohort counts

TRACERx_dge$counts <- TRACERx_dge$counts[which(rownames(TRACERx_dge$counts) %in% rownames(discovery_dge$counts)),]

dim(TRACERx_dge$counts)

discovery_dge$counts <- discovery_dge$counts[which(rownames(discovery_dge$counts) %in% rownames(TRACERx_dge$counts)),]

dim(discovery_dge$counts)


# put genes in same order as train

TRACERx_dge$counts <- TRACERx_dge$counts[rownames(discovery_dge$counts),]

TRACERx_dge <- DGEList(counts = TRACERx_dge$counts, samples = TRACERx_dge$samples)

head(TRACERx_dge$samples$lib.size)

TRACERx_dge <- calcNormFactors(TRACERx_dge) # Trimmed mean of m values (TMM) normalization


# Reference combat TRACERx count matrix with discovery cohort count matrix


TRACERx_dge$counts <- ComBat(cbind(edgeR::cpm(discovery_dge$counts,log=T), edgeR::cpm(TRACERx_dge$counts,log=T)),
                             c(rep(1,ncol(discovery_dge$counts)),rep(2,ncol(TRACERx_dge$counts))),
                             ref.batch = 1)

TRACERx_dge$counts <- TRACERx_dge$counts[,(ncol(discovery_dge$counts)+1):ncol(TRACERx_dge$counts)]

length(which(VI_predictor_genes %in% rownames(TRACERx_dge$counts))) # 48/48 genes present

# use final VI predictor model to predict

test_TRACERx <- data.frame(t(TRACERx_dge$counts[which(rownames(TRACERx_dge$counts) %in% VI_predictor_genes),]))

# convert to h2o format and add dummy class variable

test_TRACERx$response <- rep(c(0,1),length.out = ncol(TRACERx_dge$counts))

test_TRACERx <- as.h2o(test_TRACERx)

# Identify predictors and response

y <- "response"

x <- setdiff(names(test_TRACERx), y)

pred_test <- h2o.predict(aml_final_leader, test_TRACERx)  

pred_test <- data.frame(as.matrix(pred_test))

TRACERx_dge$samples$VI_predictor <- as.numeric(pred_test$p1)


TRACERx_dge$samples <- TRACERx_dge$samples %>%
  group_by(tumour_id_muttable_cruk.x) %>%
  dplyr::mutate(mean_VI_predictor = mean(VI_predictor, na.rm=TRUE)) %>%
  distinct(tumour_id_muttable_cruk.x, .keep_all = TRUE)


dim(TRACERx_dge$samples) # 137 stage I/II LUAD patients

# Converting to 5 year OS, DSS, RFS data below 

TRACERx_dge$samples <- TRACERx_dge$samples %>%
  mutate(os_5_year = ifelse(cens_os == 1 & os_time > 1825, 0 , cens_os)) %>%
  mutate(time_os_5_year = ((os_time / 365) * 12)) %>%
  mutate(time_os_5_year = ifelse(time_os_5_year >= 60, 60, time_os_5_year)) %>%
  mutate(dfs_5_year = ifelse(cens_dfs == 1 & dfs_time > 1825, 0 , cens_dfs)) %>%
  mutate(time_dfs_5_year = ((dfs_time / 365) * 12)) %>%
  mutate(time_dfs_5_year = ifelse(time_dfs_5_year >= 60, 60, time_dfs_5_year)) %>%
  mutate(rfs = ifelse(Relapse_cat_new == 'No rec',0,1)) %>%
  mutate(rfs_5_year = ifelse(rfs == 1 & lung_event_time > 1825, 0 , rfs)) %>%
  mutate(time_rfs_5_year = ((lung_event_time / 365) * 12)) %>%
  mutate(time_rfs_5_year = ifelse(time_rfs_5_year >= 60, 60, time_rfs_5_year))


efig7G <- TRACERx_dge$samples %>% 
  mutate(Stage = factor(pStage..TNM.v7.)) %>%
  select(time_rfs_5_year, rfs_5_year, Stage,
         VI_predictor)

efig7G$Stage <- factor(efig7G$Stage, levels = c('IA','IB','IIA','IIB'))

addWorksheet(source_data, "E. Figure 7G")

writeData(source_data, sheet = "E. Figure 7G", 
          x = "E. Figure 7G", 
          startCol =  1, startRow = 1)

writeData(source_data, sheet = 'E. Figure 7G',
          x = efig7G, startCol = 1, startRow = 3)

TRACERx_cox_multi <- coxph(formula = Surv(time_rfs_5_year, rfs_5_year) ~ 
                             scale(VI_predictor) + Stage,
                           data = efig7G)

summary(TRACERx_cox_multi)

efig7G <- data.frame(
  HR       = summary(TRACERx_cox_multi)$coefficients[, "exp(coef)"],
  lower_ci = summary(TRACERx_cox_multi)$conf.int[, "lower .95"],
  upper_ci = summary(TRACERx_cox_multi)$conf.int[, "upper .95"],
  P        = as.character(summary(TRACERx_cox_multi)$coefficients[, "Pr(>|z|)"])
) %>% 
  rownames_to_column("Predictor") %>% 
  mutate(Predictor = case_when(
    Predictor == "scale(VI_predictor)" ~ "VI predictor",
    Predictor == "StageIB"     ~ "Stage IB",
    Predictor == "StageIIA"     ~ "Stage IIA",
    Predictor == "StageIIB"     ~ "Stage IIB",
    TRUE ~ Predictor
  ))

efig7G <- efig7G %>%
  mutate(Predictor = factor(Predictor, 
                            levels = c("VI predictor", "Stage IB", "Stage IIA",
                                       "Stage IIB"))) %>%
  arrange(Predictor)

efig7G <- efig7G %>%
  mutate(across(c(HR, lower_ci, upper_ci), as.numeric)) %>%
  mutate(across(c(HR, lower_ci, upper_ci), ~ signif(., digits = 3))) %>%
  mutate(P = ifelse(P != "Reference", as.character(signif(as.numeric(P), digits = 3)), P)) #%>%

hr_max <- ceiling(max(efig7G$upper_ci, na.rm = TRUE) / 10) * 10

efig7G_gg <- ggplot(efig7G, aes(x = Predictor, y = HR)) +
  coord_flip() +
  geom_hline(yintercept = 1, linetype = "dotted") +
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), width = 0.5, color = "#FF0054") +
  geom_point(pch = 15, size = 3) +
  theme_minimal() +
  theme(
    legend.position   = "none",
    plot.title        = element_text(size = 20, hjust=0.5),
    aspect.ratio      = 0.8,
    panel.grid        = element_blank(),
    axis.line.x       = element_line(),
    axis.ticks.x      = element_line(),
    axis.text.x       = element_text(size = 15),
    axis.text.y       = element_text(size = 15),
    axis.title.x      = element_text(size = 20),
    axis.title.y      = element_text(size = 20),
    text              = element_text(family = 'Calibri')
  ) +
  scale_y_continuous(trans = "log", expand = c(0,-1), breaks = c(0.05,2,8,32), 
                     limits = c(0.05, 1200)) + 
  scale_x_discrete(limits=c(rev(as.character(efig7G$Predictor)),' ')) + 
  xlab("Predictor\n") + ylab("\nHazard Ratio (95% CI)") + 
  geom_text(aes(x=Predictor, y=hr_max + 100, label=P), 
            size=5, family = 'Calibri', parse = TRUE) + 
  ggtitle("TRACERx stage I/II LUAD (n=137) 5-Year RFS\n") + 
  annotate("text", x = 5, y = 110, label = "p value", 
           size = 5, fontface = 'bold', family = 'Calibri')

efig7G_gg

plot_dimensions <- dev.size("in")

aspect_ratio <- plot_dimensions[1] / plot_dimensions[2]

aspect_ratio <- 1.821809

desired_height <- 4  

desired_width <- desired_height * aspect_ratio

TIFF_dpi <- 600

ggsave(here('figures','efig7G.tiff'),
       plot = efig7G_gg,
       device = 'tiff',
       width = desired_width,
       height = desired_height,
       units = 'in', dpi = TIFF_dpi)




# VI predictor is independent of surgery type and adjuvant therapy

######################
# EXTENDED DATA FIG 7H
######################

efig7H <- TRACERx_dge_stageI$samples %>% 
  mutate(Surgery = factor(Surgery.type)) %>%
  mutate(Adjuvant = factor(Adjuvant.treatment)) %>%
  select(time_rfs_5_year, rfs_5_year, Surgery, Adjuvant,
         VI_predictor)

levels(efig7H$Surgery) <- c('Lobectomy','Sublobar')

addWorksheet(source_data, "E. Figure 7H")

writeData(source_data, sheet = "E. Figure 7H", 
          x = "E. Figure 7H", 
          startCol =  1, startRow = 1)

writeData(source_data, sheet = 'E. Figure 7H',
          x = efig7H, startCol = 1, startRow = 3)

TRACERx_cox_multi <- coxph(formula = Surv(time_rfs_5_year, rfs_5_year) ~ 
                             scale(VI_predictor) + Surgery + Adjuvant,
                           data = efig7H)

summary(TRACERx_cox_multi)

efig7H <- data.frame(
  HR       = summary(TRACERx_cox_multi)$coefficients[, "exp(coef)"],
  lower_ci = summary(TRACERx_cox_multi)$conf.int[, "lower .95"],
  upper_ci = summary(TRACERx_cox_multi)$conf.int[, "upper .95"],
  P        = as.character(summary(TRACERx_cox_multi)$coefficients[, "Pr(>|z|)"])
) %>% 
  rownames_to_column("Predictor") %>% 
  mutate(Predictor = case_when(
    Predictor == "scale(VI_predictor)" ~ "VI predictor",
    Predictor == "SurgerySublobar"     ~ "Sublobar",
    Predictor == "AdjuvantAdjuvant"     ~ "Adjuvant",
    TRUE ~ Predictor
  ))

efig7H <- efig7H %>%
  bind_rows(tibble(Predictor = c("Lobectomy", "No adjuvant"),
                   HR       = 1,
                   lower_ci = NA,
                   upper_ci = NA,
                   P        = "Reference"))

efig7H <- efig7H %>%
  mutate(Predictor = factor(Predictor, 
                            levels = c("VI predictor", "Sublobar", "Lobectomy",
                                       "Adjuvant", "No adjuvant"))) %>%
  arrange(Predictor)

efig7H <- efig7H %>%
  mutate(across(c(HR, lower_ci, upper_ci), as.numeric)) %>%
  mutate(across(c(HR, lower_ci, upper_ci), ~ signif(., digits = 3))) %>%
  mutate(P = ifelse(P != "Reference", as.character(signif(as.numeric(P), digits = 3)), P)) #%>%

hr_max <- ceiling(max(efig7H$upper_ci, na.rm = TRUE) / 10) * 10

efig7H_gg <- ggplot(efig7H, aes(x = Predictor, y = HR)) +
  coord_flip() +
  geom_hline(yintercept = 1, linetype = "dotted") +
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), width = 0.5, color = "#FF0054") +
  geom_point(pch = 15, size = 3) +
  theme_minimal() +
  theme(
    legend.position   = "none",
    plot.title        = element_text(size = 20, hjust=0.5),
    aspect.ratio      = 0.8,
    panel.grid        = element_blank(),
    axis.line.x       = element_line(),
    axis.ticks.x      = element_line(),
    axis.text.x       = element_text(size = 15),
    axis.text.y       = element_text(size = 15),
    axis.title.x      = element_text(size = 20),
    axis.title.y      = element_text(size = 20),
    text              = element_text(family = 'Calibri')
  ) +
  scale_y_continuous(trans = "log", expand = c(0,-1), breaks = c(0.05,2,8,32), 
                     limits = c(0.05, 1200)) + 
  scale_x_discrete(limits=c(rev(as.character(efig7H$Predictor)),' ')) + 
  xlab("Predictor\n") + ylab("\nHazard Ratio (95% CI)") + 
  geom_text(aes(x=Predictor, y=hr_max + 100, label=P), 
            size=5, family = 'Calibri', parse = TRUE) + 
  ggtitle("TRACERx stage I LUAD (n=82) 5-Year RFS\n") + 
  annotate("text", x = 6, y = 110, label = "p value", 
           size = 5, fontface = 'bold', family = 'Calibri')

efig7H_gg

plot_dimensions <- dev.size("in")

aspect_ratio <- plot_dimensions[1] / plot_dimensions[2]

aspect_ratio <- 1.821809

desired_height <- 4  

desired_width <- desired_height * aspect_ratio

TIFF_dpi <- 600

ggsave(here('figures','efig7H.tiff'),
       plot = efig7H_gg,
       device = 'tiff',
       width = desired_width,
       height = desired_height,
       units = 'in', dpi = TIFF_dpi)





fig4F_5 <- TRACERx_dge_stageI$samples %>%
  select(time_rfs_5_year, rfs_5_year, VI_predictor)

addWorksheet(source_data, "Figure 4F-5")

writeData(source_data, sheet = "Figure 4F-5", 
          x = "Figure 4F-5", 
          startCol =  1, startRow = 1)

writeData(source_data, sheet = 'Figure 4F-5',
          x = fig4F_5, startCol = 1, startRow = 3)

# independence from TNM stage I (not significant)

TRACERx_cox <- coxph(formula = Surv(time_rfs_5_year, rfs_5_year) ~ scale(VI_predictor) + pStage..TNM.v7.,
                     data = TRACERx_dge_stageI$samples)

summary(TRACERx_cox)

# univariate coxph analyis of VI predictor association with outcome

TRACERx_cox <- coxph(formula = Surv(time_rfs_5_year, rfs_5_year) ~ scale(VI_predictor),
                     data = TRACERx_dge_stageI$samples)

names(TRACERx_cox$coefficients) <- 'TRACERx (n=82)'

summary(TRACERx_cox)

h2o.removeAll(timeout_secs = 0, retained_elements = c())






################
# LOAD TCGA DATA
################

# Query platform Illumina HiSeq with a list of barcode

query <- GDCquery(project = "TCGA-LUAD", 
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  experimental.strategy = "RNA-Seq",
                  legacy = FALSE)

query$results[[1]] <- query$results[[1]] %>% 
  dplyr::filter(sample_type == "Primary Tumor")

# Download a list of barcodes with platform IlluminaHiSeq_RNASeqV2

GDCdownload(query)

# Prepare expression matrix with geneID in the rows and samples in the columns
# rsem.genes.results as values

TCGA_LUAD_RnaseqSE <- GDCprepare(query)

TCGA_LUAD_RNAseq_counts <- assay(TCGA_LUAD_RnaseqSE,"unstranded")

TCGA_LUAD_RNAseq_pheno <- as.data.frame(colData(TCGA_LUAD_RnaseqSE))

table(TCGA_LUAD_RNAseq_pheno$ajcc_pathologic_stage)

# subset to stage IA

TCGA_LUAD_RNAseq_pheno$ajcc_pathologic_stage <- as.factor(TCGA_LUAD_RNAseq_pheno$ajcc_pathologic_stage)

TCGA_LUAD_RNAseq_pheno <- TCGA_LUAD_RNAseq_pheno %>% filter(ajcc_pathologic_stage == "Stage IA")

TCGA_LUAD_RNAseq_pheno$ajcc_pathologic_stage <- droplevels(TCGA_LUAD_RNAseq_pheno$ajcc_pathologic_stage)

dim(TCGA_LUAD_RNAseq_pheno) # 137 Stage IA LUAD samples

TCGA_LUAD_RNAseq_counts <- TCGA_LUAD_RNAseq_counts[,which(colnames(TCGA_LUAD_RNAseq_counts) %in% rownames(TCGA_LUAD_RNAseq_pheno))]


# create DGE list object

TCGA_LUAD_dge <- DGEList(counts = TCGA_LUAD_RNAseq_counts, samples = TCGA_LUAD_RNAseq_pheno)

# load saved dge for exact reproducibility (needed survival data missing otherwise)

TCGA_LUAD_dge <- readRDS(here('data','TCGA_LUAD_stage_IA_dge'))

# convert gene names to hgnc

TCGA_LUAD_dge <- ensemble_to_hgnc(TCGA_LUAD_dge)


################
# GENE FILTERING
################

# Filter TCGA counts to filtered discovery cohort counts

TCGA_LUAD_dge$counts <- TCGA_LUAD_dge$counts[which(rownames(TCGA_LUAD_dge$counts) %in% rownames(discovery_dge$counts)),]

dim(TCGA_LUAD_dge$counts)


discovery_dge$counts <- discovery_dge$counts[which(rownames(discovery_dge$counts) %in% rownames(TCGA_LUAD_dge$counts)),]

dim(discovery_dge$counts)


# put genes in same order as train

TCGA_LUAD_dge$counts <- TCGA_LUAD_dge$counts[rownames(discovery_dge$counts),]

TCGA_LUAD_dge <- DGEList(counts = TCGA_LUAD_dge$counts, samples = TCGA_LUAD_dge$samples)

head(TCGA_LUAD_dge$samples$lib.size)

# Trimmed mean of m values (TMM) normalization

TCGA_LUAD_dge <- calcNormFactors(TCGA_LUAD_dge)




#######################
# CORRECT BATCH EFFECTS
#######################

# Reference combat TCGA count matrix with discovery cohort count matrix


TCGA_LUAD_dge$counts <- ComBat(cbind(cpm(discovery_dge$counts,log=T), cpm(TCGA_LUAD_dge$counts,log=T)),
                               c(rep(1,ncol(discovery_dge$counts)),rep(2,ncol(TCGA_LUAD_dge$counts))), 
                               ref.batch = 1)

TCGA_LUAD_dge$counts <- TCGA_LUAD_dge$counts[,(ncol(discovery_dge$counts)+1):ncol(TCGA_LUAD_dge$counts)]

length(which(VI_predictor_genes %in% rownames(TCGA_LUAD_dge$counts))) # 43/48 all genes present



############################
# GENERATE MODEL PREDICTIONS
############################

# use final VI predictor model to predict on TCGA LUAD data

test_TCGA <- data.frame(t(TCGA_LUAD_dge$counts[which(rownames(TCGA_LUAD_dge$counts) %in% VI_predictor_genes),]))

# convert to h2o format and add dummy class variable

test_TCGA$response <- rep(c(1,0),length.out = ncol(TCGA_LUAD_dge$counts))

test_TCGA <- as.h2o(test_TCGA)

# Identify predictors and response

y <- "response"

x <- setdiff(names(test_TCGA), y)

# For binary classification, response should be a factor

aml_final_leader <- h2o.loadModel(here('data','GLM_1_AutoML_1_20230802_104952'))

pred_test <- h2o.predict(aml_final_leader, test_TCGA)  

pred_test <- data.frame(as.matrix(pred_test))


TCGA_LUAD_dge$samples$VI_predictor <- as.numeric(pred_test$p1)



###################
# OUTCOME ANALYSIS
###################

# Filter to patients with outcome data
# Determine alive/Dead at last follow up and maximum follow up time to 5 years

fig4F_3 <- TCGA_LUAD_dge$samples %>%
  filter(!is.na(days_to_last_follow_up)) %>%
  mutate(
    months_to_last_follow_up = (days_to_last_follow_up / 365) * 12,
    vital_status_bin = ifelse(vital_status == "Dead", 1, 0), 
    follow_up_5yr = pmin(months_to_last_follow_up, 60)) %>%
  select(follow_up_5yr, vital_status_bin, VI_predictor)

addWorksheet(source_data, "Figure 4F-3")

writeData(source_data, sheet = "Figure 4F-3", 
          x = "Figure 4F-3", 
          startCol =  1, startRow = 1)

writeData(source_data, sheet = 'Figure 4F-3',
          x = fig4F_3, startCol = 1, startRow = 3,
          rowNames = TRUE)

# univariate coxph analyis of VI predictor association with outcome

TCGA_cox <- coxph(formula = Surv(follow_up_5yr, vital_status_bin) ~ 
                    scale(VI_predictor),
                  data=fig4F_3)

summary(TCGA_cox)

names(TCGA_cox$coefficients) <- 'TCGA (n=127)'

h2o.removeAll(timeout_secs = 0, retained_elements = c())

###################
# LOAD UPPSALA DATA
###################

GSE81089 <- getGEO("GSE81089", GSEMatrix = TRUE)

getGEOSuppFiles("GSE81089")

uppsala_counts <- read.table(gzfile(here('GSE81089','GSE81089_readcounts_featurecounts.tsv.gz')), 
                           header = TRUE, sep = "\t", row.names = 1)

uppsala_pheno <- pData(GSE81089[[1]])

rownames(uppsala_pheno) <- uppsala_pheno$title

# remove matched normal samples from pheno data

uppsala_pheno <- uppsala_pheno %>%
  filter(!grepl("^matched.sample_", title))

# clean up the count data column names

colnames(uppsala_counts) <- sapply(strsplit(colnames(uppsala_counts), "_", fixed = TRUE), `[`, 1)

uppsala_counts <- uppsala_counts %>% select(all_of(rownames(uppsala_pheno)))

# remove unwanted TC% row from count data

uppsala_counts <- uppsala_counts %>%
  filter(rownames(.) != "TC%")

# clean up phenotype data columns

colnames(uppsala_pheno) <- gsub(":ch1", "", colnames(uppsala_pheno))

uppsala_pheno <- uppsala_pheno %>%
  mutate(
    smoking = as.factor(smoking),
    histology = as.factor(histology),
    `stage tnm` = as.factor(`stage tnm`),
    age = as.numeric(age),
    gender = as.factor(gender),
    `ps who` = as.factor(`ps who`)
  )

levels(uppsala_pheno$smoking) <- c("Current", "Former", "Never")

levels(uppsala_pheno$histology) <- c("LUSC", "LUAD", "LC/NOS")

levels(uppsala_pheno$`stage tnm`) <- c("IA", "IB", "IIA", "IIB", "IIIA", "IIIB", "IV")

# Subset to LUAD samples only

uppsala_pheno <- uppsala_pheno %>% 
  dplyr::filter(histology == "LUAD") %>% 
  dplyr::filter(`stage tnm` == "IA")

uppsala_counts <- uppsala_counts[,which(colnames(uppsala_counts) %in% rownames(uppsala_pheno))]


# create DGE list object

uppsala_dge <- DGEList(counts = uppsala_counts, samples = uppsala_pheno) # create DGE list

# convert gene names to hgnc

uppsala_dge <- ensemble_to_hgnc(uppsala_dge)



################
# GENE FILTERING
################

# Filter uppsala counts to filtered discovery set counts

uppsala_dge$counts <- uppsala_dge$counts[which(rownames(uppsala_dge$counts) %in% rownames(discovery_dge$counts)),]

dim(uppsala_dge$counts)

discovery_dge$counts <- discovery_dge$counts[which(rownames(discovery_dge$counts) %in% rownames(uppsala_dge$counts)),]

dim(discovery_dge$counts)


# put genes in same order as train

uppsala_dge$counts <- uppsala_dge$counts[rownames(discovery_dge$counts),]


uppsala_dge <- DGEList(counts = uppsala_dge$counts, samples = uppsala_dge$samples)

head(uppsala_dge$samples$lib.size)

uppsala_dge <- calcNormFactors(uppsala_dge) # Trimmed mean of m values (TMM) normalization


#######################
# CORRECT BATCH EFFECTS
#######################

# Reference combat uppsala count matrix with discovery count matrix

uppsala_dge$counts <- ComBat(cbind(cpm(discovery_dge$counts,log=T), cpm(uppsala_dge$counts,log=T)),
                             c(rep(1,ncol(discovery_dge$counts)),rep(2,ncol(uppsala_dge$counts))),
                             ref.batch = 1)

uppsala_dge$counts <- uppsala_dge$counts[,(ncol(discovery_dge$counts)+1):ncol(uppsala_dge$counts)]

length(which(VI_predictor_genes %in% rownames(uppsala_dge$counts))) # 48/48 genes present


############################
# GENERATE MODEL PREDICTIONS
############################

# use final VI predictor model to predict on uppsala LUAD data

test_uppsala <- data.frame(t(uppsala_dge$counts[which(rownames(uppsala_dge$counts) %in% VI_predictor_genes),]))

# convert to h2o format and add dummy class variable

test_uppsala$response <- rep(c(1,0),length.out = ncol(uppsala_dge$counts))

test_uppsala <- as.h2o(test_uppsala)

# Identify predictors and response

y <- "response"

x <- setdiff(names(test_uppsala), y)

aml_final_leader <- h2o.loadModel(here('data','GLM_1_AutoML_1_20230802_104952'))

pred_test <- h2o.predict(aml_final_leader, test_uppsala)  

pred_test <- data.frame(as.matrix(pred_test))

uppsala_dge$samples$VI_predictor <- as.numeric(pred_test$p1)


###################
# OUTCOME ANALYSIS
###################

# Calculate months to last contact or death andd maximum follow up time to 5 years

fig4F_4 <- uppsala_dge$samples %>%
  mutate(
    months_to_last_contact_or_death = as.vector(difftime(vital.date, surgery.date, units = "days") / 365 * 12),
    dead = as.numeric(dead),
    time_by_5yr = pmin(months_to_last_contact_or_death, 60)) %>%
  select(time_by_5yr, dead, VI_predictor)

addWorksheet(source_data, "Figure 4F-4")

writeData(source_data, sheet = "Figure 4F-4", 
          x = "Figure 4F-4", 
          startCol =  1, startRow = 1)

writeData(source_data, sheet = 'Figure 4F-4',
          x = fig4F_4, startCol = 1, startRow = 3,
          rowNames = TRUE)

# univariate coxph analyis of VI predictor association with outcome

uppsala_cox <- coxph(formula = Surv(time_by_5yr, dead) ~ scale(VI_predictor),
                     data = fig4F_4)

summary(uppsala_cox)

names(uppsala_cox$coefficients) <- 'Uppsala (n=44)'

h2o.removeAll(timeout_secs = 0, retained_elements = c())




# plotting all cox model results from all datasets together

###########
# FIGURE 4F
###########

fig4F <- list('Uppsala' = uppsala_cox,
              'TCGA' = TCGA_cox,
              'TRACERx' = TRACERx_cox,
              'Validation, VI-' = validation_cox_VI_neg,
              'Validation' = validation_cox,
              'Discovery' = discovery_cox)


pval <- grobTree(textGrob(round(summary(uppsala_cox)$coefficients[,5],3),
                          x=0.85,  y=0.1, 
                          gp=gpar(fontsize=17, 
                                  fontfamily = 'Calibri')),
                 textGrob(round(summary(TCGA_cox)$coefficients[,5],3),
                          x=0.85,  y=0.26, 
                          gp=gpar(fontsize=17, 
                                  fontfamily = 'Calibri')),
                 textGrob(round(summary(TRACERx_cox)$coefficients[,5],3),
                          x=0.85,  y=0.42, 
                          gp=gpar(fontsize=17, 
                                  fontfamily = 'Calibri')),
                 textGrob(round(summary(validation_cox_VI_neg)$coefficients[,5],3),
                          x=0.85,  y=0.58, 
                          gp=gpar(fontsize=17, 
                                  fontfamily = 'Calibri')),
                 textGrob(round(summary(validation_cox)$coefficients[,5],3),
                          x=0.85,  y=0.74, 
                          gp=gpar(fontsize=17, 
                                  fontfamily = 'Calibri')),
                 textGrob(eval(bquote(scientific_10(round(summary(discovery_cox)$coefficients[,5],7)))),
                          x=0.85,  y=0.91, 
                          gp=gpar(fontsize=17, 
                                  fontfamily = 'Calibri')))

fig4F_gg <- modelplot(fig4F, exponentiate = TRUE,
                      size = 1.5, linewidth = 1.5,
                      conf_level = 0.95) + 
  theme_classic() +
  ggtitle("VI predictor association with stage I LUAD outcome\n") +
  scale_color_manual(values = c('black','black','black','black','black','black')) +
  scale_shape_manual(values = c(15)) +
  aes(shape = 'test') +
  labs(x = '\nHazard Ratio (95% CI)\n',
       y = 'Patient cohort\n',
       color = 'Dataset') +
  annotation_custom(pval) +
  theme(plot.title = element_text(size=22, hjust = 0.5),
        aspect.ratio = 1,
        axis.title.x = element_text(size=22),
        axis.title.y = element_text(size=22),
        axis.text.x = element_text(size=22),
        axis.text.y = element_text(size=22),
        legend.title = element_text(size=22),
        legend.key.size = unit(0.6, "cm"),
        legend.text = element_text(size=22),
        legend.position = 'none',
        text = element_text(family = 'Calibri'),
        plot.margin = margin(1,1,1.5,1.5, "cm")) +
  scale_x_continuous(breaks=c(0,1,2,3,4,5,6),
                     limits=c(0,9)) +
  scale_y_discrete() +
  geom_vline(xintercept = 1, color = 'grey', linetype = 'dashed')


fig4F_gg

plot_dimensions <- dev.size("in")

aspect_ratio <- plot_dimensions[1] / plot_dimensions[2]

aspect_ratio <- 1.40866

desired_height <- 7 

desired_width <- desired_height * aspect_ratio

TIFF_dpi <- 600

ggsave(here('figures','fig4F.tiff'),
       plot = fig4F_gg,
       device = 'tiff',
       width = desired_width,
       height = desired_height,
       units = 'in', dpi = TIFF_dpi)




# association of VI predictor with pathology features in VI- samples

fig4G <- validation_dge$samples %>% 
  filter(LMPVI != 'VI (G3)') %>%
  select(Neutrophils,Plasma,Invasive_size,
         Total_size,Lepidic,Acinar,Papillary,
         Solid,Micropapillary,Cribriform,Mitosis_per_10_hpf,
         LI,Necrosis,STAS,VI_predictor) %>%
  rename(`Total size` = Total_size,
         Mitosis = Mitosis_per_10_hpf,
         `Invasive size` = Invasive_size)

addWorksheet(source_data, "Figure 4G")

writeData(source_data, sheet = "Figure 4G", 
          x = "Figure 4G", 
          startCol =  1, startRow = 1)

writeData(source_data, sheet = "Figure 4G",
          x = fig4G, startCol = 1, startRow = 3,
          rowNames = TRUE)

path_features <- fig4G

fig4G <- list()

for (i in 1:(ncol(path_features)-1)){
  
  m <- lm(path_features$VI_predictor ~ path_features[,i])
  
  fig4G[[i]] <- c(summary(m)$coefficients[2,3], summary(m)$coefficients[2,4])
  
}

names(fig4G) <- names(path_features)[names(path_features) %!in% 'VI_predictor']

fig4G <- data.frame(do.call(rbind, fig4G))

colnames(fig4G) <- c('t_value','FDR')

fig4G$`Pathology features` <- rownames(fig4G)

fig4G$FDR <- -log10(p.adjust(fig4G$FDR, method = 'bonferroni', n = nrow(fig4G)))

fig4G$sig <- ifelse(fig4G$FDR > -log10(0.05),'p.adj < 0.05','p.adj > 0.05')

fig4G <- fig4G[order(fig4G$FDR),]

fig4G <- fig4G %>%
  select(FDR, `Pathology features`, sig)

###########
# FIGURE 4G
###########

fig4G_gg <- ggplot(fig4G, aes(x=FDR, y=reorder(`Pathology features`, FDR), fill = sig)) + 
  geom_dotplot(binaxis='y', stackdir='center') +
  geom_vline(xintercept = -log10(0.05), linetype="dotted") +
  theme_classic() +
  scale_fill_manual(values = c("#FF0054","#999999")) +
  ggtitle("Association with predictor score\n in VI- Tumors (n=42)\n") +
  labs(x = "\n-log10 p.adj", y = "Pathology feature\n") +
  theme(plot.title = element_text(size=22, hjust = 0.5),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        legend.position="right", axis.text=element_text(size=20),
        legend.title = element_blank(),
        legend.key.size = unit(0.6, "cm"),
        legend.text = element_text(size=20),
        text = element_text(family = 'Calibri')) +
  scale_x_continuous(breaks = -1:2, labels = c(1,0,1,2))

fig4G_gg

plot_dimensions <- dev.size("in")

aspect_ratio <- plot_dimensions[1] / plot_dimensions[2]

aspect_ratio <- 1.142202

desired_height <- 7

desired_width <- desired_height * aspect_ratio

TIFF_dpi <- 600

ggsave(here('figures','fig4G.tiff'),
       plot = fig4G_gg,
       device = 'tiff',
       width = desired_width,
       height = desired_height,
       units = 'in', dpi = TIFF_dpi)




# Define a function to calculate AUROCs

########
# FIG 4H
########

perform_roc_analysis <- function(data, filter_condition, 
                                 outcome_variable = "VI", title) {
  
  tmp <- data %>% filter(eval(parse(text = filter_condition)))
  gene_roc <- roc(as.factor(tmp[[outcome_variable]]), 
                  tmp$VI_predictor, levels = c(0, 1), ci = TRUE)
  auc <- round(pROC::auc(gene_roc), 4)
  ci_min <- gene_roc$ci[1]
  ci_max <- gene_roc$ci[3]
  n <- nrow(tmp)
  list(auc = auc, ci_min = ci_min, ci_max = ci_max, n = n)
  
}

# VI predictions by LUAD histologic pattern

conditions_histology <- c("TRUE",
                          "LMPVI != 'LMP (G1)'",
                          "Solid < 1",
                          "Solid > 1",
                          "Micropapillary < 1",
                          "Micropapillary > 1",
                          "Cribriform < 1",
                          "Cribriform > 1",
                          "Acinar < 1",
                          "Acinar > 1",
                          "Papillary < 1",
                          "Papillary > 1",
                          "Lepidic < 1",
                          "Lepidic > 1")

titles_histology <- c("Predicting VI in Validation Set",
                      "Predicting VI in non-LMP subset",
                      "Predicting VI, No Solid",
                      "Predicting VI, Any Solid",
                      "Predicting VI, No Micropapillary",
                      "Predicting VI, Any Micropapillary",
                      "Predicting VI, No Cribriform",
                      "Predicting VI, Any Cribriform",
                      "Predicting VI, No Acinar",
                      "Predicting VI, Any Acinar",
                      "Predicting VI, No Papillary",
                      "Predicting VI, Any Papillary",
                      "Predicting VI, No Lepidic",
                      "Predicting VI, Any Lepidic")

# Perform the analysis for each condition

results_histology <- lapply(seq_along(conditions_histology), 
                            function(i) perform_roc_analysis(validation_dge$samples, 
                                                             conditions_histology[i], 
                                                             "VI", 
                                                             titles_histology[i]))

# Extract the results

aucs_histology <- sapply(results_histology, function(res) res$auc)

ci_mins_histology <- sapply(results_histology, function(res) res$ci_min)

ci_maxs_histology <- sapply(results_histology, function(res) res$ci_max)

ns_histology <- sapply(results_histology, function(res) res$n)

fig4H <- data.frame(AUCs = aucs_histology, ci_min = ci_mins_histology, ci_max = ci_maxs_histology)

fig4H$Histology <- c('All samples (n=42 VI-,17 VI+)', 'non-LMP (n=35 VI-,17 VI+)', 
                     '0% Solid (n=30 VI-,4 VI+)', '1-100% Solid (n=12 VI-,13 VI+)', 
                     '0% Micropapillary (n=29 VI-,12 VI+)', '1-100% Micropapillary (n=13 VI-,5 VI+)', 
                     '0% Cribriform (n=32 VI-,9 VI+)', '1-100% Cribriform (n=10 VI-,8 VI+)', 
                     '0% Acinar (n=7 VI-,6 VI+)', '1-100% Acinar (n=35 VI-,11 VI+)', 
                     '0% Papillary (n=18 VI-,7 VI+)', '1-100% Papillary (n=24 VI-,10 VI+)', 
                     '0% Lepidic (n=13 VI-,8 VI+)', '1-100% Lepidic (n=29 VI-,9 VI+)')


addWorksheet(source_data, "Figure 4H")

writeData(source_data, sheet = "Figure 4H", 
          x = "Figure 4H", 
          startCol =  1, startRow = 1)

writeData(source_data, sheet = 'Figure 4H',
          x = fig4H, startCol = 1, startRow = 3)

fig4H_gg <- ggplot(fig4H, aes(x = reorder(Histology, AUCs), y = AUCs)) + 
  geom_bar(stat = "identity", width = 0.8, alpha = 0.5) +
  coord_flip() +
  geom_abline(slope = 0, intercept = 0.50, linetype = "dashed", alpha = 0.7, color = "#FF0054") + 
  geom_errorbar(aes(x = Histology, ymin = ci_min, ymax = ci_max), width = 0.2, colour = "black", alpha = 0.9, linewidth = 1) +
  theme_classic() +
  labs(title = "VI predictor performance \nby LUAD histopathology\n", 
       x = "Subgroup\n", y = "\nAUROC") +
  theme(plot.title = element_text(size = 22, hjust = 0.5),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        legend.position = "top", 
        axis.text = element_text(size = 20),
        text = element_text(family = 'Calibri'))

fig4H_gg

plot_dimensions <- dev.size("in")

aspect_ratio <- plot_dimensions[1] / plot_dimensions[2]

aspect_ratio <- 1.426606

desired_height <- 7  

desired_width <- desired_height * aspect_ratio

TIFF_dpi <- 600

ggsave(here('figures','fig4H.tiff'),
       plot = fig4H_gg,
       device = 'tiff',
       width = desired_width,
       height = desired_height,
       units = 'in', dpi = TIFF_dpi)





# VI predictions of other invasion types

########
# FIG 4I
########

conditions_invasion <- c("TRUE",
                         "TRUE",
                         "LMPVI != 'VI (G3)'",
                         "TRUE",
                         "LMPVI != 'VI (G3)'",
                         "TRUE",
                         "LMPVI != 'VI (G3)'",
                         "TRUE")

outcome_variables_invasion <- c("VI",
                                "VPI",
                                "VPI",
                                "LI",
                                "LI",
                                "STAS",
                                "STAS",
                                "AnyInvasion")

titles_invasion <- c("Predicting VI",
                     "Predicting VPI",
                     "Predicting VPI (in absence of VI)",
                     "Predicting LI",
                     "Predicting LI (in absence of VI)",
                     "Predicting STAS",
                     "Predicting STAS (in absence of VI)",
                     "Predicting Any Invasion")

# Perform the analysis for each condition

results_invasion <- lapply(seq_along(conditions_invasion), 
                           function(i) perform_roc_analysis(validation_dge$samples, 
                                                            conditions_invasion[i], 
                                                            outcome_variables_invasion[i], 
                                                            titles_invasion[i]))

# Extract the results

aucs_invasion <- sapply(results_invasion, function(res) res$auc)

ci_mins_invasion <- sapply(results_invasion, function(res) res$ci_min)

ci_maxs_invasion <- sapply(results_invasion, function(res) res$ci_max)

ns_invasion <- sapply(results_invasion, function(res) res$n)

fig4I <- data.frame(AUCs = aucs_invasion, ci_min = ci_mins_invasion, ci_max = ci_maxs_invasion)

fig4I$Invasion <- c('VI (n=42 VI-,17 VI+)',
                    'VPI (n=47 VPI-,12 VPI+)',
                    'VPI (VI-) (n=37 VPI-,5 VPI+)',
                    'LI (n=33 LI-,26 LI+)',
                    'LI (VI-) (n=25 LI-,17 LI+)',
                    'STAS (n=36 STAS-,23 STAS+)',
                    'STAS (VI-) (n=29 STAS-,13 STAS+)',
                    'Any Invasion (n=20-,39+)')

addWorksheet(source_data, "Figure 4I")

writeData(source_data, sheet = "Figure 4I", 
          x = "Figure 4I", 
          startCol =  1, startRow = 1)

writeData(source_data, sheet = 'Figure 4I',
          x = fig4I, startCol = 1, startRow = 3)

fig4I_gg <- ggplot(fig4I, aes(x = reorder(Invasion, AUCs), y = AUCs)) + 
  geom_bar(stat = "identity", width = 0.8, alpha = 0.5) +
  coord_flip() +
  geom_errorbar(aes(x = Invasion, ymin = ci_min, ymax = ci_max), width = 0.2, 
                colour = "black", alpha = 0.9, linewidth = 1) +
  theme_classic() +
  geom_abline(slope = 0, intercept = 0.50, linetype = "dashed", 
              alpha = 0.7, color = "red") + 
  labs(title = "VI predictor performance\n by LUAD invasion type\n", 
       x = "Subgroup\n", 
       y = "\nAUROC") +
  theme(plot.title = element_text(size = 22, hjust = 0.5),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        legend.position = "top", 
        axis.text = element_text(size = 20),
        text = element_text(family = 'Calibri')) +
  scale_y_continuous(limits = c(0, 1))

fig4I_gg

ggsave(here('figures','fig4I.tiff'),
       plot = fig4I_gg,
       device = 'tiff',
       width = desired_width,
       height = desired_height,
       units = 'in', dpi = TIFF_dpi)



h2o.shutdown(prompt = FALSE)

# save source data file

saveWorkbook(source_data, 
             here('data','steiner_source_data_bulkRNAseq.xlsx'), 
             overwrite = TRUE)


sessionInfo()


