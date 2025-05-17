# Steiner_NatComm


#### Step 1: clone this repository to a local folder
#### Step 2: use reviewer access token to download GSE273377_exp_count.txt.gz file from subseries GSE273377 within superseries GSE273528, unzip and place in Steiner_NatComm/data
#### Step 3: download GSE273378_RAW.tar	from subseries GSE273378 within GSE273528, unzip and place GSE273378_RAW folder in Steiner_NatComm/data
#### Step 4: ensure that R version and installed package versions match versions below for each script prior to running
#### Step 5: run Steiner_NatComm_bulkRNAseq_analysis.R
#### Step 6: run Steiner_NatComm_bulkRNAseq_biopsies_analysis.R (pending GEO data submission)
#### Step 7: run Steiner_NatComm_salcher_conversion.R
#### Step 8: Install miniconda/4.11.0 and follow the remianing cytospace installation instructions. To install v1.0.6 specifically, clone cytospace v1.0.6 into Steiner_NatComm repository folder: ```git clone --branch v1.0.6 --single-branch https://github.com/digitalcytometry/cytospace.git```. 
#### Step 9: run cytospace.sh
#### Step 10: run Steiner_NatComm_stRNAseq_analysis.R

Session info for Steiner_NatComm_bulkRNAseq_analysis.R
```
sessionInfo()

R version 4.2.1 (2022-06-23)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: AlmaLinux 8.10 (Cerulean Leopard)

Matrix products: default
BLAS:   /share/pkg.7/r/4.2.1/install/lib64/R/lib/libRblas.so
LAPACK: /share/pkg.7/r/4.2.1/install/lib64/R/lib/libRlapack.so

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] grid      stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] openxlsx_4.2.5.2            GEOquery_2.66.0             stringr_1.5.0              
 [4] fst_0.9.8                   TCGAbiolinks_2.25.3         ggprism_1.0.5              
 [7] rstatix_0.7.2               scales_1.2.1                enrichplot_1.18.3          
[10] forcats_1.0.0               modelsummary_1.4.1          cmprsk_2.2-11              
[13] clusterProfiler_4.6.0       ComplexUpset_1.3.3          rcompanion_2.4.26          
[16] corrplot_0.92               RColorBrewer_1.1-3          tableone_0.13.2            
[19] ggrepel_0.9.3               CePa_0.8.0                  estimate_1.0.13            
[22] tibble_3.1.8                phenoTest_1.46.0            BMA_3.18.17                
[25] rrcov_1.7-2                 inline_0.3.19               leaps_3.1                  
[28] Heatplus_3.6.0              annotate_1.76.0             XML_3.99-0.10              
[31] AnnotationDbi_1.60.0        enrichR_3.1                 tidyr_1.3.0                
[34] plyr_1.8.8                  robustbase_0.95-1           data.table_1.14.6          
[37] ggfortify_0.4.15            ggcorrplot_0.1.4            msigdbr_7.5.1              
[40] broom_1.0.3                 h2o_3.40.0.4                caret_6.0-93               
[43] SummarizedExperiment_1.28.0 Biobase_2.58.0              GenomicRanges_1.50.2       
[46] GenomeInfoDb_1.34.9         IRanges_2.32.0              S4Vectors_0.36.2           
[49] BiocGenerics_0.44.0         MatrixGenerics_1.10.0       matrixStats_0.63.0         
[52] rms_6.5-0                   SparseM_1.81                Hmisc_4.8-0                
[55] Formula_1.2-4               lattice_0.20-45             ggsci_3.0.0                
[58] plotROC_2.3.0               pROC_1.18.0                 GSVA_1.46.0                
[61] edgeR_3.40.2                limma_3.54.1                pheatmap_1.0.12            
[64] sva_3.46.0                  BiocParallel_1.32.5         genefilter_1.80.3          
[67] mgcv_1.8-41                 nlme_3.1-162                survminer_0.4.9            
[70] ggpubr_0.6.0                ggplot2_3.4.2               survival_3.5-3             
[73] dplyr_1.1.0                 readxl_1.4.2               

loaded via a namespace (and not attached):
  [1] rsvd_1.0.5                  svglite_2.1.1               class_7.3-21               
  [4] foreach_1.5.2               lmtest_0.9-40               crayon_1.5.2               
  [7] MASS_7.3-58.2               rhdf5filters_1.10.0         backports_1.4.1            
 [10] ellipse_0.4.3               GOSemSim_2.24.0             rlang_1.1.1                
 [13] XVector_0.38.0              HDO.db_0.99.1               irlba_2.3.5.1              
 [16] filelock_1.0.2              rjson_0.2.21                bit64_4.0.5                
 [19] glue_1.6.2                  parallel_4.2.1              DOSE_3.24.2                
 [22] tidyselect_1.2.0            km.ci_0.5-6                 zoo_1.8-11                 
 [25] xtable_1.8-4                MatrixModels_0.5-1          magrittr_2.0.3             
 [28] evaluate_0.21               cli_3.6.0                   zlibbioc_1.44.0            
 [31] rstudioapi_0.14             rpart_4.1.19                fastmatch_1.1-3            
 [34] treeio_1.22.0               BiocSingular_1.14.0         xfun_0.39                  
 [37] gson_0.0.9                  cluster_2.1.4               caTools_1.18.2             
 [40] tidygraph_1.2.3             KEGGREST_1.38.0             expm_0.999-7               
 [43] quantreg_5.94               ape_5.7                     listenv_0.9.0              
 [46] Biostrings_2.66.0           png_0.1-7                   future_1.31.0              
 [49] ipred_0.9-13                withr_2.5.0                 ggforce_0.4.1              
 [52] bitops_1.0-7                RBGL_1.74.0                 cellranger_1.1.0           
 [55] GSEABase_1.60.0             pcaPP_2.0-3                 hardhat_1.2.0              
 [58] e1071_1.7-13                survey_4.1-1                pillar_1.9.0               
 [61] gplots_3.1.3                cachem_1.0.8                multcomp_1.4-22            
 [64] hopach_2.58.0               DelayedMatrixStats_1.20.0   vctrs_0.5.2                
 [67] ellipsis_0.3.2              generics_0.1.3              nortest_1.0-4              
 [70] lava_1.7.1                  tools_4.2.1                 foreign_0.8-84             
 [73] tweenr_2.0.2                munsell_0.5.0               fgsea_1.24.0               
 [76] proxy_0.4-27                DelayedArray_0.24.0         fastmap_1.1.1              
 [79] compiler_4.2.1              abind_1.4-5                 gt_0.10.1                  
 [82] DescTools_0.99.47           GenomeInfoDbData_1.2.9      prodlim_2019.11.13         
 [85] gridExtra_2.3               deldir_1.0-6                utf8_1.2.3                 
 [88] BiocFileCache_2.6.0         recipes_1.0.4               jsonlite_1.8.4             
 [91] gld_2.6.6                   graph_1.76.0                ScaledMatrix_1.6.0         
 [94] tidytree_0.4.2              carData_3.0-5               sparseMatrixStats_1.10.0   
 [97] lazyeval_0.2.2              car_3.1-1                   latticeExtra_0.6-30        
[100] fstcore_0.9.14              checkmate_2.1.0             rmarkdown_2.22             
[103] sandwich_3.0-2              cowplot_1.1.1               webshot_0.5.4              
[106] downloader_0.4              igraph_1.3.2                HDF5Array_1.26.0           
[109] systemfonts_1.0.4           htmltools_0.5.4             memoise_2.0.1              
[112] modeltools_0.2-23           locfit_1.5-9.7              graphlayouts_0.8.4         
[115] viridisLite_0.4.1           digest_0.6.31               assertthat_0.2.1           
[118] rappdirs_0.3.3              KMsurv_0.1-5                RSQLite_2.2.20             
[121] TCGAbiolinksGUI.data_1.18.0 yulab.utils_0.0.6           future.apply_1.10.0        
[124] Exact_3.2                   blob_1.2.3                  survMisc_0.5.6             
[127] splines_4.2.1               Rhdf5lib_1.20.0             RCurl_1.98-1.7             
[130] hms_1.1.2                   rhdf5_2.42.0                colorspace_2.1-0           
[133] base64enc_0.1-3             aplot_0.1.9                 libcoin_1.0-9              
[136] nnet_7.3-19                 Rcpp_1.0.10                 coin_1.4-2                 
[139] mvtnorm_1.1-3               multcompView_0.1-8          fansi_1.0.4                
[142] tzdb_0.3.0                  tables_0.9.17               parallelly_1.34.0          
[145] ModelMetrics_1.2.2.2        R6_2.5.1                    lifecycle_1.0.3            
[148] polspline_1.1.22            rootSolve_1.8.2.3           zip_2.2.2                  
[151] curl_5.0.0                  ggsignif_0.6.4              Matrix_1.5-3               
[154] qvalue_2.30.0               TH.data_1.1-1               iterators_1.0.14           
[157] gower_1.0.1                 htmlwidgets_1.6.2           polyclip_1.10-4            
[160] beachmat_2.14.0             biomaRt_2.54.0              purrr_1.0.1                
[163] shadowtext_0.1.2            gridGraphics_0.5-1          timechange_0.2.0           
[166] rvest_1.0.3                 globals_0.16.2              insight_0.19.3             
[169] lmom_2.9                    htmlTable_2.4.1             patchwork_1.1.2            
[172] codetools_0.2-19            lubridate_1.9.2             GO.db_3.16.0               
[175] gtools_3.9.4                prettyunits_1.1.1           SingleCellExperiment_1.20.1
[178] dbplyr_2.3.0                gtable_0.3.1                DBI_1.1.3                  
[181] ggfun_0.1.1                 httr_1.4.4                  KernSmooth_2.23-20         
[184] stringi_1.7.12              progress_1.2.2              farver_2.1.1               
[187] reshape2_1.4.4              viridis_0.6.2               Rgraphviz_2.42.0           
[190] ggtree_3.6.2                timeDate_4022.108           DT_0.27                    
[193] xml2_1.3.3                  boot_1.3-28.1               kableExtra_1.3.4           
[196] readr_2.1.4                 interp_1.1-3                ggplotify_0.1.1            
[199] Category_2.64.0             DEoptimR_1.0-11             bit_4.0.5                  
[202] scatterpie_0.1.8            jpeg_0.1-10                 hgu133a.db_3.13.0          
[205] ggraph_2.1.0                pkgconfig_2.0.3             babelgene_22.9             
[208] mitools_2.4                 knitr_1.43

```
Session info for Steiner_NatComm_bulkRNAseq_biopsies_analysis.R
```
sessionInfo()
R version 4.2.1 (2022-06-23)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: AlmaLinux 8.10 (Cerulean Leopard)

Matrix products: default
BLAS:   /share/pkg.7/r/4.2.1/install/lib64/R/lib/libRblas.so
LAPACK: /share/pkg.7/r/4.2.1/install/lib64/R/lib/libRlapack.so

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8   
 [6] LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] grid      stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] RColorBrewer_1.1-3  ggpubr_0.6.0        tidyr_1.3.0         stringr_1.5.0       ggplot2_3.4.2       pROC_1.18.0         sva_3.46.0          BiocParallel_1.32.5
 [9] genefilter_1.80.3   mgcv_1.8-41         nlme_3.1-162        tibble_3.1.8        tableone_0.13.2     dplyr_1.1.0         readxl_1.4.2        h2o_3.40.0.4       
[17] openxlsx_4.2.5.2    here_1.0.1          lmerTest_3.1-3      lme4_1.1-31         Matrix_1.5-3        edgeR_3.40.2        limma_3.54.1       

loaded via a namespace (and not attached):
  [1] utf8_1.2.3                  spatstat.explore_3.0-6      reticulate_1.28             tidyselect_1.2.1            RSQLite_2.2.20             
  [6] AnnotationDbi_1.60.0        htmlwidgets_1.6.2           Rtsne_0.16                  munsell_0.5.0               ragg_1.2.5                 
 [11] codetools_0.2-19            ica_1.0-3                   units_0.8-0                 future_1.31.0               miniUI_0.1.1.1             
 [16] withr_2.5.0                 spatstat.random_3.1-3       colorspace_2.1-0            progressr_0.13.0            Biobase_2.58.0             
 [21] rstudioapi_0.14             Seurat_4.3.0.1              stats4_4.2.1                SingleCellExperiment_1.20.1 ROCR_1.0-11                
 [26] ggsignif_0.6.4              tensor_1.5                  listenv_0.9.0               labeling_0.4.2              MatrixGenerics_1.10.0      
 [31] GenomeInfoDbData_1.2.9      polyclip_1.10-4             farver_2.1.1                bit64_4.0.5                 rprojroot_2.0.3            
 [36] parallelly_1.34.0           vctrs_0.5.2                 generics_0.1.3              R6_2.5.1                    GenomeInfoDb_1.34.9        
 [41] locfit_1.5-9.7              bitops_1.0-7                spatstat.utils_3.0-1        cachem_1.0.8                DelayedArray_0.24.0        
 [46] promises_1.2.0.1            scales_1.3.0                gtable_0.3.1                globals_0.16.2              goftest_1.2-3              
 [51] rlang_1.1.1                 systemfonts_1.0.4           splines_4.2.1               rstatix_0.7.2               lazyeval_0.2.2             
 [56] spatstat.geom_3.0-6         broom_1.0.3                 reshape2_1.4.4              abind_1.4-5                 backports_1.4.1            
 [61] httpuv_1.6.9                tools_4.2.1                 ellipsis_0.3.2              proxy_0.4-27                BiocGenerics_0.44.0        
 [66] ggridges_0.5.4              Rcpp_1.0.10                 plyr_1.8.8                  zlibbioc_1.44.0             purrr_1.0.1                
 [71] RCurl_1.98-1.7              deldir_1.0-6                pbapply_1.7-0               cowplot_1.1.1               S4Vectors_0.36.2           
 [76] zoo_1.8-11                  haven_2.5.1                 SeuratObject_4.1.3          SummarizedExperiment_1.28.0 ggrepel_0.9.3              
 [81] cluster_2.1.4               survey_4.1-1                magrittr_2.0.3              data.table_1.14.6           scattermore_0.8            
 [86] lmtest_0.9-40               RANN_2.6.1                  ggnewscale_0.4.8            fitdistrplus_1.1-8          matrixStats_0.63.0         
 [91] SPATA2_3.1.4                hms_1.1.2                   patchwork_1.1.2             mime_0.12                   fftwtools_0.9-11           
 [96] xtable_1.8-4                XML_3.99-0.10               jpeg_0.1-10                 IRanges_2.32.0              gridExtra_2.3              
[101] compiler_4.2.1              KernSmooth_2.23-20          crayon_1.5.2                minqa_1.2.5                 SPATAData_1.0.0            
[106] htmltools_0.5.4             later_1.3.1                 tzdb_0.3.0                  tiff_0.1-11                 ggprism_1.0.5              
[111] DBI_1.1.3                   MASS_7.3-58.2               boot_1.3-28.1               car_3.1-1                   readr_2.1.4                
[116] mitools_2.4                 cli_3.6.4                   parallel_4.2.1              igraph_1.3.2                forcats_1.0.0              
[121] GenomicRanges_1.50.2        pkgconfig_2.0.3             numDeriv_2016.8-1.1         sp_1.6-0                    plotly_4.10.1              
[126] spatstat.sparse_3.0-0       annotate_1.76.0             XVector_0.38.0              rematch_1.0.1               digest_0.6.31              
[131] sctransform_0.3.5           RcppAnnoy_0.0.20            spatstat.data_3.0-0         Biostrings_2.66.0           cellranger_1.1.0           
[136] leiden_0.4.3                uwot_0.1.14                 curl_5.0.0                  shiny_1.7.4                 EBImage_4.40.1             
[141] nloptr_2.0.3                lifecycle_1.0.3             jsonlite_1.8.4              carData_3.0-5               viridisLite_0.4.1          
[146] fansi_1.0.4                 labelled_2.10.0             pillar_1.9.0                lattice_0.20-45             KEGGREST_1.38.0            
[151] fastmap_1.1.1               httr_1.4.4                  survival_3.5-3              glue_1.8.0                  zip_2.2.2                  
[156] png_0.1-7                   bit_4.0.5                   class_7.3-21                stringi_1.7.12              blob_1.2.3                 
[161] textshaping_0.3.6           memoise_2.0.1               e1071_1.7-13                irlba_2.3.5.1               future.apply_1.10.0  


```
Session info for Steiner_NatComm_stRNAseq_analysis.R
```
sessionInfo()

R version 4.2.1 (2022-06-23)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: AlmaLinux 8.10 (Cerulean Leopard)

Matrix products: default
BLAS:   /share/pkg.7/r/4.2.1/install/lib64/R/lib/libRblas.so
LAPACK: /share/pkg.7/r/4.2.1/install/lib64/R/lib/libRlapack.so

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8   
 [6] LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
 [1] parallel  grid      stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] psych_2.2.9                 FNN_1.1.3.1                 R.utils_2.12.2              R.oo_1.25.0                 R.methodsS3_1.8.2          
 [6] here_1.0.1                  openxlsx_4.2.5.2            SeuratDisk_0.0.0.9021       stxBrain.SeuratData_0.1.1   SeuratData_0.2.2           
[11] ggprism_1.0.5               colorspace_2.1-0            emmeans_1.8.4-1             GWmodel_2.3-1               robustbase_0.95-1          
[16] enrichR_3.1                 scMC_1.0.0                  tibble_3.1.8                sf_1.0-7                    spgwr_0.6-36               
[21] spData_2.2.2                sp_1.6-0                    purrr_1.0.1                 SPATA2_3.1.4                clusterProfiler_4.6.0      
[26] CePa_0.8.0                  estimate_1.0.13             forcats_1.0.0               ggnewscale_0.4.8            lmerTest_3.1-3             
[31] ggfortify_0.4.15            lme4_1.1-31                 corrplot_0.92               plotROC_2.3.0               pROC_1.18.0                
[36] scMerge_1.14.0              textshape_1.7.3             ggcorrplot_0.1.4            rstatix_0.7.2               tidyr_1.3.0                
[41] HDF5Array_1.26.0            rhdf5_2.42.0                DelayedArray_0.24.0         Matrix_1.5-3                sva_3.46.0                 
[46] BiocParallel_1.32.5         genefilter_1.80.3           mgcv_1.8-41                 nlme_3.1-162                cowplot_1.1.1              
[51] data.table_1.14.6           stringr_1.5.0               ape_5.7                     reshape2_1.4.4              viridis_0.6.2              
[56] viridisLite_0.4.1           mclust_6.0.0                ggpubr_0.6.0                msigdbr_7.5.1               GSVA_1.46.0                
[61] pheatmap_1.0.12             proxy_0.4-27                scran_1.26.2                scater_1.26.1               scuttle_1.8.4              
[66] patchwork_1.1.2             harmony_0.1.1               Rcpp_1.0.10                 dittoSeq_1.10.0             rmarkdown_2.22             
[71] SingleCellExperiment_1.20.1 readxl_1.4.2                hdf5r_1.3.8                 readbitmap_0.1.5            RColorBrewer_1.1-3         
[76] ggplot2_3.4.2               future_1.31.0               SeuratObject_4.1.3          Seurat_4.3.0.1              dplyr_1.1.0                
[81] plyr_1.8.8                  SummarizedExperiment_1.28.0 Biobase_2.58.0              GenomicRanges_1.50.2        GenomeInfoDb_1.34.9        
[86] IRanges_2.32.0              S4Vectors_0.36.2            BiocGenerics_0.44.0         MatrixGenerics_1.10.0       matrixStats_0.63.0         
[91] edgeR_3.40.2                limma_3.54.1                h2o_3.40.0.4               

loaded via a namespace (and not attached):
  [1] KEGGREST_1.38.0           locfit_1.5-9.7            gtsummary_1.6.3           lattice_0.20-45           spatstat.utils_3.0-1      vctrs_0.5.2              
  [7] utf8_1.2.3                blob_1.2.3                withr_2.5.0               foreign_0.8-84            lifecycle_1.0.3           cellranger_1.1.0         
 [13] munsell_0.5.0             ragg_1.2.5                ScaledMatrix_1.6.0        codetools_0.2-19          kableExtra_1.3.4          lmtest_0.9-40            
 [19] M3Drop_1.24.0             spatialreg_1.3-1          annotate_1.76.0           parallelly_1.34.0         fastmatch_1.1-3           tidycmprsk_0.2.0         
 [25] metapod_1.6.0             Rtsne_0.16                stringi_1.7.12            sctransform_0.3.5         polyclip_1.10-4           rhdf5filters_1.10.0      
 [31] yulab.utils_0.0.6         sandwich_3.0-2            goftest_1.2-3             cluster_2.1.4             ggraph_2.1.0              pkgconfig_2.0.3          
 [37] prettyunits_1.1.1         ruv_0.9.7.1               sparseMatrixStats_1.10.0  ggridges_0.5.4            lubridate_1.9.2           timechange_0.2.0         
 [43] estimability_1.4.1        httr_1.4.4                gt_0.10.1                 fftwtools_0.9-11          igraph_1.3.2              progress_1.2.2           
 [49] treeio_1.22.0             beachmat_2.14.0           graphlayouts_0.8.4        ggfun_0.1.1               gson_0.0.9                V8_4.2.2                 
 [55] htmltools_0.5.4           miniUI_0.1.1.1            pillar_1.9.0              later_1.3.1               fitdistrplus_1.1-8        glue_1.8.0               
 [61] DBI_1.1.3                 gtable_0.3.1              GOSemSim_2.24.0           rsvd_1.0.5                caTools_1.18.2            latticeExtra_0.6-30      
 [67] fastmap_1.1.1             spdep_1.3-1               AnnotationDbi_1.60.0      reldist_1.7-1             broom_1.0.3               checkmate_2.1.0          
 [73] promises_1.2.0.1          webshot_0.5.4             confuns_1.0.3             modelsummary_1.4.1        textshaping_0.3.6         mnormt_2.1.1             
 [79] ggforce_0.4.1             hms_1.1.2                 png_0.1-7                 ggtree_3.6.2              spatstat.explore_3.0-6    lazyeval_0.2.2           
 [85] Formula_1.2-4             crayon_1.5.2              svglite_2.1.1             boot_1.3-28.1             rstan_2.21.8              tidyselect_1.2.1         
 [91] xfun_0.52                 BiocSingular_1.14.0       splines_4.2.1             spacetime_1.2-8           loo_2.5.1                 survival_3.5-3           
 [97] EBImage_4.40.1            rappdirs_0.3.3            bit64_4.0.5               jpeg_0.1-10               ggsignif_0.6.4            htmlTable_2.4.1          
[103] xtable_1.8-4              DT_0.27                   cachem_1.0.8              StanHeaders_2.21.0-7      abind_1.4-5               systemfonts_1.0.4        
[109] mime_0.12                 rjson_0.2.21              aplot_0.1.9               ggrepel_0.9.3             processx_3.8.0            insight_0.19.3           
[115] spatstat.sparse_3.0-0     numDeriv_2016.8-1.1       tools_4.2.1               cli_3.6.4                 densEstBayes_1.0-2.1      magrittr_2.0.3           
[121] future.apply_1.10.0       ggplotify_0.1.1           DelayedMatrixStats_1.20.0 ggbeeswarm_0.7.1          rematch_1.0.1             qvalue_2.30.0            
[127] fgsea_1.24.0              ica_1.0-3                 pbapply_1.7-0             ggrastr_1.0.1             s2_1.1.2                  plotly_4.10.1            
[133] survMisc_0.5.6            tweenr_2.0.2              Rgraphviz_2.42.0          multcomp_1.4-22           zlibbioc_1.44.0           zip_2.2.2                
[139] inline_0.3.19             shadowtext_0.1.2          tzdb_0.3.0                ps_1.7.2                  fansi_1.0.4               tidygraph_1.2.3          
[145] GSEABase_1.60.0           xts_0.13.2                TH.data_1.1-1             tensor_1.5                ROCR_1.0-11               KernSmooth_2.23-20       
[151] readr_2.1.4               backports_1.4.1           XVector_0.38.0            interp_1.1-3              farver_2.1.1              bit_4.0.5                
[157] gplots_3.1.3              RANN_2.6.1                shiny_1.7.4               KMsurv_0.1-5              scattermore_0.8           DOSE_3.24.2              
[163] scatterpie_0.1.8          RcppAnnoy_0.0.20          downloader_0.4            SPATAData_1.0.0           rstudioapi_0.14           minqa_1.2.5              
[169] Rhdf5lib_1.20.0           spatstat.geom_3.0-6       intervals_0.15.2          tables_0.9.17             sfsmisc_1.1-14            gtools_3.9.4             
[175] beeswarm_0.4.0            listenv_0.9.0             generics_0.1.3            base64enc_0.1-3           XML_3.99-0.10             pkgbuild_1.4.0           
[181] e1071_1.7-13              spatstat.data_3.0-0       dqrng_0.3.0               GenomeInfoDbData_1.2.9    Biostrings_2.66.0         progressr_0.13.0         
[187] evaluate_0.21             memoise_2.0.1             coda_0.19-4               knitr_1.43                vipor_0.4.5               httpuv_1.6.9             
[193] class_7.3-21              irlba_2.3.5.1             classInt_0.4-8            dbscan_1.1-11             jsonlite_1.8.4            Hmisc_4.8-0              
[199] km.ci_0.5-6               RSpectra_0.16-1           babelgene_22.9            bmp_0.3                   digest_0.6.31             rprojroot_2.0.3          
[205] bitops_1.0-7              LearnBayes_2.15.1         RSQLite_2.2.20            globals_0.16.2            compiler_4.2.1            nnet_7.3-19              
[211] reticulate_1.28           statmod_1.5.0             zoo_1.8-11                carData_3.0-5             gridGraphics_0.5-1        rlang_1.1.1              
[217] nloptr_2.0.3              uwot_0.1.14               wk_0.7.1                  rvest_1.0.3               bdsmatrix_1.3-6           mvtnorm_1.1-3            
[223] htmlwidgets_1.6.2         labeling_0.4.2            callr_3.7.3               leiden_0.4.3              ComplexUpset_1.3.3        curl_5.0.0               
[229] concaveman_1.1.0          DEoptimR_1.0-11           BiocNeighbors_1.16.0      scales_1.3.0              RcppParallel_5.1.6        graph_1.76.0             
[235] startupmsg_0.9.6          enrichplot_1.18.3         HDO.db_0.99.1             deldir_1.0-6              gridExtra_2.3             tiff_0.1-11              
[241] bluster_1.8.0             distr_2.9.1               bbmle_1.0.25              RCurl_1.98-1.7            car_3.1-1                 GO.db_3.16.0             
[247] MASS_7.3-58.2             broom.helpers_1.13.0      ellipsis_0.3.2            tidytree_0.4.2            spatstat.random_3.1-3     xml2_1.3.3               
[253] survminer_0.4.9           rpart_4.1.19              R6_2.5.1                  units_0.8-0

```
