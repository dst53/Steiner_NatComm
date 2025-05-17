
#######
# SETUP
#######

# load R packages (all installed under R 4.2.1)

library(dplyr)
library(zellkonverter)
library(basilisk)
library(here)

# download Salcher atlas

url <- 'https://zenodo.org/records/6411868/files/extended_atlas.h5ad.gz?download=1'

destfile <- here('data','extended_atlas.h5ad.gz')

options(timeout = 3600)  # Set timeout to 1 hour

download.file(url, destfile, mode = "wb")

gunzip(here('data','extended_atlas.h5ad.gz'), overwrite = TRUE, remove = FALSE)


# convert anndata object to sce object

Salcher_lung_cancer_atlas <- readH5AD(here('data','extended_atlas.h5ad'))

saveRDS(Salcher_lung_cancer_atlas, file = here('data','extended_atlas.rds'))



