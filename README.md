# Steiner_NatComm

#### Prereqs:
##### Homebrew for system libraries if not already installed
    /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
##### R 4.2.1
##### A few R packages have C/C++ dependencies; install these libraries via Homebrew:
    brew update
    brew install \
        openssl@3 \
        hdf5 \
        freetype \
        libpng \
        pkg-config
##### add the following to ~/.zshrc or ~/.bash_profile:
    export PKG_CONFIG_PATH="/opt/homebrew/opt/freetype/lib/pkgconfig:/opt/homebrew/opt/libpng/lib/pkgconfig${PKG_CONFIG_PATH:+:}$PKG_CONFIG_PATH"
    export LDFLAGS="-L/opt/homebrew/opt/openssl@3/lib -L/opt/homebrew/opt/hdf5/lib"
    export CPPFLAGS="-I/opt/homebrew/opt/openssl@3/include -I/opt/homebrew/opt/hdf5/include"
#### Step 1: clone this repository to a local folder
#### Step 2: use reviewer access token to download GSE273377_exp_count.txt.gz file from subseries GSE273377 within superseries GSE273528, unzip and place in Steiner_NatComm/data
#### Step 3: download GSE273378_RAW.tar	from subseries GSE273378 within GSE273528, unzip and place GSE273378_RAW folder in /data
#### Step 4: This project uses renv for package management; do renv::restore() to reproduce exactly
#### Step 5: run Steiner_NatComm_bulkRNAseq_analysis.R
#### Step 6: run Steiner_NatComm_salcher_conversion.R
#### Step 7: Install miniconda/4.11.0 and follow the remianing cytospace installation instructions. To install v1.0.6 specifically, clone cytospace v1.0.6 into Steiner_NatComm repository folder: ```git clone --branch v1.0.6 --single-branch https://github.com/digitalcytometry/cytospace.git```. 
#### Step 8: run cytospace.sh
#### Step 9: run Steiner_NatComm_stRNAseq_analysis.R
#### Step 10: download GSE298714_biopsies_exp_count.txt.gz file from subseries GSE298714 within superseries GSE273528, unzip and place in Steiner_NatComm/data
#### Step 11: run Steiner_NatComm_bulkRNAseq_biopsies_analysis.R (pending GEO data submission)

