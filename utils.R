# utils.R

library(RColorBrewer)

# function for data subsetting

'%!in%' <- function(x,y)!('%in%'(x,y))

# function for reformatting p values into scientific notation for plotting

scientific_10 <- function(x) {
  parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
}


# function to convert ensemble gene names to hgnc on edgeR DGEList class

ensemble_to_hgnc <- function(dge_list) {
  
  dge_list$counts <- as.data.frame(dge_list$counts)
  
  has_version <- grepl("\\.", rownames(dge_list$counts)[1])
  
  if (has_version) {
    dge_list$counts <- dge_list$counts %>%
      mutate(ensembl_gene_id_version = rownames(.))
    
    dge_list$counts <- dge_list$counts %>%
      left_join(gene_annotations, 
                by = c("ensembl_gene_id_version" = "gene_id_version"))
    
  } else {
    dge_list$counts <- dge_list$counts %>%
      mutate(ensembl_gene_id = rownames(.))
    
    dge_list$counts <- dge_list$counts %>%
      left_join(gene_annotations, by = c("ensembl_gene_id" = "gene_id"))
  }
  
  dge_list$counts <- dge_list$counts %>%
    filter(!is.na(gene_name)) %>%
    mutate(gene_name = make.names(gene_name, unique = TRUE)) %>%
    column_to_rownames("gene_name") %>%
    select(-one_of("ensembl_gene_id", "ensembl_gene_id_version", 
                   names(gene_annotations)))
  
  colnames(dge_list$counts) <- rownames(dge_list$samples)
  
  return(dge_list)
}

# function to z score genes

zscore.rows <- function(x){
  
  if(is.null(nrow(x)) == FALSE){
    
    return(t(apply(x, 1, function(x) (x - mean(na.omit(x)))/sd(na.omit(x)))))
  }
  
  else{
    
    return(x)
    
  }
}

# define color palettes

VI_cluster_colors <- c("1" = "#cab2d6","2" = "#b2df8a",
                       "3" = "#E31A1C","4" = "#A6CEE3")

colorss <- colorRampPalette(brewer.pal(name="Paired", n = 12))(16)

bluered <- colorRampPalette(c("#007FFF", "white", "red"))(256)

