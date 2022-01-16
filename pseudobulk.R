pseudobulk <- function(seurat, info, clust.var, clust){
  counts <- seurat[["RNA"]]@counts[ , seurat@meta.data[[clust.var]] == clust]
  counts.meta <- seurat@meta.data[colnames(counts),]
  counts.sample <- lapply(rownames(info), function(x) counts[, rownames(counts.meta)[counts.meta$SampleID==x]])
  counts.sample.sum <- lapply(counts.sample, function(x) apply(as.data.frame(x), 1, sum))
  counts.df <- data.frame(matrix(unlist(counts.sample.sum), ncol=length(counts.sample.sum), byrow=F))
  colnames(counts.df) <- rownames(info)
  rownames(counts.df) <- rownames(counts)
  return(counts.df)
}
