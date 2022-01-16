# Test for differential expression within each cell type (or cluster) across experimental groups using mixed/hierarchical model in limma with dream extension.
# This version is designed for simple model: ~ group + 1|sample
# Can be edited as needed for more complex models

# "tests" are contrast names
# "gene.info" is metadata for genes, such as ENTREZ ID, transcript length, gene type, etc.
# "contrasts" define differential expression contrasts among groups. Naming must conform to design:
# e.g., "AvBinF = groupA.F-groupB.F, AvBinM = groupALS.M-groupWT.M, Interaction = (groupA.F-groupB.F)-(groupA.M-groupB.M), AvB = (groupA.F + groupA.M)/2 - (groupB.F + groupB.M)/2"

library(limma)
library(edgeR)
library(variancePartition)
library(BiocParallel)

param = SnowParam(4, "SOCK", progressbar=TRUE)
register(param)

runLimmaDream <- function(seurat, clust.var, clust, contrasts){
  # Subset data for cells of interest
  keep.cells <- seurat@meta.data[[clust.var]] == clust
  info <- seurat@meta.data[keep.cells, c("SampleID", "Group")]
  temp <- data.frame(seurat[["RNA"]]@counts[ , keep.cells])

  # Set up limma object
  cnt <- DGEList(temp)
  cnt$samples$cell <- cell <- rownames(info)
  cnt$samples$sample <- sample <- info$SampleID
  cnt$samples$group <- group <- info$Group # combination of genotype and sex

  cnt$genes <- gene.info

  # Remove genes with no or low expression - adjust filtering as desired
  #stricter
  #keep.exprs <- filterByExpr(cnt, group=group, min.prop = 0.3, min.count=3)

  #looser
  keep.exprs <- rowMeans(cpm(cnt)>1) >= 0.3 #looser
  cnt <- cnt[keep.exprs, , keep.lib.sizes=FALSE]

  # Normalize
  cnt <- calcNormFactors(cnt, method = "TMMwsp")

  # estimate weights using linear mixed model of dream
  # voomWithDreamWeights() replaces voom() to estimate precision weights
  form <- ~ 0 + group + (1|sample) # model formula

  pdf(paste0("meanvar.", clust, ".pdf"))
  v.dream = voomWithDreamWeights(cnt, form, cnt$samples, plot=TRUE)
  dev.off()

  # define contrasts, edit according to study
  design <- model.matrix(~0+group) # just to get colnames for contrasts in this case, not the actual model
  contrasts.command <- paste0("makeContrasts(", contrasts, ", levels = colnames(design))")
  contr.matrix <- eval(parse(text=contrasts.command))

  # fit dream model with contrasts
  # dream() replaces lmFit() to estimate regression coefficients
  fit = dream(v.dream, form, cnt$samples, contr.matrix) #, ddf ="Kenward-Roger"

  # Results tables
  tables <- lapply(colnames(contr.matrix), function(x) as.data.frame(topTable(fit, coef=x, number=Inf, sort.by = "P", p.value = 1 )))
names(tables) <- colnames(contr.matrix)

  lapply(names(tables), function(x) write.table(tables[[x]], file=paste("limma", clust, x, "txt", sep="."), sep="\t", quote=F, row.names=F))

  results <- list(cnt, v.dream, fit)
  results <- c(results, tables)
  names(results) <- c("cnt", "v.dream", "fit", names(tables))

  return(results)
}

# Example to run on HPC with SLURM
# for celltype in `cat cell.types`; do sbatch -p cpu_short --mem-per-cpu=24G --cpus-per-task 5 -t 00-12:00:00 --wrap="Rscript --vanilla dream.R ${celltype}"; done
