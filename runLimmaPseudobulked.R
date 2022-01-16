# Test for differential expression within each cell type (or cluster) across experimental groups using pseudobulk data in limma.
# Function also produces diagnostic plots for normalization, mean-variance relationship, as well as interactive results plots from Glimma
# This version is designed for simple model: ~ group + batch
# Can be edited as needed for other models

# "clust" is cell type or cluster of interest (must match name of item in "pseudobulk")
# "pseudobulk" is a list of gene by sample count matrix, named by cell type or cluster usually
# "info" is sample metadata (not cell metadata). n rows = n biological samples
# "contrasts" define differential expression contrasts among groups. Naming must conform to design:
# e.g., "AvBinF = groupA.F-groupB.F, AvBinM = groupALS.M-groupWT.M, Interaction = (groupA.F-groupB.F)-(groupA.M-groupB.M), AvB = (groupA.F + groupA.M)/2 - (groupB.F + groupB.M)/2"

library(limma)
library(edgeR)
library(Glimma)

runLimmaPseudobulked <- function(clust, pseudobulk, info, contrasts){
  setwd(wd2)
  keep.samples <- colSums(pseudobulk[[clust]])>0
  info <- info[keep.samples, ]
  temp <- pseudobulk[[clust]]
  temp <- temp[ , keep.samples]

  cnt <- DGEList(temp)
  cnt$samples$sample <- sample <- rownames(info)
  cnt$samples$batch <- batch <- info$Batch
  cnt$samples$group <- group <- info$Group

  cnt$genes <- gene.info

  design <- model.matrix(~0+group+batch)

  ## Diagnostic plots
  lcpm <- cpm(cnt, log=TRUE)
  L <- mean(cnt$samples$lib.size) * 1e-6
  M <- median(cnt$samples$lib.size) * 1e-6
  lcpm.cutoff <- log2(10/M + 2/L)
  library(RColorBrewer)
  col.group <- as.factor(group)
  cols <- brewer.pal(nlevels(col.group), "Set1")
  col.group <- cols[col.group]
  nsamples <- ncol(cnt)

  # Plot data before and after gene filtering
  pdf(paste0("./diagnostics/filtdens.",clust,".pdf"), width=14)
  par(mfrow=c(1,2))
  plot(density(lcpm[ , 1]), col= col.group[1], lwd=2, las=2, main = "Raw data", xlab="log(CPM)")
  abline(v = lcpm.cutoff, lty=3)
  for (i in 2:nsamples){
    den <- density(lcpm[ , i])
    lines(den$x, den$y, col = col.group[i], lwd=2)
  }
  legend("topright", rownames(cnt$samples), text.col = col.group, bty = "n")

  # Remove genes with no or low expression
  keep.exprs <- filterByExpr(cnt, group = group)
  cnt <- cnt[keep.exprs, , keep.lib.sizes = FALSE]

  lcpm <- cpm(cnt, log=TRUE)
  plot(density(lcpm[,1]), col=col.group[1], lwd=2, las=2, main="Filtered data", xlab="log(CPM)")
  abline(v=lcpm.cutoff, lty=3)
  for (i in 2:nsamples){
    den <- density(lcpm[,i])
    lines(den$x, den$y, col=col.group[i], lwd=2)
  }
  legend("topright", rownames(cnt$samples), text.col=col.group, bty="n")
  dev.off()

  # Boxplots of normalized and unnormalized data
  x <- cnt
  x$samples$norm.factors <- 1
  x$counts[,1] <- ceiling(x$counts[,1]*0.05)
  x$counts[,2] <- x$counts[,2]*5
  lcpm2 <- cpm(x, log=TRUE)
  pdf(paste0("./diagnostics/normbox.",clust,".pdf"), width=14)
  par(mfrow=c(1,2))
  par(mar=c(8, 4.1, 4.1, 2.1))
  boxplot(lcpm2, las=2, col=col.group, main="Unnormalized data", ylab="log(CPM)")

  # Normalize
  cnt <- calcNormFactors(cnt, method = "TMMwsp")
  lcpm <- cpm(cnt, log=TRUE)
  par(mar=c(8, 4.1, 4.1, 2.1))
  boxplot(lcpm, las=2, col=col.group, main="Normalized data", ylab="log(CPM)")
  dev.off()

  # Scale data with voom transformation, generate limma-voom object, plot mean-variance trend
  pdf(paste0("./diagnostics/meanvar.",clust,".pdf"), width=14)
  par(mfrow=c(1,2), mar=c(5.1, 4.1, 4.1, 2.1))

  v <- voom(cnt, design, plot=TRUE)

  contrasts.command <- paste0("makeContrasts(", contrasts, ", levels = colnames(design))")
  contr.matrix <- eval(parse(text=contrasts.command))

  vfit <- lmFit(v, design)
  vfit <- contrasts.fit(vfit, contr.matrix)
  efit <- eBayes(vfit, robust=TRUE)
  eres <- decideTests(efit, p.value = 0.05)

  plotSA(efit, main="Final model: Mean-variance trend")
  dev.off()

  # Create and save results tables
  tables <- lapply(tests, function(a) as.data.frame(topTable(efit, coef=a, p.value = 1, n=Inf, sort.by = "P")))
  names(tables) <- tests

  lapply(names(tables), function(q) write.table(tables[[q]], file=paste0("./limma/", q, ".limma.", clust, ".txt"), sep="\t", quote=F, row.names=F))

  results <- list(info, cnt, contr.matrix, v, vfit, efit, eres)
  results <- c(results, tables)
  names(results) <- c("info", "cnt", "contr.matrix", "v", "vfit", "efit", "eres", tests)

  # Run Glimma for interactive plots
  dir.create(paste0("./interactive.plots/MDS/", clust))
  setwd(paste0("./interactive.plots/MDS/", clust))
  glMDSPlot(v$E, labels=sample, groups=info, launch=FALSE)

  c=1
  for (t in tests){
    dir.create(paste0(wd2, "/interactive.plots/", t, "/", clust))
    setwd(paste0(wd2, "/interactive.plots/", t, "/", clust))
    glMDPlot(efit, coef=t, status=abs(eres[ , t]), main=t, side.main="SYMBOL", counts=v$E, groups=group, launch=FALSE)
    c=c+1
  }

  return(results)
}
