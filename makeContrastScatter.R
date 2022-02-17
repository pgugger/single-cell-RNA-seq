# Function to make scatterplots of logFC across two contrasts, coloring points that are significant in both

library(ggplot2)
library(ggrepel)
library(scales)

makeContrastScatter <- function(contrast1, contrast2, lim.res, p.threshold=0.05, max.labels=25, plot.title = "", out.dir="."){

  dir.create(out.dir)
  setwd(out.dir)

  # Merge limma results table for relevant contrasts
  lim.merge <- merge(x=lim.res[[contrast1]], y=lim.res[[contrast2]],  all=T, by="row.names")
rownames(lim.merge) <- lim.merge$Row.names

  # Correlation of logFC across contrasts
  cor.res <- cor(lim.merge$logFC.x, lim.merge$logFC.y, use="pairwise.complete.obs")

  # Define colors, alpha, and size based on significance and quadrant
  col.vec <- ifelse(lim.merge$adj.P.Val.x < p.threshold & lim.merge$adj.P.Val.y < p.threshold & lim.merge$logFC.x > 0 & lim.merge$logFC.y > 0, "red2", ifelse(lim.merge$adj.P.Val.x < p.threshold & lim.merge$adj.P.Val.y < p.threshold & lim.merge$logFC.x < 0 & lim.merge$logFC.y < 0, "royalblue", ifelse(lim.merge$adj.P.Val.x < p.threshold & lim.merge$adj.P.Val.y < p.threshold & lim.merge$logFC.x > 0 & lim.merge$logFC.y < 0, "black", ifelse(lim.merge$adj.P.Val.x < p.threshold & lim.merge$adj.P.Val.y < p.threshold & lim.merge$logFC.x < 0 & lim.merge$logFC.y > 0, "black", "gray70"))))

  alpha.vec <- ifelse(lim.merge$adj.P.Val.x < p.threshold & lim.merge$adj.P.Val.y < p.threshold & lim.merge$logFC.x > 0 & lim.merge$logFC.y > 0, 1, ifelse(lim.merge$adj.P.Val.x < p.threshold & lim.merge$adj.P.Val.y < p.threshold & lim.merge$logFC.x < 0 & lim.merge$logFC.y < 0, 1, ifelse(lim.merge$adj.P.Val.x < p.threshold & lim.merge$adj.P.Val.y < p.threshold & lim.merge$logFC.x > 0 & lim.merge$logFC.y < 0, 1, ifelse(lim.merge$adj.P.Val.x < p.threshold & lim.merge$adj.P.Val.y < p.threshold & lim.merge$logFC.x < 0 & lim.merge$logFC.y > 0, 1, 0.5))))

  size.vec <- ifelse(lim.merge$adj.P.Val.x < p.threshold & lim.merge$adj.P.Val.y < p.threshold & lim.merge$logFC.x > 0 & lim.merge$logFC.y > 0, 1, ifelse(lim.merge$adj.P.Val.x < p.threshold & lim.merge$adj.P.Val.y < p.threshold & lim.merge$logFC.x < 0 & lim.merge$logFC.y < 0, 1, ifelse(lim.merge$adj.P.Val.x < p.threshold & lim.merge$adj.P.Val.y < p.threshold & lim.merge$logFC.x > 0 & lim.merge$logFC.y < 0, 1, ifelse(lim.merge$adj.P.Val.x < p.threshold & lim.merge$adj.P.Val.y < p.threshold & lim.merge$logFC.x < 0 & lim.merge$logFC.y > 0, 1, 0.5))))

  pca <- prcomp(~logFC.x+logFC.y, data=lim.merge, center = TRUE, scale = TRUE, na.action = na.omit)

  label.vec <- ifelse((rownames(lim.merge) %in% c(names(head(sort(pca$x[,1]), max.labels)), names(tail(sort(pca$x[,1]), max.labels)), names(head(sort(pca$x[,2]), max.labels)), names(tail(sort(pca$x[,2]), max.labels)))) & lim.merge$adj.P.Val.x < p.threshold & lim.merge$adj.P.Val.y < p.threshold, rownames(lim.merge), NA)

  out.plot <- ggplot(lim.merge, aes(x=logFC.x, y=logFC.y)) + geom_hline(yintercept=0, size=0.1) + geom_vline(xintercept=0, size=0.1) + geom_point(color=col.vec, alpha=alpha.vec, size=size.vec) + theme(panel.background = element_rect(fill = "white"), panel.grid = element_blank(), axis.line = element_line(color="black", size = 0.5)) + ggtitle(plot.title) + xlab(paste0("logFC ", contrast1)) + ylab(paste0("logFC ", contrast2)) + geom_text_repel(size=2, aes(label = label.vec, fontface = "italic"), box.padding = unit(0.1, "lines"), max.overlaps=25)

  ggsave(filename = paste0("scatterplot.", plot.title, ".", contrast1, ".", contrast2, ".pdf"), plot = out.plot, device = "pdf")

  res <- list(cor.res, lim.merge, out.plot)
  names(res) <- c("pearson.correlation", "limma.merged", "ggplot.output")

  return(res)
}

# mapply to iterate over limma results for each cell type
# scatter.out <- mapply(x=limma.results, y=names(limma.results), function(x,y) makeContrastScatter(contrast1 = "KO_v_Control", contrast2 = "doubleKO_v_KO", lim.res = x, p.threshold = 0.1, max.labels = 100, plot.title = y, out.dir = "."), SIMPLIFY = F, USE.NAMES = T)
