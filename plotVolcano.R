# Function to make volcano plots from results of runLimmaPseudobulked() or runLimmaDream()

library(ggplot2)
library(ggrepel)

plotVolcano <- function(lim.res, contr, p.thresh=0.05, p.type="adj.P.Val", n.lab){
  x <- lim.res[[contr]][ , "logFC"]
  y <- -log10(lim.res[[contr]][ , p.type])
  pt.col <- ifelse(lim.res[[contr]][ , p.type] < p.thresh & lim.res[[contr]][ , "logFC"] > 0, "red2", ifelse(lim.res[[contr]][ , p.type] < p.thresh & lim.res[[contr]][ , "logFC"] < 0, "royalblue", "gray30"))
  a <- ifelse(lim.res[[contr]][ , p.type], 0.8, 0.5)

  pt.lab <- ifelse(lim.res[[contr]][ , p.type] < p.thresh, rownames(lim.res[[contr]]), NA)
  pt.lab <- c(pt.lab[1:n.lab], rep(NA, length(pt.lab)-n.lab))

  dat <- data.frame(x=x, y=y, pt.col=pt.col, a=a, pt.lab=pt.lab)

  ggplot(dat, aes(x=x, y=y)) + geom_point(color=pt.col) + theme(axis.line = element_line(colour = "black"), panel.background = element_rect(fill = "white"), panel.grid = element_blank()) + ggtitle(contr) + xlab("logFC") + ylab("-log(FDR-adjusted P)") + geom_text_repel(aes(label = pt.lab, fontface = "italic"), box.padding = unit(0.25, "lines"), max.overlaps=50)
}

# Example looping across cell types for specific contrast
# lapply(limma.results, function(x) plotVolcano(lim.res=x, contr="KO_v_Control", p.thresh=0.05, p.type="adj.P.Val", n.lab=20))
