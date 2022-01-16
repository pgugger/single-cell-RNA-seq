# Wrapper functions to run ROAST and CAMERA using results object from runLimmaPseudobulked()
# Outputs results tables

# Example of annotation data
#library(gage)
#data(kegg.gs)
#kegg.mmu <- kegg.gsets("mmu")
#kegg.mmu <- kegg.mmu$kg.sets

# CAMERA
runCamera <- function(lim.res, celltype, annotations){
  gene.set.indices <- ids2indices(annotations, lim.res$cnt$genes$ENTREZID)
  gene.set.indices <- gene.set.indices[sapply(gene.set.indices, function(x) length(x) > 4)]

  cam.res <- lapply(1:ncol(lim.res$contr.matrix), function(x) camera(lim.res$v, gene.set.indices, lim.res$v$design, contrast=lim.res$contr.matrix[ , x]))

  names(cam.res) <- tests

  lapply(names(cam.res), function(x) write.table(cam.res[[x]], paste0("camera.kegg.", celltype, ".", x, ".txt"), sep="\t", quote=F, col.names=NA))

  return(cam.res)
}

# camera.kegg.by.cell.type <- mapply(lim.res=limma.results.by.cell.type, celltype=names(limma.results.by.cell.type), runCamera, SIMPLIFY = FALSE, USE.NAMES = TRUE)


# ROAST
runRoast <- function(lim.res, celltype, annotations){
  gene.set.indices <- ids2indices(annotations, lim.res$cnt$genes$ENTREZID)
  gene.set.indices <- gene.set.indices[sapply(gene.set.indices, function(x) length(x) > 4)]

  roa.res <- lapply(1:ncol(lim.res$contr.matrix), function(x) mroast(lim.res$v, gene.set.indices, lim.res$v$design, contrast=lim.res$contr.matrix[ , x], nrot=9999, sort = "mixed"))

  names(roa.res) <- tests

  lapply(names(roa.res), function(x) write.table(roa.res[[x]], paste0("roast.kegg.", celltype, ".", x, ".txt"), sep="\t", quote=F, col.names=NA))

  return(roa.res)
}

# roast.kegg.by.cell.type <- mapply(lim.res=limma.results.by.cell.type, celltype=names(limma.results.by.cell.type), runRoastKEGG, SIMPLIFY = FALSE, USE.NAMES = TRUE)
