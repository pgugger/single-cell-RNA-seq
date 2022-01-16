# single-cell-RNA-seq
Functions for single-cell and single-nucleus RNA-Seq analyses.

pseudobulk: Generate pseudobulk data per cell type or cluster from a Seurat object

runLimmaPseudobulked: Differential expression by cell type or cluster using limma with pseudobulk data 

runLimmaDream: Differential expression by cell type or cluster using limma with dream to enable mixed model (1|sample)

runRoast: wrapper for ROAST pathway analysis using output from runLimmaPseudobulked

runCamera: wrapper for Camera pathway analysis using output from runLimmaPseudobulked


Coming soon: 
1) wrapper function to import CellRanger output and organize for Seurat
2) Example workflow
