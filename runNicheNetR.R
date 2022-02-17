# Run NicheNetR for each cell type with DEGs as receiver, all others in object as senders
# Produces numerous summary plots and tables
# Modified from NicheNetR developers' tutorial
# Assumes mouse
# Assumes Seurat object has "integrated" assay

library(nichenetr)
library(tidyverse)
library(circlize)
library(ggplot2)

# ligand, receptor, target database
ligand_target_matrix <- readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))

lr_network <- readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))

weighted_networks <- readRDS(url("https://zenodo.org/record/3260758/files/weighted_networks.rds"))

weighted_networks_lr <- weighted_networks$lr_sig %>% inner_join(lr_network %>% distinct(from,to), by = c("from","to"))

# Human to mouse by orthology
lr_network <- lr_network %>% mutate(from = convert_human_to_mouse_symbols(from), to = convert_human_to_mouse_symbols(to)) %>% drop_na()
colnames(ligand_target_matrix) = ligand_target_matrix %>% colnames() %>% convert_human_to_mouse_symbols()
rownames(ligand_target_matrix) = ligand_target_matrix %>% rownames() %>% convert_human_to_mouse_symbols()

ligand_target_matrix <- ligand_target_matrix %>% .[!is.na(rownames(ligand_target_matrix)), !is.na(colnames(ligand_target_matrix))]

weighted_networks_lr <- weighted_networks_lr %>% mutate(from = convert_human_to_mouse_symbols(from), to = convert_human_to_mouse_symbols(to)) %>% drop_na()


# Core function

runNicheNetR <- function(receiver.cell, sender.cells, lim.res, lfc, contr, seurat) {

  expressed_genes_receiver <- get_expressed_genes(receiver.cell, seurat, pct = 0.10, assay_oi="integrated")

  background_expressed_genes <- expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]

  ## sender
  sender_celltypes <- sender.cells

  list_expressed_genes_sender <- sender_celltypes %>% unique() %>% lapply(get_expressed_genes, seurat, 0.10, assay_oi="integrated") # lapply to get the expressed genes of every sender cell type separately here
  names(list_expressed_genes_sender) <- sender_celltypes
  expressed_genes_sender <- list_expressed_genes_sender %>% unlist() %>% unique()

  # Genes of interest (DEGs)
  geneset_oi <- rownames(lim.res[[receiver.cell]][[contr]])[lim.res[[receiver.cell]][[contr]]$adj.P.Val < 0.1]

  # Define expressed ligands and receptors
  ligands = lr_network %>% pull(from) %>% unique()
  receptors = lr_network %>% pull(to) %>% unique()

  #expressed_ligands_by_cell <- lapply(list_expressed_genes_sender, function(x) intersect(ligands, x))

  expressed_ligands = intersect(ligands, expressed_genes_sender)
  expressed_receptors = intersect(receptors, expressed_genes_receiver)

  potential_ligands = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()

  # Perform analysis
  ligand_activities = predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)

  ligand_activities = ligand_activities %>% dplyr::arrange(-pearson) %>% dplyr::mutate(rank = rank(desc(pearson)))
  #ligand_activities
  write.table(ligand_activities, paste0("ligand.activities.", receiver.cell, ".txt"), quote=F, row.names = F, sep="\t")

  # Table of logFC for genes by cell type...
  DE_table_all <- lapply(sender_celltypes, function(x) data.frame(gene = rownames(lfc[[x]]), logFC = lfc[[x]]$logFC) )

  names(DE_table_all) <- sender_celltypes

  DE_table_all <- Reduce(function(x, y) merge(x, y, by = "gene", all = TRUE), DE_table_all)

  rownames(DE_table_all) <- DE_table_all$gene
  colnames(DE_table_all) <- c("gene", sender_celltypes)
  DE_table_all[is.na(DE_table_all)] = 0

  DE_table_all_sort <- DE_table_all[ligand_activities$test_ligand, sender_celltypes]

  best_upstream_ligands = ligand_activities %>% top_n(20, pearson) %>% arrange(-pearson) %>% pull(test_ligand) %>% unique()

  dotplot = DotPlot(seurat, features = rev(best_upstream_ligands), assay = "integrated", cols = c("gray90", "darkblue")) + coord_flip() + theme(legend.text = element_text(size = 10), legend.title = element_text(size = 12))
  dotplot + theme(legend.position = "none", axis.ticks = element_blank(), axis.title.x = element_text(size = 12), axis.text.y = element_text(face = "italic", size = 9), axis.text.x = element_text(size = 9,  angle = 90,hjust = 0)) + ylab("Expression in Sender") + xlab("") + scale_y_discrete(position = "right")
  ggsave(paste0("best.upstream.ligands.dotplot.", receiver.cell, ".pdf"), plot = last_plot(), device = "pdf",  width = 8, height = 6)

  # Infer receptor and target of top ligands
  # Targets
  active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links, geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 250) %>% bind_rows() %>% drop_na()

  active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.33)

  order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev() %>% make.names()
  order_targets = active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links)) %>% make.names()
  rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
  colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23

  write.table(active_ligand_target_links, paste0("ligand.target.links.", receiver.cell, ".txt"), quote=F, col.names = NA, sep="\t")

  vis_ligand_target = active_ligand_target_links[order_targets, order_ligands,  drop=F] %>% t()
  p_ligand_target_network = vis_ligand_target %>% make_heatmap_ggplot("Prioritized ligands","Predicted target genes", color = "purple", legend_position = "top", x_axis_position = "top", legend_title = "Regulatory potential") + theme(axis.text.x = element_text(face = "italic"), axis.text.y = element_text(face = "italic")) + scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.0045,0.0090))
  p_ligand_target_network

  ggsave(paste0("ligand.target.heatmap.", receiver.cell, ".pdf"), plot = last_plot(), device = "pdf",  width = 10, height = 5)


  # Receptors
  lr_network_top = lr_network %>% filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% distinct(from,to)
  best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()

  lr_network_top_df_large = weighted_networks_lr %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)

  lr_network_top_df = lr_network_top_df_large %>% spread("from","weight",fill = 0)
  lr_network_top_matrix = lr_network_top_df %>% dplyr::select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)

  write.table(lr_network_top_df, paste0("ligand.receptor.", receiver.cell, ".txt"), quote=F, row.names = F, sep="\t")

  dist_receptors = dist(lr_network_top_matrix, method = "binary")
  hclust_receptors = hclust(dist_receptors, method = "ward.D2")
  order_receptors = hclust_receptors$labels[hclust_receptors$order]

  dist_ligands = dist(lr_network_top_matrix %>% t(), method = "binary")
  hclust_ligands = hclust(dist_ligands, method = "ward.D2")
  order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]

  order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix))
  order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix))

  vis_ligand_receptor_network = lr_network_top_matrix[order_receptors, order_ligands_receptor]
  rownames(vis_ligand_receptor_network) = order_receptors %>% make.names()
  colnames(vis_ligand_receptor_network) = order_ligands_receptor %>% make.names()
  p_ligand_receptor_network = vis_ligand_receptor_network %>% t() %>% make_heatmap_ggplot("Ligands", "Receptors", color = "mediumvioletred", x_axis_position = "top", legend_title = "Prior interaction potential") + theme(axis.text.x = element_text(face = "italic"), axis.text.y = element_text(face = "italic"))
  p_ligand_receptor_network

  ggsave(paste0("ligand.receptor.heatmap.", receiver.cell, ".pdf"), plot = last_plot(), device = "pdf",  width = 20, height = 5)

  # Receptors of top-ranked ligands, but after considering only bona fide ligand-receptor interactions documented in literature and publicly available databases

  lr_network_strict = lr_network %>% filter(database != "ppi_prediction_go" & database != "ppi_prediction")
  ligands_bona_fide = lr_network_strict %>% pull(from) %>% unique()
  receptors_bona_fide = lr_network_strict %>% pull(to) %>% unique()

  lr_network_top_df_large_strict = lr_network_top_df_large %>% distinct(from,to) %>% inner_join(lr_network_strict, by = c("from","to")) %>% distinct(from,to)
  lr_network_top_df_large_strict = lr_network_top_df_large_strict %>% inner_join(lr_network_top_df_large, by = c("from","to"))

  lr_network_top_df_strict = lr_network_top_df_large_strict %>% spread("from","weight",fill = 0)
  lr_network_top_matrix_strict = lr_network_top_df_strict %>% dplyr::select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df_strict$to)

  write.table(lr_network_top_df_strict, paste0("ligand.receptor.bona.fide.", receiver.cell, ".txt"), quote=F, row.names = F, sep="\t")

  dist_receptors = dist(lr_network_top_matrix_strict, method = "binary")
  hclust_receptors = hclust(dist_receptors, method = "ward.D2")
  order_receptors = hclust_receptors$labels[hclust_receptors$order]

  dist_ligands = dist(lr_network_top_matrix_strict %>% t(), method = "binary")
  hclust_ligands = hclust(dist_ligands, method = "ward.D2")
  order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]

  order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix_strict))
  order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix_strict))

  vis_ligand_receptor_network_strict = lr_network_top_matrix_strict[order_receptors, order_ligands_receptor]
  rownames(vis_ligand_receptor_network_strict) = order_receptors %>% make.names()
  colnames(vis_ligand_receptor_network_strict) = order_ligands_receptor %>% make.names()
  p_ligand_receptor_network_strict = vis_ligand_receptor_network_strict %>% t() %>% make_heatmap_ggplot("Ligands","Receptors", color = "mediumvioletred", x_axis_position = "top", legend_title = "Prior interaction potential\n(bona fide)") + theme(axis.text.x = element_text(face = "italic"), axis.text.y = element_text(face = "italic"))
  p_ligand_receptor_network_strict

  ggsave(paste0("ligand.receptor.bona.fide.", receiver.cell, ".pdf"), plot = last_plot(), device = "pdf",  width = 6, height = 5)

  # Show logFC
  # Combine ligand activities with DE information
  ligand_activities_de = ligand_activities %>% dplyr::select(test_ligand, pearson) %>% dplyr::rename(ligand = test_ligand) %>% left_join(DE_table_all %>% dplyr::rename(ligand = gene))
  ligand_activities_de[is.na(ligand_activities_de)] = 0

  write.table(ligand_activities_de, paste0("ligand.activity.logfc.", receiver.cell, ".txt"), quote=F, row.names = F, sep="\t")

  # make LFC heatmap
    lfc_matrix = ligand_activities_de  %>% dplyr::select(-ligand, -pearson) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities_de$ligand)
  rownames(lfc_matrix) = rownames(lfc_matrix) %>% make.names()

  #order_ligands = order_ligands[order_ligands %in% rownames(lfc_matrix)]
  vis_ligand_lfc = lfc_matrix[order_ligands, ]

  colnames(vis_ligand_lfc) = vis_ligand_lfc %>% colnames() %>% make.names()

  p_ligand_lfc = vis_ligand_lfc %>% make_threecolor_heatmap_ggplot("Prioritized ligands", paste0("LogFC ", contr, " Expression in Sender"), low_color = "midnightblue", mid_color = "white", mid = median(vis_ligand_lfc), high_color = "red",legend_position = "top", x_axis_position = "top", legend_title = "logFC") + theme(axis.text.y = element_text(face = "italic"))
  p_ligand_lfc

  ggsave(paste0("ligand.logfc.", receiver.cell, ".pdf"), plot = last_plot(), device = "pdf",  width = 5, height = 5)

  # Most useful summary plots together
  # ligand activity heatmap
  ligand_pearson_matrix = ligand_activities %>% dplyr::select(pearson) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities$test_ligand)

  rownames(ligand_pearson_matrix) = rownames(ligand_pearson_matrix) %>% make.names()
  colnames(ligand_pearson_matrix) = colnames(ligand_pearson_matrix) %>% make.names()

  vis_ligand_pearson = ligand_pearson_matrix[order_ligands, ] %>% as.matrix(ncol = 1) %>% magrittr::set_colnames("Pearson")
  p_ligand_pearson = vis_ligand_pearson %>% make_heatmap_ggplot("Prioritized ligands","Ligand activity", color = "darkorange",legend_position = "top", x_axis_position = "top", legend_title = "Pearson correlation coefficient\ntarget gene prediction ability)") + theme(axis.text.y = element_text(face = "italic"), legend.text = element_text(size = 9))

  # ligand expression Seurat dotplot
  order_ligands_adapted = order_ligands
  order_ligands_adapted[order_ligands_adapted == "H2.M3"] = "H2-M3"
  order_ligands_adapted[order_ligands_adapted == "H2.T23"] = "H2-T23"

  rotated_dotplot = DotPlot(seurat %>% subset(cell.type.refined.final %in% sender_celltypes), features = order_ligands_adapted , cols = c("gray90", "darkblue")) + coord_flip() + theme(legend.text = element_text(size = 10), legend.title = element_text(size = 12)) # flip of coordinates necessary because we want to show ligands in the rows when combining all plots

  vis_ligand_receptor_network_for_summary = lr_network_top_matrix[order_receptors, order_ligands]
  rownames(vis_ligand_receptor_network_for_summary) = order_receptors %>% make.names()
  colnames(vis_ligand_receptor_network_for_summary) = order_ligands %>% make.names()
  p_ligand_receptor_network_for_summary = vis_ligand_receptor_network_for_summary %>% t() %>% make_heatmap_ggplot("Ligands", "Receptors", color = "mediumvioletred", x_axis_position = "top", legend_title = "Prior interaction potential") + theme(axis.text.x = element_text(face = "italic"), axis.text.y = element_text(face = "italic"))

  figures_without_legend = cowplot::plot_grid(
    p_ligand_pearson + theme(legend.position = "none", axis.ticks = element_blank()) + theme(axis.title.x = element_text()),
    rotated_dotplot + theme(legend.position = "none", axis.ticks = element_blank(), axis.title.x = element_text(size = 12), axis.text.y = element_text(face = "italic", size = 9), axis.text.x = element_text(size = 9,  angle = 90,hjust = 0)) + ylab("Expression in Sender") + xlab("") + scale_y_discrete(position = "right"),
    p_ligand_lfc + theme(legend.position = "none", axis.ticks = element_blank()) + theme(axis.title.x = element_text()) + ylab(""),
    p_ligand_target_network + theme(legend.position = "none", axis.ticks = element_blank()) + ylab(""),
    p_ligand_receptor_network_for_summary + theme(legend.position = "none", axis.ticks = element_blank()) + ylab(""),
    align = "hv",
    nrow = 1,
    rel_widths = c(ncol(vis_ligand_pearson) + 3, ncol(vis_ligand_lfc), ncol(vis_ligand_lfc),  ncol(vis_ligand_target) + 1, nrow(vis_ligand_receptor_network_for_summary))
    )
  legends = cowplot::plot_grid(
    ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_pearson)),
    ggpubr::as_ggplot(ggpubr::get_legend(rotated_dotplot)),
    ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_lfc)),
    ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_target_network)),
    ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_receptor_network_for_summary)),
    nrow = 1,
    align = "h", rel_widths = c(1.5, 1, 1, 1, 1))

  combined_plot = cowplot::plot_grid(figures_without_legend, legends, rel_heights = c(10,5), nrow = 2, align = "hv")
  combined_plot

  ggsave(paste0("summary.", receiver.cell, ".pdf"), plot = last_plot(), device = "pdf",  width = 35, height = 7)

  return(list(sender_celltypes=sender_celltypes, geneset_oi=geneset_oi, background_expressed_genes=background_expressed_genes, ligand_activities=ligand_activities, best_upstream_ligands=best_upstream_ligands, best_upstream_receptors=best_upstream_receptors))
}

# Example
# Organize data
# degs.KO_v_Control <- lapply(limma.results, function(x) rownames(x[["KO_v_Control"]])[which(x[["KO_v_Control"]]["adj.P.Val"] < 0.1)])
# receivers <- names(degs.KO_v_Control)[which(lengths(degs.KO_v_Control) > 0)]
# senders <- unique(seurat.object$cell.type)
# logfc is simply a table with logFC for all genes, even those not tested in limma

# Run
# nichenetr.KO_v_Control <- lapply(receivers, function(x) tryCatch(runNicheNetR(receiver.cell=x, sender.cells=senders, lim.res=limma.results[senders], lfc=logfc, contr="KO_v_Control", seurat=seurat.object), error=function(e) NULL))
# names(nichenetr.LysMCre_v_Control) <- receivers
