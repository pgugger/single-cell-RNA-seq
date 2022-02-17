# Function to make circos plots from results of runNicheNetR()
# First pass visualization; usually further customization needed for publication

makeCircos <- function(receiver.cell, best_upstream_ligands, geneset_oi, sender_celltypes, best_upstream_receptors){
  # Assign ligands to max expressing cell type
  # Alternative would be to use logFC or cell proportion or more complex metric used in tutorial
  ligand_type_indication_df <- tibble(
    ligand_type = colnames(cluster.averages)[apply(cluster.averages[best_upstream_ligands,], 1, which.max)],
    ligand = best_upstream_ligands)

  active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links, geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 250) %>% bind_rows()

  active_ligand_target_links_df = active_ligand_target_links_df %>% mutate(target_type = "deg") %>% inner_join(ligand_type_indication_df) # if you want ot make circos plots for multiple gene sets, combine the different data frames and differentiate which target belongs to which gene set via the target type

  ligand_remove <- active_ligand_target_links_df$ligand[is.na(active_ligand_target_links_df$weight)]

  ligand_type_indication_df <- ligand_type_indication_df[ifelse(ligand_type_indication_df$ligand %in% ligand_remove, FALSE, TRUE), ]
  active_ligand_target_links_df <- active_ligand_target_links_df[ifelse(active_ligand_target_links_df$ligand %in% ligand_remove, FALSE, TRUE), ]

  cutoff_include_all_ligands = active_ligand_target_links_df$weight %>% quantile(0.66, na.rm=T)

  active_ligand_target_links_df_circos = active_ligand_target_links_df %>% filter(weight > cutoff_include_all_ligands)

  #ligands_to_remove = setdiff(active_ligand_target_links_df$ligand %>% unique(), active_ligand_target_links_df_circos$ligand %>% unique())
  #targets_to_remove = setdiff(active_ligand_target_links_df$target %>% unique(), active_ligand_target_links_df_circos$target %>% unique())

  circos_links = active_ligand_target_links_df # %>% filter(!target %in% targets_to_remove &!ligand %in% ligands_to_remove)

  grid_col_ligand <- rainbow(length(unique(circos_links$ligand_type))) #sender_celltypes
  names(grid_col_ligand) <- sort(unique(circos_links$ligand_type)) #sender_celltypes

  grid_col_target =c(
    "deg" = "tomato")

  grid_col_tbl_ligand = tibble(ligand_type = grid_col_ligand %>% names(), color_ligand_type = grid_col_ligand)
  grid_col_tbl_target = tibble(target_type = grid_col_target %>% names(), color_target_type = grid_col_target)

  circos_links = circos_links %>% mutate(ligand = paste(ligand," ")) # extra space: make a difference between a gene as ligand and a gene as target!
  circos_links = circos_links %>% inner_join(grid_col_tbl_ligand) %>% inner_join(grid_col_tbl_target)
  #circos_links = circos_links[complete.cases(circos_links), ]
  links_circle = circos_links %>% dplyr::select(ligand,target, weight)


  ligand_color = circos_links %>% distinct(ligand,color_ligand_type)
  grid_ligand_color = ligand_color$color_ligand_type %>% set_names(ligand_color$ligand)
  target_color = circos_links %>% distinct(target,color_target_type)
  grid_target_color = target_color$color_target_type %>% set_names(target_color$target)

  grid_col =c(grid_ligand_color,grid_target_color)

  # give the option that links in the circos plot will be transparant ~ ligand-target potential score
  transparency = circos_links %>% mutate(weight =(weight-min(weight))/(max(weight)-min(weight))) %>% mutate(transparency = 1-weight) %>% .$transparency


  target_order = circos_links$target %>% unique()
  ligand_order = ligand_type_indication_df[order(ligand_type_indication_df$ligand_type), "ligand"] %>% pull(ligand) %>% c(paste(.," "))  %>% intersect(circos_links$ligand)
  order = c(ligand_order,target_order)


  width_same_cell_same_ligand_type = 0.5
  width_different_cell = 6
  width_ligand_target = 15
  width_same_cell_same_target_type = 0.5

  ligand_type_table <- table(ligand_type_indication_df$ligand_type) - 1

  gaps = c(
    head(as.vector(unlist(lapply(ligand_type_table, function(x) c(rep( width_same_cell_same_ligand_type, x ), width_different_cell)))), -1),
    width_ligand_target,
    rep(width_same_cell_same_target_type, times = (circos_links %>% filter(target_type == "deg") %>% distinct(target) %>% nrow() -1)),
    width_ligand_target
  )

  sender.legend <- grid_col_tbl_ligand
  receiver.legend <- c(paste0("Target in ", receiver.cell), "tomato")


  #[ligand_type_table != 0]
  pdf(paste0("circos.ligand.target.", receiver.cell, ".pdf"), height = 15, width = 15)
  circos.par(gap.degree = gaps)
  chordDiagram(links_circle, directional = 1, order=order, link.sort = TRUE, link.decreasing = FALSE, grid.col = grid_col, transparency = transparency, diffHeight = 0.000, direction.type = c("diffHeight", "arrows"), link.arr.type = "big.arrow", link.visible = links_circle$weight >= cutoff_include_all_ligands, annotationTrack = "grid",
               preAllocateTracks = list(track.height = 0.075))
  # we go back to the first track and customize sector labels
  circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
                facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 1)
  }, bg.border = NA) #


  library(ComplexHeatmap)
  sender.leg <- Legend(at = sender.legend$ligand_type, type = "points",
         legend_gp = gpar(col = sender.legend$color_ligand_type), background=sender.legend$color_ligand_type, title_position = "topleft",
         title = "Sender", labels_gp = gpar(fontsize = 14), title_gp = gpar(fontsize = 16, fontface = "bold"))

  receiver.leg <- Legend(at = receiver.legend[[1]], type = "points",
         legend_gp = gpar(col = receiver.legend[[2]]), background=receiver.legend[[2]], title_position = "topleft",
         title = "Receiver", labels_gp = gpar(fontsize = 14), title_gp = gpar(fontsize = 16, fontface = "bold"))

  lgd_list_vertical = packLegend(sender.leg, receiver.leg)

  draw(lgd_list_vertical, x = unit(4, "mm"), y = unit(4, "mm"), just = c("left", "bottom"))

  circos.clear()
  dev.off()



  #grid_col_ligand <- rainbow(length(sender_celltypes))
  #names(grid_col_ligand) <- sender_celltypes
  lr_network_top_df = weighted_networks_lr %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors) %>% dplyr::rename(ligand = from, receptor = to)

  lr_network_top_df = lr_network_top_df %>% mutate(receptor_type = "deg_receptor") %>% inner_join(ligand_type_indication_df)

  grid_col_receptor =c(
    "deg_receptor" = "darkred")

  grid_col_tbl_ligand = tibble(ligand_type = grid_col_ligand %>% names(), color_ligand_type = grid_col_ligand)
  grid_col_tbl_receptor = tibble(receptor_type = grid_col_receptor %>% names(), color_receptor_type = grid_col_receptor)

  circos_links = lr_network_top_df %>% mutate(ligand = paste(ligand," ")) # extra space: make a difference between a gene as ligand and a gene as receptor!
  circos_links = circos_links %>% inner_join(grid_col_tbl_ligand) %>% inner_join(grid_col_tbl_receptor)
  links_circle = circos_links %>% dplyr::select(ligand,receptor, weight)

  ligand_color = circos_links %>% distinct(ligand,color_ligand_type)
  grid_ligand_color = ligand_color$color_ligand_type %>% set_names(ligand_color$ligand)
  receptor_color = circos_links %>% distinct(receptor,color_receptor_type)
  grid_receptor_color = receptor_color$color_receptor_type %>% set_names(receptor_color$receptor)

  grid_col = c(grid_ligand_color,grid_receptor_color)

  # give the option that links in the circos plot will be transparent ~ ligand-receptor potential score
  transparency = circos_links %>% mutate(weight =(weight-min(weight))/(max(weight)-min(weight))) %>% mutate(transparency = 1-weight) %>% .$transparency
  #Prepare the circos visualization: order ligands and receptors

  receptor_order = circos_links$receptor %>% unique()
  #ligand_order = c(CAF_specific_ligands,general_ligands,endothelial_specific_ligands) %>% c(paste(.," ")) %>% intersect(circos_links$ligand)
  order = c(ligand_order, receptor_order)
  #Prepare the circos visualization: define the gaps between the different segments

  width_same_cell_same_ligand_type = 0.5
  width_different_cell = 6
  width_ligand_receptor = 15
  width_same_cell_same_receptor_type = 0.5

  gaps = c(
    head(as.vector(unlist(lapply(ligand_type_table, function(x) c(rep( width_same_cell_same_ligand_type, x ), width_different_cell)))), -1),
    width_ligand_receptor,
    rep(width_same_cell_same_receptor_type, times = (circos_links %>% filter(receptor_type == "deg_receptor") %>% distinct(receptor) %>% nrow() -1)),
    width_ligand_receptor
  )

  sender.legend <- grid_col_tbl_ligand
  receiver.legend <- c(paste0("Receptor in ", receiver.cell), "darkred")


  pdf(paste0("circos.ligand.receptor.", receiver.cell, ".pdf"), height = 15, width = 15)
  circos.par(gap.degree = gaps)
  chordDiagram(links_circle, directional = 1, order=order, link.sort = TRUE, link.decreasing = FALSE, grid.col = grid_col, transparency = transparency, diffHeight = 0.000, direction.type = c("diffHeight", "arrows"),link.arr.type = "big.arrow", link.visible = links_circle$weight >= cutoff_include_all_ligands,annotationTrack = "grid",
               preAllocateTracks = list(track.height = 0.075))
  # we go back to the first track and customize sector labels
  circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
                facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 0.8)
  }, bg.border = NA)

  library(ComplexHeatmap)
  sender.leg <- Legend(at = sender.legend$ligand_type, type = "points",
                       legend_gp = gpar(col = sender.legend$color_ligand_type), background=sender.legend$color_ligand_type, title_position = "topleft",
                       title = "Sender", labels_gp = gpar(fontsize = 14), title_gp = gpar(fontsize = 16, fontface = "bold"))

  receiver.leg <- Legend(at = receiver.legend[[1]], type = "points",
                         legend_gp = gpar(col = receiver.legend[[2]]), background=receiver.legend[[2]], title_position = "topleft",
                         title = "Receiver", labels_gp = gpar(fontsize = 14), title_gp = gpar(fontsize = 16, fontface = "bold"))

  lgd_list_vertical = packLegend(sender.leg, receiver.leg)

  draw(lgd_list_vertical, x = unit(4, "mm"), y = unit(4, "mm"), just = c("left", "bottom"))

  circos.clear()
  dev.off()
}

# Example
# Generate cluster expression averages
cluster.averages <- AverageExpression(seurat.object, assays = "integrated")
cluster.averages <- cluster.averages[["integrated"]][, senders]

# Remove NULL runs from runNicheNetR() object
nichenetr.KO_v_Control <- nichenetr.KO_v_Control[unlist(lapply(nichenetr.KO_v_Control, function(x) !is.null(x)))]

# Make Circos plots
# lapply(names(nichenetr.KO_v_Control), function(x) tryCatch(makeCircos(x, nichenetr.KO_v_Control[[x]]$best_upstream_ligands, nichenetr.KO_v_Control[[x]]$geneset_oi, nichenetr.KO_v_Control[[x]]$sender_celltypes, nichenetr.KO_v_Control[[x]]$best_upstream_receptors), error=function(e) NULL) )
