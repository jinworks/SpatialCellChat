#' @title Spatial Cell-Cell Communication Distance Plot
#' @description
#' SpatialCellChat uses cell-cell distances as communication constraint. Use this
#' function to visualize reasonable cell-cell communication distances in the
#' spatial transcriptomics data for both "Secreted" signalings and "Contact" signalings.
#'
#' @param object SpatialCellChat object.
#' @param signaling.type Character. One of c("both","Secreted","Contact").
#' @param interaction.range Numeric. By default 250um. The maximum interaction/diffusion length of ligands (Unit: microns).
#' This hard threshold is used to filter out the connections between spatially distant individual cells.
#' @param contact.range Numeric. By default 10um. The interaction range (Unit: microns) to restrict the contact-dependent signaling.
#' For spatial transcriptomics in a single-cell resolution, `contact.range` is approximately equal to the estimated
#' cell diameter (i.e., the cell center-to-center distance), which means that contact-dependent and juxtacrine
#' signaling can only happens when the two cells are contact to each other.
#' @param tol NULL or Numeric. set a distance tolerance when computing cell-cell distances and contact adjacent matrix.
#' By default `tol = NULL` means distance tolerance equals to spot.size/2 or cell.diameter/2.
#' @param density.alpha Numeric. Control the transparency of the plot.
#'
#' @return ggplot
#' @export
spatialCCCDistPlot <- function(
    object,
    signaling.type=c("both","Secreted","Contact"),
    interaction.range = 250,
    contact.range = 10,
    tol = NULL, # will be removed in the future
    density.alpha=0.5
){
  signaling.type <- match.arg(signaling.type)
  if(is.null(object@images$result.computeCellDistance)){
    data.spatial <- object@images$coordinates %>% BiocGenerics::as.data.frame()
    ratio <- object@images$spatial.factors[["ratio"]]
    if(is.null(tol)) tol <- object@images$spatial.factors[["tol"]] else NULL

    # compute the cell-to-cell distances
    # res is a list object, containing `d.spatial` matrix and `adj.contact` matrix!
    res <- computeCellDistance(coordinates = data.spatial,
                               ratio = ratio,
                               interaction.range = interaction.range,
                               contact.range = contact.range,
                               tol = tol)
  } else {
    res <- object@images$result.computeCellDistance
  }

  # long-range distance
  d.spatial <- res$d.spatial
  # short-range distance adjacent matrix for contact-dependent and juxtacrine signaling
  adj.contact <- d.spatial*res$adj.contact

  df.Secreted <- data.frame(
    x = d.spatial@x
  )
  df.Contact <- data.frame(
    x = adj.contact@x
  )

  my.theme <- theme_bw() +theme(
    # axis.ticks = element_blank(),
    # axis.text = element_blank(),
    axis.title.y = element_blank(),
    # plot.background = element_rect(linetype = "transparent")
  )

  p.Secreted <- ggplot(df.Secreted) +
    geom_density( aes(x = x, y = ..density..), fill="#69b3a2",alpha= density.alpha)+
    scale_y_continuous(expand = expansion(c(0, 0)), breaks = NULL, labels = NULL) +
    scale_x_continuous(limits = c(0, (interaction.range+10)),breaks = seq(0, (interaction.range+10), by = 20))+
    xlab("interaction range (um)")+my.theme

  p.Contact <- ggplot(df.Contact) +
    geom_density( aes(x = x, y = ..density..),fill= "#404080",alpha=density.alpha)+
    scale_y_continuous(expand = expansion(c(0, 0)), breaks = NULL, labels = NULL) +
    scale_x_continuous(limits = c(0, (interaction.range+10)),breaks = seq(0, (interaction.range+10), by = 20))+
    xlab("contact range (um)")+my.theme

  if(signaling.type == "both"){
    p <- patchwork::wrap_plots(p.Secreted,p.Contact,nrow = 2,ncol = 1)
    return(p)
  } else if(signaling.type == "Secreted"){
    return(p.Secreted)
  } else if(signaling.type == "Contact"){
    return(p.Contact)
  }
}

#' @title Visualize spatial cell groups such as signaling sources and targets
#'
#' @description
#' This function takes a CellChat object as input, and then plot cell groups of interest.
#'
#' @param object cellchat object
#' @param color.use defining the color for each cell group
#' @param group.by Name of one metadata columns to group (color) cells. Default is the defined cell groups in CellChat object
# #' @param label whether labeling cell groups
# #' @param label.size Sets the size of the labels
#' @param sources.use a vector giving the index or the name of source cell groups
#' @param targets.use a vector giving the index or the name of target cell groups
#' @param idents.use a vector giving the index or the name of cell groups of interest
#'
#' @param proportion data frame of cell type proportion (rows are cells and columns are cell types)
#' @param radius dumeric value about the radius of each pie chart, if NULL, we will calculate it inside the function
#' @param radius.size numeric scaling factor for the pie chart radius
#' @param alpha the transparency of individual spot
#' @param shape.by the shape of individual spot
#' @param title.name title name
#' @param point.size the size of spots
#' @param legend.size the size of legend
#' @param legend.text.size the text size on the legend
#' @param legend.position,ncol,byrow parameters for configurating the plot
#' @param legend.spacing a two-elements vector respectively specifying legend.key.spacing.x and legend.key.spacing.y for spacing apart legend key-label pairs
#' @return
#' @export
spatialDimPlot <- function(object, color.use = NULL, group.by = NULL, sources.use = NULL, targets.use = NULL, idents.use = NULL,
                           proportion = NULL, radius = NULL, radius.size = 0.55,
                           alpha = 1, shape.by = 16, title.name = NULL, point.size = 2.4,
                           legend.size = 3, legend.text.size = 8, legend.position = "right", legend.spacing = c(0, -8), ncol = 1, byrow = FALSE){
  coordinates <- object@images$coordinates
  if (NCOL(coordinates) == 2) {
    colnames(coordinates) <- c("x_cent","y_cent")
    temp_coordinates = coordinates
    coordinates[,1] = temp_coordinates[,2]
    coordinates[,2] = temp_coordinates[,1]
  } else {
    stop("Please check the input 'coordinates' and make sure it is a two column matrix.")
  }

  if (is.null(group.by)) {
    labels <- object@idents
  } else {
    labels = object@meta[,group.by]
    # avoid loss of factor levels
    if(!is.factor(labels)) labels <- factor(labels)
  }
  cells.level <- levels(labels)

  if (!is.null(idents.use)) {
    if (is.numeric(idents.use)) {
      idents.use <- cells.level[idents.use]
    }
    group <- rep("Others", length(labels))
    group[(labels %in% idents.use)] <- idents.use
    group = factor(group, levels = c(idents.use, "Others"))

    if (is.null(color.use)) {
      color.use.all <- scPalette(nlevels(labels))

      # get the first 2 colors, then convert the 3th color to "grey90"
      color.use <- color.use.all[match(c(idents.use), levels(labels))]

      color.use[nlevels(group)] <- "grey90"
      names(color.use) <- c(idents.use, "Others")
    }
    labels <- group
  }

  if (is.null(sources.use) & is.null(targets.use)){
    if (is.null(color.use)) {
      color.use <- scPalette(nlevels(labels))
      # Assign Colors by Factor in ggplot2
      # https://statisticalpoint.com/color-by-factor-ggplot2/#:~:text=How%20to%20Assign%20Colors%20by%20Factor%20in%20ggplot2,the%20following%20syntax%3A%20ggplot%28df%2C%20aes%28x%3Dx_variable%2C%20y%3Dy_variable%2C%20color%3Dcolor_variable%29%29%20%2B
      names(color.use) <- levels(labels)
    }
  } else {
    if (is.numeric(sources.use)) {
      sources.use <- cells.level[sources.use]
    }
    if (is.numeric(targets.use)) {
      targets.use <- cells.level[targets.use]
    }

    group <- rep("Others", length(labels))
    group[(labels %in% sources.use)] <- sources.use
    group[(labels %in% targets.use)] <- targets.use
    group = factor(group, levels = c(sources.use, targets.use, "Others"))

    if (is.null(color.use)) {
      color.use.all <- scPalette(nlevels(labels))
      color.use <- color.use.all[match(c(sources.use, targets.use), levels(labels))]
      color.use[nlevels(group)] <- "grey90"
      names(color.use) <- c(sources.use, targets.use, "Others")
    }
    labels <- group
  }



  if (is.null(proportion)) {
    gg <- ggplot(data = coordinates,aes(x=x_cent,y=y_cent,colour = labels))+
      geom_point(alpha = alpha, size = point.size, shape=shape.by) +
      scale_color_manual(values = color.use, na.value = "grey90") +
      guides(color = guide_legend(override.aes = list(size=legend.size), ncol = ncol, byrow = byrow))
  } else {
    proportion = as.data.frame(proportion)
    #proportion = proportion[,mixedsort(colnames(proportion))]
    df <- cbind(proportion,coordinates)
    ct.select = levels(labels) # colnames(proportion)
    if(is.null(radius)){
      radius = (max(coordinates$x_cent) - min(coordinates$x_cent)) * (max(coordinates$y_cent) - min(coordinates$y_cent))
      radius = radius / nrow(coordinates)
      radius = radius / pi
      radius = sqrt(radius) * radius.size
    }else{
      radius = radius
    }
    gg <- ggplot() + scatterpie::geom_scatterpie(aes(x=x_cent, y=y_cent,r = radius), data = df,
                                                 cols = ct.select, color=NA) + #+ coord_fixed(ratio = 1*max(data$x)/max(data$y)) + scale_fill_manual(values =  colors)
      scale_fill_manual(values = color.use, na.value = "grey90") +
      guides(fill = guide_legend(override.aes = list(size=legend.size), ncol = ncol, byrow = byrow))
  }
  gg <- gg + theme(legend.position = legend.position, legend.key.spacing.y = unit(legend.spacing[2], 'pt'), legend.key.spacing.x = unit(legend.spacing[1], 'pt')) +
    theme(legend.title = element_blank(), legend.text = element_text(size = legend.text.size, margin = margin(l = 0)))  + # , legend.key.size = unit(0.4, "inches")
    theme(panel.background = element_blank(),axis.ticks = element_blank(), axis.text = element_blank()) + xlab(NULL) + ylab(NULL) + theme(legend.key = element_blank()) +
    coord_fixed()+scale_y_reverse()
  #coord_fixed(ratio = 1*diff(range(coordinates$x_cent))/diff(range(coordinates$y_cent))) + scale_y_reverse()

  if (!is.null(title.name)){
    gg <- gg + ggtitle(title.name) + theme(plot.title = element_text(hjust = 0.5, vjust = 0, size = 10))
  }
  return(gg)

}


#' A spatially feature plots
#'
#' This function takes a CellChat object as input, and then plot gene expression distribution.
#'
#' @param object cellchat object
#' @param features a char vector containing features to visualize. `features` can be genes or column names of `object@meta`.
#' @param signaling signalling names to visualize
#' @param pairLR.use a data frame consisting of one column named "interaction_name", defining the L-R pairs of interest
#' @param slot.name Character. "data.signaling" or "data". By default slot.name = "data.signaling".
#' @param enriched.only  whether only return the identified enriched signaling genes in the database. Default = TRUE, returning the significantly enriched signaling interactions
#' @param do.group set `do.group = TRUE` when only showing enriched signaling based on cell group-level communication; set `do.group = FALSE` when only showing enriched signaling based on individual cell-level communication
#' @param thresh threshold of the p-value for determining significant interaction when visualizing links at the level of ligands/receptors;
#' @param color.heatmap A character string or vector indicating the colormap option to use. It can be the avaibale color palette in brewer.pal() or viridis_pal() (e.g., "Spectral","viridis")
#' @param n.colors,direction n.colors: number of basic colors to generate from color palette; direction: Sets the order of colors in the scale. If 1, the default colors are used. If -1, the order of colors is reversed.
#' @param do.binary,cutoff whether binarizing the expression using a given cutoff
#' @param color.use defining the color for cells/spots expressing ligand only, expressing receptor only, expressing both ligand & receptor and cells/spots without expression of given ligands and receptors
#' @param alpha the transparency of individual spot
#' @param point.size the size of cell slot
#' @param shape.by the shape of individual spot
#' @param legend.size the size of legend
#' @param legend.text.size the text size on the legend
#' @param ncol number of columns if plotting multiple plots
#' @param show.legend whether show each figure legend
#' @param show.legend.combined whether show the figure legend for the last plot
#' @return
#' @export
#'
spatialFeaturePlot <- function(object, features = NULL, signaling = NULL, pairLR.use = NULL,slot.name=c("data.signaling","data"),
                               enriched.only = TRUE,thresh = 0.05, do.group = TRUE,
                               color.heatmap = NULL, n.colors = 8, direction = -1,
                               do.binary = FALSE, cutoff = NULL, color.use = NULL, alpha = 1,
                               point.size = 1, legend.size = 3, legend.text.size = 8, shape.by = 16, ncol = NULL,
                               show.legend = TRUE, show.legend.combined = FALSE){
  coords <- object@images$coordinates
  if (ncol(coords) == 2) {
    colnames(coords) <- c("x_cent","y_cent")
    temp_coord = coords
    coords[,1] = temp_coord[,2]
    coords[,2] = temp_coord[,1]
  } else {
    stop("Please check the input 'coordinates' and make sure it is a two column matrix.")
  }
  slot.name <- match.arg(slot.name)
  data <- methods::slot(object,slot.name)
  meta <- object@meta

  if (is.null(color.heatmap)) {
    colormap <- grDevices::colorRampPalette(c("#FFFFEF","#FFFF3E", "#FF9F00", "#FF1A00", "#930000","#050000"))(64)
    colormap[1] <- "#E5E5E5"
  } else if (length(color.heatmap) == 1) {
    colormap <- tryCatch({
      RColorBrewer::brewer.pal(n = n.colors, name = color.heatmap)
    }, error = function(e) {
      scales::viridis_pal(option = color.heatmap, direction = -1)(n.colors)
    })
    if (direction == -1) {
      colormap <- rev(colormap)
    }
    colormap <- colorRampPalette(colormap)(99)
    colormap[1] <- "#E5E5E5"
  } else {
    colormap <- color.heatmap
  }

  # View(colormap)
  if (is.null(features) & is.null(signaling) & is.null(pairLR.use)){
    stop("Please input either features, signaling or pairLR.use.")
  }
  if (!is.null(features) & !is.null(signaling)){
    stop("Please don't input features or signaling simultaneously.")
  }
  if (!is.null(features) & !is.null(pairLR.use)){
    stop("Please don't input features or pairLR.use simultaneously.")
  }
  if (!is.null(signaling) & !is.null(pairLR.use)){
    stop("Please don't input signaling or pairLR.use simultaneously.")
  }

  df <- data.frame(x = coords[, 1], y = coords[, 2])
  if (!do.binary) {
    if (!is.null(signaling)) {
      res <- extractEnrichedLR(object, signaling = signaling, geneLR.return = TRUE, enriched.only = enriched.only, thresh = thresh, do.group = do.group)
      feature.use <- res$geneLR
    } else if (!is.null(pairLR.use)) {
      if (is.character(pairLR.use)) {
        pairLR.use <- data.frame(interaction_name = pairLR.use)
      }
      if (enriched.only) {
        if (do.group) {
          object@net$prob[object@net$pval > thresh] <- 0
          pairLR.use.name <- pairLR.use$interaction_name[pairLR.use$interaction_name %in% dimnames(object@net$prob)[[3]]]
          prob <- object@net$prob[,,pairLR.use.name, drop = FALSE]
          prob.sum <- apply(prob > 0, 3, sum)
          names(prob.sum) <- pairLR.use.name
          signaling.includes <- names(prob.sum)[prob.sum > 0]
          pairLR.use <- pairLR.use[pairLR.use$interaction_name %in% signaling.includes, , drop = FALSE]
        } else {
          pairLR.use.name <- pairLR.use$interaction_name[pairLR.use$interaction_name %in% dimnames(object@net$prob.cell)[[3]]]
          prob.cell <- object@net$prob.cell[,,pairLR.use.name, drop = FALSE]
          prob.sum <- apply(prob.cell > 0, 3, sum)
          names(prob.sum) <- pairLR.use.name
          signaling.includes <- names(prob.sum)[prob.sum > 0]
          pairLR.use <- pairLR.use[pairLR.use$interaction_name %in% signaling.includes, , drop = FALSE]
        }
        if (length(pairLR.use$interaction_name) == 0) {
          stop(paste0('There is no significant communication related with the input `pairLR.use`. Set `enriched.only = FALSE` to show non-significant signaling.'))
        }
      }
      LR.pair <- object@LR$LRsig[pairLR.use$interaction_name, c("ligand","receptor")]
      geneL <- unique(LR.pair$ligand)
      geneR <- unique(LR.pair$receptor)
      geneL <- extractGeneSubset(geneL, object@DB$complex, object@DB$geneInfo)
      geneR <- extractGeneSubset(geneR, object@DB$complex, object@DB$geneInfo)
      feature.use <- c(geneL, geneR)
    } else {
      feature.use <- features
    }
    if (length(intersect(feature.use, rownames(data))) > 0) {
      feature.use <- feature.use[feature.use %in% rownames(data)]
      data.use <- data[feature.use, , drop = FALSE]
      data.use <- as.matrix(data.use)# sparse->dense coercion

    } else if (length(intersect(feature.use, colnames(meta))) > 0) {
      feature.use <- feature.use[feature.use %in% colnames(meta)]
      data.use <- t(meta[ ,feature.use, drop = FALSE])
    } else {
      stop("Please check your input! ")
    }
    if (!is.null(cutoff)) {
      cat("Applying a cutoff of ",cutoff,"to the values...", '\n')
      data.use[data.use <= cutoff] <- NA
    }


    if (is.null(ncol)) {
      if (length(feature.use) > 9) {
        ncol <- 5
      } else {
        ncol <- min(length(feature.use), 5)
      }
    }
    numFeature = length(feature.use)
    gg <- vector("list", numFeature)
    for (i in seq_len(numFeature)) {
      feature.name <- feature.use[i]
      df$feature.data <- data.use[i, ]
      g <- ggplot(data = df, aes(x, y)) +
        geom_point(aes(colour = feature.data), alpha = alpha, size=point.size, shape=shape.by) +
        scale_colour_gradientn(colours = colormap, guide = guide_colorbar(title = NULL, ticks = T, label = T, barwidth = 0.5), na.value = "grey90") +
        theme(legend.position = "right") +
        theme(legend.title = element_blank(), legend.text = element_text(size = legend.text.size), legend.key.size = unit(0.15, "inches"))  + # , legend.key.size = unit(0.4, "inches")
        ggtitle(feature.name) + theme(plot.title = element_text(hjust = 0.5, vjust = 0, size = 10))+
        theme(panel.background = element_blank(),axis.ticks = element_blank(), axis.text = element_blank()) + xlab(NULL) + ylab(NULL) +
        theme(legend.key = element_blank())
      # g <- g + coord_fixed(ratio = 1*diff(range(coords$x_cent))/diff(range(coords$y_cent))) + scale_y_reverse()
      g <- g + coord_fixed() + scale_y_reverse()
      if (!show.legend) {
        g <- g + theme(legend.position = "none")
      }
      if (show.legend.combined & i == numFeature) {
        g <- g + theme(legend.position = "right", legend.key.height = grid::unit(0.15, "in"), legend.key.width = grid::unit(0.5, "in"), legend.title = element_blank(),legend.key = element_blank())
      }
      gg[[i]] <- g
    }
    if (ncol > 1) {
      gg <- patchwork::wrap_plots(gg, ncol = ncol)
    } else {
      gg <- gg[[1]]
    }

  } else {

    if (is.null(color.use)) {
      color.use <- ggPalette(4)
      color.use[4] <- "grey90"
    }
    color.use1 = color.use
    if (!is.null(signaling)) {
      res <- extractEnrichedLR(object, signaling = signaling, enriched.only = enriched.only, thresh = thresh, do.group = do.group)
      # gene.pair = searchPair(signaling = signaling, pairLR.use = object@LR$LRsig, key = "pathway_name", matching.exact = T, pair.only = T)
      # LR.pair <- gene.pair[res$interaction_name, c("ligand","receptor")]
      LR.pair <- object@LR$LRsig[res$interaction_name, c("ligand","receptor")]
    } else if (!is.null(pairLR.use)) {
      if (is.character(pairLR.use)) {
        pairLR.use <- data.frame(interaction_name = pairLR.use)
      }
      if (enriched.only) {
        if (do.group) {
          object@net$prob[object@net$pval > thresh] <- 0
          pairLR.use.name <- pairLR.use$interaction_name[pairLR.use$interaction_name %in% dimnames(object@net$prob)[[3]]]
          prob <- object@net$prob[,,pairLR.use.name, drop = FALSE]
          prob.sum <- apply(prob > 0, 3, sum)
          names(prob.sum) <- pairLR.use.name
          signaling.includes <- names(prob.sum)[prob.sum > 0]
          pairLR.use <- pairLR.use[pairLR.use$interaction_name %in% signaling.includes, , drop = FALSE]
        } else {
          pairLR.use.name <- pairLR.use$interaction_name[pairLR.use$interaction_name %in% dimnames(object@net$prob.cell)[[3]]]
          prob.cell <- object@net$prob.cell[,,pairLR.use.name, drop = FALSE]
          prob.sum <- apply(prob.cell > 0, 3, sum)
          names(prob.sum) <- pairLR.use.name
          signaling.includes <- names(prob.sum)[prob.sum > 0]
          pairLR.use <- pairLR.use[pairLR.use$interaction_name %in% signaling.includes, , drop = FALSE]
        }
        if (length(pairLR.use$interaction_name) == 0) {
          stop(paste0('There is no significant communication related with the input `pairLR.use`. Set `enriched.only = FALSE` to show non-significant signaling.'))
        }
      }
      LR.pair <- object@LR$LRsig[pairLR.use$interaction_name, c("ligand","receptor")]
    } else {
      stop("Please input either `pairLR.use` or `signaling` for `binary` mode!")
    }
    geneL <- as.character(LR.pair$ligand)
    geneR <- as.character(LR.pair$receptor)
    # compute the expression of ligand or receptor
    complex_input <- object@DB$complex
    data <- as.matrix(data)
    dataL <- computeExpr_LR(geneL, data, complex_input)
    dataR <- computeExpr_LR(geneR, data, complex_input)
    rownames(dataL) <- geneL; rownames(dataR) <- geneR;
    # data.use <- matrix(0, nrow = nrow(dataL)*2, ncol = ncol(dataL))
    # data.use[seq_len(nrow(data.use)) %% 2 == 1, ] <- dataL
    # data.use[seq_len(nrow(data.use)) %% 2 == 0, ] <- dataR
    # rownames(data.use)[seq_len(nrow(data.use)) %% 2 == 1] <- geneL
    # rownames(data.use)[seq_len(nrow(data.use)) %% 2 == 0] <- geneR

    feature.use <- rownames(LR.pair)
    numFeature = nrow(LR.pair)
    if (is.null(ncol)) {
      if (length(feature.use) > 9) {
        ncol <- 4
      } else {
        ncol <- min(length(feature.use), 4)
      }
    }
    if (is.null(cutoff)) {
      stop("A `cutoff` must be provided when plotting expression in binary mode! " )
    }
    gg <- vector("list", numFeature)
    for (i in seq_len(numFeature)) {
      feature.name <- feature.use[i]
      idx1 = dataL[i, ] > cutoff
      idx2 = dataR[i, ] > cutoff
      idx3 = idx1 & idx2
      group = rep("None",ncol(dataL))
      group[idx1] = geneL[i]
      group[idx2] = geneR[i]
      group[idx3] = "Both"
      group = factor(group, levels = c(geneL[i],geneR[i],"Both","None"))
      color.use <- color.use1
      names(color.use) <- c(geneL[i],geneR[i],"Both","None")

      if (length(setdiff(levels(group), unique(group))) > 0) {
        color.use <- color.use[names(color.use) %in% unique(group)]
        group = droplevels(group, exclude = setdiff(levels(group), unique(group)))
      }

      df$feature.data <- group
      g <- ggplot(data = df, aes(x, y)) +
        geom_point(aes(colour = feature.data), alpha = alpha, size=point.size, shape=shape.by) +
        scale_color_manual(values = color.use, na.value = "grey90") +
        theme(legend.position = "right") +
        theme(legend.title = element_blank(), legend.text = element_text(size = legend.text.size), legend.key.size = unit(0.15, "inches"))  + # , legend.key.size = unit(0.4, "inches")
        guides(color = guide_legend(override.aes = list(size=legend.size))) +
        ggtitle(feature.name) + theme(plot.title = element_text(hjust = 0.5, vjust = 0, size = 10))+
        theme(panel.background = element_blank(),axis.ticks = element_blank(), axis.text = element_blank()) + xlab(NULL) + ylab(NULL) +
        theme(legend.key = element_blank())
      g <- g + coord_fixed() + scale_y_reverse()
      # g <- g + coord_fixed(ratio = 1*diff(range(coords$x_cent))/diff(range(coords$y_cent))) + scale_y_reverse()
      if (!show.legend) {
        g <- g + theme(legend.position = "none")
      }
      if (show.legend.combined & i == numFeature) {
        g <- g + theme(legend.position = "right", legend.key.height = grid::unit(0.15, "in"), legend.key.width = grid::unit(0.5, "in"), legend.title = element_blank(),legend.key = element_blank())
      }
      gg[[i]] <- g
    }
    if (ncol > 1) {
      gg <- patchwork::wrap_plots(gg, ncol = ncol)
    } else {
      gg <- gg[[1]]
    }

  }
  return(gg)
}


#' Spatially visualizing any computed scores
#'
#' @param object cellchat object
#' @param signaling a signaling pathway or ligand-receptor pair to visualize
#' @param method the method for scoring
#' @param measure measurements to show. When method = "centrality", centality measurements can be c("outdeg_unweighted", "indeg_unweighted","outdeg","indeg","hub","authority","eigen","page_rank","betweenness","flowbet","info")
#' @param measure.name the measure names to show, which should correspond to the defined `measure`
#' @param merge whether to merge two scoring into one single plot
#' @param do.group set `do.group = TRUE` when computing scores of each cell group; set `do.group = FALSE` when computing scores of each individual cell
#' @param slot.name the slot name of object that is used to compute scores of signaling networks
#' @param color.heatmap A character string or vector indicating the colormap option to use. It can be the avaibale color palette in viridis_pal() or brewer.pal()
#' @param n.colors,direction n.colors: number of basic colors to generate from color palette; direction: Sets the order of colors in the scale. If 1, the default colors are used. If -1, the order of colors is reversed.
#' @param alpha the transparency of individual spot
#' @param do.binary whether to plot the scores in a binary mode
#' @param color.use colors to use when `do.binary = TRUE`
#' @param point.size the size of cell slot
#' @param shape.by the shape of individual spot
#' @param legend.size the size of legend
#' @param legend.text.size the text size on the legend
#' @param legend.spacing a two-elements vector respectively specifying legend.key.spacing.x and legend.key.spacing.y for spacing apart legend key-label pairs
#' @param ncol number of columns when plotting multiple plots
#' @param combined whether outputting an combined plot using patchwork
#' @param show.legend whether show each figure legend
#' @param show.legend.combined whether show the figure legend for the last plot
#' @return
#' @export

spatialVisual_scoring <- function(object, signaling, method = "centrality", measure = c("outdeg","indeg"), measure.name = c("outgoing", "incoming"), merge = FALSE,
                                  slot.name = "netP", do.group = FALSE,
                                  color.heatmap = "BuGn", n.colors = 9, direction = 1, alpha = 1,
                                  do.binary = FALSE, color.use = NULL,
                                  point.size = 0.8, legend.size = 3, legend.text.size = 8, legend.spacing = c(0, -8), shape.by = 16, ncol = NULL,
                                  combined = TRUE, show.legend = TRUE, show.legend.combined = FALSE){
  coords <- object@images$coordinates
  for (i in 1:length(measure)) {
    if (measure[i] == "outdeg") {
      measure.name[i] = "Outgoing"
    } else if (measure[i] == "indeg") {
      measure.name[i] = "Incoming"
    }
  }
  if (ncol(coords) == 2) {
    colnames(coords) <- c("x_cent","y_cent")
    temp_coord = coords
    coords[,1] = temp_coord[,2]
    coords[,2] = temp_coord[,1]
  } else {
    stop("Please check the input 'coordinates' and make sure it is a two column matrix.")
  }
  if (length(signaling) > 1) {
    stop("Only one signaling is allowed as an input.")
  }
  if (method == "centrality") {

    if (do.group) {
      if (length(slot(object, slot.name)$centr) == 0) {
        stop("Please run `netAnalysis_computeCentrality` to compute the network centrality scores! ")
      }
      centr <- slot(object, slot.name)$centr[measure,,signaling, drop = FALSE]
      label <- object@idents
      cell.levels <- levels(label)
      value <- matrix(NA, nrow = length(measure), ncol = nrow(coords))
      for (i in 1:length(cell.levels)) {
        value[, label == cell.levels[i]] <- matrix(rep(centr[,i,1], sum(label == cell.levels[i])), byrow = FALSE, nrow = length(measure))
      }

    } else {
      if (length(slot(object, slot.name)$centr.cell) == 0) {
        stop("Please run `netAnalysis_computeCentrality` to compute the network centrality scores! ")
      }
      centr <- slot(object, slot.name)$centr.cell[measure,,signaling]
      value <- matrix(centr, nrow = length(measure), ncol = nrow(coords))
    }
  }

  if (!merge) {
    if (!do.binary) {
      if (length(color.heatmap) == 1) {
        colormap <- tryCatch({
          RColorBrewer::brewer.pal(n = n.colors, name = color.heatmap)
        }, error = function(e) {
          scales::viridis_pal(option = color.heatmap, direction = -1)(n.colors)
        })
        if (direction == -1) {
          colormap <- rev(colormap)
        }
        colormap <- grDevices::colorRampPalette(colormap)(99)
        colormap[1] <- "#E5E5E5"
      } else {
        colormap <- color.heatmap
      }
    } else {
      if (is.null(color.use)) {
        color.use <- c("#2171b5")
      }
      color.use <- c("grey90",color.use)
    }

    numFeature = length(measure)
    if (is.null(ncol)) {
      if (numFeature > 9) {
        ncol <- 4
      } else {
        ncol <- min(numFeature, 2)
      }
    }
    df <- data.frame(x = coords[, 1], y = coords[, 2])
    gg <- vector("list", numFeature)
    for (i in seq_len(numFeature)) {
      feature.name <- paste0(measure.name[i], " scores of ", signaling, " signaling")
      if (!do.binary) {
        df$feature.data <- value[i, ]
        g <- ggplot(data = df, aes(x, y)) +
          geom_point(aes(colour = feature.data), alpha = alpha, size=point.size, shape=shape.by) +
          scale_colour_gradientn(colours = colormap, guide = guide_colorbar(title = NULL, ticks = T, label = T, barwidth = 0.5), na.value = "grey90")
        # metR::geom_contour_fill(aes(z = feature.data))

      } else {
        df$feature.data <- (value[i, ] > 0) * 1
        df$feature.data <- factor(df$feature.data, levels = c("0","1"))
        g <- ggplot(data = df, aes(x, y)) +
          geom_point(aes(colour = feature.data), alpha = alpha, size=point.size, shape=shape.by) +
          scale_color_manual(values = color.use, na.value = "grey90")
        show.legend <- FALSE; show.legend.combined <- FALSE
      }
      g <- g + theme(legend.position = "right") + theme(legend.title = element_blank(), legend.text = element_text(size = legend.text.size), legend.key.size = unit(0.15, "inches"))  + # , legend.key.size = unit(0.4, "inches")
        ggtitle(feature.name) + theme(plot.title = element_text(hjust = 0.5, vjust = 0, size = 10))+
        theme(panel.background = element_blank(),axis.ticks = element_blank(), axis.text = element_blank()) + xlab(NULL) + ylab(NULL) +
        theme(legend.key = element_blank())
      g <- g + coord_fixed() + scale_y_reverse()
      # g <- g + coord_fixed(ratio = 1*diff(range(coords$x_cent))/diff(range(coords$y_cent))) + scale_y_reverse()
      if (!show.legend) {
        g <- g + theme(legend.position = "none")
      }
      if (show.legend.combined & i == numFeature) {
        g <- g + theme(legend.position = "right", legend.key.height = grid::unit(0.15, "in"), legend.key.width = grid::unit(0.3, "in"), legend.title = element_blank(),legend.key = element_blank())
      }
      gg[[i]] <- g
    }
  } else {
    if (!do.binary) {
      if (length(color.heatmap) != 3) {
        warning("Users should provide three color schemes when merging the two measurements. Default(\"Blues\",\"Reds\",\"Purples\") will be used when no input.")
        color.heatmap <- c("Blues","Reds","Purples")
      }
      colormap.all <- list()
      for (i in 1:length(color.heatmap)) {
        colormap <- tryCatch({
          RColorBrewer::brewer.pal(n = n.colors, name = color.heatmap[i])
        }, error = function(e) {
          scales::viridis_pal(option = color.heatmap[i], direction = -1)(n.colors)
        })
        if (direction == -1) {
          colormap <- rev(colormap)
        }
        colormap <- colormap[2:length(colormap)]
        colormap <- grDevices::colorRampPalette(colormap)(99)
        # if (i != 3) {
        #   colormap[1] <- "#E5E5E5"
        # }

        colormap.all[[i]] <- colormap
      }
    } else {
      if (is.null(color.use)) {
        color.use <- c("#2171b5", "#cb181d", "#6a51a3")
        # color.use <- c("#7631F6", "#d92061","#d94801")
      }
      color.use <- c("grey90",color.use)
    }

    numFeature = length(measure)
    if (numFeature != 2) {
      stop("`merge` only works when visualizing two measurements!")
    }
    ncol <- numFeature+1
    measure <- c(measure, 'Merged')
    feature.name <- paste0(measure[3], " scores of ", signaling, " signaling")
    df <- data.frame(x = coords[, 1], y = coords[, 2], value1 = value[1, ], value2 = value[2, ])
    df0 <- subset(df, value1 == 0 & value2 == 0)
    df1 <- subset(df, (value1 > 0 & value2 == 0))
    df2 <- subset(df, (value1 == 0 & value2 > 0) )
    df3 <- subset(df, value1 > 0 & value2 > 0)
    df3$value3 = apply(df3[,c("value1", "value2")], 1, mean)
    if (!do.binary) {
      g <- ggplot(mapping = aes(x, y)) + geom_point(data = df0, colour = "#E5E5E5", size=point.size, shape=shape.by)+
        ggnewscale::new_scale_color() +
        geom_point(data = df1, aes(colour = value1), alpha = alpha, size=point.size, shape=shape.by) +
        scale_colour_gradientn(colours = colormap.all[[1]], guide = guide_colorbar(title = measure.name[1], ticks = T, label = T, barwidth = 0.3, order = 1), na.value = "grey90") +
        ggnewscale::new_scale_color() +
        geom_point(data = df2, aes(color = value2), alpha = alpha, size=point.size, shape=shape.by) +
        scale_colour_gradientn(colours = colormap.all[[2]], guide = guide_colorbar(title = measure.name[2], ticks = T, label = T, barwidth = 0.3, order = 2), na.value = "grey90") +
        ggnewscale::new_scale_color() +
        geom_point(data = df3, aes(color = value3), alpha = alpha, size=point.size, shape=shape.by) +
        scale_colour_gradientn(colours = colormap.all[[3]], guide = guide_colorbar(title = measure[3], ticks = T, label = T, barwidth = 0.3, order = 3), na.value = "grey90")

    } else {
      idx1 = df$value1 > 0 & df$value2 == 0
      idx2 = df$value1 == 0 & df$value2 > 0
      idx3 = df$value1 > 0 & df$value2 > 0
      group = rep("None",nrow(df))
      group[idx1] = measure.name[1]
      group[idx2] = measure.name[2]
      group[idx3] = "Both"
      group = factor(group, levels = c("None", measure.name[1],measure.name[2],"Both"))
      names(color.use) <- c("None", measure.name[1],measure.name[2],"Both")

      if (length(setdiff(levels(group), unique(group))) > 0) {
        color.use <- color.use[names(color.use) %in% unique(group)]
        group = droplevels(group, exclude = setdiff(levels(group), unique(group)))
      }
      df$feature.data <- group
      g <- ggplot(data = df, aes(x, y)) +
        geom_point(aes(colour = feature.data), alpha = alpha, size=point.size, shape=shape.by) +
        scale_color_manual(values = color.use, na.value = "grey90")
    }

    g <- g + theme(legend.position = "right") +
      theme(legend.text = element_text(size = legend.text.size), legend.key.size = unit(0.1, "inches"), legend.title = element_text(size = legend.text.size))  + # , legend.key.size = unit(0.4, "inches")
      ggtitle(feature.name) + theme(plot.title = element_text(hjust = 0.5, vjust = 0, size = 10))+
      theme(panel.background = element_blank(),axis.ticks = element_blank(), axis.text = element_blank()) + xlab(NULL) + ylab(NULL) +
      theme(legend.key = element_blank())
    g <- g + coord_fixed() + scale_y_reverse()
    # g <- g + coord_fixed(ratio = 1*diff(range(coords$x_cent))/diff(range(coords$y_cent))) + scale_y_reverse()
    if (!show.legend) {
      g <- g + theme(legend.position = "none")
    }
    gg <- g
    combined <- FALSE
  }


  if (combined) {
    gg <- patchwork::wrap_plots(gg, ncol = ncol)
  }
  return(gg)
}



#' ggplot theme in CellChat
#'
#' @return
#' @export
#'
#' @examples
#' @importFrom ggplot2 theme_classic element_rect theme element_blank element_line element_text
CellChat_theme_opts <- function() {
  theme(strip.background = element_rect(colour = "white", fill = "white")) +
    theme_classic() +
    theme(panel.border = element_blank()) +
    theme(axis.line.x = element_line(color = "black")) +
    theme(axis.line.y = element_line(color = "black")) +
    theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) +
    theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) +
    theme(panel.background = element_rect(fill = "white")) +
    theme(legend.key = element_blank()) + theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5))
}


#' Generate ggplot2 colors
#'
#' @param n number of colors to generate
#' @importFrom grDevices hcl
#' @export
#'
ggPalette <- function(n) {
  hues = seq(15, 375, length = n + 1)
  grDevices::hcl(h = hues, l = 65, c = 100)[1:n]
}

#' Generate colors from a customed color palette for sc visualization
#'
#' @param n number of colors
#'
#' @return A color palette for plotting
#' @importFrom grDevices colorRampPalette
#'
#' @export
#'
scPalette <- function(n) {
  colorSpace <- c('#E41A1C','#377EB8','#4DAF4A','#984EA3','#F29403','#F781BF','#BC9DCC','#A65628','#54B0E4','#222F75','#1B9E77','#B2DF8A',
                  '#E3BE00','#FB9A99','#E7298A','#910241','#00CDD1','#A6CEE3','#CE1261','#5E4FA2','#8CA77B','#00441B','#DEDC00','#B3DE69','#8DD3C7','#999999')
  if (n <= length(colorSpace)) {
    colors <- colorSpace[1:n]
  } else {
    colors <- grDevices::colorRampPalette(colorSpace)(n)
  }
  return(colors)
}

#' Visualize the inferred cell-cell communication network
#'
#' Automatically save plots in the current working directory.
#'
#' @param object CellChat object
#' @param signaling a signaling pathway name
#' @param signaling.name alternative signaling pathway name to show on the plot
#' @param color.use the character vector defining the color of each cell group
#' @param vertex.receiver a numeric vector giving the index of the cell groups as targets in the first hierarchy plot
#' @param top the fraction of interactions to show (0 < top <= 1)
#' @param sources.use a vector giving the index or the name of source cell groups
#' @param targets.use a vector giving the index or the name of target cell groups.
#' @param remove.isolate whether remove the isolate nodes in the communication network
#' @param weight.scale whether scale the edge weight
#' @param vertex.weight The weight of vertex: either a scale value or a vector

#'
#' Default is a scale value being 1, indicating all vertex is plotted in the same size;
#'
#' Set `vertex.weight` as a vector to plot vertex in different size; setting `vertex.weight = NULL` will have vertex with different size that are portional to the number of cells in each cell group.
#'
#' @param vertex.weight.max the maximum weight of vertex; defualt = max(vertex.weight)
#' @param vertex.size.max the maximum vertex size for visualization
#' @param edge.weight.max.individual the maximum weight of edge when plotting the individual L-R netwrok; defualt = max(net)
#' @param edge.weight.max.aggregate the maximum weight of edge when plotting the aggregated signaling pathway network
#' @param edge.width.max The maximum edge width for visualization
#' @param layout "hierarchy", "circle" or "chord"
#' @param height height of plot
#' @param thresh threshold of the p-value for determining significant interaction
#' @param pt.title font size of the text
#' @param title.space the space between the title and plot
#' @param vertex.label.cex The label size of vertex in the network
#' @param out.format the format of output figures: svg, png and pdf
#'
#' Parameters below are set for "spatial" diagram. Please also check the function `netVisual_spatial` for more parameters.
#' @param alpha.image the transparency of individual spots
#' @param point.size the size of spots
#'
#' @param from,to,bidirection Deprecated. Use `sources.use`,`targets.use`
#' @param vertex.size Deprecated. Use `vertex.weight`
#'
#' Parameters below are set for "chord" diagram. Please also check the function `netVisual_chord_cell` for more parameters.
#' @param group A named group labels for making multiple-group Chord diagrams. The sector names should be used as the names in the vector.
#' The order of group controls the sector orders and if group is set as a factor, the order of levels controls the order of groups.
#' @param cell.order a char vector defining the cell type orders (sector orders)
#' @param small.gap Small gap between sectors.
#' @param big.gap Gap between the different sets of sectors, which are defined in the `group` parameter
#' @param scale scale each sector to same width; default = FALSE; however, it is set to be TRUE when remove.isolate = TRUE
#' @param reduce if the ratio of the width of certain grid compared to the whole circle is less than this value, the grid is removed on the plot. Set it to value less than zero if you want to keep all tiny grid.
#' @param show.legend whether show the figure legend
#' @param legend.pos.x,legend.pos.y adjust the legend position
#' @param nCol number of columns when displaying the network mediated by ligand-receptor using "circle" or "chord"
#'
#' @param ... other parameters (e.g.,vertex.label.cex, vertex.label.color, alpha.edge, label.edge, edge.label.color, edge.label.cex, edge.curved)
#'  passing to `netVisual_hierarchy1`,`netVisual_hierarchy2`,`netVisual_circle`. NB: some parameters might be not supported
#' @importFrom svglite svglite
#' @importFrom grDevices dev.off pdf
#'
#' @return
#' @export
#'
#' @examples
#'
netVisual <- function(object, signaling, signaling.name = NULL, color.use = NULL, vertex.receiver = NULL, sources.use = NULL, targets.use = NULL, top = 1, remove.isolate = FALSE,
                      vertex.weight = 1, vertex.weight.max = NULL, vertex.size.max = NULL,
                      weight.scale = TRUE, edge.weight.max.individual = NULL, edge.weight.max.aggregate = NULL, edge.width.max=8,
                      layout = c("circle","hierarchy","chord","spatial"), height = 5, thresh = 0.05, pt.title = 12, title.space = 6, vertex.label.cex = 0.8,from = NULL, to = NULL, bidirection = NULL,vertex.size = NULL,
                      out.format = c("svg","png"),
                      alpha.image = 0.15, point.size = 1.5,
                      group = NULL,cell.order = NULL,small.gap = 1, big.gap = 10, scale = FALSE, reduce = -1, show.legend = FALSE, legend.pos.x = 20,legend.pos.y = 20, nCol = NULL,
                      ...) {
  layout <- match.arg(layout)
  if (!is.null(vertex.size)) {
    warning("'vertex.size' is deprecated. Use `vertex.weight`")
  }
  if (is.null(vertex.weight)) {
    vertex.weight <- as.numeric(table(object@idents))
  }
  if (is.null(vertex.size.max)) {
    if (length(unique(vertex.weight)) == 1) {
      vertex.size.max <- 5
    } else {
      vertex.size.max <- 15
    }
  }
  pairLR <- searchPair(signaling = signaling, pairLR.use = object@LR$LRsig, key = "pathway_name", matching.exact = T, pair.only = F)

  if (is.null(signaling.name)) {
    signaling.name <- signaling
  }
  net <- object@net

  pairLR.use.name <- dimnames(net$prob)[[3]]
  pairLR.name <- intersect(rownames(pairLR), pairLR.use.name)
  pairLR <- pairLR[pairLR.name, ]
  prob <- net$prob
  pval <- net$pval

  prob[pval > thresh] <- 0
  if (length(pairLR.name) > 1) {
    pairLR.name.use <- pairLR.name[apply(prob[,,pairLR.name], 3, sum) != 0]
  } else {
    pairLR.name.use <- pairLR.name[sum(prob[,,pairLR.name]) != 0]
  }


  if (length(pairLR.name.use) == 0) {
    stop(paste0('There is no significant communication of ', signaling.name))
  } else {
    pairLR <- pairLR[pairLR.name.use,]
  }
  nRow <- length(pairLR.name.use)

  prob <- prob[,,pairLR.name.use]
  pval <- pval[,,pairLR.name.use]

  if (is.null(nCol)) {
    nCol <- min(length(pairLR.name.use), 2)
  }

  if (length(dim(prob)) == 2) {
    prob <- replicate(1, prob, simplify="array")
    pval <- replicate(1, pval, simplify="array")
  }
  #  prob <-(prob-min(prob))/(max(prob)-min(prob))
  if (is.null(edge.weight.max.individual)) {
    edge.weight.max.individual = max(prob)
  }
  prob.sum <- apply(prob, c(1,2), sum)
  #  prob.sum <-(prob.sum-min(prob.sum))/(max(prob.sum)-min(prob.sum))
  if (is.null(edge.weight.max.aggregate)) {
    edge.weight.max.aggregate = max(prob.sum)
  }

  if (layout == "hierarchy") {
    if (is.element("svg", out.format)) {
      svglite::svglite(file = paste0(signaling.name, "_hierarchy_individual.svg"), width = 8, height = nRow*height)
      par(mfrow=c(nRow,2), mar = c(5, 4, 4, 2) +0.1)
      for (i in 1:length(pairLR.name.use)) {
        #signalName_i <- paste0(pairLR$ligand[i], "-",pairLR$receptor[i], sep = "")
        signalName_i <- pairLR$interaction_name_2[i]
        prob.i <- prob[,,i]
        netVisual_hierarchy1(prob.i, vertex.receiver = vertex.receiver, sources.use = sources.use, targets.use = targets.use, remove.isolate = remove.isolate, top = top, color.use = color.use, vertex.weight = vertex.weight, vertex.weight.max = vertex.weight.max, vertex.size.max = vertex.size.max, weight.scale = weight.scale, edge.weight.max = edge.weight.max.individual, edge.width.max=edge.width.max, title.name = signalName_i, vertex.label.cex = vertex.label.cex,...)
        netVisual_hierarchy2(prob.i, vertex.receiver = setdiff(1:nrow(prob.i),vertex.receiver), sources.use = sources.use, targets.use = targets.use, remove.isolate = remove.isolate, top = top, color.use = color.use, vertex.weight = vertex.weight, vertex.weight.max = vertex.weight.max, vertex.size.max = vertex.size.max, weight.scale = weight.scale, edge.weight.max = edge.weight.max.individual, edge.width.max=edge.width.max, title.name = signalName_i, vertex.label.cex = vertex.label.cex,...)
      }
      dev.off()
    }
    if (is.element("png", out.format)) {
      grDevices::png(paste0(signaling.name, "_hierarchy_individual.png"), width = 8, height = nRow*height, units = "in",res = 300)
      par(mfrow=c(nRow,2), mar = c(5, 4, 4, 2) +0.1)
      for (i in 1:length(pairLR.name.use)) {
        signalName_i <- pairLR$interaction_name_2[i]
        prob.i <- prob[,,i]
        netVisual_hierarchy1(prob.i, vertex.receiver = vertex.receiver, sources.use = sources.use, targets.use = targets.use, remove.isolate = remove.isolate, top = top, color.use = color.use, vertex.weight = vertex.weight, vertex.weight.max = vertex.weight.max, vertex.size.max = vertex.size.max, weight.scale = weight.scale, edge.weight.max = edge.weight.max.individual, edge.width.max=edge.width.max, title.name = signalName_i, vertex.label.cex = vertex.label.cex,...)
        netVisual_hierarchy2(prob.i, vertex.receiver = setdiff(1:nrow(prob.i),vertex.receiver), sources.use = sources.use, targets.use = targets.use, remove.isolate = remove.isolate, top = top, color.use = color.use, vertex.weight = vertex.weight, vertex.weight.max = vertex.weight.max, vertex.size.max = vertex.size.max, weight.scale = weight.scale, edge.weight.max =edge.weight.max.individual, edge.width.max=edge.width.max, title.name = signalName_i, vertex.label.cex = vertex.label.cex,...)
      }
      dev.off()
    }
    if (is.element("pdf", out.format)) {
      # grDevices::pdf(paste0(signaling.name, "_hierarchy_individual.pdf"), width = 8, height = nRow*height)
      grDevices::cairo_pdf(paste0(signaling.name, "_hierarchy_individual.pdf"), width = 8, height = nRow*height)
      par(mfrow=c(nRow,2), mar = c(5, 4, 4, 2) +0.1)
      for (i in 1:length(pairLR.name.use)) {
        signalName_i <- pairLR$interaction_name_2[i]
        prob.i <- prob[,,i]
        netVisual_hierarchy1(prob.i, vertex.receiver = vertex.receiver, sources.use = sources.use, targets.use = targets.use, remove.isolate = remove.isolate, top = top, color.use = color.use, vertex.weight = vertex.weight, vertex.weight.max = vertex.weight.max, vertex.size.max = vertex.size.max, weight.scale = weight.scale, edge.weight.max = edge.weight.max.individual, edge.width.max=edge.width.max, title.name = signalName_i, vertex.label.cex = vertex.label.cex,...)
        netVisual_hierarchy2(prob.i, vertex.receiver = setdiff(1:nrow(prob.i),vertex.receiver), sources.use = sources.use, targets.use = targets.use, remove.isolate = remove.isolate, top = top, color.use = color.use, vertex.weight = vertex.weight, vertex.weight.max = vertex.weight.max, vertex.size.max = vertex.size.max, weight.scale = weight.scale, edge.weight.max =edge.weight.max.individual, edge.width.max=edge.width.max, title.name = signalName_i, vertex.label.cex = vertex.label.cex,...)
      }
      dev.off()
    }


    if (is.element("svg", out.format)) {
      svglite::svglite(file = paste0(signaling.name, "_hierarchy_aggregate.svg"), width = 7, height = 1*height)
      par(mfrow=c(1,2), ps = pt.title)
      netVisual_hierarchy1(prob.sum, vertex.receiver = vertex.receiver, sources.use = sources.use, targets.use = targets.use, remove.isolate = remove.isolate, top = top, color.use = color.use, vertex.weight = vertex.weight, vertex.weight.max = vertex.weight.max, vertex.size.max = vertex.size.max, weight.scale = weight.scale, edge.weight.max = edge.weight.max.aggregate, edge.width.max=edge.width.max,title.name = NULL, vertex.label.cex = vertex.label.cex,...)
      netVisual_hierarchy2(prob.sum, vertex.receiver = setdiff(1:nrow(prob.sum),vertex.receiver), sources.use = sources.use, targets.use = targets.use, remove.isolate = remove.isolate, top = top, color.use = color.use, vertex.weight = vertex.weight, vertex.weight.max = vertex.weight.max, vertex.size.max = vertex.size.max, weight.scale = weight.scale, edge.weight.max = edge.weight.max.aggregate, edge.width.max=edge.width.max,title.name = NULL, vertex.label.cex = vertex.label.cex,...)
      graphics::mtext(paste0(signaling.name, " signaling pathway network"), side = 3, outer = TRUE, cex = 1, line = -title.space)
      dev.off()
    }
    if (is.element("png", out.format)) {
      grDevices::png(paste0(signaling.name, "_hierarchy_aggregate.png"), width = 7, height = 1*height, units = "in",res = 300)
      par(mfrow=c(1,2), ps = pt.title)
      netVisual_hierarchy1(prob.sum, vertex.receiver = vertex.receiver, sources.use = sources.use, targets.use = targets.use, remove.isolate = remove.isolate, top = top, color.use = color.use, vertex.weight = vertex.weight, vertex.weight.max = vertex.weight.max, vertex.size.max = vertex.size.max, weight.scale = weight.scale, edge.weight.max = edge.weight.max.aggregate, edge.width.max=edge.width.max, title.name = NULL, vertex.label.cex = vertex.label.cex,...)
      netVisual_hierarchy2(prob.sum, vertex.receiver = setdiff(1:nrow(prob.sum),vertex.receiver), sources.use = sources.use, targets.use = targets.use, remove.isolate = remove.isolate, top = top, color.use = color.use, vertex.weight = vertex.weight, vertex.weight.max = vertex.weight.max, vertex.size.max = vertex.size.max, weight.scale = weight.scale, edge.weight.max = edge.weight.max.aggregate, edge.width.max=edge.width.max,title.name = NULL, vertex.label.cex = vertex.label.cex,...)
      graphics::mtext(paste0(signaling.name, " signaling pathway network"), side = 3, outer = TRUE, cex = 1, line = -title.space)
      dev.off()
    }
    if (is.element("pdf", out.format)) {
      # grDevices::pdf(paste0(signaling.name, "_hierarchy_aggregate.pdf"), width = 7, height = 1*height)
      grDevices::cairo_pdf(paste0(signaling.name, "_hierarchy_aggregate.pdf"), width = 7, height = 1*height)
      par(mfrow=c(1,2), ps = pt.title)
      netVisual_hierarchy1(prob.sum, vertex.receiver = vertex.receiver, sources.use = sources.use, targets.use = targets.use, remove.isolate = remove.isolate, top = top, color.use = color.use, vertex.weight = vertex.weight, vertex.weight.max = vertex.weight.max, vertex.size.max = vertex.size.max, weight.scale = weight.scale, edge.weight.max = edge.weight.max.aggregate, edge.width.max=edge.width.max, title.name = NULL, vertex.label.cex = vertex.label.cex,...)
      netVisual_hierarchy2(prob.sum, vertex.receiver = setdiff(1:nrow(prob.sum),vertex.receiver), sources.use = sources.use, targets.use = targets.use, remove.isolate = remove.isolate, top = top, color.use = color.use, vertex.weight = vertex.weight, vertex.weight.max = vertex.weight.max, vertex.size.max = vertex.size.max, weight.scale = weight.scale, edge.weight.max = edge.weight.max.aggregate, edge.width.max=edge.width.max, title.name = NULL, vertex.label.cex = vertex.label.cex,...)
      graphics::mtext(paste0(signaling.name, " signaling pathway network"), side = 3, outer = TRUE, cex = 1, line = -title.space)
      dev.off()
    }

  } else if (layout == "circle") {
    if (is.element("svg", out.format)) {
      svglite::svglite(file = paste0(signaling.name,"_", layout, "_individual.svg"), width = height, height = nRow*height)
      # par(mfrow=c(nRow,1))
      par(mfrow = c(ceiling(length(pairLR.name.use)/nCol), nCol), xpd=TRUE)
      for (i in 1:length(pairLR.name.use)) {
        #signalName_i <- paste0(pairLR$ligand[i], "-",pairLR$receptor[i], sep = "")
        signalName_i <- pairLR$interaction_name_2[i]
        prob.i <- prob[,,i]
        netVisual_circle(prob.i, sources.use = sources.use, targets.use = targets.use, remove.isolate = remove.isolate, top = top, color.use = color.use, vertex.weight = vertex.weight, vertex.weight.max = vertex.weight.max, vertex.size.max = vertex.size.max, weight.scale = weight.scale, edge.weight.max = edge.weight.max.individual, edge.width.max=edge.width.max, title.name = signalName_i, vertex.label.cex = vertex.label.cex,...)
      }
      dev.off()
    }
    if (is.element("png", out.format)) {
      grDevices::png(paste0(signaling.name,"_", layout, "_individual.png"), width = height, height = nRow*height, units = "in",res = 300)
      # par(mfrow=c(nRow,1))
      par(mfrow = c(ceiling(length(pairLR.name.use)/nCol), nCol), xpd=TRUE)
      for (i in 1:length(pairLR.name.use)) {
        #signalName_i <- paste0(pairLR$ligand[i], "-",pairLR$receptor[i], sep = "")
        signalName_i <- pairLR$interaction_name_2[i]
        prob.i <- prob[,,i]
        netVisual_circle(prob.i, sources.use = sources.use, targets.use = targets.use, remove.isolate = remove.isolate, top = top, color.use = color.use, vertex.weight = vertex.weight, vertex.weight.max = vertex.weight.max, vertex.size.max = vertex.size.max, weight.scale = weight.scale, edge.weight.max = edge.weight.max.individual, edge.width.max=edge.width.max, title.name = signalName_i, vertex.label.cex = vertex.label.cex,...)
      }
      dev.off()
    }
    if (is.element("pdf", out.format)) {
      # grDevices::pdf(paste0(signaling.name,"_", layout, "_individual.pdf"), width = height, height = nRow*height)
      grDevices::cairo_pdf(paste0(signaling.name,"_", layout, "_individual.pdf"), width = height, height = nRow*height)
      # par(mfrow=c(nRow,1))
      par(mfrow = c(ceiling(length(pairLR.name.use)/nCol), nCol), xpd=TRUE)
      for (i in 1:length(pairLR.name.use)) {
        #signalName_i <- paste0(pairLR$ligand[i], "-",pairLR$receptor[i], sep = "")
        signalName_i <- pairLR$interaction_name_2[i]
        prob.i <- prob[,,i]
        netVisual_circle(prob.i, sources.use = sources.use, targets.use = targets.use, remove.isolate = remove.isolate, top = top, color.use = color.use, vertex.weight = vertex.weight, vertex.weight.max = vertex.weight.max, vertex.size.max = vertex.size.max, weight.scale = weight.scale, edge.weight.max = edge.weight.max.individual, edge.width.max=edge.width.max,title.name = signalName_i, vertex.label.cex = vertex.label.cex,...)
      }
      dev.off()
    }

    #  prob.sum <- apply(prob, c(1,2), sum)
    #  prob.sum <-(prob.sum-min(prob.sum))/(max(prob.sum)-min(prob.sum))
    if (is.element("svg", out.format)) {
      svglite(file = paste0(signaling.name,"_", layout,  "_aggregate.svg"), width = height, height = 1*height)
      netVisual_circle(prob.sum, sources.use = sources.use, targets.use = targets.use, remove.isolate = remove.isolate, top = top, color.use = color.use, vertex.weight = vertex.weight, vertex.weight.max = vertex.weight.max, vertex.size.max = vertex.size.max, weight.scale = weight.scale, edge.weight.max = edge.weight.max.aggregate, edge.width.max=edge.width.max,title.name = paste0(signaling.name, " signaling pathway network"), vertex.label.cex = vertex.label.cex,...)
      dev.off()
    }
    if (is.element("png", out.format)) {
      grDevices::png(paste0(signaling.name,"_", layout,  "_aggregate.png"), width = height, height = 1*height, units = "in",res = 300)
      netVisual_circle(prob.sum, sources.use = sources.use, targets.use = targets.use, remove.isolate = remove.isolate, top = top, color.use = color.use, vertex.weight = vertex.weight, vertex.weight.max = vertex.weight.max, vertex.size.max = vertex.size.max, weight.scale = weight.scale, edge.weight.max = edge.weight.max.aggregate, edge.width.max=edge.width.max,title.name = paste0(signaling.name, " signaling pathway network"), vertex.label.cex = vertex.label.cex,...)
      dev.off()
    }
    if (is.element("pdf", out.format)) {
      # grDevices::pdf(paste0(signaling.name,"_", layout,  "_aggregate.pdf"), width = height, height = 1*height)
      grDevices::cairo_pdf(paste0(signaling.name,"_", layout,  "_aggregate.pdf"), width = height, height = 1*height)
      netVisual_circle(prob.sum, sources.use = sources.use, targets.use = targets.use, remove.isolate = remove.isolate, top = top, color.use = color.use, vertex.weight = vertex.weight, vertex.weight.max = vertex.weight.max, vertex.size.max = vertex.size.max, weight.scale = weight.scale, edge.weight.max = edge.weight.max.aggregate, edge.width.max=edge.width.max, title.name = paste0(signaling.name, " signaling pathway network"), vertex.label.cex = vertex.label.cex,...)
      dev.off()
    }
  } else if (layout == "spatial") {
    coordinates <- object@images$coordinates
    if (is.element("svg", out.format)) {
      svglite::svglite(file = paste0(signaling.name,"_", layout, "_individual.svg"), width = height, height = nRow*height)
      # par(mfrow=c(nRow,1))
      par(mfrow = c(ceiling(length(pairLR.name.use)/nCol), nCol), xpd=TRUE)
      for (i in 1:length(pairLR.name.use)) {
        #signalName_i <- paste0(pairLR$ligand[i], "-",pairLR$receptor[i], sep = "")
        signalName_i <- pairLR$interaction_name_2[i]
        prob.i <- prob[,,i]

        netVisual_spatial(prob.i, coordinates = coordinates, labels = labels, alpha.image = alpha.image, point.size = point.size, sources.use = sources.use, targets.use = targets.use, idents.use = idents.use, remove.isolate = remove.isolate, top = top, color.use = color.use, vertex.weight = vertex.weight, vertex.weight.max = vertex.weight.max, vertex.size.max = vertex.size.max, weight.scale = weight.scale, edge.weight.max = edge.weight.max, edge.width.max=edge.width.max,title.name = signalName_i, vertex.label.cex = vertex.label.cex,...)

      }
      dev.off()
    }
    if (is.element("png", out.format)) {
      grDevices::png(paste0(signaling.name,"_", layout, "_individual.png"), width = height, height = nRow*height, units = "in",res = 300)
      # par(mfrow=c(nRow,1))
      par(mfrow = c(ceiling(length(pairLR.name.use)/nCol), nCol), xpd=TRUE)
      for (i in 1:length(pairLR.name.use)) {
        #signalName_i <- paste0(pairLR$ligand[i], "-",pairLR$receptor[i], sep = "")
        signalName_i <- pairLR$interaction_name_2[i]
        prob.i <- prob[,,i]
        netVisual_spatial(prob.i, coordinates = coordinates, labels = labels, alpha.image = alpha.image, point.size = point.size, sources.use = sources.use, targets.use = targets.use, idents.use = idents.use, remove.isolate = remove.isolate, top = top, color.use = color.use, vertex.weight = vertex.weight, vertex.weight.max = vertex.weight.max, vertex.size.max = vertex.size.max, weight.scale = weight.scale, edge.weight.max = edge.weight.max, edge.width.max=edge.width.max,title.name = signalName_i, vertex.label.cex = vertex.label.cex,...)

      }
      dev.off()
    }
    if (is.element("pdf", out.format)) {
      # grDevices::pdf(paste0(signaling.name,"_", layout, "_individual.pdf"), width = height, height = nRow*height)
      grDevices::cairo_pdf(paste0(signaling.name,"_", layout, "_individual.pdf"), width = height, height = nRow*height)
      # par(mfrow=c(nRow,1))
      par(mfrow = c(ceiling(length(pairLR.name.use)/nCol), nCol), xpd=TRUE)
      for (i in 1:length(pairLR.name.use)) {
        #signalName_i <- paste0(pairLR$ligand[i], "-",pairLR$receptor[i], sep = "")
        signalName_i <- pairLR$interaction_name_2[i]
        prob.i <- prob[,,i]
        netVisual_spatial(prob.i, coordinates = coordinates, labels = labels, alpha.image = alpha.image, point.size = point.size, sources.use = sources.use, targets.use = targets.use, idents.use = idents.use, remove.isolate = remove.isolate, top = top, color.use = color.use, vertex.weight = vertex.weight, vertex.weight.max = vertex.weight.max, vertex.size.max = vertex.size.max, weight.scale = weight.scale, edge.weight.max = edge.weight.max, edge.width.max=edge.width.max,title.name = signalName_i, vertex.label.cex = vertex.label.cex,...)

      }
      dev.off()
    }

    #  prob.sum <- apply(prob, c(1,2), sum)
    #  prob.sum <-(prob.sum-min(prob.sum))/(max(prob.sum)-min(prob.sum))
    if (is.element("svg", out.format)) {
      svglite(file = paste0(signaling.name,"_", layout,  "_aggregate.svg"), width = height, height = 1*height)
      netVisual_spatial(prob.sum, coordinates = coordinates, labels = labels, alpha.image = alpha.image, point.size = point.size, sources.use = sources.use, targets.use = targets.use, idents.use = idents.use, remove.isolate = remove.isolate, top = top, color.use = color.use, vertex.weight = vertex.weight, vertex.weight.max = vertex.weight.max, vertex.size.max = vertex.size.max, weight.scale = weight.scale, edge.weight.max = edge.weight.max, edge.width.max=edge.width.max,title.name = paste0(signaling.name, " signaling pathway network"), vertex.label.cex = vertex.label.cex,...)
      dev.off()
    }
    if (is.element("png", out.format)) {
      grDevices::png(paste0(signaling.name,"_", layout,  "_aggregate.png"), width = height, height = 1*height, units = "in",res = 300)
      netVisual_spatial(prob.sum, coordinates = coordinates, labels = labels, alpha.image = alpha.image, point.size = point.size, sources.use = sources.use, targets.use = targets.use, idents.use = idents.use, remove.isolate = remove.isolate, top = top, color.use = color.use, vertex.weight = vertex.weight, vertex.weight.max = vertex.weight.max, vertex.size.max = vertex.size.max, weight.scale = weight.scale, edge.weight.max = edge.weight.max, edge.width.max=edge.width.max,title.name = paste0(signaling.name, " signaling pathway network"), vertex.label.cex = vertex.label.cex,...)
      dev.off()
    }
    if (is.element("pdf", out.format)) {
      # grDevices::pdf(paste0(signaling.name,"_", layout,  "_aggregate.pdf"), width = height, height = 1*height)
      grDevices::cairo_pdf(paste0(signaling.name,"_", layout,  "_aggregate.pdf"), width = height, height = 1*height)
      netVisual_spatial(prob.sum, coordinates = coordinates, labels = labels, alpha.image = alpha.image, point.size = point.size, sources.use = sources.use, targets.use = targets.use, idents.use = idents.use, remove.isolate = remove.isolate, top = top, color.use = color.use, vertex.weight = vertex.weight, vertex.weight.max = vertex.weight.max, vertex.size.max = vertex.size.max, weight.scale = weight.scale, edge.weight.max = edge.weight.max, edge.width.max=edge.width.max,title.name = paste0(signaling.name, " signaling pathway network"), vertex.label.cex = vertex.label.cex,...)
      dev.off()
    }
  } else if (layout == "chord") {
    if (is.element("svg", out.format)) {

      svglite::svglite(file = paste0(signaling.name,"_", layout, "_individual.svg"), width = height, height = nRow*height)
      par(mfrow = c(ceiling(length(pairLR.name.use)/nCol), nCol), xpd=TRUE)
      #  gg <- vector("list", length(pairLR.name.use))
      for (i in 1:length(pairLR.name.use)) {
        title.name <- pairLR$interaction_name_2[i]
        net <- prob[,,i]
        netVisual_chord_cell_internal(net, color.use = color.use, sources.use = sources.use, targets.use = targets.use, remove.isolate = remove.isolate,
                                      group = group, cell.order = cell.order,
                                      lab.cex = vertex.label.cex,small.gap = small.gap, big.gap = big.gap,
                                      scale = scale, reduce = reduce,
                                      title.name = title.name, show.legend = show.legend, legend.pos.x = legend.pos.x,legend.pos.y=legend.pos.y)
      }
      dev.off()
    }
    if (is.element("png", out.format)) {
      grDevices::png(paste0(signaling.name,"_", layout, "_individual.png"), width = height, height = nRow*height, units = "in",res = 300)
      par(mfrow = c(ceiling(length(pairLR.name.use)/nCol), nCol), xpd=TRUE)
      #  gg <- vector("list", length(pairLR.name.use))
      for (i in 1:length(pairLR.name.use)) {
        title.name <- pairLR$interaction_name_2[i]
        net <- prob[,,i]
        netVisual_chord_cell_internal(net, color.use = color.use, sources.use = sources.use, targets.use = targets.use, remove.isolate = remove.isolate,
                                      group = group, cell.order = cell.order,
                                      lab.cex = vertex.label.cex,small.gap = small.gap, big.gap = big.gap,
                                      scale = scale, reduce = reduce,
                                      title.name = title.name, show.legend = show.legend, legend.pos.x = legend.pos.x,legend.pos.y=legend.pos.y)
      }
      dev.off()
    }
    if (is.element("pdf", out.format)) {
      # grDevices::pdf(paste0(signaling.name,"_", layout, "_individual.pdf"), width = height, height = nRow*height)
      grDevices::cairo_pdf(paste0(signaling.name,"_", layout, "_individual.pdf"), width = height, height = nRow*height)
      par(mfrow = c(ceiling(length(pairLR.name.use)/nCol), nCol), xpd=TRUE)
      #  gg <- vector("list", length(pairLR.name.use))
      for (i in 1:length(pairLR.name.use)) {
        title.name <- pairLR$interaction_name_2[i]
        net <- prob[,,i]
        netVisual_chord_cell_internal(net, color.use = color.use, sources.use = sources.use, targets.use = targets.use, remove.isolate = remove.isolate,
                                      group = group, cell.order = cell.order,
                                      lab.cex = vertex.label.cex,small.gap = small.gap, big.gap = big.gap,
                                      scale = scale, reduce = reduce,
                                      title.name = title.name, show.legend = show.legend, legend.pos.x = legend.pos.x,legend.pos.y=legend.pos.y)
      }
      dev.off()
    }

    #  prob.sum <- apply(prob, c(1,2), sum)
    if (is.element("svg", out.format)) {
      svglite(file = paste0(signaling.name,"_", layout,  "_aggregate.svg"), width = height, height = 1*height)
      netVisual_chord_cell_internal(prob.sum, color.use = color.use, sources.use = sources.use, targets.use = targets.use, remove.isolate = remove.isolate,
                                    group = group, cell.order = cell.order,
                                    lab.cex = vertex.label.cex,small.gap = small.gap, big.gap = big.gap,
                                    scale = scale, reduce = reduce,
                                    title.name = paste0(signaling.name, " signaling pathway network"), show.legend = show.legend, legend.pos.x = legend.pos.x,legend.pos.y=legend.pos.y)
      dev.off()
    }
    if (is.element("png", out.format)) {
      grDevices::png(paste0(signaling.name,"_", layout,  "_aggregate.png"), width = height, height = 1*height, units = "in",res = 300)
      netVisual_chord_cell_internal(prob.sum, color.use = color.use, sources.use = sources.use, targets.use = targets.use, remove.isolate = remove.isolate,
                                    group = group, cell.order = cell.order,
                                    lab.cex = vertex.label.cex,small.gap = small.gap, big.gap = big.gap,
                                    scale = scale, reduce = reduce,
                                    title.name = paste0(signaling.name, " signaling pathway network"), show.legend = show.legend, legend.pos.x = legend.pos.x,legend.pos.y=legend.pos.y)
      dev.off()
    }
    if (is.element("pdf", out.format)) {
      # grDevices::pdf(paste0(signaling.name,"_", layout,  "_aggregate.pdf"), width = height, height = 1*height)
      grDevices::cairo_pdf(paste0(signaling.name,"_", layout,  "_aggregate.pdf"), width = height, height = 1*height)
      netVisual_chord_cell_internal(prob.sum, color.use = color.use, sources.use = sources.use, targets.use = targets.use, remove.isolate = remove.isolate,
                                    group = group, cell.order = cell.order,
                                    lab.cex = vertex.label.cex,small.gap = small.gap, big.gap = big.gap,
                                    scale = scale, reduce = reduce,
                                    title.name = paste0(signaling.name, " signaling pathway network"), show.legend = show.legend, legend.pos.x = legend.pos.x,legend.pos.y=legend.pos.y)
      dev.off()
    }
  }

}


#' Visualize the inferred signaling network of signaling pathways by aggregating all L-R pairs
#'
#' @param object CellChat object
#' @param signaling a signaling pathway name
#' @param signaling.name alternative signaling pathway name to show on the plot
#' @param color.use the character vector defining the color of each cell group
#' @param vertex.receiver a numeric vector giving the index of the cell groups as targets in the first hierarchy plot
#' @param sources.use a vector giving the index or the name of source cell groups
#' @param targets.use a vector giving the index or the name of target cell groups.
#' @param idents.use a vector giving the index or the name of cell groups of interest.
#' @param remove.isolate whether remove the isolate nodes in the communication network
#' @param top the fraction of interactions to show
#' @param weight.scale whether scale the edge weight
#' @param vertex.weight The weight of vertex: either a scale value or a vector
#'
#' Default is a scale value being 1, indicating all vertex is plotted in the same size;
#'
#' Set `vertex.weight` as a vector to plot vertex in different size; setting `vertex.weight = NULL` will have vertex with different size that are portional to the number of cells in each cell group.
#'
#' @param vertex.weight.max the maximum weight of vertex; defualt = max(vertex.weight)
#' @param vertex.size.max the maximum vertex size for visualization
#' @param edge.weight.max the maximum weight of edge; defualt = max(net)
#' @param edge.width.max The maximum edge width for visualization
#' @param layout "hierarchy", "circle", "chord" or "spatial"
#' @param thresh threshold of the p-value for determining significant interaction
#' @param pt.title font size of the text
#' @param title.space the space between the title and plot
#' @param vertex.label.cex The label size of vertex in the network
#'
#' Parameters below are set for "spatial" diagram. Please also check the function `netVisual_spatial` for more parameters.
#' @param alpha.image the transparency of individual spots
#' @param point.size the size of spots
#'
#' Parameters below are set for "chord" diagram. Please also check the function `netVisual_chord_cell` for more parameters.
#' @param group A named group labels for making multiple-group Chord diagrams. The sector names should be used as the names in the vector.
#' The order of group controls the sector orders and if group is set as a factor, the order of levels controls the order of groups.
#' @param cell.order a char vector defining the cell type orders (sector orders)
#' @param small.gap Small gap between sectors.
#' @param big.gap Gap between the different sets of sectors, which are defined in the `group` parameter
#' @param scale scale each sector to same width; default = FALSE; however, it is set to be TRUE when remove.isolate = TRUE
#' @param reduce if the ratio of the width of certain grid compared to the whole circle is less than this value, the grid is removed on the plot. Set it to value less than zero if you want to keep all tiny grid.
#' @param show.legend whether show the figure legend
#' @param legend.pos.x,legend.pos.y adjust the legend position
#'
#' @param ... other parameters (e.g.,vertex.label.cex, vertex.label.color, alpha.edge, label.edge, edge.label.color, edge.label.cex, edge.curved)
#'  passing to `netVisual_hierarchy1`,`netVisual_hierarchy2`,`netVisual_circle`,`netVisual_spatial`. NB: some parameters might be not supported
#' @importFrom grDevices recordPlot
#'
#' @return  an object of class "recordedplot" or ggplot
#' @export
#'
#'
netVisual_aggregate <- function(object, signaling, signaling.name = NULL, color.use = NULL, thresh = 0.05, vertex.receiver = NULL, sources.use = NULL, targets.use = NULL, idents.use = NULL, top = 1, remove.isolate = FALSE,
                                vertex.weight = 1, vertex.weight.max = NULL, vertex.size.max = NULL,
                                weight.scale = TRUE, edge.weight.max = NULL, edge.width.max=8,
                                layout = c("circle","hierarchy","chord","spatial"),
                                pt.title = 12, title.space = 6, vertex.label.cex = 0.8,
                                alpha.image = 0.15, point.size = 1.5,
                                group = NULL,cell.order = NULL,small.gap = 1, big.gap = 10, scale = FALSE, reduce = -1, show.legend = FALSE, legend.pos.x = 20,legend.pos.y = 20,
                                ...) {
  layout <- match.arg(layout)
  if (is.null(vertex.weight)) {
    vertex.weight <- as.numeric(table(object@idents))
  }
  if (is.null(vertex.size.max)) {
    if (length(unique(vertex.weight)) == 1) {
      vertex.size.max <- 5
    } else {
      vertex.size.max <- 15
    }
  }
  pairLR <- searchPair(signaling = signaling, pairLR.use = object@LR$LRsig, key = "pathway_name", matching.exact = T, pair.only = T)

  if (is.null(signaling.name)) {
    signaling.name <- signaling
  }
  net <- object@net

  pairLR.use.name <- dimnames(net$prob)[[3]]
  pairLR.name <- intersect(rownames(pairLR), pairLR.use.name)
  pairLR <- pairLR[pairLR.name, ]
  prob <- net$prob
  pval <- net$pval

  prob[pval > thresh] <- 0
  if (length(pairLR.name) > 1) {
    pairLR.name.use <- pairLR.name[apply(prob[,,pairLR.name], 3, sum) != 0]
  } else {
    pairLR.name.use <- pairLR.name[sum(prob[,,pairLR.name]) != 0]
  }


  if (length(pairLR.name.use) == 0) {
    stop(paste0('There is no significant communication of ', signaling.name))
  } else {
    pairLR <- pairLR[pairLR.name.use,]
  }
  nRow <- length(pairLR.name.use)

  prob <- prob[,,pairLR.name.use]
  pval <- pval[,,pairLR.name.use]

  if (length(dim(prob)) == 2) {
    prob <- replicate(1, prob, simplify="array")
    pval <- replicate(1, pval, simplify="array")
  }
  # prob <-(prob-min(prob))/(max(prob)-min(prob))

  if (layout == "hierarchy") {
    prob.sum <- apply(prob, c(1,2), sum)
    # prob.sum <-(prob.sum-min(prob.sum))/(max(prob.sum)-min(prob.sum))
    if (is.null(edge.weight.max)) {
      edge.weight.max = max(prob.sum)
    }
    par(mfrow=c(1,2), ps = pt.title)
    netVisual_hierarchy1(prob.sum, vertex.receiver = vertex.receiver, sources.use = sources.use, targets.use = targets.use, remove.isolate = remove.isolate, top = top, color.use = color.use, vertex.weight = vertex.weight, vertex.weight.max = vertex.weight.max, vertex.size.max = vertex.size.max, weight.scale = weight.scale, edge.weight.max = edge.weight.max, edge.width.max=edge.width.max, title.name = NULL, vertex.label.cex = vertex.label.cex,...)
    netVisual_hierarchy2(prob.sum, vertex.receiver = setdiff(1:nrow(prob.sum),vertex.receiver), sources.use = sources.use, targets.use = targets.use, remove.isolate = remove.isolate, top = top, color.use = color.use, vertex.weight = vertex.weight, vertex.weight.max = vertex.weight.max, vertex.size.max = vertex.size.max, weight.scale = weight.scale, edge.weight.max = edge.weight.max, edge.width.max=edge.width.max, title.name = NULL, vertex.label.cex = vertex.label.cex,...)
    graphics::mtext(paste0(signaling.name, " signaling pathway network"), side = 3, outer = TRUE, cex = 1, line = -title.space)
    # https://www.andrewheiss.com/blog/2016/12/08/save-base-graphics-as-pseudo-objects-in-r/
    # grid.echo()
    # gg <-  grid.grab()
    gg <- recordPlot()
  } else if (layout == "circle") {
    prob.sum <- apply(prob, c(1,2), sum)
    # prob.sum <-(prob.sum-min(prob.sum))/(max(prob.sum)-min(prob.sum))
    gg <- netVisual_circle(prob.sum, sources.use = sources.use, targets.use = targets.use, idents.use = idents.use, remove.isolate = remove.isolate, top = top, color.use = color.use, vertex.weight = vertex.weight, vertex.weight.max = vertex.weight.max, vertex.size.max = vertex.size.max, weight.scale = weight.scale, edge.weight.max = edge.weight.max, edge.width.max=edge.width.max,title.name = paste0(signaling.name, " signaling pathway network"), vertex.label.cex = vertex.label.cex,...)
  }  else if (layout == "spatial") {
    prob.sum <- apply(prob, c(1,2), sum)
    if (vertex.weight == "incoming"){
      if (length(slot(object, "netP")$centr) == 0) {
        stop("Please run `netAnalysis_computeCentrality` to compute the network centrality scores! ")
      }
      # vertex.weight = object@netP$centr[[signaling]]$indeg
      vertex.weight = object@netP$centr["indeg", , signaling]
    }
    if (vertex.weight == "outgoing"){
      if (length(slot(object, "netP")$centr) == 0) {
        stop("Please run `netAnalysis_computeCentrality` to compute the network centrality scores! ")
      }
      # vertex.weight = object@netP$centr[[signaling]]$outdeg
      vertex.weight = object@netP$centr["outdeg", , signaling]
    }
    coordinates <- object@images$coordinates
    labels <- object@idents
    gg <- netVisual_spatial(prob.sum, coordinates = coordinates, labels = labels, alpha.image = alpha.image, point.size = point.size, sources.use = sources.use, targets.use = targets.use, idents.use = idents.use, remove.isolate = remove.isolate, top = top, color.use = color.use, vertex.weight = vertex.weight, vertex.weight.max = vertex.weight.max, vertex.size.max = vertex.size.max, weight.scale = weight.scale, edge.weight.max = edge.weight.max, edge.width.max=edge.width.max,title.name = paste0(signaling.name, " signaling pathway network"), vertex.label.cex = vertex.label.cex,...)

  } else if (layout == "chord") {
    prob.sum <- apply(prob, c(1,2), sum)
    gg <- netVisual_chord_cell_internal(prob.sum, color.use = color.use, sources.use = sources.use, targets.use = targets.use, remove.isolate = remove.isolate,
                                        group = group, cell.order = cell.order,
                                        lab.cex = vertex.label.cex,small.gap = small.gap, big.gap = big.gap,
                                        scale = scale, reduce = reduce,
                                        title.name = paste0(signaling.name, " signaling pathway network"), show.legend = show.legend, legend.pos.x = legend.pos.x, legend.pos.y= legend.pos.y)
  }

  return(gg)

}



#' Visualize the inferred signaling network of individual L-R pairs
#'
#' @param object CellChat object
#' @param signaling a signaling pathway name
#' @param signaling.name alternative signaling pathway name to show on the plot
#' @param pairLR.use a char vector or a data frame consisting of one column named "interaction_name", defining the L-R pairs of interest
#' @param color.use the character vector defining the color of each cell group
#' @param vertex.receiver a numeric vector giving the index of the cell groups as targets in the first hierarchy plot
#' @param sources.use a vector giving the index or the name of source cell groups
#' @param targets.use a vector giving the index or the name of target cell groups.
#' @param remove.isolate whether remove the isolate nodes in the communication network
#' @param top the fraction of interactions to show
#' @param weight.scale whether scale the edge weight
#' @param vertex.weight The weight of vertex: either a scale value or a vector.
#'
#' Default is a scale value being 1, indicating all vertex is plotted in the same size;
#'
#' Set `vertex.weight` as a vector to plot vertex in different size; setting `vertex.weight = NULL` will have vertex with different size that are portional to the number of cells in each cell group.
#'
#' @param vertex.weight.max the maximum weight of vertex; defualt = max(vertex.weight)
#' @param vertex.size.max the maximum vertex size for visualization
#' @param vertex.label.cex The label size of vertex in the network
#' @param edge.weight.max the maximum weight of edge; defualt = max(net)
#' @param edge.width.max The maximum edge width for visualization
#' @param graphics.init whether do graphics initiation using par(...). If graphics.init=FALSE, USERS can use par() in a more fexible way
#' @param layout "hierarchy", "circle" or "chord"
#' @param height height of plot
#' @param thresh threshold of the p-value for determining significant interaction
# #' @param from,to,bidirection Deprecated. Use `sources.use`,`targets.use`
# #' @param vertex.size Deprecated. Use `vertex.weight`

#' Parameters below are set for "spatial" diagram. Please also check the function `netVisual_spatial` for more parameters.
#' @param alpha.image the transparency of individual spots
#' @param point.size the size of spots
#'
#' Parameters below are set for "chord" diagram. Please also check the function `netVisual_chord_cell` for more parameters.
#' @param group A named group labels for making multiple-group Chord diagrams. The sector names should be used as the names in the vector.
#' The order of group controls the sector orders and if group is set as a factor, the order of levels controls the order of groups.
#' @param cell.order a char vector defining the cell type orders (sector orders)
#' @param small.gap Small gap between sectors.
#' @param big.gap Gap between the different sets of sectors, which are defined in the `group` parameter
#' @param scale scale each sector to same width; default = FALSE; however, it is set to be TRUE when remove.isolate = TRUE
#' @param reduce if the ratio of the width of certain grid compared to the whole circle is less than this value, the grid is removed on the plot. Set it to value less than zero if you want to keep all tiny grid.
#' @param show.legend whether show the figure legend
#' @param legend.pos.x,legend.pos.y adjust the legend position
#' @param nCol number of columns when displaying the figures using "circle" or "chord"
#'
#' @param ... other parameters (e.g.,vertex.label.cex, vertex.label.color, alpha.edge, label.edge, edge.label.color, edge.label.cex, edge.curved)
#'  passing to `netVisual_hierarchy1`,`netVisual_hierarchy2`,`netVisual_circle`. NB: some parameters might be not supported
#' @importFrom grDevices dev.off pdf
#'
#' @return  an object of class "recordedplot"
#' @export
#'
#'
netVisual_individual <- function(object, signaling, signaling.name = NULL, pairLR.use = NULL, color.use = NULL, vertex.receiver = NULL, sources.use = NULL, targets.use = NULL, top = 1, remove.isolate = FALSE,
                                 vertex.weight = 1, vertex.weight.max = NULL, vertex.size.max = NULL, vertex.label.cex = 0.8,
                                 weight.scale = TRUE, edge.weight.max = NULL, edge.width.max=8, graphics.init = TRUE,
                                 layout = c("circle","hierarchy","chord","spatial"), height = 5, thresh = 0.05, #from = NULL, to = NULL, bidirection = NULL,vertex.size = NULL,
                                 alpha.image = 0.15, point.size = 1.5,
                                 group = NULL,cell.order = NULL,small.gap = 1, big.gap = 10, scale = FALSE, reduce = -1, show.legend = FALSE, legend.pos.x = 20, legend.pos.y = 20, nCol = NULL,
                                 ...) {
  layout <- match.arg(layout)
  # if (!is.null(vertex.size)) {
  #   warning("'vertex.size' is deprecated. Use `vertex.weight`")
  # }
  if (is.null(vertex.weight)) {
    vertex.weight <- as.numeric(table(object@idents))
  }
  if (is.null(vertex.size.max)) {
    if (length(unique(vertex.weight)) == 1) {
      vertex.size.max <- 5
    } else {
      vertex.size.max <- 15
    }
  }

  pairLR <- searchPair(signaling = signaling, pairLR.use = object@LR$LRsig, key = "pathway_name", matching.exact = T, pair.only = F)

  if (is.null(signaling.name)) {
    signaling.name <- signaling
  }
  net <- object@net

  pairLR.use.name <- dimnames(net$prob)[[3]]
  pairLR.name <- intersect(rownames(pairLR), pairLR.use.name)
  if (!is.null(pairLR.use)) {
    if (is.data.frame(pairLR.use)) {
      pairLR.name <- intersect(pairLR.name, as.character(pairLR.use$interaction_name))
    } else {
      pairLR.name <- intersect(pairLR.name, as.character(pairLR.use))
    }

    if (length(pairLR.name) == 0) {
      stop("There is no significant communication for the input L-R pairs!")
    }
  }

  pairLR <- pairLR[pairLR.name, ]
  prob <- net$prob
  pval <- net$pval

  prob[pval > thresh] <- 0
  if (length(pairLR.name) > 1) {
    pairLR.name.use <- pairLR.name[apply(prob[,,pairLR.name], 3, sum) != 0]
  } else {
    pairLR.name.use <- pairLR.name[sum(prob[,,pairLR.name]) != 0]
  }

  if (length(pairLR.name.use) == 0) {
    stop(paste0('There is no significant communication of ', signaling.name))
  } else {
    pairLR <- pairLR[pairLR.name.use,]
  }

  nRow <- length(pairLR.name.use)

  prob <- prob[,,pairLR.name.use]
  pval <- pval[,,pairLR.name.use]

  if (is.null(nCol)) {
    nCol <- min(length(pairLR.name.use), 2)
  }

  if (length(dim(prob)) == 2) {
    prob <- replicate(1, prob, simplify="array")
    pval <- replicate(1, pval, simplify="array")
  }

  # prob <-(prob-min(prob))/(max(prob)-min(prob))
  if (is.null(edge.weight.max)) {
    edge.weight.max = max(prob)
  }

  if (layout == "hierarchy") {
    if (graphics.init) {
      par(mfrow=c(nRow,2), mar = c(5, 4, 4, 2) +0.1)
    }

    for (i in 1:length(pairLR.name.use)) {
      signalName_i <- pairLR$interaction_name_2[i]
      prob.i <- prob[,,i]
      netVisual_hierarchy1(prob.i, vertex.receiver = vertex.receiver, sources.use = sources.use, targets.use = targets.use, remove.isolate = remove.isolate, top = top, color.use = color.use, vertex.weight = vertex.weight, vertex.weight.max = vertex.weight.max, vertex.size.max = vertex.size.max, weight.scale = weight.scale, edge.weight.max = edge.weight.max, edge.width.max=edge.width.max, title.name = signalName_i,...)
      netVisual_hierarchy2(prob.i, vertex.receiver = setdiff(1:nrow(prob.i),vertex.receiver), sources.use = sources.use, targets.use = targets.use, remove.isolate = remove.isolate, top = top, color.use = color.use, vertex.weight = vertex.weight, vertex.weight.max = vertex.weight.max, vertex.size.max = vertex.size.max, weight.scale = weight.scale, edge.weight.max = edge.weight.max, edge.width.max=edge.width.max, title.name = signalName_i,...)
    }
    # grid.echo()
    # gg <-  grid.grab()
    gg <- recordPlot()

  } else if (layout == "circle") {
    # par(mfrow=c(nRow,1))
    if (graphics.init) {
      par(mfrow = c(ceiling(length(pairLR.name.use)/nCol), nCol), xpd=TRUE)
    }
    gg <- vector("list", length(pairLR.name.use))
    for (i in 1:length(pairLR.name.use)) {
      signalName_i <- pairLR$interaction_name_2[i]
      prob.i <- prob[,,i]
      gg[[i]] <- netVisual_circle(prob.i, sources.use = sources.use, targets.use = targets.use, remove.isolate = remove.isolate, top = top, color.use = color.use, vertex.weight = vertex.weight, vertex.weight.max = vertex.weight.max, vertex.size.max = vertex.size.max, weight.scale = weight.scale, edge.weight.max = edge.weight.max, edge.width.max=edge.width.max, title.name = signalName_i, vertex.label.cex = vertex.label.cex,...)
    }
  } else if (layout == "spatial") {
    # par(mfrow=c(nRow,1))
    if (graphics.init) {
      par(mfrow = c(ceiling(length(pairLR.name.use)/nCol), nCol), xpd=TRUE)
    }
    coordinates <- object@images$coordinates
    gg <- vector("list", length(pairLR.name.use))
    for (i in 1:length(pairLR.name.use)) {
      signalName_i <- pairLR$interaction_name_2[i]
      prob.i <- prob[,,i]
      gg[[i]] <- netVisual_spatial(prob.i, coordinates = coordinates, labels = labels, alpha.image = alpha.image, point.size = point.size, sources.use = sources.use, targets.use = targets.use, idents.use = idents.use, remove.isolate = remove.isolate, top = top, color.use = color.use, vertex.weight = vertex.weight, vertex.weight.max = vertex.weight.max, vertex.size.max = vertex.size.max, weight.scale = weight.scale, edge.weight.max = edge.weight.max, edge.width.max=edge.width.max,title.name = signalName_i, vertex.label.cex = vertex.label.cex,...)
    }
  } else if (layout == "chord") {
    if (graphics.init) {
      par(mfrow = c(ceiling(length(pairLR.name.use)/nCol), nCol), xpd=TRUE)
    }

    gg <- vector("list", length(pairLR.name.use))
    for (i in 1:length(pairLR.name.use)) {
      title.name <- pairLR$interaction_name_2[i]
      net <- prob[,,i]
      gg[[i]] <- netVisual_chord_cell_internal(net, color.use = color.use, sources.use = sources.use, targets.use = targets.use, remove.isolate = remove.isolate,
                                               group = group, cell.order = cell.order,
                                               lab.cex = vertex.label.cex,small.gap = small.gap, big.gap = big.gap,
                                               scale = scale, reduce = reduce,
                                               title.name = title.name, show.legend = show.legend, legend.pos.x = legend.pos.x, legend.pos.y = legend.pos.y)
    }
  }
  return(gg)
}



#' Hierarchy plot of cell-cell communications sending to cell groups in vertex.receiver
#'
#' The width of edges represent the strength of the communication.
#'
#' @param net a weighted matrix defining the signaling network
#' @param vertex.receiver  a numeric vector giving the index of the cell groups as targets in the first hierarchy plot
#' @param color.use the character vector defining the color of each cell group
#' @param title.name alternative signaling pathway name to show on the plot
#' @param sources.use a vector giving the index or the name of source cell groups
#' @param targets.use a vector giving the index or the name of target cell groups.
#' @param remove.isolate whether remove the isolate nodes in the communication network
#' @param top the fraction of interactions to show
#' @param weight.scale whether rescale the edge weights
#' @param vertex.weight The weight of vertex: either a scale value or a vector
#' @param vertex.weight.max the maximum weight of vertex; defualt = max(vertex.weight)
#' @param vertex.size.max the maximum vertex size for visualization
#' @param edge.weight.max the maximum weight of edge; defualt = max(net)
#' @param edge.width.max The maximum edge width for visualization
#' @param label.dist the distance between labels and dot position
#' @param space.v the space between different columns in the plot
#' @param space.h the space between different rows in the plot
#' @param edge.curved Specifies whether to draw curved edges, or not.
#' This can be a logical or a numeric vector or scalar.
#' First the vector is replicated to have the same length as the number of
#' edges in the graph. Then it is interpreted for each edge separately.
#' A numeric value specifies the curvature of the edge; zero curvature means
#' straight edges, negative values means the edge bends clockwise, positive
#' values the opposite. TRUE means curvature 0.5, FALSE means curvature zero
#' @param shape The shape of the vertex, currently "circle", "square",
#' "csquare", "rectangle", "crectangle", "vrectangle", "pie" (see
#' vertex.shape.pie), "sphere", and "none" are supported, and only by the
#' plot.igraph command. "none" does not draw the vertices at all, although
#' vertex label are plotted (if given). See shapes for details about vertex
#' shapes and vertex.shape.pie for using pie charts as vertices.
#' @param margin The amount of empty space below, over, at the left and right
#'  of the plot, it is a numeric vector of length four. Usually values between
#'  0 and 0.5 are meaningful, but negative values are also possible, that will
#'  make the plot zoom in to a part of the graph. If it is shorter than four
#'  then it is recycled.
#' @param vertex.label.cex The label size of vertex
#' @param vertex.label.color The color of label for vertex
#' @param arrow.width The width of arrows
#' @param arrow.size the size of arrow
#' @param alpha.edge the transprency of edge
#' @param label.edge whether label edge
#' @param edge.label.color The color for single arrow
#' @param edge.label.cex The size of label for arrows
#' @param vertex.size Deprecated. Use `vertex.weight`
#' @importFrom igraph graph_from_adjacency_matrix ends E V layout_
#' @importFrom grDevices adjustcolor recordPlot
#' @importFrom shape Arrows
#' @return  an object of class "recordedplot"
#' @export
netVisual_hierarchy1 <- function(net, vertex.receiver, color.use = NULL, title.name = NULL,  sources.use = NULL, targets.use = NULL, remove.isolate = FALSE, top = 1,
                                 weight.scale = FALSE, vertex.weight=20, vertex.weight.max = NULL, vertex.size.max = NULL,
                                 edge.weight.max = NULL, edge.width.max=8, alpha.edge = 0.6,
                                 label.dist = 2.8, space.v = 1.5, space.h = 1.6, shape= NULL, label.edge=FALSE,edge.curved=0, margin=0.2,
                                 vertex.label.cex=0.6,vertex.label.color= "black",arrow.width=1,arrow.size = 0.2,edge.label.color='black',edge.label.cex=0.5, vertex.size = NULL){
  if (!is.null(vertex.size)) {
    warning("'vertex.size' is deprecated. Use `vertex.weight`")
  }
  if (is.null(vertex.size.max)) {
    if (length(unique(vertex.weight)) == 1) {
      vertex.size.max <- 5
    } else {
      vertex.size.max <- 15
    }
  }
  options(warn = -1)
  thresh <- stats::quantile(net, probs = 1-top)
  net[net < thresh] <- 0
  cells.level <- rownames(net)

  if ((!is.null(sources.use)) | (!is.null(targets.use))) {
    df.net <- reshape2::melt(net, value.name = "value")
    colnames(df.net)[1:2] <- c("source","target")
    # keep the interactions associated with sources and targets of interest
    if (!is.null(sources.use)){
      if (is.numeric(sources.use)) {
        sources.use <- cells.level[sources.use]
      }
      df.net <- subset(df.net, source %in% sources.use)
    }
    if (!is.null(targets.use)){
      if (is.numeric(targets.use)) {
        targets.use <- cells.level[targets.use]
      }
      df.net <- subset(df.net, target %in% targets.use)
    }
    df.net$source <- factor(df.net$source, levels = cells.level)
    df.net$target <- factor(df.net$target, levels = cells.level)
    df.net$value[is.na(df.net$value)] <- 0
    net <- tapply(df.net[["value"]], list(df.net[["source"]], df.net[["target"]]), sum)
  }
  net[is.na(net)] <- 0

  if (remove.isolate) {
    idx1 <- which(Matrix::rowSums(net) == 0)
    idx2 <- which(Matrix::colSums(net) == 0)
    idx <- intersect(idx1, idx2)
    net <- net[-idx, ]
    net <- net[, -idx]
  }

  if (is.null(color.use)) {
    color.use <- scPalette(nrow(net))
  }

  if (is.null(vertex.weight.max)) {
    vertex.weight.max <- max(vertex.weight)
  }
  vertex.weight <- vertex.weight/vertex.weight.max*vertex.size.max+6

  m <- length(vertex.receiver)
  net2 <- net
  reorder.row <- c(vertex.receiver, setdiff(1:nrow(net),vertex.receiver))
  net2 <- net2[reorder.row,vertex.receiver]
  # Expand out to symmetric (M+N)x(M+N) matrix
  m1 <- nrow(net2); n1 <- ncol(net2)
  net3 <- rbind(cbind(matrix(0, m1, m1), net2), matrix(0, n1, m1+n1))

  row.names(net3) <- c(row.names(net)[vertex.receiver], row.names(net)[setdiff(1:m1,vertex.receiver)], rep("",m))
  colnames(net3) <- row.names(net3)
  color.use3 <- c(color.use[vertex.receiver], color.use[setdiff(1:m1,vertex.receiver)], rep("#FFFFFF",length(vertex.receiver)))
  color.use3.frame <- c(color.use[vertex.receiver], color.use[setdiff(1:m1,vertex.receiver)], color.use[vertex.receiver])

  if (length(vertex.weight) != 1) {
    vertex.weight = c(vertex.weight[vertex.receiver], vertex.weight[setdiff(1:m1,vertex.receiver)],vertex.weight[vertex.receiver])
  }
  if (is.null(shape)) {
    shape <- c(rep("circle",m), rep("circle", m1-m), rep("circle",m))
  }

  g <- graph_from_adjacency_matrix(net3, mode = "directed", weighted = T)
  edge.start <- ends(g, es=E(g), names=FALSE)
  coords <- matrix(NA, nrow(net3), 2)
  coords[1:m,1] <- 0; coords[(m+1):m1,1] <- space.h; coords[(m1+1):nrow(net3),1] <- space.h/2;
  coords[1:m,2] <- seq(space.v, 0, by = -space.v/(m-1)); coords[(m+1):m1,2] <- seq(space.v, 0, by = -space.v/(m1-m-1));coords[(m1+1):nrow(net3),2] <- seq(space.v, 0, by = -space.v/(n1-1));
  coords_scale<-coords

  igraph::V(g)$size<-vertex.weight
  igraph::V(g)$color<-color.use3[igraph::V(g)]
  igraph::V(g)$frame.color <- color.use3.frame[igraph::V(g)]
  igraph::V(g)$label.color <- vertex.label.color
  igraph::V(g)$label.cex<-vertex.label.cex
  if(label.edge){
    E(g)$label<-E(g)$weight
    igraph::E(g)$label <- round(igraph::E(g)$label, digits = 1)
  }
  if (is.null(edge.weight.max)) {
    edge.weight.max <- max(igraph::E(g)$weight)
  }
  if (weight.scale == TRUE) {
    # E(g)$width<-0.3+edge.max.width/(max(E(g)$weight)-min(E(g)$weight))*(E(g)$weight-min(E(g)$weight))
    E(g)$width<- 0.3+E(g)$weight/edge.weight.max*edge.width.max
  }else{
    E(g)$width<-0.3+edge.width.max*E(g)$weight
  }

  E(g)$arrow.width<-arrow.width
  E(g)$arrow.size<-arrow.size
  E(g)$label.color<-edge.label.color
  E(g)$label.cex<-edge.label.cex
  E(g)$color<-adjustcolor(igraph::V(g)$color[edge.start[,1]],alpha.edge)

  label.dist <- c(rep(space.h*label.dist,m), rep(space.h*label.dist, m1-m),rep(0, nrow(net3)-m1))
  label.locs <- c(rep(-pi, m), rep(0, m1-m),rep(-pi, nrow(net3)-m1))
  # text.pos <- cbind(c(-space.h/1.5, space.h/10, space.h/1.2), space.v-space.v/10)
  text.pos <- cbind(c(-space.h/1.5, space.h/22, space.h/1.5), space.v-space.v/7)
  igraph::add.vertex.shape("fcircle", clip=igraph::igraph.shape.noclip,plot=mycircle, parameters=list(vertex.frame.color=1, vertex.frame.width=1))
  plot(g,edge.curved=edge.curved,layout=coords_scale,margin=margin,rescale=T,vertex.shape="fcircle", vertex.frame.width = c(rep(1,m1), rep(2,nrow(net3)-m1)),
       vertex.label.degree=label.locs, vertex.label.dist=label.dist, vertex.label.family="Helvetica")
  text(text.pos, c("Source","Target","Source"), cex = 0.8, col = c("#c51b7d","#c51b7d","#2f6661"))
  arrow.pos1 <- c(-space.h/1.5, space.v-space.v/4, space.h/100000, space.v-space.v/4)
  arrow.pos2 <- c(space.h/1.5, space.v-space.v/4, space.h/20, space.v-space.v/4)
  shape::Arrows(arrow.pos1[1], arrow.pos1[2], arrow.pos1[3], arrow.pos1[4], col = "#c51b7d",arr.lwd = 0.0001,arr.length = 0.2, lwd = 0.8,arr.type="triangle")
  shape::Arrows(arrow.pos2[1], arrow.pos2[2], arrow.pos2[3], arrow.pos2[4], col = "#2f6661",arr.lwd = 0.0001,arr.length = 0.2, lwd = 0.8,arr.type="triangle")
  if (!is.null(title.name)) {
    title.pos = c(space.h/8, space.v)
    text(title.pos[1],title.pos[2],paste0(title.name, " signaling network"), cex = 1)
  }
  # https://www.andrewheiss.com/blog/2016/12/08/save-base-graphics-as-pseudo-objects-in-r/
  # grid.echo()
  # gg <-  grid.grab()
  gg <- recordPlot()
  return(gg)
}


#' Hierarchy plot of cell-cell communication sending to cell groups not in vertex.receiver
#'
#' This function loads the significant interactions as a weighted matrix, and colors
#' represent different types of cells as a structure. The width of edges represent the strength of the communication.
#'
#' @param net a weighted matrix defining the signaling network
#' @param vertex.receiver  a numeric vector giving the index of the cell groups as targets in the first hierarchy plot
#' @param color.use the character vector defining the color of each cell group
#' @param title.name alternative signaling pathway name to show on the plot
#' @param sources.use a vector giving the index or the name of source cell groups
#' @param targets.use a vector giving the index or the name of target cell groups.
#' @param remove.isolate whether remove the isolate nodes in the communication network
#' @param top the fraction of interactions to show
#' @param weight.scale whether rescale the edge weights
#' @param vertex.weight The weight of vertex: either a scale value or a vector
#' @param vertex.weight.max the maximum weight of vertex; defualt = max(vertex.weight)
#' @param vertex.size.max the maximum vertex size for visualization
#' @param edge.weight.max the maximum weight of edge; defualt = max(net)
#' @param edge.width.max The maximum edge width for visualization
#' @param label.dist the distance between labels and dot position
#' @param space.v the space between different columns in the plot
#' @param space.h the space between different rows in the plot
#' @param label.edge Whether or not shows the label of edges (number of connections between different cell types)
#' @param edge.curved Specifies whether to draw curved edges, or not.
#' This can be a logical or a numeric vector or scalar.
#' First the vector is replicated to have the same length as the number of
#' edges in the graph. Then it is interpreted for each edge separately.
#' A numeric value specifies the curvature of the edge; zero curvature means
#' straight edges, negative values means the edge bends clockwise, positive
#' values the opposite. TRUE means curvature 0.5, FALSE means curvature zero
#' @param shape The shape of the vertex, currently "circle", "square",
#' "csquare", "rectangle", "crectangle", "vrectangle", "pie" (see
#' vertex.shape.pie), "sphere", and "none" are supported, and only by the
#' plot.igraph command. "none" does not draw the vertices at all, although
#' vertex label are plotted (if given). See shapes for details about vertex
#' shapes and vertex.shape.pie for using pie charts as vertices.
#' @param margin The amount of empty space below, over, at the left and right
#'  of the plot, it is a numeric vector of length four. Usually values between
#'  0 and 0.5 are meaningful, but negative values are also possible, that will
#'  make the plot zoom in to a part of the graph. If it is shorter than four
#'  then it is recycled.
#' @param vertex.label.cex The label size of vertex
#' @param vertex.label.color The color of label for vertex
#' @param arrow.width The width of arrows
#' @param arrow.size the size of arrow
#' @param alpha.edge the transprency of edge
#' @param edge.label.color The color for single arrow
#' @param edge.label.cex The size of label for arrows
#' @param vertex.size Deprecated. Use `vertex.weight`
#' @importFrom igraph graph_from_adjacency_matrix ends E V layout_
#' @importFrom grDevices adjustcolor recordPlot
#' @importFrom shape Arrows
#' @return  an object of class "recordedplot"
#' @export
netVisual_hierarchy2 <-function(net, vertex.receiver, color.use = NULL, title.name = NULL, sources.use = NULL, targets.use = NULL, remove.isolate = FALSE, top = 1,
                                weight.scale = FALSE, vertex.weight=20, vertex.weight.max = NULL, vertex.size.max = NULL,
                                edge.weight.max = NULL, edge.width.max=8,alpha.edge = 0.6,
                                label.dist = 2.8, space.v = 1.5, space.h = 1.6, shape= NULL, label.edge=FALSE,edge.curved=0, margin=0.2,
                                vertex.label.cex=0.6,vertex.label.color= "black",arrow.width=1,arrow.size = 0.2,edge.label.color='black',edge.label.cex=0.5, vertex.size = NULL){
  if (!is.null(vertex.size)) {
    warning("'vertex.size' is deprecated. Use `vertex.weight`")
  }
  if (is.null(vertex.size.max)) {
    if (length(unique(vertex.weight)) == 1) {
      vertex.size.max <- 5
    } else {
      vertex.size.max <- 15
    }
  }
  options(warn = -1)
  thresh <- stats::quantile(net, probs = 1-top)
  net[net < thresh] <- 0

  if ((!is.null(sources.use)) | (!is.null(targets.use))) {
    df.net <- reshape2::melt(net, value.name = "value")
    colnames(df.net)[1:2] <- c("source","target")
    # keep the interactions associated with sources and targets of interest
    if (!is.null(sources.use)){
      if (is.numeric(sources.use)) {
        sources.use <- levels(object@idents)[sources.use]
      }
      df.net <- subset(df.net, source %in% sources.use)
    }
    if (!is.null(targets.use)){
      if (is.numeric(targets.use)) {
        targets.use <- levels(object@idents)[targets.use]
      }
      df.net <- subset(df.net, target %in% targets.use)
    }
    cells.level <- levels(object@idents)
    df.net$source <- factor(df.net$source, levels = cells.level)
    df.net$target <- factor(df.net$target, levels = cells.level)
    df.net$value[is.na(df.net$value)] <- 0
    net <- tapply(df.net[["value"]], list(df.net[["source"]], df.net[["target"]]), sum)
  }
  net[is.na(net)] <- 0

  if (remove.isolate) {
    idx1 <- which(Matrix::rowSums(net) == 0)
    idx2 <- which(Matrix::colSums(net) == 0)
    idx <- intersect(idx1, idx2)
    net <- net[-idx, ]
    net <- net[, -idx]
  }


  if (is.null(color.use)) {
    color.use <- scPalette(nrow(net))
  }

  if (is.null(vertex.weight.max)) {
    vertex.weight.max <- max(vertex.weight)
  }
  vertex.weight <- vertex.weight/vertex.weight.max*vertex.size.max+6

  m <- length(vertex.receiver)
  m0 <- nrow(net)-length(vertex.receiver)
  net2 <- net
  reorder.row <- c(setdiff(1:nrow(net),vertex.receiver), vertex.receiver)
  net2 <- net2[reorder.row,vertex.receiver]
  # Expand out to symmetric (M+N)x(M+N) matrix
  m1 <- nrow(net2); n1 <- ncol(net2)
  net3 <- rbind(cbind(matrix(0, m1, m1), net2), matrix(0, n1, m1+n1))
  row.names(net3) <- c(row.names(net)[setdiff(1:m1,vertex.receiver)],row.names(net)[vertex.receiver],  rep("",m))
  colnames(net3) <- row.names(net3)
  color.use3 <- c(color.use[setdiff(1:m1,vertex.receiver)],color.use[vertex.receiver],  rep("#FFFFFF",length(vertex.receiver)))
  color.use3.frame <- c(color.use[setdiff(1:m1,vertex.receiver)], color.use[vertex.receiver], color.use[vertex.receiver])


  if (length(vertex.weight) != 1) {
    vertex.weight = c(vertex.weight[setdiff(1:m1,vertex.receiver)], vertex.weight[vertex.receiver], vertex.weight[vertex.receiver])
  }
  if (is.null(shape)) {
    shape <- rep("circle",nrow(net3))
  }

  g <- graph_from_adjacency_matrix(net3, mode = "directed", weighted = T)
  edge.start <- ends(g, es=igraph::E(g), names=FALSE)
  coords <- matrix(NA, nrow(net3), 2)
  coords[1:m0,1] <- 0; coords[(m0+1):m1,1] <- space.h; coords[(m1+1):nrow(net3),1] <- space.h/2;
  coords[1:m0,2] <- seq(space.v, 0, by = -space.v/(m0-1)); coords[(m0+1):m1,2] <- seq(space.v, 0, by = -space.v/(m1-m0-1));coords[(m1+1):nrow(net3),2] <- seq(space.v, 0, by = -space.v/(n1-1));
  coords_scale<-coords

  igraph::V(g)$size<-vertex.weight
  igraph::V(g)$color<-color.use3[igraph::V(g)]
  igraph::V(g)$frame.color <- color.use3.frame[igraph::V(g)]
  igraph::V(g)$label.color <- vertex.label.color
  igraph::V(g)$label.cex<-vertex.label.cex
  if(label.edge){
    igraph::E(g)$label<-igraph::E(g)$weight
    igraph::E(g)$label <- round(igraph::E(g)$label, digits = 1)
  }
  if (is.null(edge.weight.max)) {
    edge.weight.max <- max(igraph::E(g)$weight)
  }
  if (weight.scale == TRUE) {
    # E(g)$width<-0.3+edge.max.width/(max(E(g)$weight)-min(E(g)$weight))*(E(g)$weight-min(E(g)$weight))
    igraph::E(g)$width<- 0.3+igraph::E(g)$weight/edge.weight.max*edge.width.max
  }else{
    igraph::E(g)$width<-0.3+edge.width.max*igraph::E(g)$weight
  }
  igraph::E(g)$arrow.width<-arrow.width
  igraph::E(g)$arrow.size<-arrow.size
  igraph::E(g)$label.color<-edge.label.color
  igraph::E(g)$label.cex<-edge.label.cex
  igraph::E(g)$color<-adjustcolor(igraph::V(g)$color[edge.start[,1]],alpha.edge)

  label.dist <- c(rep(space.h*label.dist,m), rep(space.h*label.dist, m1-m),rep(0, nrow(net3)-m1))
  label.locs <- c(rep(-pi, m0), rep(0, m1-m0),rep(-pi, nrow(net3)-m1))
  #text.pos <- cbind(c(-space.h/1.5, space.h/10, space.h/1.2), space.v-space.v/10)
  text.pos <- cbind(c(-space.h/1.5, space.h/22, space.h/1.5), space.v-space.v/7)
  igraph::add.vertex.shape("fcircle", clip=igraph::igraph.shape.noclip,plot=mycircle, parameters=list(vertex.frame.color=1, vertex.frame.width=1))
  plot(g,edge.curved=edge.curved,layout=coords_scale,margin=margin,rescale=T,vertex.shape="fcircle", vertex.frame.width = c(rep(1,m1), rep(2,nrow(net3)-m1)),
       vertex.label.degree=label.locs, vertex.label.dist=label.dist, vertex.label.family="Helvetica")
  text(text.pos, c("Source","Target","Source"), cex = 0.8, col = c("#c51b7d","#2f6661","#2f6661"))

  arrow.pos1 <- c(-space.h/1.5, space.v-space.v/4, space.h/100000, space.v-space.v/4)
  arrow.pos2 <- c(space.h/1.5, space.v-space.v/4, space.h/20, space.v-space.v/4)
  shape::Arrows(arrow.pos1[1], arrow.pos1[2], arrow.pos1[3], arrow.pos1[4], col = "#c51b7d",arr.lwd = 0.0001,arr.length = 0.2, lwd = 0.8,arr.type="triangle")
  shape::Arrows(arrow.pos2[1], arrow.pos2[2], arrow.pos2[3], arrow.pos2[4], col = "#2f6661",arr.lwd = 0.0001,arr.length = 0.2, lwd = 0.8,arr.type="triangle")

  if (!is.null(title.name)) {
    title.pos = c(space.h/8, space.v)
    text(title.pos[1],title.pos[2],paste0(title.name, " signaling network"), cex = 1)
  }
  # https://www.andrewheiss.com/blog/2016/12/08/save-base-graphics-as-pseudo-objects-in-r/
  # grid.echo()
  # gg <-  grid.grab()
  gg <- recordPlot()
  return(gg)
}


#' Circle plot of cell-cell communication network
#'
#' The width of edges represent the strength of the communication.
#'
#' @param net A weighted matrix representing the connections
#' @param color.use Colors represent different cell groups
#' @param title.name the name of the title
#' @param sources.use a vector giving the index or the name of source cell groups
#' @param targets.use a vector giving the index or the name of target cell groups.
#' @param idents.use a vector giving the index or the name of cell groups of interest.
#' @param remove.isolate whether remove the isolate nodes in the communication network
#' @param top the fraction of interactions to show
#' @param weight.scale whether scale the weight
#' @param vertex.weight The weight of vertex: either a scale value or a vector
#' @param vertex.weight.max the maximum weight of vertex; defualt = max(vertex.weight)
#' @param vertex.size.max the maximum vertex size for visualization
#' @param vertex.label.cex The label size of vertex
#' @param vertex.label.color The color of label for vertex
#' @param edge.weight.max the maximum weight of edge; defualt = max(net)
#' @param edge.width.max The maximum edge width for visualization
#' @param label.edge Whether or not shows the label of edges
#' @param alpha.edge the transprency of edge
#' @param edge.label.color The color for single arrow
#' @param edge.label.cex The size of label for arrows
#' @param edge.curved Specifies whether to draw curved edges, or not.
#' This can be a logical or a numeric vector or scalar.
#' First the vector is replicated to have the same length as the number of
#' edges in the graph. Then it is interpreted for each edge separately.
#' A numeric value specifies the curvature of the edge; zero curvature means
#' straight edges, negative values means the edge bends clockwise, positive
#' values the opposite. TRUE means curvature 0.5, FALSE means curvature zero
#' @param shape The shape of the vertex, currently "circle", "square",
#' "csquare", "rectangle", "crectangle", "vrectangle", "pie" (see
#' vertex.shape.pie), "sphere", and "none" are supported, and only by the
#' plot.igraph command. "none" does not draw the vertices at all, although
#' vertex label are plotted (if given). See shapes for details about vertex
#' shapes and vertex.shape.pie for using pie charts as vertices.
#' @param layout The layout specification. It must be a call to a layout
#' specification function.
#' @param margin The amount of empty space below, over, at the left and right
#'  of the plot, it is a numeric vector of length four. Usually values between
#'  0 and 0.5 are meaningful, but negative values are also possible, that will
#'  make the plot zoom in to a part of the graph. If it is shorter than four
#'  then it is recycled.
#' @param arrow.width The width of arrows
#' @param arrow.size the size of arrow
# #' @param from,to,bidirection Deprecated. Use `sources.use`,`targets.use`
#' @param vertex.size Deprecated. Use `vertex.weight`
#' @importFrom igraph graph_from_adjacency_matrix ends E V layout_ in_circle
#' @importFrom grDevices recordPlot
#' @return  an object of class "recordedplot"
#' @export
netVisual_circle <-function(net, color.use = NULL,title.name = NULL, sources.use = NULL, targets.use = NULL, idents.use = NULL, remove.isolate = FALSE, top = 1,
                            weight.scale = FALSE, vertex.weight = 20, vertex.weight.max = NULL, vertex.size.max = NULL, vertex.label.cex=1,vertex.label.color= "black",
                            edge.weight.max = NULL, edge.width.max=8, alpha.edge = 0.6, label.edge = FALSE,edge.label.color='black',edge.label.cex=0.8,
                            edge.curved=0.2,shape='circle',layout=in_circle(), margin=0.2, vertex.size = NULL,
                            arrow.width=1,arrow.size = 0.2){
  if (!is.null(vertex.size)) {
    warning("'vertex.size' is deprecated. Use `vertex.weight`")
  }
  if (is.null(vertex.size.max)) {
    if (length(unique(vertex.weight)) == 1) {
      vertex.size.max <- 5
    } else {
      vertex.size.max <- 15
    }
  }
  options(warn = -1)
  thresh <- stats::quantile(net, probs = 1-top)
  net[net < thresh] <- 0

  if ((!is.null(sources.use)) | (!is.null(targets.use)) | (!is.null(idents.use)) ) {
    if (is.null(rownames(net))) {
      stop("The input weighted matrix should have rownames!")
    }
    cells.level <- rownames(net)
    df.net <- reshape2::melt(net, value.name = "value")
    colnames(df.net)[1:2] <- c("source","target")
    # keep the interactions associated with sources and targets of interest
    if (!is.null(sources.use)){
      if (is.numeric(sources.use)) {
        sources.use <- cells.level[sources.use]
      }
      df.net <- subset(df.net, source %in% sources.use)
    }
    if (!is.null(targets.use)){
      if (is.numeric(targets.use)) {
        targets.use <- cells.level[targets.use]
      }
      df.net <- subset(df.net, target %in% targets.use)
    }
    if (!is.null(idents.use)) {
      if (is.numeric(idents.use)) {
        idents.use <- cells.level[idents.use]
      }
      df.net <- filter(df.net, (source %in% idents.use) | (target %in% idents.use))
    }
    df.net$source <- factor(df.net$source, levels = cells.level)
    df.net$target <- factor(df.net$target, levels = cells.level)
    df.net$value[is.na(df.net$value)] <- 0
    net <- tapply(df.net[["value"]], list(df.net[["source"]], df.net[["target"]]), sum)
  }
  net[is.na(net)] <- 0


  if (remove.isolate) {
    idx1 <- which(Matrix::rowSums(net) == 0)
    idx2 <- which(Matrix::colSums(net) == 0)
    idx <- intersect(idx1, idx2)
    net <- net[-idx, ]
    net <- net[, -idx]
  }

  g <- graph_from_adjacency_matrix(net, mode = "directed", weighted = T)
  edge.start <- igraph::ends(g, es=igraph::E(g), names=FALSE)
  coords<-layout_(g,layout)
  if(nrow(coords)!=1){
    coords_scale=scale(coords)
  }else{
    coords_scale<-coords
  }
  if (is.null(color.use)) {
    color.use = scPalette(length(igraph::V(g)))
  }
  if (is.null(vertex.weight.max)) {
    vertex.weight.max <- max(vertex.weight)
  }
  vertex.weight <- vertex.weight/vertex.weight.max*vertex.size.max+5

  loop.angle<-ifelse(coords_scale[igraph::V(g),1]>0,atan(coords_scale[igraph::V(g),2]/coords_scale[igraph::V(g),1]),-pi+atan(coords_scale[igraph::V(g),2]/coords_scale[igraph::V(g),1]))
  igraph::V(g)$size<-vertex.weight
  igraph::V(g)$color<-color.use[igraph::V(g)]
  igraph::V(g)$frame.color <- color.use[igraph::V(g)]
  igraph::V(g)$label.color <- vertex.label.color
  igraph::V(g)$label.cex<-vertex.label.cex
  if(label.edge){
    igraph::E(g)$label<-igraph::E(g)$weight
    igraph::E(g)$label <- round(igraph::E(g)$label, digits = 1)
  }
  if (is.null(edge.weight.max)) {
    edge.weight.max <- max(igraph::E(g)$weight)
  }
  if (weight.scale == TRUE) {
    #E(g)$width<-0.3+edge.width.max/(max(E(g)$weight)-min(E(g)$weight))*(E(g)$weight-min(E(g)$weight))
    igraph::E(g)$width<- 0.3+igraph::E(g)$weight/edge.weight.max*edge.width.max
  }else{
    igraph::E(g)$width<-0.3+edge.width.max*igraph::E(g)$weight
  }

  igraph::E(g)$arrow.width<-arrow.width
  igraph::E(g)$arrow.size<-arrow.size
  igraph::E(g)$label.color<-edge.label.color
  igraph::E(g)$label.cex<-edge.label.cex
  igraph::E(g)$color<- grDevices::adjustcolor(igraph::V(g)$color[edge.start[,1]],alpha.edge)

  igraph::E(g)$loop.angle <- 0
  if(sum(edge.start[,2]==edge.start[,1])!=0){
    igraph::E(g)$loop.angle[which(edge.start[,2]==edge.start[,1])]<-loop.angle[edge.start[which(edge.start[,2]==edge.start[,1]),1]]
  }
  radian.rescale <- function(x, start=0, direction=1) {
    c.rotate <- function(x) (x + start) %% (2 * pi) * direction
    c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
  }
  label.locs <- radian.rescale(x=1:length(igraph::V(g)), direction=-1, start=0)
  label.dist <- vertex.weight/max(vertex.weight)+2
  plot(g,edge.curved=edge.curved,vertex.shape=shape,layout=coords_scale,margin=margin, vertex.label.dist=label.dist,
       vertex.label.degree=label.locs, vertex.label.family="Helvetica", edge.label.family="Helvetica") # "sans"
  if (!is.null(title.name)) {
    text(0,1.5,title.name, cex = 1.1)
  }
  # https://www.andrewheiss.com/blog/2016/12/08/save-base-graphics-as-pseudo-objects-in-r/
  # grid.echo()
  # gg <-  grid.grab()
  gg <- recordPlot()
  return(gg)
}



#' generate circle symbol
#'
#' @param coords coordinates of points
#' @param v vetex
#' @param params parameters
#' @importFrom graphics symbols
#' @return
mycircle <- function(coords, v=NULL, params) {
  vertex.color <- params("vertex", "color")
  if (length(vertex.color) != 1 && !is.null(v)) {
    vertex.color <- vertex.color[v]
  }
  vertex.size  <- 1/200 * params("vertex", "size")
  if (length(vertex.size) != 1 && !is.null(v)) {
    vertex.size <- vertex.size[v]
  }
  vertex.frame.color <- params("vertex", "frame.color")
  if (length(vertex.frame.color) != 1 && !is.null(v)) {
    vertex.frame.color <- vertex.frame.color[v]
  }
  vertex.frame.width <- params("vertex", "frame.width")
  if (length(vertex.frame.width) != 1 && !is.null(v)) {
    vertex.frame.width <- vertex.frame.width[v]
  }

  mapply(coords[,1], coords[,2], vertex.color, vertex.frame.color,
         vertex.size, vertex.frame.width,
         FUN=function(x, y, bg, fg, size, lwd) {
           symbols(x=x, y=y, bg=bg, fg=fg, lwd=lwd,
                   circles=size, add=TRUE, inches=FALSE)
         })
}


#' Spatial plot of cell-cell communication network
#'
#' Autocrine interactions are omitted on this plot. Group centroids may be not accurate for some data due to complex geometry.
#' The width of edges represent the strength of the communication.
#'
#' @param net A weighted matrix representing the connections
#' @param coordinates a data matrix in which each row gives the spatial locations/coordinates of each cell/spot
#' @param labels a vector giving the group label of each cell/spot. The length should be the same as the number of rows in `coordinates`
#' @param color.use Colors represent different cell groups
#' @param title.name the name of the title
#' @param sources.use a vector giving the index or the name of source cell groups
#' @param targets.use a vector giving the index or the name of target cell groups
#' @param idents.use a vector giving the index or the name of cell groups of interest
#' @param remove.isolate whether remove the isolate nodes in the communication network
#' @param remove.loop whether remove the self-loop in the communication network. Default: TRUE
#' @param top the fraction of interactions to show
#' @param weight.scale whether scale the weight
#' @param vertex.weight The weight of vertex: either a scale value or a vector
#' @param vertex.weight.max the maximum weight of vertex; defualt = max(vertex.weight)
#' @param vertex.size.max the maximum vertex size for visualization
#' @param vertex.label.cex The label size of vertex
#' @param vertex.label.color The color of label for vertex
#' @param edge.weight.max the maximum weight of edge; defualt = max(net)
#' @param edge.width.max The maximum edge width for visualization
#' @param alpha.edge the transprency of edge
#' @param edge.curved Specifies whether to draw curved edges, or not.
#' This can be a logical or a numeric vector or scalar.
#' First the vector is replicated to have the same length as the number of
#' edges in the graph. Then it is interpreted for each edge separately.
#' A numeric value specifies the curvature of the edge; zero curvature means
#' straight edges, negative values means the edge bends clockwise, positive
#' values the opposite. TRUE means curvature 0.5, FALSE means curvature zero
#' @param arrow.angle the angle of arrows
#' @param alpha.image the transparency of individual spots
# #' @param arrow.width The width of arrows
#' @param arrow.size the size of arrow
#' @param point.size the size of spots
#' @param legend.size the size of legend
#' @importFrom igraph graph_from_adjacency_matrix get.edgelist ends E V
#' @import ggplot2
#' @importFrom ggnetwork geom_nodetext_repel
#' @return  an object of ggplot
#' @export
netVisual_spatial <-function(net, coordinates, labels, color.use = NULL,title.name = NULL, sources.use = NULL, targets.use = NULL, idents.use = NULL, remove.isolate = FALSE, remove.loop = TRUE, top = 1,
                             weight.scale = FALSE, vertex.weight = 20, vertex.weight.max = NULL, vertex.size.max = NULL, vertex.label.cex = 3,vertex.label.color= "black",
                             edge.weight.max = NULL, edge.width.max=8, edge.curved=0.2, alpha.edge = 0.6, arrow.angle = 5, arrow.size = 0.2, alpha.image = 0.15, point.size = 1, legend.size = 3){
  cells.level <- rownames(net)
  if (ncol(coordinates) == 2) {
    colnames(coordinates) <- c("x_cent","y_cent")
    temp_coordinates = coordinates
    coordinates[,1] = temp_coordinates[,2]
    coordinates[,2] = temp_coordinates[,1]
  } else {
    stop("Please check the input 'coordinates' and make sure it is a two column matrix.")
  }
  num_cluster <- length(cells.level)
  node_coords <- matrix(0, nrow = num_cluster, ncol = 2)
  for (i in c(1:num_cluster)) {
    node_coords[i,1] <- median(coordinates[as.character(labels) == cells.level[i], 1])
    node_coords[i,2] <- median(coordinates[as.character(labels) == cells.level[i], 2])
  }
  rownames(node_coords) <- cells.level

  if (is.null(vertex.size.max)) {
    if (length(unique(vertex.weight)) == 1) {
      vertex.size.max <- 5
    } else {
      vertex.size.max <- 15
    }
  }
  options(warn = -1)
  thresh <- stats::quantile(net, probs = 1-top)
  net[net < thresh] <- 0

  if ((!is.null(sources.use)) | (!is.null(targets.use)) | (!is.null(idents.use)) ) {
    if (is.null(rownames(net))) {
      stop("The input weighted matrix should have rownames!")
    }
    df.net <- reshape2::melt(net, value.name = "value")
    colnames(df.net)[1:2] <- c("source","target")
    # keep the interactions associated with sources and targets of interest
    if (!is.null(sources.use)){
      if (is.numeric(sources.use)) {
        sources.use <- cells.level[sources.use]
      }
      df.net <- subset(df.net, source %in% sources.use)
    }
    if (!is.null(targets.use)){
      if (is.numeric(targets.use)) {
        targets.use <- cells.level[targets.use]
      }
      df.net <- subset(df.net, target %in% targets.use)
    }
    if (!is.null(idents.use)) {
      if (is.numeric(idents.use)) {
        idents.use <- cells.level[idents.use]
      }
      df.net <- filter(df.net, (source %in% idents.use) | (target %in% idents.use))
    }
    df.net$source <- factor(df.net$source, levels = cells.level)
    df.net$target <- factor(df.net$target, levels = cells.level)
    df.net$value[is.na(df.net$value)] <- 0
    net <- tapply(df.net[["value"]], list(df.net[["source"]], df.net[["target"]]), sum)
  }
  net[is.na(net)] <- 0


  if (remove.loop) {
    diag(net) <- 0
  }
  if (remove.isolate) {
    idx1 <- which(Matrix::rowSums(net) == 0)
    idx2 <- which(Matrix::colSums(net) == 0)
    idx <- intersect(idx1, idx2)
    net <- net[-idx, ]
    net <- net[, -idx]
    node_coords <- node_coords[-idx, ]
    cells.level <- cells.level[-idx]
  }

  g <- graph_from_adjacency_matrix(net, mode = "directed", weighted = T)
  edgelist <- get.edgelist(g)
  # loop_curve = c()
  # for (i in c(1:nrow(edgelist))) {
  #   if (edgelist[i,1] == edgelist[i,2]){
  #     loop_curve  = c(loop_curve ,i)
  #   }
  # }
  # edgelist <- edgelist[-loop_curve,]

  edges <- data.frame(node_coords[edgelist[,1],,drop = FALSE], node_coords[edgelist[,2],,drop = FALSE])
  colnames(edges) <- c("X1","Y1","X2","Y2")
  node_coords = data.frame(node_coords)
  node_idents = factor(cells.level, levels = cells.level)
  node_family = data.frame(node_coords,node_idents)
  if (is.null(color.use)) {
    color.use = scPalette(length(igraph::V(g)))
    names(color.use) <- cells.level
  }
  if (is.null(vertex.weight.max)) {
    vertex.weight.max <- max(vertex.weight)
  }
  vertex.weight <- vertex.weight/vertex.weight.max*vertex.size.max+1
  if (is.null(edge.weight.max)) {
    edge.weight.max <- max(igraph::E(g)$weight)
  }
  # width of edge
  if (weight.scale == TRUE) {
    igraph::E(g)$width<- 0.3+igraph::E(g)$weight/edge.weight.max*edge.width.max
  }else{
    igraph::E(g)$width<-0.3+edge.width.max*igraph::E(g)$weight
  }

  gg <- ggplot(data=node_family,aes(X1, X2)) +
    geom_curve(aes(x=X1, y=Y1, xend = X2, yend = Y2), data=edges, size = igraph::E(g)$width, curvature = edge.curved, alpha = alpha.edge, arrow = arrow(angle = arrow.angle, type = "closed",length = unit(arrow.size, "inches")),colour=color.use[edgelist[,1]]) +
    geom_point(aes(X1, X2,colour = node_idents), data=node_family, size = vertex.weight,show.legend = TRUE) +scale_color_manual(values = color.use) +
    guides(color = guide_legend(override.aes = list(size=legend.size))) +
    xlab(NULL) + ylab(NULL)  +
    coord_fixed() + theme(aspect.ratio = 1) +
    theme(panel.background = element_blank(),panel.border = element_blank(),axis.text=element_blank(),axis.ticks = element_blank(),legend.title = element_blank())+ theme(legend.key = element_blank())
  gg <- gg + geom_point(aes(x_cent, y_cent), data = coordinates,colour = color.use[labels],alpha = alpha.image, size = point.size, show.legend = FALSE)
  gg <- gg + scale_y_reverse()
  if (vertex.label.cex > 0){
    gg <- gg + ggnetwork::geom_nodetext_repel(aes(label = node_idents), color="black", size = vertex.label.cex)
  }
  if (!is.null(title.name)){
    gg <- gg + ggtitle(title.name) + theme(plot.title = element_text(hjust = 0.5, vjust = 0))
  }

  gg
  return(gg)

}







#' Circle plot showing differential cell-cell communication network between two datasets
#'
#' The width of edges represent the relative number of interactions or interaction strength.
#' Red (or blue) colored edges represent increased (or decreased) signaling in the second dataset compared to the first one.
#'
#' @param object A merged CellChat objects
#' @param comparison a numerical vector giving the datasets for comparison in object.list; e.g., comparison = c(1,2)
#' @param measure "count" or "weight". "count": comparing the number of interactions; "weight": comparing the total interaction weights (strength)
#' @param color.use Colors represent different cell groups
#' @param color.edge Colors for indicating whether the signaling is increased (`color.edge[1]`) or decreased (`color.edge[2]`)
#' @param title.name the name of the title
#' @param sources.use a vector giving the index or the name of source cell groups
#' @param targets.use a vector giving the index or the name of target cell groups.
#' @param remove.isolate whether remove the isolate nodes in the communication network
#' @param top the fraction of interactions to show
#' @param weight.scale whether scale the weight
#' @param vertex.weight The weight of vertex: either a scale value or a vector
#' @param vertex.weight.max the maximum weight of vertex; defualt = max(vertex.weight)
#' @param vertex.size.max the maximum vertex size for visualization
#' @param vertex.label.cex The label size of vertex
#' @param vertex.label.color The color of label for vertex
#' @param edge.weight.max the maximum weight of edge; defualt = max(net)
#' @param edge.width.max The maximum edge width for visualization
#' @param label.edge Whether or not shows the label of edges
#' @param alpha.edge the transprency of edge
#' @param edge.label.color The color for single arrow
#' @param edge.label.cex The size of label for arrows
#' @param edge.curved Specifies whether to draw curved edges, or not.
#' This can be a logical or a numeric vector or scalar.
#' First the vector is replicated to have the same length as the number of
#' edges in the graph. Then it is interpreted for each edge separately.
#' A numeric value specifies the curvature of the edge; zero curvature means
#' straight edges, negative values means the edge bends clockwise, positive
#' values the opposite. TRUE means curvature 0.5, FALSE means curvature zero
#' @param shape The shape of the vertex, currently "circle", "square",
#' "csquare", "rectangle", "crectangle", "vrectangle", "pie" (see
#' vertex.shape.pie), "sphere", and "none" are supported, and only by the
#' plot.igraph command. "none" does not draw the vertices at all, although
#' vertex label are plotted (if given). See shapes for details about vertex
#' shapes and vertex.shape.pie for using pie charts as vertices.
#' @param layout The layout specification. It must be a call to a layout
#' specification function.
#' @param margin The amount of empty space below, over, at the left and right
#'  of the plot, it is a numeric vector of length four. Usually values between
#'  0 and 0.5 are meaningful, but negative values are also possible, that will
#'  make the plot zoom in to a part of the graph. If it is shorter than four
#'  then it is recycled.
#' @param arrow.width The width of arrows
#' @param arrow.size the size of arrow
# #' @param from,to,bidirection Deprecated. Use `sources.use`,`targets.use`
# #' @param vertex.size Deprecated. Use `vertex.weight`
#' @importFrom igraph graph_from_adjacency_matrix ends E V layout_ in_circle
#' @importFrom grDevices recordPlot
#' @return  an object of class "recordedplot"
#' @export
netVisual_diffInteraction <- function(object, comparison = c(1,2), measure = c("count", "weight", "count.merged", "weight.merged"), color.use = NULL, color.edge = c('#b2182b','#2166ac'), title.name = NULL, sources.use = NULL, targets.use = NULL, remove.isolate = FALSE, top = 1,
                                      weight.scale = FALSE, vertex.weight = 20, vertex.weight.max = NULL, vertex.size.max = 15, vertex.label.cex=1,vertex.label.color= "black",
                                      edge.weight.max = NULL, edge.width.max=8, alpha.edge = 0.6, label.edge = FALSE,edge.label.color='black',edge.label.cex=0.8,
                                      edge.curved=0.2,shape='circle',layout=in_circle(), margin=0.2,
                                      arrow.width=1,arrow.size = 0.2){
  options(warn = -1)
  measure <- match.arg(measure)
  obj1 <- object@net[[comparison[1]]][[measure]]
  obj2 <- object@net[[comparison[2]]][[measure]]
  net.diff <- obj2 - obj1
  if (measure %in% c("count", "count.merged")) {
    if (is.null(title.name)) {
      title.name = "Differential number of interactions"
    }
  } else if (measure %in% c("weight", "weight.merged")) {
    if (is.null(title.name)) {
      title.name = "Differential interaction strength"
    }
  }
  net <- net.diff
  if ((!is.null(sources.use)) | (!is.null(targets.use))) {
    df.net <- reshape2::melt(net, value.name = "value")
    colnames(df.net)[1:2] <- c("source","target")
    # keep the interactions associated with sources and targets of interest
    if (!is.null(sources.use)){
      if (is.numeric(sources.use)) {
        sources.use <- rownames(net.diff)[sources.use]
      }
      df.net <- subset(df.net, source %in% sources.use)
    }
    if (!is.null(targets.use)){
      if (is.numeric(targets.use)) {
        targets.use <- rownames(net.diff)[targets.use]
      }
      df.net <- subset(df.net, target %in% targets.use)
    }
    cells.level <- rownames(net.diff)
    df.net$source <- factor(df.net$source, levels = cells.level)
    df.net$target <- factor(df.net$target, levels = cells.level)
    df.net$value[is.na(df.net$value)] <- 0
    net <- tapply(df.net[["value"]], list(df.net[["source"]], df.net[["target"]]), sum)
    net[is.na(net)] <- 0
  }

  if (remove.isolate) {
    idx1 <- which(Matrix::rowSums(net) == 0)
    idx2 <- which(Matrix::colSums(net) == 0)
    idx <- intersect(idx1, idx2)
    net <- net[-idx, ]
    net <- net[, -idx]
  }

  net[abs(net) < stats::quantile(abs(net), probs = 1-top)] <- 0

  g <- graph_from_adjacency_matrix(net, mode = "directed", weighted = T)
  edge.start <- igraph::ends(g, es=igraph::E(g), names=FALSE)
  coords<-layout_(g,layout)
  if(nrow(coords)!=1){
    coords_scale=scale(coords)
  }else{
    coords_scale<-coords
  }
  if (is.null(color.use)) {
    color.use = scPalette(length(igraph::V(g)))
  }
  if (is.null(vertex.weight.max)) {
    vertex.weight.max <- max(vertex.weight)
  }
  vertex.weight <- vertex.weight/vertex.weight.max*vertex.size.max+5

  loop.angle<-ifelse(coords_scale[igraph::V(g),1]>0,-atan(coords_scale[igraph::V(g),2]/coords_scale[igraph::V(g),1]),pi-atan(coords_scale[igraph::V(g),2]/coords_scale[igraph::V(g),1]))
  igraph::V(g)$size<-vertex.weight
  igraph::V(g)$color<-color.use[igraph::V(g)]
  igraph::V(g)$frame.color <- color.use[igraph::V(g)]
  igraph::V(g)$label.color <- vertex.label.color
  igraph::V(g)$label.cex<-vertex.label.cex
  if(label.edge){
    igraph::E(g)$label<-igraph::E(g)$weight
    igraph::E(g)$label <- round(igraph::E(g)$label, digits = 1)
  }
  igraph::E(g)$arrow.width<-arrow.width
  igraph::E(g)$arrow.size<-arrow.size
  igraph::E(g)$label.color<-edge.label.color
  igraph::E(g)$label.cex<-edge.label.cex
  #igraph::E(g)$color<- grDevices::adjustcolor(igraph::V(g)$color[edge.start[,1]],alpha.edge)
  igraph::E(g)$color <- ifelse(igraph::E(g)$weight > 0, color.edge[1],color.edge[2])
  igraph::E(g)$color <- grDevices::adjustcolor(igraph::E(g)$color, alpha.edge)

  igraph::E(g)$weight <- abs(igraph::E(g)$weight)

  if (is.null(edge.weight.max)) {
    edge.weight.max <- max(igraph::E(g)$weight)
  }
  if (weight.scale == TRUE) {
    #E(g)$width<-0.3+edge.width.max/(max(E(g)$weight)-min(E(g)$weight))*(E(g)$weight-min(E(g)$weight))
    igraph::E(g)$width<- 0.3+igraph::E(g)$weight/edge.weight.max*edge.width.max
  }else{
    igraph::E(g)$width<-0.3+edge.width.max*igraph::E(g)$weight
  }


  if(sum(edge.start[,2]==edge.start[,1])!=0){
    igraph::E(g)$loop.angle[which(edge.start[,2]==edge.start[,1])]<-loop.angle[edge.start[which(edge.start[,2]==edge.start[,1]),1]]
  }
  radian.rescale <- function(x, start=0, direction=1) {
    c.rotate <- function(x) (x + start) %% (2 * pi) * direction
    c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
  }
  label.locs <- radian.rescale(x=1:length(igraph::V(g)), direction=-1, start=0)
  label.dist <- vertex.weight/max(vertex.weight)+2
  plot(g,edge.curved=edge.curved,vertex.shape=shape,layout=coords_scale,margin=margin, vertex.label.dist=label.dist,
       vertex.label.degree=label.locs, vertex.label.family="Helvetica", edge.label.family="Helvetica") # "sans"
  if (!is.null(title.name)) {
    text(0,1.5,title.name, cex = 1.1)
  }
  # https://www.andrewheiss.com/blog/2016/12/08/save-base-graphics-as-pseudo-objects-in-r/
  # grid.echo()
  # gg <-  grid.grab()
  gg <- recordPlot()
  return(gg)
}


#' Visualization of network using heatmap
#'
#' This heatmap can be used to show differential number of interactions or interaction strength in the cell-cell communication network between two datasets;
#' the number of interactions or interaction strength in a single dataset
#' the inferred cell-cell communication network in single dataset, defined by `signaling`
#'
#' When show differential number of interactions or interaction strength in the cell-cell communication network between two datasets, the width of edges represent the relative number of interactions or interaction strength.
#' Red (or blue) colored edges represent increased (or decreased) signaling in the second dataset compared to the first one.
#'
#' The top colored bar plot represents the sum of column of values displayed in the heatmap. The right colored bar plot represents the sum of row of values.
#'
#'
#' @param object A merged CellChat object or a single CellChat object
#' @param comparison a numerical vector giving the datasets for comparison in object.list; e.g., comparison = c(1,2)
#' @param measure "count" or "weight". "count": comparing the number of interactions; "weight": comparing the total interaction weights (strength)
#' @param signaling a character vector giving the name of signaling networks in a single CellChat object
#' @param slot.name the slot name of object. Set is to be "netP" if input signaling is a pathway name; Set is to be "net" if input signaling is a ligand-receptor pair
#' @param thresh threshold of the p-value for determining significant interaction
#' @param color.use the character vector defining the color of each cell group
#' @param color.heatmap A vector of two colors corresponding to max/min values, or a color name in brewer.pal only when the data in the heatmap do not contain negative values
#' @param title.name the name of the title
#' @param width width of heatmap
#' @param height height of heatmap
#' @param font.size fontsize in heatmap
#' @param font.size.title font size of the title
#' @param x.lab.rot do rotation for the x-tick labels
#' @param cluster.rows whether cluster rows
#' @param cluster.cols whether cluster columns
#' @param sources.use a vector giving the index or the name of source cell groups
#' @param targets.use a vector giving the index or the name of target cell groups.
#' @param remove.isolate whether to remove the isolate nodes in the communication network
#' @param remove.diagonal whether to remove diagonal elements, i.e., ignore the communication from a cell to itself.
#' @param row.show,col.show a vector giving the index or the name of row or columns to show in the heatmap
#' @importFrom methods slot
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#' @importFrom ComplexHeatmap Heatmap HeatmapAnnotation anno_barplot rowAnnotation
#' @return  an object of ComplexHeatmap
#' @export
netVisual_heatmap <- function(
    object,
    comparison = c(1, 2),
    measure = c("count", "weight"),
    signaling = NULL,
    slot.name = c("netP", "net"),
    thresh = 0.05,
    color.use = NULL,
    color.heatmap = c("#2166ac", "#b2182b"),
    title.name = NULL,
    width = NULL,
    height = NULL,
    font.size = 8,
    font.size.title = 10,
    x.lab.rot = 90,
    cluster.rows = FALSE,
    cluster.cols = FALSE,
    sources.use = NULL,
    targets.use = NULL,
    remove.isolate = FALSE,
    remove.diagonal = FALSE,
    row.show = NULL,
    col.show = NULL
) {
  # obj1 <- object.list[[comparison[1]]]
  # obj2 <- object.list[[comparison[2]]]
  if (!is.null(measure)) {
    measure <- match.arg(measure)
  }
  slot.name <- match.arg(slot.name)
  if (object@options$mode=="merged") {
    message("Do heatmap based on a merged object \n")
    obj1 <- object@net[[comparison[1]]][[measure]]
    obj2 <- object@net[[comparison[2]]][[measure]]
    net.diff <- obj2 - obj1

    if (measure == "count") {
      if (is.null(title.name)) {
        title.name = "Differential number of interactions"
      }
    } else if (measure == "weight") {
      if (is.null(title.name)) {
        title.name = "Differential interaction strength"
      }
    }
    legend.name = "Relative values"
  } else {
    message("Do heatmap based on a single object \n")
    if (!is.null(signaling)) {
      prob <- slot(object, slot.name)$prob
      if (slot.name == "net") {
        prob[object@net$pval > thresh] <- 0
      }
      net.diff <- prob[,,signaling]
      if (is.null(title.name)) {
        title.name = paste0(signaling, " signaling network")
      }
      legend.name <- "Communication Prob."
    } else if (!is.null(measure)) {
      net.diff <- object@net[[measure]]
      if (measure == "count") {
        if (is.null(title.name)) {
          title.name = "Number of interactions"
        }
      } else if (measure == "weight") {
        if (is.null(title.name)) {
          title.name = "Interaction strength"
        }
      }
      legend.name <- title.name
    }
  }

  net <- net.diff


  if ((!is.null(sources.use)) | (!is.null(targets.use))) {
    df.net <- reshape2::melt(net, value.name = "value")
    colnames(df.net)[1:2] <- c("source","target")
    # keep the interactions associated with sources and targets of interest
    if (!is.null(sources.use)){
      if (is.numeric(sources.use)) {
        sources.use <- rownames(net.diff)[sources.use]
      }
      df.net <- subset(df.net, source %in% sources.use)
    }
    if (!is.null(targets.use)){
      if (is.numeric(targets.use)) {
        targets.use <- rownames(net.diff)[targets.use]
      }
      df.net <- subset(df.net, target %in% targets.use)
    }
    cells.level <- rownames(net.diff)
    df.net$source <- factor(df.net$source, levels = cells.level)
    df.net$target <- factor(df.net$target, levels = cells.level)
    df.net$value[is.na(df.net$value)] <- 0
    net <- tapply(df.net[["value"]], list(df.net[["source"]], df.net[["target"]]), sum)
  }
  net[is.na(net)] <- 0

  if (remove.isolate) {
    idx1 <- which(Matrix::rowSums(net) == 0)
    idx2 <- which(Matrix::colSums(net) == 0)
    idx <- intersect(idx1, idx2)
    net <- net[-idx, ]
    net <- net[, -idx]
    if (length(idx) > 0) {
      net <- net[-idx, ]
      net <- net[, -idx]
    }
  }
  if (remove.diagonal) {
    diag(net) <- 0
  }

  mat <- net
  if (is.null(color.use)) {
    color.use <- scPalette(ncol(mat))
  }
  names(color.use) <- colnames(mat)

  if (!is.null(row.show)) {
    mat <- mat[row.show, ]
  }
  if (!is.null(col.show)) {
    mat <- mat[ ,col.show]
    color.use <- color.use[col.show]
  }


  if (min(mat) < 0) {
    color.heatmap.use = colorRamp3(c(min(mat), 0, max(mat)), c(color.heatmap[1], "#f7f7f7", color.heatmap[2]))
    colorbar.break <- c(round(min(mat, na.rm = T), digits = nchar(sub(".*\\.(0*).*","\\1",min(mat, na.rm = T)))+1), 0, round(max(mat, na.rm = T), digits = nchar(sub(".*\\.(0*).*","\\1",max(mat, na.rm = T)))+1))
    # color.heatmap.use = colorRamp3(c(seq(min(mat), -(max(mat)-min(max(mat)))/9, length.out = 4), 0, seq((max(mat)-min(max(mat)))/9, max(mat), length.out = 4)), RColorBrewer::brewer.pal(n = 9, name = color.heatmap))
  } else {
    if (length(color.heatmap) == 3) {
      color.heatmap.use = colorRamp3(c(0, min(mat), max(mat)), color.heatmap)
    } else if (length(color.heatmap) == 2) {
      color.heatmap.use = colorRamp3(c(min(mat), max(mat)), color.heatmap)
    } else if (length(color.heatmap) == 1) {
      color.heatmap.use = grDevices::colorRampPalette((RColorBrewer::brewer.pal(n = 9, name = color.heatmap)))(100)
    }
    colorbar.break <- c(round(min(mat, na.rm = T), digits = nchar(sub(".*\\.(0*).*","\\1",min(mat, na.rm = T)))+1), round(max(mat, na.rm = T), digits = nchar(sub(".*\\.(0*).*","\\1",max(mat, na.rm = T)))+1))
  }
  # col_fun(as.vector(mat))

  df<- data.frame(group = colnames(mat)); rownames(df) <- colnames(mat)
  col_annotation <- HeatmapAnnotation(df = df, col = list(group = color.use),which = "column",
                                      show_legend = FALSE, show_annotation_name = FALSE,
                                      simple_anno_size = grid::unit(0.2, "cm"))
  row_annotation <- HeatmapAnnotation(df = df, col = list(group = color.use), which = "row",
                                      show_legend = FALSE, show_annotation_name = FALSE,
                                      simple_anno_size = grid::unit(0.2, "cm"))

  ha1 = rowAnnotation(Strength = anno_barplot(rowSums(abs(mat)), border = FALSE,gp = gpar(fill = color.use, col=color.use)), show_annotation_name = FALSE)
  ha2 = HeatmapAnnotation(Strength = anno_barplot(colSums(abs(mat)), border = FALSE,gp = gpar(fill = color.use, col=color.use)), show_annotation_name = FALSE)

  if (sum(mat) == 0) {
    stop("No significant cell-cell communication is inferred!",'\n')
  } else if (sum(abs(unique(as.vector(mat))) > 0) == 1) {
    color.heatmap.use = c("white", color.heatmap.use)
  } else {
    mat[mat == 0] <- NA
  }
  ht1 = Heatmap(mat, col = color.heatmap.use, na_col = "white", name = legend.name,
                bottom_annotation = col_annotation, left_annotation =row_annotation, top_annotation = ha2, right_annotation = ha1,
                cluster_rows = cluster.rows,cluster_columns = cluster.rows,
                row_names_side = "left",row_names_rot = 0,row_names_gp = gpar(fontsize = font.size),column_names_gp = gpar(fontsize = font.size),
                # width = unit(width, "cm"), height = unit(height, "cm"),
                column_title = title.name,column_title_gp = gpar(fontsize = font.size.title),column_names_rot = x.lab.rot,
                row_title = "Sources (Sender)",row_title_gp = gpar(fontsize = font.size.title),row_title_rot = 90,
                heatmap_legend_param = list(title_gp = gpar(fontsize = 8, fontface = "plain"),title_position = "leftcenter-rot",
                                            border = NA, #at = colorbar.break,
                                            legend_height = unit(20, "mm"),labels_gp = gpar(fontsize = 8),grid_width = unit(2, "mm"))
  )
  draw(ht1)
  # return(ht1)
}


#' Visualization of (differential) number of interactions
#'
#' @param object A merged CellChat object or a single CellChat object
#' @param comparison a numerical vector giving the datasets for comparison in object.list; e.g., comparison = c(1,2)
#' @param measure "count" or "weight". "count": comparing the number of interactions; "weight": comparing the total interaction weights (strength)
#' @param sources.use a vector giving the index or the name of source cell groups
#' @param targets.use a vector giving the index or the name of target cell groups.
#' @param invert.source,invert.target retain the complementary set
#' @param signaling a character vector giving the name of signaling networks in a single CellChat object
#' @param slot.name the slot name of object. Set is to be "netP" if input signaling is a pathway name; Set is to be "net" if input signaling is a ligand-receptor pair
#' @param color.use the character vector defining the color of each cell group
#' @param title.name the name of the title
#' @param x.lab.rot do rotation for the x-tick labels
#' @param ... Parameters passing to `barplot_internal`
#' @importFrom methods slot
#' @return  an object of ggplot
#' @export
netVisual_barplot <- function(object, comparison = c(1,2), measure = c("count", "weight"), sources.use = NULL, targets.use = NULL, invert.source = FALSE, invert.target = FALSE,signaling = NULL, slot.name = c("netP", "net"), color.use = NULL,
                              title.name = NULL,x.lab.rot = FALSE,...){
  if (!is.null(measure)) {
    measure <- match.arg(measure)
  }
  slot.name <- match.arg(slot.name)
  if (is.list(object@net[[1]])) {
    message("Show differential number of interactions based on a merged object \n")
    obj1 <- object@net[[comparison[1]]][[measure]]
    obj2 <- object@net[[comparison[2]]][[measure]]
    net.diff <- obj2 - obj1

    if (measure == "count") {
      if (is.null(title.name)) {
        title.name = "Differential number of interactions"
      }
    } else if (measure == "weight") {
      if (is.null(title.name)) {
        title.name = "Differential interaction strength"
      }
    }
  } else {
    message("Show number of interactions based on a single object \n")
    if (!is.null(signaling)) {
      net.diff <- slot(object, slot.name)$prob[,,signaling]
      if (is.null(title.name)) {
        title.name = paste0(signaling, " signaling network")
      }
    } else if (!is.null(measure)) {
      net.diff <- object@net[[measure]]
      if (measure == "count") {
        if (is.null(title.name)) {
          title.name = "Number of interactions"
        }
      } else if (measure == "weight") {
        if (is.null(title.name)) {
          title.name = "Interaction strength"
        }
      }
    }
  }

  net <- net.diff
  cells.level <- rownames(net.diff)

  if ((!is.null(sources.use)) | (!is.null(targets.use))) {
    df.net <- reshape2::melt(net, value.name = "value")
    colnames(df.net)[1:2] <- c("source","target")
    # keep the interactions associated with sources and targets of interest
    if (!is.null(sources.use)){
      if (is.numeric(sources.use)) {
        sources.use <- rownames(net.diff)[sources.use]
      }
      if (invert.source == TRUE) {
        sources.use <- setdiff(rownames(net.diff), sources.use)
      }
      df.net <- subset(df.net, source %in% sources.use)
    }
    if (!is.null(targets.use)){
      if (is.numeric(targets.use)) {
        targets.use <- rownames(net.diff)[targets.use]
      }
      if (invert.target == TRUE) {
        targets.use <- setdiff(rownames(net.diff), targets.use)
      }
      df.net <- subset(df.net, target %in% targets.use)
    }

    df.net$source <- factor(df.net$source, levels = cells.level[cells.level %in% unique(df.net$source)])
    df.net$target <- factor(df.net$target, levels = cells.level[cells.level %in% unique(df.net$target)])
  }

  if (is.null(color.use)) {
    color.use <- scPalette(length(cells.level))
  }
  names(color.use) <- cells.level
  color.use <- color.use[cells.level %in% unique(df.net$target)]

  gg <- barplot_internal(df.net, x = "target", y = "value", fill = "target", color.use = color.use, title.name = title.name,x.lab.rot = x.lab.rot,...)

  return(gg)

}


#' Show all the significant interactions (L-R pairs) from some cell groups to other cell groups
#'
#' The dot color and size represent the calculated communication probability and p-values.
#'
#' @param object CellChat object
#' @param sources.use a vector giving the index or the name of source cell groups
#' @param targets.use a vector giving the index or the name of target cell groups.
#' @param signaling a character vector giving the name of signaling pathways of interest
#' @param pairLR.use a data frame consisting of one column named either "interaction_name" or "pathway_name", defining the interactions of interest and the order of L-R on y-axis
#' @param sort.by.source,sort.by.target,sort.by.source.priority set the order of interacting cell pairs on x-axis; please check examples for details
#' @param color.heatmap A character string or vector indicating the colormap option to use. It can be the avaibale color palette in viridis_pal() or brewer.pal()
#' @param direction Sets the order of colors in the scale. If 1, the default colors are used. If -1, the order of colors is reversed.
#' @param n.colors number of basic colors to generate from color palette
#' @param thresh threshold of the p-value for determining significant interaction
#' @param comparison a numerical vector giving the datasets for comparison in the merged object; e.g., comparison = c(1,2)
#' @param group a numerical vector giving the group information of different datasets; e.g., group = c(1,2,2)
#' @param remove.isolate whether remove the entire empty column, i.e., communication between certain cell groups
#' @param max.dataset a scale, keep the communications with highest probability in max.dataset (i.e., certrain condition)
#' @param min.dataset a scale, keep the communications with lowest probability in min.dataset (i.e., certrain condition)
#' @param min.quantile,max.quantile minimum and maximum quantile cutoff values for the colorbar, may specify quantile in [0,1]
#' @param line.on whether add vertical line when doing comparison analysis for the merged object
#' @param line.size size of vertical line if added
#' @param color.text.use whether color the xtick labels according to the dataset origin when doing comparison analysis
#' @param color.text the colors for xtick labels according to the dataset origin when doing comparison analysis
#' @param title.name main title of the plot
#' @param font.size,font.size.title font size of all the text and the title name
#' @param show.legend whether show legend
#' @param grid.on,color.grid whether add grid
#' @param angle.x,vjust.x,hjust.x parameters for adjusting the rotation of xtick labels
#' @param return.data whether return the data.frame for replotting
#'
#' @return
#' @export
#'
#' @examples
#'\dontrun{
#' # show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
#' netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11), remove.isolate = FALSE)
#'
#' # show all the significant interactions (L-R pairs) associated with certain signaling pathways
#' netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11), signaling = c("CCL","CXCL"))
#'
#' # show all the significant interactions (L-R pairs) based on user's input (defined by `pairLR.use`; the order of L-R is also based on user's input)
#' pairLR.use <- extractEnrichedLR(cellchat, signaling = c("CCL","CXCL","FGF"))
#' netVisual_bubble(cellchat, sources.use = c(3,4), targets.use = c(5:8), pairLR.use = pairLR.use, remove.isolate = TRUE)
#'
#' # set the order of interacting cell pairs on x-axis
#' # (1) Default: first sort cell pairs based on the appearance of sources in levels(object@idents), and then based on the appearance of targets in levels(object@idents)
#' # (2) sort cell pairs based on the targets.use defined by users
#' netVisual_bubble(cellchat, targets.use = c("LC","Inflam. DC","cDC2","CD40LG+ TC"), pairLR.use = pairLR.use, remove.isolate = TRUE, sort.by.target = T)
#' # (3) sort cell pairs based on the sources.use defined by users
#' netVisual_bubble(cellchat, sources.use = c("FBN1+ FIB","APOE+ FIB","Inflam. FIB"), pairLR.use = pairLR.use, remove.isolate = TRUE, sort.by.source = T)
#' # (4) sort cell pairs based on the sources.use and then targets.use defined by users
#' netVisual_bubble(cellchat, sources.use = c("FBN1+ FIB","APOE+ FIB","Inflam. FIB"), targets.use = c("LC","Inflam. DC","cDC2","CD40LG+ TC"), pairLR.use = pairLR.use, remove.isolate = TRUE, sort.by.source = T, sort.by.target = T)
#' # (5) sort cell pairs based on the targets.use and then sources.use defined by users
#' netVisual_bubble(cellchat, sources.use = c("FBN1+ FIB","APOE+ FIB","Inflam. FIB"), targets.use = c("LC","Inflam. DC","cDC2","CD40LG+ TC"), pairLR.use = pairLR.use, remove.isolate = TRUE, sort.by.source = T, sort.by.target = T, sort.by.source.priority = FALSE)
#'
#'# show all the increased interactions in the second dataset compared to the first dataset
#' netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:8), remove.isolate = TRUE, max.dataset = 2)
#'
#'# show all the decreased interactions in the second dataset compared to the first dataset
#' netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:8), remove.isolate = TRUE, max.dataset = 1)
#'}
netVisual_bubble <- function(object, sources.use = NULL, targets.use = NULL, signaling = NULL, pairLR.use = NULL, sort.by.source = FALSE, sort.by.target = FALSE, sort.by.source.priority = TRUE, color.heatmap = c("Spectral","viridis"), n.colors = 10, direction = -1, thresh = 0.05,
                             comparison = NULL, group = NULL, remove.isolate = FALSE, max.dataset = NULL, min.dataset = NULL,
                             min.quantile = 0, max.quantile = 1, line.on = TRUE, line.size = 0.2, color.text.use = TRUE, color.text = NULL,
                             title.name = NULL, font.size = 10, font.size.title = 10, show.legend = TRUE,
                             grid.on = TRUE, color.grid = "grey90", angle.x = 90, vjust.x = NULL, hjust.x = NULL,
                             return.data = FALSE){
  color.heatmap <- match.arg(color.heatmap)
  if (is.list(object@net[[1]])) {
    message("Comparing communications on a merged object \n")
  } else {
    message("Comparing communications on a single object \n")
  }
  if (is.null(vjust.x) | is.null(hjust.x)) {
    angle=c(0, 45, 90)
    hjust=c(0, 1, 1)
    vjust=c(0, 1, 0.5)
    vjust.x = vjust[angle == angle.x]
    hjust.x = hjust[angle == angle.x]
  }
  if (length(color.heatmap) == 1) {
    color.use <- tryCatch({
      RColorBrewer::brewer.pal(n = n.colors, name = color.heatmap)
    }, error = function(e) {
      scales::viridis_pal(option = color.heatmap, direction = -1)(n.colors)
    })
  } else {
    color.use <- color.heatmap
  }
  if (direction == -1) {
    color.use <- rev(color.use)
  }

  if (is.null(comparison)) {
    cells.level <- levels(object@idents)
    if (is.numeric(sources.use)) {
      sources.use <- cells.level[sources.use]
    }
    if (is.numeric(targets.use)) {
      targets.use <- cells.level[targets.use]
    }
    df.net <- subsetCommunication(object, slot.name = "net",
                                  sources.use = sources.use, targets.use = targets.use,
                                  signaling = signaling,
                                  pairLR.use = pairLR.use,
                                  thresh = thresh)
    df.net$source.target <- paste(df.net$source, df.net$target, sep = " -> ")
    source.target <- paste(rep(sources.use, each = length(targets.use)), targets.use, sep = " -> ")
    source.target.isolate <- setdiff(source.target, unique(df.net$source.target))
    if (length(source.target.isolate) > 0) {
      df.net.isolate <- as.data.frame(matrix(NA, nrow = length(source.target.isolate), ncol = ncol(df.net)))
      colnames(df.net.isolate) <- colnames(df.net)
      df.net.isolate$source.target <- source.target.isolate
      df.net.isolate$interaction_name_2 <- df.net$interaction_name_2[1]
      df.net.isolate$pval <- 1
      a <- stringr::str_split(df.net.isolate$source.target, " -> ", simplify = T)
      df.net.isolate$source <- as.character(a[, 1])
      df.net.isolate$target <- as.character(a[, 2])
      df.net <- rbind(df.net, df.net.isolate)
    }

    df.net$pval[df.net$pval > 0.05] = 1
    df.net$pval[df.net$pval > 0.01 & df.net$pval <= 0.05] = 2
    df.net$pval[df.net$pval <= 0.01] = 3
    df.net$prob[df.net$prob == 0] <- NA
    df.net$prob.original <- df.net$prob
    df.net$prob <- -1/log(df.net$prob)

    idx1 <- which(is.infinite(df.net$prob) | df.net$prob < 0)
    if (sum(idx1) > 0) {
      values.assign <- seq(max(df.net$prob, na.rm = T)*1.1, max(df.net$prob, na.rm = T)*1.5, length.out = length(idx1))
      position <- sort(prob.original[idx1], index.return = TRUE)$ix
      df.net$prob[idx1] <- values.assign[match(1:length(idx1), position)]
    }
    # rownames(df.net) <- df.net$interaction_name_2

    df.net$source <- factor(df.net$source, levels = cells.level[cells.level %in% unique(df.net$source)])
    df.net$target <- factor(df.net$target, levels = cells.level[cells.level %in% unique(df.net$target)])
    group.names <- paste(rep(levels(df.net$source), each = length(levels(df.net$target))), levels(df.net$target), sep = " -> ")

    df.net$interaction_name_2 <- as.character(df.net$interaction_name_2)
    df.net <- with(df.net, df.net[order(interaction_name_2),])
    df.net$interaction_name_2 <- factor(df.net$interaction_name_2, levels = unique(df.net$interaction_name_2))
    cells.order <- group.names
    df.net$source.target <- factor(df.net$source.target, levels = cells.order)
    df <- df.net
  } else {
    dataset.name <- names(object@net)
    df.net.all <- subsetCommunication(object, slot.name = "net",
                                      sources.use = sources.use, targets.use = targets.use,
                                      signaling = signaling,
                                      pairLR.use = pairLR.use,
                                      thresh = thresh)
    df.all <- data.frame()
    for (ii in 1:length(comparison)) {
      cells.level <- levels(object@idents[[comparison[ii]]])
      if (is.numeric(sources.use)) {
        sources.use <- cells.level[sources.use]
      }
      if (is.numeric(targets.use)) {
        targets.use <- cells.level[targets.use]
      }

      df.net <- df.net.all[[comparison[ii]]]
      df.net$interaction_name_2 <- as.character(df.net$interaction_name_2)
      df.net$source.target <- paste(df.net$source, df.net$target, sep = " -> ")
      source.target <- paste(rep(sources.use, each = length(targets.use)), targets.use, sep = " -> ")
      source.target.isolate <- setdiff(source.target, unique(df.net$source.target))
      if (length(source.target.isolate) > 0) {
        df.net.isolate <- as.data.frame(matrix(NA, nrow = length(source.target.isolate), ncol = ncol(df.net)))
        colnames(df.net.isolate) <- colnames(df.net)
        df.net.isolate$source.target <- source.target.isolate
        df.net.isolate$interaction_name_2 <- df.net$interaction_name_2[1]
        df.net.isolate$pval <- 1
        a <- stringr::str_split(df.net.isolate$source.target, " -> ", simplify = T)
        df.net.isolate$source <- as.character(a[, 1])
        df.net.isolate$target <- as.character(a[, 2])
        df.net <- rbind(df.net, df.net.isolate)
      }

      df.net$source <- factor(df.net$source, levels = cells.level[cells.level %in% unique(df.net$source)])
      df.net$target <- factor(df.net$target, levels = cells.level[cells.level %in% unique(df.net$target)])
      group.names <- paste(rep(levels(df.net$source), each = length(levels(df.net$target))), levels(df.net$target), sep = " -> ")
      group.names0 <- group.names
      group.names <- paste0(group.names0, " (", dataset.name[comparison[ii]], ")")

      if (nrow(df.net) > 0) {
        df.net$pval[df.net$pval > 0.05] = 1
        df.net$pval[df.net$pval > 0.01 & df.net$pval <= 0.05] = 2
        df.net$pval[df.net$pval <= 0.01] = 3
        df.net$prob[df.net$prob == 0] <- NA
        df.net$prob.original <- df.net$prob
        df.net$prob <- -1/log(df.net$prob)
      } else {
        df.net <- as.data.frame(matrix(NA, nrow = length(group.names), ncol = 5))
        colnames(df.net) <- c("interaction_name_2","source.target","prob","pval","prob.original")
        df.net$source.target <- group.names0
      }
      # df.net$group.names <- sub(paste0(' \\(',dataset.name[comparison[ii]],'\\)'),'',as.character(df.net$source.target))
      df.net$group.names <- as.character(df.net$source.target)
      df.net$source.target <- paste0(df.net$source.target, " (", dataset.name[comparison[ii]], ")")
      df.net$dataset <- dataset.name[comparison[ii]]
      df.all <- rbind(df.all, df.net)
    }
    if (nrow(df.all) == 0) {
      stop("No interactions are detected. Please consider changing the cell groups for analysis. ")
    }

    idx1 <- which(is.infinite(df.all$prob) | df.all$prob < 0)
    if (sum(idx1) > 0) {
      values.assign <- seq(max(df.all$prob, na.rm = T)*1.1, max(df.all$prob, na.rm = T)*1.5, length.out = length(idx1))
      position <- sort(df.all$prob.original[idx1], index.return = TRUE)$ix
      df.all$prob[idx1] <- values.assign[match(1:length(idx1), position)]
    }

    df.all$interaction_name_2[is.na(df.all$interaction_name_2)] <- df.all$interaction_name_2[!is.na(df.all$interaction_name_2)][1]

    df <- df.all
    df <- with(df, df[order(interaction_name_2),])
    df$interaction_name_2 <- factor(df$interaction_name_2, levels = unique(df$interaction_name_2))

    cells.order <- c()
    dataset.name.order <- c()
    for (i in 1:length(group.names0)) {
      for (j in 1:length(comparison)) {
        cells.order <- c(cells.order, paste0(group.names0[i], " (", dataset.name[comparison[j]], ")"))
        dataset.name.order <- c(dataset.name.order, dataset.name[comparison[j]])
      }
    }
    df$source.target <- factor(df$source.target, levels = cells.order)
  }

  min.cutoff <- quantile(df$prob, min.quantile,na.rm= T)
  max.cutoff <- quantile(df$prob, max.quantile,na.rm= T)
  df$prob[df$prob < min.cutoff] <- min.cutoff
  df$prob[df$prob > max.cutoff] <- max.cutoff


  if (remove.isolate) {
    df <- df[!is.na(df$prob), ]
    line.on <- FALSE
  }
  if (!is.null(max.dataset)) {
    # line.on <- FALSE
    # df <- df[!is.na(df$prob),]
    signaling <- as.character(unique(df$interaction_name_2))
    for (i in signaling) {
      df.i <- df[df$interaction_name_2 == i, ,drop = FALSE]
      cell <- as.character(unique(df.i$group.names))
      for (j in cell) {
        df.i.j <- df.i[df.i$group.names == j, , drop = FALSE]
        values <- df.i.j$prob
        idx.max <- which(values == max(values, na.rm = T))
        idx.min <- which(values == min(values, na.rm = T))
        #idx.na <- c(which(is.na(values)), which(!(dataset.name[comparison] %in% df.i.j$dataset)))
        dataset.na <- c(df.i.j$dataset[is.na(values)], setdiff(dataset.name[comparison], df.i.j$dataset))
        if (length(idx.max) > 0) {
          if (!(df.i.j$dataset[idx.max] %in% dataset.name[max.dataset])) {
            df.i.j$prob <- NA
          } else if ((idx.max != idx.min) & !is.null(min.dataset)) {
            if (!(df.i.j$dataset[idx.min] %in% dataset.name[min.dataset])) {
              df.i.j$prob <- NA
            } else if (length(dataset.na) > 0 & sum(!(dataset.name[min.dataset] %in% dataset.na)) > 0) {
              df.i.j$prob <- NA
            }
          }
        }
        df.i[df.i$group.names == j, "prob"] <- df.i.j$prob
      }
      df[df$interaction_name_2 == i, "prob"] <- df.i$prob
    }
    #df <- df[!is.na(df$prob), ]
  }
  if (remove.isolate) {
    df <- df[!is.na(df$prob), ]
    line.on <- FALSE
  }
  if (nrow(df) == 0) {
    stop("No interactions are detected. Please consider changing the cell groups for analysis. ")
  }
  # Re-order y-axis
  interaction_name_2.order <- intersect(object@DB$interaction[pairLR.use$interaction_name, ]$interaction_name_2, unique(df$interaction_name_2))
  df$interaction_name_2 <- factor(df$interaction_name_2, levels = interaction_name_2.order)

  # Re-order x-axis
  df$source.target = droplevels(df$source.target, exclude = setdiff(levels(df$source.target),unique(df$source.target)))
  if (sort.by.target & !sort.by.source) {
    if (!is.null(targets.use)) {
      df$target <- factor(df$target, levels = intersect(targets.use, df$target))
      df <- with(df, df[order(target, source),])
      source.target.order <- unique(as.character(df$source.target))
      df$source.target <- factor(df$source.target, levels = source.target.order)
    }
  }
  if (sort.by.source & !sort.by.target) {
    if (!is.null(sources.use)) {
      df$source <- factor(df$source, levels = intersect(sources.use, df$source))
      df <- with(df, df[order(source, target),])
      source.target.order <- unique(as.character(df$source.target))
      df$source.target <- factor(df$source.target, levels = source.target.order)
    }
  }
  if (sort.by.source & sort.by.target) {
    if (!is.null(sources.use)) {
      df$source <- factor(df$source, levels = intersect(sources.use, df$source))
      if (!is.null(targets.use)) {
        df$target <- factor(df$target, levels = intersect(targets.use, df$target))
      }
      if (sort.by.source.priority) {
        df <- with(df, df[order(source, target),])
      } else {
        df <- with(df, df[order(target, source),])
      }

      source.target.order <- unique(as.character(df$source.target))
      df$source.target <- factor(df$source.target, levels = source.target.order)
    }
  }

  g <- ggplot(df, aes(x = source.target, y = interaction_name_2, color = prob, size = pval)) +
    geom_point(pch = 16) +
    theme_linedraw() + theme(panel.grid.major = element_blank()) +
    theme(axis.text.x = element_text(angle = angle.x, hjust= hjust.x, vjust = vjust.x),
          axis.title.x = element_blank(),
          axis.title.y = element_blank()) +
    scale_x_discrete(position = "bottom")

  values <- c(1,2,3); names(values) <- c("p > 0.05", "0.01 < p < 0.05","p < 0.01")
  g <- g + scale_radius(range = c(min(df$pval), max(df$pval)), breaks = sort(unique(df$pval)),labels = names(values)[values %in% sort(unique(df$pval))], name = "p-value")
  #g <- g + scale_radius(range = c(1,3), breaks = values,labels = names(values), name = "p-value")
  if (min(df$prob, na.rm = T) != max(df$prob, na.rm = T)) {
    g <- g + scale_colour_gradientn(colors = colorRampPalette(color.use)(99), na.value = "white", limits=c(quantile(df$prob, 0,na.rm= T), quantile(df$prob, 1,na.rm= T)),
                                    breaks = c(quantile(df$prob, 0,na.rm= T), quantile(df$prob, 1,na.rm= T)), labels = c("min","max")) +
      guides(color = guide_colourbar(barwidth = 0.5, title = "Commun. Prob."))
  } else {
    g <- g + scale_colour_gradientn(colors = colorRampPalette(color.use)(99), na.value = "white") +
      guides(color = guide_colourbar(barwidth = 0.5, title = "Commun. Prob."))
  }

  g <- g + theme(text = element_text(size = font.size),plot.title = element_text(size=font.size.title)) +
    theme(legend.title = element_text(size = 8), legend.text = element_text(size = 6))

  if (grid.on) {
    if (length(unique(df$source.target)) > 1) {
      g <- g + geom_vline(xintercept=seq(1.5, length(unique(df$source.target))-0.5, 1),lwd=0.1,colour=color.grid)
    }
    if (length(unique(df$interaction_name_2)) > 1) {
      g <- g + geom_hline(yintercept=seq(1.5, length(unique(df$interaction_name_2))-0.5, 1),lwd=0.1,colour=color.grid)
    }
  }
  if (!is.null(title.name)) {
    g <- g + ggtitle(title.name) + theme(plot.title = element_text(hjust = 0.5))
  }

  if (!is.null(comparison)) {
    if (line.on) {
      xintercept = seq(0.5+length(dataset.name[comparison]), length(group.names0)*length(dataset.name[comparison]), by = length(dataset.name[comparison]))
      g <- g + geom_vline(xintercept = xintercept, linetype="dashed", color = "grey60", size = line.size)
    }
    if (color.text.use) {
      if (is.null(group)) {
        group <- 1:length(comparison)
        names(group) <- dataset.name[comparison]
      }
      if (is.null(color.text)) {
        color <- ggPalette(length(unique(group)))
      } else {
        color <- color.text
      }
      names(color) <- names(group[!duplicated(group)])
      color <- color[group]
      #names(color) <- dataset.name[comparison]
      dataset.name.order <- levels(df$source.target)
      dataset.name.order <- stringr::str_match(dataset.name.order, "\\(.*\\)")
      dataset.name.order <- stringr::str_sub(dataset.name.order, 2, stringr::str_length(dataset.name.order)-1)
      xtick.color <- color[dataset.name.order]
      g <- g + theme(axis.text.x = element_text(colour = xtick.color))
    }
  }
  if (!show.legend) {
    g <- g + theme(legend.position = "none")
  }
  if (return.data) {
    return(list(communication = df, gg.obj = g))
  } else {
    return(g)
  }

}




#' Chord diagram for visualizing cell-cell communication for a signaling pathway
#'
#' Names of cell states will be displayed in this chord diagram
#'
#' @param object CellChat object
#' @param signaling a character vector giving the name of signaling networks
#' @param net a weighted matrix or a data frame with three columns defining the cell-cell communication network
#' @param slot.name the slot name of object: slot.name = "net" when visualizing cell-cell communication network per each ligand-receptor pair associated with a given signaling pathway;
#' slot.name = "netP" when visualizing cell-cell communication network at the level of signaling pathways
#' @param color.use colors for the cell groups
#' @param group A named group labels for making multiple-group Chord diagrams. The sector names should be used as the names in the vector.
#' The order of group controls the sector orders and if group is set as a factor, the order of levels controls the order of groups.
#' @param cell.order a char vector defining the cell type orders (sector orders)
#' @param sources.use a vector giving the index or the name of source cell groups
#' @param targets.use a vector giving the index or the name of target cell groups.
#' @param lab.cex font size for the text
#' @param small.gap Small gap between sectors.
#' @param big.gap Gap between the different sets of sectors, which are defined in the `group` parameter
#' @param annotationTrackHeight annotationTrack Height
#' @param remove.isolate whether remove sectors without any links
#' @param link.visible whether plot the link. The value is logical, if it is set to FALSE, the corresponding link will not plotted, but the space is still ocuppied. The format is a matrix with names or a data frame with three columns
#' @param scale scale each sector to same width; default = FALSE; however, it is set to be TRUE when remove.isolate = TRUE
#' @param link.target.prop If the Chord diagram is directional, for each source sector, whether to draw bars that shows the proportion of target sectors.
#' @param reduce if the ratio of the width of certain grid compared to the whole circle is less than this value, the grid is removed on the plot. Set it to value less than zero if you want to keep all tiny grid.
#' @param directional Whether links have directions. 1 means the direction is from the first column in df to the second column, -1 is the reverse, 0 is no direction, and 2 for two directional.
#' @param transparency Transparency of link colors
#' @param link.border border for links, single scalar or a matrix with names or a data frame with three columns
#' @param title.name title name
#' @param show.legend whether show the figure legend
#' @param legend.pos.x,legend.pos.y adjust the legend position
#' @param nCol number of columns when displaying the figures
#' @param thresh threshold of the p-value for determining significant interaction when visualizing links at the level of ligands/receptors;
#' @param ... other parameters passing to chordDiagram
#' @return an object of class "recordedplot"
#' @export

netVisual_chord_cell <- function(object, signaling = NULL, net = NULL, slot.name = "netP",
                                 color.use = NULL,group = NULL,cell.order = NULL,
                                 sources.use = NULL, targets.use = NULL,
                                 lab.cex = 0.8,small.gap = 1, big.gap = 10, annotationTrackHeight = c(0.03),
                                 remove.isolate = FALSE, link.visible = TRUE, scale = FALSE, directional = 1,link.target.prop = TRUE, reduce = -1,
                                 transparency = 0.4, link.border = NA,
                                 title.name = NULL, show.legend = FALSE, legend.pos.x = 20, legend.pos.y = 20, nCol = NULL,
                                 thresh = 0.05,...){

  if (!is.null(signaling)) {
    pairLR <- searchPair(signaling = signaling, pairLR.use = object@LR$LRsig, key = "pathway_name", matching.exact = T, pair.only = F)
    net <- object@net

    pairLR.use.name <- dimnames(net$prob)[[3]]
    pairLR.name <- intersect(rownames(pairLR), pairLR.use.name)
    pairLR <- pairLR[pairLR.name, ]
    prob <- net$prob
    pval <- net$pval

    prob[pval > thresh] <- 0
    if (length(pairLR.name) > 1) {
      pairLR.name.use <- pairLR.name[apply(prob[,,pairLR.name], 3, sum) != 0]
    } else {
      pairLR.name.use <- pairLR.name[sum(prob[,,pairLR.name]) != 0]
    }

    if (length(pairLR.name.use) == 0) {
      stop(paste0('There is no significant communication of ', signaling))
    } else {
      pairLR <- pairLR[pairLR.name.use,]
    }
    nRow <- length(pairLR.name.use)

    prob <- prob[,,pairLR.name.use]

    if (length(dim(prob)) == 2) {
      prob <- replicate(1, prob, simplify="array")
    }

    if (slot.name == "netP") {
      message("Plot the aggregated cell-cell communication network at the signaling pathway level")
      net <- apply(prob, c(1,2), sum)
      if (is.null(title.name)) {
        title.name <- paste0(signaling, " signaling pathway network")
      }
      # par(mfrow = c(1,1), xpd=TRUE)
      # par(mar = c(5, 4, 4, 2))
      gg <- netVisual_chord_cell_internal(net, color.use = color.use, group = group, cell.order = cell.order, sources.use = sources.use, targets.use = targets.use,
                                          lab.cex = lab.cex,small.gap = small.gap, big.gap = big.gap,annotationTrackHeight = annotationTrackHeight,
                                          remove.isolate = remove.isolate, link.visible = link.visible, scale = scale, directional = directional,link.target.prop = link.target.prop, reduce = reduce,
                                          transparency = transparency, link.border = link.border,
                                          title.name = title.name, show.legend = show.legend, legend.pos.x = legend.pos.x, legend.pos.y = legend.pos.y, ...)
    } else if (slot.name == "net") {
      message("Plot the cell-cell communication network per each ligand-receptor pair associated with a given signaling pathway")
      if (is.null(nCol)) {
        nCol <- min(length(pairLR.name.use), 2)
      }
      #   layout(matrix(1:length(pairLR.name.use), ncol = nCol))
      # par(xpd=TRUE)
      # par(mfrow = c(ceiling(length(pairLR.name.use)/nCol), nCol), xpd=TRUE, mar = c(5, 4, 4, 2) +0.1)
      par(mfrow = c(ceiling(length(pairLR.name.use)/nCol), nCol), xpd=TRUE)
      gg <- vector("list", length(pairLR.name.use))
      for (i in 1:length(pairLR.name.use)) {
        #par(mar = c(5, 4, 4, 2))
        title.name <- pairLR$interaction_name_2[i]
        net <- prob[,,i]
        gg[[i]] <- netVisual_chord_cell_internal(net, color.use = color.use, group = group,cell.order = cell.order,sources.use = sources.use, targets.use = targets.use,
                                                 lab.cex = lab.cex,small.gap = small.gap,big.gap = big.gap, annotationTrackHeight = annotationTrackHeight,
                                                 remove.isolate = remove.isolate, link.visible = link.visible, scale = scale, directional = directional,link.target.prop = link.target.prop, reduce = reduce,
                                                 transparency = transparency, link.border = link.border,
                                                 title.name = title.name, show.legend = show.legend, legend.pos.x = legend.pos.x, legend.pos.y = legend.pos.y, ...)
      }
    }

  } else if (!is.null(net)) {
    gg <- netVisual_chord_cell_internal(net, color.use = color.use, group = group,cell.order = cell.order,sources.use = sources.use, targets.use = targets.use,
                                        lab.cex = lab.cex,small.gap = small.gap, big.gap = big.gap,annotationTrackHeight = annotationTrackHeight,
                                        remove.isolate = remove.isolate, link.visible = link.visible, scale = scale, directional = directional,link.target.prop = link.target.prop, reduce = reduce,
                                        transparency = transparency, link.border = link.border,
                                        title.name = title.name, show.legend = show.legend, legend.pos.x = legend.pos.x,legend.pos.y=legend.pos.y, ...)
  } else {
    stop("Please assign values to either `signaling` or `net`")
  }

  return(gg)
}


#' Chord diagram for visualizing cell-cell communication from a weighted adjacency matrix or a data frame
#'
#' Names of cell states/groups will be displayed in this chord diagram
#'
#' @param net a weighted matrix or a data frame with three columns defining the cell-cell communication network
#' @param color.use colors for the cell groups
#' @param group A named group labels for making multiple-group Chord diagrams. The sector names should be used as the names in the vector.
#' The order of group controls the sector orders and if group is set as a factor, the order of levels controls the order of groups.
#' @param cell.order a char vector defining the cell type orders (sector orders)
#' @param sources.use a vector giving the index or the name of source cell groups
#' @param targets.use a vector giving the index or the name of target cell groups.
#' @param lab.cex font size for the text
#' @param small.gap Small gap between sectors.
#' @param big.gap Gap between the different sets of sectors, which are defined in the `group` parameter
#' @param annotationTrackHeight annotationTrack Height
#' @param remove.isolate whether remove sectors without any links
#' @param link.visible whether plot the link. The value is logical, if it is set to FALSE, the corresponding link will not plotted, but the space is still ocuppied. The format is a matrix with names or a data frame with three columns
#' @param scale scale each sector to same width; default = FALSE; however, it is set to be TRUE when remove.isolate = TRUE
#' @param link.target.prop If the Chord diagram is directional, for each source sector, whether to draw bars that shows the proportion of target sectors.
#' @param reduce if the ratio of the width of certain grid compared to the whole circle is less than this value, the grid is removed on the plot. Set it to value less than zero if you want to keep all tiny grid.
#' @param directional Whether links have directions. 1 means the direction is from the first column in df to the second column, -1 is the reverse, 0 is no direction, and 2 for two directional.
#' @param transparency Transparency of link colors
#' @param link.border border for links, single scalar or a matrix with names or a data frame with three columns
#' @param title.name title name of the plot
#' @param show.legend whether show the figure legend
#' @param legend.pos.x,legend.pos.y adjust the legend position
#' @param ... other parameters passing to chordDiagram
#' @importFrom circlize circos.clear chordDiagram circos.track circos.text get.cell.meta.data
#' @importFrom grDevices recordPlot
#' @return an object of class "recordedplot"
#' @export

netVisual_chord_cell_internal <- function(net, color.use = NULL, group = NULL, cell.order = NULL,
                                          sources.use = NULL, targets.use = NULL,
                                          lab.cex = 0.8,small.gap = 1, big.gap = 10, annotationTrackHeight = c(0.03),
                                          remove.isolate = FALSE, link.visible = TRUE, scale = FALSE, directional = 1, link.target.prop = TRUE, reduce = -1,
                                          transparency = 0.4, link.border = NA,
                                          title.name = NULL, show.legend = FALSE, legend.pos.x = 20, legend.pos.y = 20,...){
  if (inherits(x = net, what = c("matrix", "Matrix"))) {
    cell.levels <- union(rownames(net), colnames(net))
    net <- reshape2::melt(net, value.name = "prob")
    colnames(net)[1:2] <- c("source","target")
  } else if (is.data.frame(net)) {
    if (all(c("source","target", "prob") %in% colnames(net)) == FALSE) {
      stop("The input data frame must contain three columns named as source, target, prob")
    }
    cell.levels <- as.character(union(net$source,net$target))
  }
  if (!is.null(cell.order)) {
    cell.levels <- cell.order
  }
  net$source <- as.character(net$source)
  net$target <- as.character(net$target)

  # keep the interactions associated with sources and targets of interest
  if (!is.null(sources.use)){
    if (is.numeric(sources.use)) {
      sources.use <- cell.levels[sources.use]
    }
    net <- subset(net, source %in% sources.use)
  }
  if (!is.null(targets.use)){
    if (is.numeric(targets.use)) {
      targets.use <- cell.levels[targets.use]
    }
    net <- subset(net, target %in% targets.use)
  }
  # remove the interactions with zero values
  net <- subset(net, prob > 0)
  if(dim(net)[1]<=0){message("No interaction between those cells")}
  # create a fake data if keeping the cell types (i.e., sectors) without any interactions
  if (!remove.isolate) {
    cells.removed <- setdiff(cell.levels, as.character(union(net$source,net$target)))
    if (length(cells.removed) > 0) {
      net.fake <- data.frame(cells.removed, cells.removed, 1e-10*sample(length(cells.removed), length(cells.removed)))
      colnames(net.fake) <- colnames(net)
      net <- rbind(net, net.fake)
      link.visible <- net[, 1:2]
      link.visible$plot <- FALSE
      if(nrow(net) > nrow(net.fake)){
        link.visible$plot[1:(nrow(net) - nrow(net.fake))] <- TRUE
      }
      # directional <- net[, 1:2]
      # directional$plot <- 0
      # directional$plot[1:(nrow(net) - nrow(net.fake))] <- 1
      # link.arr.type = "big.arrow"
      # message("Set scale = TRUE when remove.isolate = FALSE")
      scale = TRUE
    }
  }

  df <- net
  cells.use <- union(df$source,df$target)

  # define grid order
  order.sector <- cell.levels[cell.levels %in% cells.use]

  # define grid color
  if (is.null(color.use)){
    color.use = scPalette(length(cell.levels))
    names(color.use) <- cell.levels
  } else if (is.null(names(color.use))) {
    names(color.use) <- cell.levels
  }
  grid.col <- color.use[order.sector]
  names(grid.col) <- order.sector

  # set grouping information
  if (!is.null(group)) {
    group <- group[names(group) %in% order.sector]
  }

  # define edge color
  edge.color <- color.use[as.character(df$source)]

  if (directional == 0 | directional == 2) {
    link.arr.type = "triangle"
  } else {
    link.arr.type = "big.arrow"
  }

  circos.clear()
  chordDiagram(df,
               order = order.sector,
               col = edge.color,
               grid.col = grid.col,
               transparency = transparency,
               link.border = link.border,
               directional = directional,
               direction.type = c("diffHeight","arrows"),
               link.arr.type = link.arr.type, # link.border = "white",
               annotationTrack = "grid",
               annotationTrackHeight = annotationTrackHeight,
               preAllocateTracks = list(track.height = max(strwidth(order.sector))),
               small.gap = small.gap,
               big.gap = big.gap,
               link.visible = link.visible,
               scale = scale,
               group = group,
               link.target.prop = link.target.prop,
               reduce = reduce,
               ...)
  circos.track(track.index = 1, panel.fun = function(x, y) {
    xlim = get.cell.meta.data("xlim")
    xplot = get.cell.meta.data("xplot")
    ylim = get.cell.meta.data("ylim")
    sector.name = get.cell.meta.data("sector.index")
    circos.text(mean(xlim), ylim[1], sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5),cex = lab.cex)
  }, bg.border = NA)

  # https://jokergoo.github.io/circlize_book/book/legends.html
  if (show.legend) {
    lgd <- ComplexHeatmap::Legend(at = names(grid.col), type = "grid", legend_gp = grid::gpar(fill = grid.col), title = "Cell State")
    ComplexHeatmap::draw(lgd, x = unit(1, "npc")-unit(legend.pos.x, "mm"), y = unit(legend.pos.y, "mm"), just = c("right", "bottom"))
  }

  if(!is.null(title.name)){
    # title(title.name, cex = 1)
    text(-0, 1.02, title.name, cex=1)
  }
  circos.clear()
  gg <- recordPlot()
  return(gg)
}


#' Chord diagram for visualizing cell-cell communication for a set of ligands/receptors or signaling pathways
#'
#' Names of ligands/receptors or signaling pathways will be displayed in this chord diagram
#'
#' @param object CellChat object
#' @param slot.name the slot name of object: slot.name = "net" when visualizing links at the level of ligands/receptors; slot.name = "netP" when visualizing links at the level of signaling pathways
#' @param signaling a character vector giving the name of signaling networks
#' @param pairLR.use a data frame consisting of one column named either "interaction_name" or "pathway_name", defining the interactions of interest
#' @param net A data frame consisting of the interactions of interest.
#' net should have at least three columns: "source","target" and "interaction_name" when visualizing links at the level of ligands/receptors;
#' "source","target" and "pathway_name" when visualizing links at the level of signaling pathway; "interaction_name" and "pathway_name" must be the matched names in CellChatDB$interaction.
#' @param sources.use a vector giving the index or the name of source cell groups
#' @param targets.use a vector giving the index or the name of target cell groups.
#' @param color.use colors for the cell groups
#' @param lab.cex font size for the text
#' @param small.gap Small gap between sectors.
#' @param big.gap Gap between the different sets of sectors, which are defined in the `group` parameter
#' @param annotationTrackHeight annotationTrack Height
#' @param link.visible whether plot the link. The value is logical, if it is set to FALSE, the corresponding link will not plotted, but the space is still ocuppied. The format is a matrix with names or a data frame with three columns
#' @param scale scale each sector to same width; default = FALSE; however, it is set to be TRUE when remove.isolate = TRUE
#' @param link.target.prop If the Chord diagram is directional, for each source sector, whether to draw bars that shows the proportion of target sectors.
#' @param reduce if the ratio of the width of certain grid compared to the whole circle is less than this value, the grid is removed on the plot. Set it to value less than zero if you want to keep all tiny grid.
#' @param directional Whether links have directions. 1 means the direction is from the first column in df to the second column, -1 is the reverse, 0 is no direction, and 2 for two directional.
#' @param transparency Transparency of link colors
#' @param link.border border for links, single scalar or a matrix with names or a data frame with three columns
#' @param title.name title name of the plot
#' @param show.legend whether show the figure legend
#' @param legend.pos.x,legend.pos.y adjust the legend position
#' @param thresh threshold of the p-value for determining significant interaction when visualizing links at the level of ligands/receptors;
#' @param ... other parameters to chordDiagram
#' @importFrom circlize circos.clear chordDiagram circos.track circos.text get.cell.meta.data
#' @importFrom dplyr select %>% group_by summarize
#' @importFrom grDevices recordPlot
#' @importFrom stringr str_split
#' @return an object of class "recordedplot"
#' @export

netVisual_chord_gene <- function(object, slot.name = "net", color.use = NULL,
                                 signaling = NULL, pairLR.use = NULL, net = NULL,
                                 sources.use = NULL, targets.use = NULL,
                                 lab.cex = 0.8,small.gap = 1, big.gap = 10, annotationTrackHeight = c(0.03),
                                 link.visible = TRUE, scale = FALSE, directional = 1, link.target.prop = TRUE, reduce = -1,
                                 transparency = 0.4, link.border = NA,
                                 title.name = NULL, legend.pos.x = 20, legend.pos.y = 20, show.legend = TRUE,
                                 thresh = 0.05,
                                 ...){
  if (!is.null(pairLR.use)) {
    if (!is.data.frame(pairLR.use)) {
      stop("pairLR.use should be a data frame with a signle column named either 'interaction_name' or 'pathway_name' ")
    } else if ("pathway_name" %in% colnames(pairLR.use)) {
      message("slot.name is set to be 'netP' when pairLR.use contains signaling pathways")
      slot.name = "netP"
    }
  }

  if (!is.null(pairLR.use) & !is.null(signaling)) {
    stop("Please do not assign values to 'signaling' when using 'pairLR.use'")
  }

  if (is.null(net)) {
    prob <- slot(object, "net")$prob
    pval <- slot(object, "net")$pval
    prob[pval > thresh] <- 0
    net <- reshape2::melt(prob, value.name = "prob")
    colnames(net)[1:3] <- c("source","target","interaction_name")

    pairLR = dplyr::select(object@LR$LRsig, c("interaction_name_2", "pathway_name",  "ligand",  "receptor" ,"annotation","evidence"))
    idx <- match(net$interaction_name, rownames(pairLR))
    temp <- pairLR[idx,]
    net <- cbind(net, temp)
  }

  if (!is.null(signaling)) {
    pairLR.use <- data.frame()
    for (i in 1:length(signaling)) {
      pairLR.use.i <- searchPair(signaling = signaling[i], pairLR.use = object@LR$LRsig, key = "pathway_name", matching.exact = T, pair.only = T)
      pairLR.use <- rbind(pairLR.use, pairLR.use.i)
    }
  }

  if (!is.null(pairLR.use)){
    if ("interaction_name" %in% colnames(pairLR.use)) {
      net <- subset(net,interaction_name %in% pairLR.use$interaction_name)
    } else if ("pathway_name" %in% colnames(pairLR.use)) {
      net <- subset(net, pathway_name %in% as.character(pairLR.use$pathway_name))
    }
  }

  if (slot.name == "netP") {
    net <- dplyr::select(net, c("source","target","pathway_name","prob"))
    net$source_target <- paste(net$source, net$target, sep = "sourceTotarget")
    net <- net %>% dplyr::group_by(source_target, pathway_name) %>% dplyr::summarize(prob = sum(prob))
    a <- stringr::str_split(net$source_target, "sourceTotarget", simplify = T)
    net$source <- as.character(a[, 1])
    net$target <- as.character(a[, 2])
    net$ligand <- net$pathway_name
    net$receptor <- " "
  }

  # keep the interactions associated with sources and targets of interest
  if (!is.null(sources.use)){
    if (is.numeric(sources.use)) {
      sources.use <- levels(object@idents)[sources.use]
    }
    net <- subset(net, source %in% sources.use)
  } else {
    sources.use <- levels(object@idents)
  }
  if (!is.null(targets.use)){
    if (is.numeric(targets.use)) {
      targets.use <- levels(object@idents)[targets.use]
    }
    net <- subset(net, target %in% targets.use)
  } else {
    targets.use <- levels(object@idents)
  }
  # remove the interactions with zero values
  df <- subset(net, prob > 0)

  if (nrow(df) == 0) {
    stop("No signaling links are inferred! ")
  }

  if (length(unique(net$ligand)) == 1) {
    message("You may try the function `netVisual_chord_cell` for visualizing individual signaling pathway")
  }

  df$id <- 1:nrow(df)
  # deal with duplicated sector names
  ligand.uni <- unique(df$ligand)
  for (i in 1:length(ligand.uni)) {
    df.i <- df[df$ligand == ligand.uni[i], ]
    source.uni <- unique(df.i$source)
    for (j in 1:length(source.uni)) {
      df.i.j <- df.i[df.i$source == source.uni[j], ]
      df.i.j$ligand <- paste0(df.i.j$ligand, paste(rep(' ',j-1),collapse = ''))
      df$ligand[df$id %in% df.i.j$id] <- df.i.j$ligand
    }
  }
  receptor.uni <- unique(df$receptor)
  for (i in 1:length(receptor.uni)) {
    df.i <- df[df$receptor == receptor.uni[i], ]
    target.uni <- unique(df.i$target)
    for (j in 1:length(target.uni)) {
      df.i.j <- df.i[df.i$target == target.uni[j], ]
      df.i.j$receptor <- paste0(df.i.j$receptor, paste(rep(' ',j-1),collapse = ''))
      df$receptor[df$id %in% df.i.j$id] <- df.i.j$receptor
    }
  }

  cell.order.sources <- levels(object@idents)[levels(object@idents) %in% sources.use]
  cell.order.targets <- levels(object@idents)[levels(object@idents) %in% targets.use]

  df$source <- factor(df$source, levels = cell.order.sources)
  df$target <- factor(df$target, levels = cell.order.targets)
  # df.ordered.source <- df[with(df, order(source, target, -prob)), ]
  # df.ordered.target <- df[with(df, order(target, source, -prob)), ]
  df.ordered.source <- df[with(df, order(source, -prob)), ]
  df.ordered.target <- df[with(df, order(target, -prob)), ]

  order.source <- unique(df.ordered.source[ ,c('ligand','source')])
  order.target <- unique(df.ordered.target[ ,c('receptor','target')])

  # define sector order
  order.sector <- c(order.source$ligand, order.target$receptor)

  # define cell type color
  if (is.null(color.use)){
    color.use = scPalette(nlevels(object@idents))
    names(color.use) <- levels(object@idents)
    color.use <- color.use[levels(object@idents) %in% as.character(union(df$source,df$target))]
  } else if (is.null(names(color.use))) {
    names(color.use) <- levels(object@idents)
    color.use <- color.use[levels(object@idents) %in% as.character(union(df$source,df$target))]
  }

  # define edge color
  edge.color <- color.use[as.character(df.ordered.source$source)]
  names(edge.color) <- as.character(df.ordered.source$source)

  # define grid colors
  grid.col.ligand <- color.use[as.character(order.source$source)]
  names(grid.col.ligand) <- as.character(order.source$source)
  grid.col.receptor <- color.use[as.character(order.target$target)]
  names(grid.col.receptor) <- as.character(order.target$target)
  grid.col <- c(as.character(grid.col.ligand), as.character(grid.col.receptor))
  names(grid.col) <- order.sector

  df.plot <- df.ordered.source[ ,c('ligand','receptor','prob')]

  if (directional == 2) {
    link.arr.type = "triangle"
  } else {
    link.arr.type = "big.arrow"
  }
  circos.clear()
  chordDiagram(df.plot,
               order = order.sector,
               col = edge.color,
               grid.col = grid.col,
               transparency = transparency,
               link.border = link.border,
               directional = directional,
               direction.type = c("diffHeight","arrows"),
               link.arr.type = link.arr.type,
               annotationTrack = "grid",
               annotationTrackHeight = annotationTrackHeight,
               preAllocateTracks = list(track.height = max(strwidth(order.sector))),
               small.gap = small.gap,
               big.gap = big.gap,
               link.visible = link.visible,
               scale = scale,
               link.target.prop = link.target.prop,
               reduce = reduce,
               ...)

  circos.track(track.index = 1, panel.fun = function(x, y) {
    xlim = get.cell.meta.data("xlim")
    xplot = get.cell.meta.data("xplot")
    ylim = get.cell.meta.data("ylim")
    sector.name = get.cell.meta.data("sector.index")
    circos.text(mean(xlim), ylim[1], sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5),cex = lab.cex)
  }, bg.border = NA)

  # https://jokergoo.github.io/circlize_book/book/legends.html
  if (show.legend) {
    lgd <- ComplexHeatmap::Legend(at = names(color.use), type = "grid", legend_gp = grid::gpar(fill = color.use), title = "Cell State")
    ComplexHeatmap::draw(lgd, x = unit(1, "npc")-unit(legend.pos.x, "mm"), y = unit(legend.pos.y, "mm"), just = c("right", "bottom"))
  }

  circos.clear()
  if(!is.null(title.name)){
    text(-0, 1.02, title.name, cex=1)
  }
  gg <- recordPlot()
  return(gg)
}




#' River plot showing the associations of latent patterns with cell groups and ligand-receptor pairs or signaling pathways
#'
#' River (alluvial) plot shows the correspondence between the inferred latent patterns and cell groups as well as ligand-receptor pairs or signaling pathways.
#'
#' The thickness of the flow indicates the contribution of the cell group or signaling pathway to each latent pattern. The height of each pattern is proportional to the number of its associated cell groups or signaling pathways.
#'
#' Outgoing patterns reveal how the sender cells coordinate with each other as well as how they coordinate with certain signaling pathways to drive communication.
#'
#' Incoming patterns show how the target cells coordinate with each other as well as how they coordinate with certain signaling pathways to respond to incoming signaling.
#'
#' @param object CellChat object
#' @param slot.name the slot name of object that is used to compute centrality measures of signaling networks
#' @param pattern "outgoing" or "incoming"
#' @param cutoff the threshold for filtering out weak links
#' @param sources.use a vector giving the index or the name of source cell groups of interest
#' @param targets.use a vector giving the index or the name of target cell groups of interest
#' @param signaling a character vector giving the name of signaling pathways of interest
#' @param color.use the character vector defining the color of each cell group
#' @param color.use.pattern the character vector defining the color of each pattern
#' @param color.use.signaling the character vector defining the color of each signaling
#' @param do.order whether reorder the cell groups or signaling according to their similarity
#' @param main.title the title of plot
#' @param font.size font size of the text
#' @param font.size.title font size of the title
#' @importFrom methods slot
#' @importFrom stats cutree dist hclust
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#' @import ggalluvial
# #' @importFrom ggalluvial geom_stratum geom_flow to_lodes_form
#' @importFrom ggplot2 geom_text scale_x_discrete scale_fill_manual theme ggtitle
#' @importFrom cowplot plot_grid ggdraw draw_label
#' @return
#' @export
#'
#' @examples
netAnalysis_river <- function(object, slot.name = "netP", pattern = c("outgoing","incoming"), cutoff = 0.5,
                              sources.use = NULL, targets.use = NULL, signaling = NULL,
                              color.use = NULL, color.use.pattern = NULL, color.use.signaling = "grey50",
                              do.order = FALSE, main.title = NULL,
                              font.size = 2.5, font.size.title = 12){
  message("Please make sure you have load `library(ggalluvial)` when running this function")
  requireNamespace("ggalluvial")
  #  suppressMessages(require(ggalluvial))
  res.pattern <- methods::slot(object, slot.name)$pattern[[pattern]]
  data1 = res.pattern$pattern$cell
  data2 = res.pattern$pattern$signaling
  if (is.null(color.use.pattern)) {
    nPatterns <- length(unique(data1$Pattern))
    if (pattern == "outgoing") {
      color.use.pattern = ggPalette(nPatterns*2)[seq(1,nPatterns*2, by = 2)]
    } else if (pattern == "incoming") {
      color.use.pattern = ggPalette(nPatterns*2)[seq(2,nPatterns*2, by = 2)]
    }
  }
  if (is.null(main.title)) {
    if (pattern == "outgoing") {
      main.title = "Outgoing communication patterns of secreting cells"
    } else if (pattern == "incoming") {
      main.title = "Incoming communication patterns of target cells"
    }
  }

  if (is.null(data2)) {
    data1$Contribution[data1$Contribution < cutoff] <- 0
    plot.data <- data1
    nPatterns<-length(unique(plot.data$Pattern))
    nCellGroup<-length(unique(plot.data$CellGroup))
    if (is.null(color.use)) {
      color.use <- scPalette(nCellGroup)
    }
    if (is.null(color.use.pattern)){
      color.use.pattern <- ggPalette(nPatterns)
    }

    plot.data.long <- to_lodes_form(plot.data, axes = 1:2, id = "connection")
    if (do.order) {
      mat = tapply(plot.data[["Contribution"]], list(plot.data[["CellGroup"]], plot.data[["Pattern"]]), sum)
      d <- dist(as.matrix(mat))
      hc <- hclust(d, "ave")
      k <- length(unique(grep("Pattern", plot.data.long$stratum[plot.data.long$Contribution != 0], value = T)))
      cluster <- hc %>% cutree(k)
      order.name <- order(cluster)
      plot.data.long$stratum <- factor(plot.data.long$stratum, levels = c(names(cluster)[order.name], colnames(mat)))
      color.use <- color.use[order.name]
    }
    color.use.all <- c(color.use, color.use.pattern)
    gg <- ggplot(plot.data.long,aes(x = factor(x, levels = c("CellGroup", "Pattern")),y=Contribution,
                                    stratum = stratum, alluvium = connection,
                                    fill = stratum, label = stratum)) +
      geom_flow(width = 1/3,aes.flow = "backward") +
      geom_stratum(width=1/3,size=0.1,color="black", alpha = 0.8, linetype = 1) +
      geom_text(stat = "stratum", size = font.size) +
      scale_x_discrete(limits = c(),  labels=c("Cell groups", "Patterns")) +
      scale_fill_manual(values = alpha(color.use.all, alpha = 0.8), drop = FALSE) +
      theme_bw()+
      theme(legend.position = "none",
            axis.title = element_blank(),
            axis.text.y= element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor  = element_blank(),
            panel.border = element_blank(),
            axis.ticks = element_blank(),axis.text=element_text(size=10))+
      ggtitle(main.title)

  } else {
    data1$Contribution[data1$Contribution < cutoff] <- 0
    plot.data <- data1
    nPatterns<-length(unique(plot.data$Pattern))
    nCellGroup<-length(unique(plot.data$CellGroup))
    cells.level = levels(object@idents)
    if (is.null(color.use)) {
      color.use <- scPalette(length(cells.level))[cells.level %in% unique(plot.data$CellGroup)]
    }
    if (is.null(color.use.pattern)){
      color.use.pattern <- ggPalette(nPatterns)
    }
    if (!is.null(sources.use)) {
      if (is.numeric(sources.use)) {
        sources.use <- cells.level[sources.use]
      }
      plot.data <- subset(plot.data, CellGroup %in% sources.use)
    }
    if (!is.null(targets.use)) {
      if (is.numeric(targets.use)) {
        targets.use <- cells.level[targets.use]
      }
      plot.data <- subset(plot.data, CellGroup %in% targets.use)
    }
    ## connect cell groups with patterns
    plot.data.long <- to_lodes_form(plot.data, axes = 1:2, id = "connection")
    if (do.order) {
      mat = tapply(plot.data[["Contribution"]], list(plot.data[["CellGroup"]], plot.data[["Pattern"]]), sum)
      d <- dist(as.matrix(mat))
      hc <- hclust(d, "ave")
      k <- length(unique(grep("Pattern", plot.data.long$stratum[plot.data.long$Contribution != 0], value = T)))
      cluster <- hc %>% cutree(k)
      order.name <- order(cluster)
      plot.data.long$stratum <- factor(plot.data.long$stratum, levels = c(names(cluster)[order.name], colnames(mat)))
      color.use <- color.use[order.name]
    }
    color.use.all <- c(color.use, color.use.pattern)
    StatStratum <- ggalluvial::StatStratum
    gg1 <- ggplot(plot.data.long,aes(x = factor(x, levels = c("CellGroup", "Pattern")),y=Contribution,
                                     stratum = stratum, alluvium = connection,
                                     fill = stratum, label = stratum)) +
      geom_flow(width = 1/3,aes.flow = "backward") +
      geom_stratum(width=1/3,size=0.1,color="black", alpha = 0.8, linetype = 1) +
      geom_text(stat = "stratum", size = font.size) +
      scale_x_discrete(limits = c(),  labels=c("Cell groups", "Patterns")) +
      scale_fill_manual(values = alpha(color.use.all, alpha = 0.8), drop = FALSE) +
      theme_bw()+
      theme(legend.position = "none",
            axis.title = element_blank(),
            axis.text.y= element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor  = element_blank(),
            panel.border = element_blank(),
            axis.ticks = element_blank(),axis.text=element_text(size=10)) +
      theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

    ## connect patterns with signaling
    data2$Contribution[data2$Contribution < cutoff] <- 0
    plot.data <- data2
    nPatterns<-length(unique(plot.data$Pattern))
    nSignaling<-length(unique(plot.data$Signaling))
    if (length(color.use.signaling) == 1) {
      color.use.all <- c(color.use.pattern, rep(color.use.signaling, nSignaling))
    } else {
      color.use.all <- c(color.use.pattern, color.use.signaling)
    }

    if (!is.null(signaling)) {
      plot.data <- plot.data[plot.data$Signaling %in% signaling, ]
    }

    plot.data.long <- ggalluvial::to_lodes_form(plot.data, axes = 1:2, id = "connection")
    if (do.order) {
      mat = tapply(plot.data[["Contribution"]], list(plot.data[["Signaling"]], plot.data[["Pattern"]]), sum)
      mat[is.na(mat)] <- 0; mat <- mat[-which(rowSums(mat) == 0), ]
      d <- dist(as.matrix(mat))
      hc <- hclust(d, "ave")
      k <- length(unique(grep("Pattern", plot.data.long$stratum[plot.data.long$Contribution != 0], value = T)))
      cluster <- hc %>% cutree(k)
      order.name <- order(cluster)
      plot.data.long$stratum <- factor(plot.data.long$stratum, levels = c(colnames(mat),names(cluster)[order.name]))
    }

    gg2 <- ggplot(plot.data.long,aes(x = factor(x, levels = c("Pattern", "Signaling")),y= Contribution,
                                     stratum = stratum, alluvium = connection,
                                     fill = stratum, label = stratum)) +
      geom_flow(width = 1/3,aes.flow = "forward") +
      geom_stratum(width=1/3,size=0.1,color="black", alpha = 0.8, linetype = 1) +
      geom_text(stat = "stratum", size = font.size) + # 2.5
      scale_x_discrete(limits = c(),  labels=c("Patterns", "Signaling")) +
      scale_fill_manual(values = alpha(color.use.all, alpha = 0.8), drop = FALSE) +
      theme_bw()+
      theme(legend.position = "none",
            axis.title = element_blank(),
            axis.text.y= element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor  = element_blank(),
            panel.border = element_blank(),
            axis.ticks = element_blank(),axis.text=element_text(size= 10))+
      theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

    ## connect cell groups with signaling
    # data1 = data1[data1$Contribution > 0,]
    # data2 = data2[data2$Contribution > 0,]

    # data3 = merge(data1, data2, by.x="Pattern", by.y="Pattern")
    # data3$Contribution <- data3$Contribution.x * data3$Contribution.y
    # data3 <- data3[,colnames(data3) %in% c("CellGroup","Signaling","Contribution")]

    # plot.data <- data3
    # nSignaling<-length(unique(plot.data$Signaling))
    # nCellGroup<-length(unique(plot.data$CellGroup))
    #
    # if (length(color.use.signaling) == 1) {
    #   color.use.signaling <- rep(color.use.signaling, nSignaling)
    # }
    #
    #
    # ## connect cell groups with patterns
    # plot.data.long <- to_lodes_form(plot.data, axes = 1:2, id = "connection")
    # if (do.order) {
    #   mat = tapply(plot.data[["Contribution"]], list(plot.data[["CellGroup"]], plot.data[["Signaling"]]), sum)
    #   d <- dist(as.matrix(mat))
    #   hc <- hclust(d, "ave")
    #   k <- length(unique(grep("Signaling", plot.data.long$stratum[plot.data.long$Contribution != 0], value = T)))
    #   cluster <- hc %>% cutree(k)
    #   order.name <- order(cluster)
    #   plot.data.long$stratum <- factor(plot.data.long$stratum, levels = c(names(cluster)[order.name], colnames(mat)))
    #   color.use <- color.use[order.name]
    # }
    # color.use.all <- c(color.use, color.use.signaling)

    # gg3 <- ggplot(plot.data.long, aes(x = factor(x, levels = c("CellGroup", "Signaling")),y=Contribution,
    #                                  stratum = stratum, alluvium = connection,
    #                                  fill = stratum, label = stratum)) +
    #   geom_flow(width = 1/3,aes.flow = "forward") +
    #   geom_stratum(width=1/3,size=0.1,color="black", alpha = 0.8, linetype = 1) +
    #   geom_text(stat = "stratum", size = 2.5) +
    #   scale_x_discrete(limits = c(),  labels=c("Cell groups", "Signaling")) +
    #   scale_fill_manual(values = alpha(color.use.all, alpha = 0.8), drop = FALSE) +
    #   theme_bw()+
    #   theme(legend.position = "none",
    #         axis.title = element_blank(),
    #         axis.text.y= element_blank(),
    #         panel.grid.major = element_blank(),
    #         panel.grid.minor  = element_blank(),
    #         panel.border = element_blank(),
    #         axis.ticks = element_blank(),axis.text=element_text(size=10)) +
    #   theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))


    gg <- cowplot::plot_grid(gg1, gg2,align = "h", nrow = 1)
    title <- cowplot::ggdraw() + cowplot::draw_label(main.title,size = font.size.title)
    gg <- cowplot::plot_grid(title, gg, ncol=1, rel_heights=c(0.1, 1))
  }
  return(gg)
}

#' Dot plots showing the associations of latent patterns with cell groups and ligand-receptor pairs or signaling pathways
#'
#' Using a contribution score of each cell group to each signaling pathway computed by multiplying W by H obtained from `identifyCommunicationPatterns`, we constructed a dot plot in which the dot size is proportion to the contribution score to show association between cell group and their enriched signaling pathways.
#'
#' @param object CellChat object
#' @param slot.name the slot name of object that is used to compute centrality measures of signaling networks
#' @param pattern "outgoing" or "incoming"
#' @param cutoff the threshold for filtering out weak links. Default is 1/R where R is the number of latent patterns. We set the elements in W and H to be zero if they are less than `cutoff`.
#' @param color.use the character vector defining the color of each cell group
#' @param pathway.show the character vector defining the signaling to show
#' @param group.show the character vector defining the cell group to show
#' @param shape the shape of the symbol: 21 for circle and 22 for square
#' @param dot.size a range defining the size of the symbol
#' @param dot.alpha transparency
#' @param main.title the title of plot
#' @param font.size font size of the text
#' @param font.size.title font size of the title
#' @importFrom methods slot
#' @import ggplot2
#' @importFrom dplyr group_by top_n
#' @return
#' @export
#'
#' @examples
netAnalysis_dot <- function(object, slot.name = "netP", pattern = c("outgoing","incoming"), cutoff = NULL, color.use = NULL,
                            pathway.show = NULL, group.show = NULL,
                            shape = 21, dot.size = c(1, 3), dot.alpha = 1, main.title = NULL,
                            font.size = 10, font.size.title = 12){
  pattern <- match.arg(pattern)
  patternSignaling <- methods::slot(object, slot.name)$pattern[[pattern]]
  data1 = patternSignaling$pattern$cell
  data2 = patternSignaling$pattern$signaling
  data = patternSignaling$data
  if (is.null(main.title)) {
    if (pattern == "outgoing") {
      main.title = "Outgoing communication patterns of secreting cells"
    } else if (pattern == "incoming") {
      main.title = "Incoming communication patterns of target cells"
    }
  }
  if (is.null(color.use)) {
    color.use <- scPalette(nlevels(data1$CellGroup))
  }
  if (is.null(cutoff)) {
    cutoff <- 1/length(unique(data1$Pattern))
  }
  options(warn = -1)
  data1$Contribution[data1$Contribution < cutoff] <- 0
  data2$Contribution[data2$Contribution < cutoff] <- 0
  data3 = merge(data1, data2, by.x="Pattern", by.y="Pattern")
  data3$Contribution <- data3$Contribution.x * data3$Contribution.y
  data3 <- data3[,colnames(data3) %in% c("CellGroup","Signaling","Contribution")]
  if (!is.null(pathway.show)) {
    data3 <- data3[data3$Signaling %in% pathway.show, ]
    pathway.add <- pathway.show[which(pathway.show %in% data3$Signaling == 0)]
    if (length(pathway.add) > 1) {
      data.add <- expand.grid(CellGroup = levels(data1$CellGroup), Signaling = pathway.add)
      data.add$Contribution <- 0
      data3 <- rbind(data3, data.add)
    }
    data3$Signaling <- factor(data3$Signaling, levels = pathway.show)
  }
  if (!is.null(group.show)) {
    data3$CellGroup <- as.character(data3$CellGroup)
    data3 <- data3[data3$CellGroup %in% group.show, ]
    data3$CellGroup <- factor(data3$CellGroup, levels = group.show)
  }

  data <- as.data.frame(as.table(data));
  data <- data[data[,3] != 0, ]
  data12 <- paste0(data[,1],data[,2])
  data312 <- paste0(data3[,1],data3[,2])
  idx1 <- which(match(data312, data12, nomatch = 0) ==0)
  data3$Contribution[idx1] <- 0
  data3$id <- data312
  data3 <- data3 %>% group_by(id) %>% top_n(1, Contribution)

  data3$Contribution[which(data3$Contribution == 0)] <- NA

  df <- data3
  gg <- ggplot(data = df, aes(x = Signaling, y = CellGroup)) +
    geom_point(aes(size =  Contribution, fill = CellGroup, colour = CellGroup), shape = shape) +
    scale_size_continuous(range = dot.size) +
    theme_linedraw() +
    scale_x_discrete(position = "bottom") +
    ggtitle(main.title) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(text = element_text(size = font.size),plot.title = element_text(size=font.size.title, face="plain"),
          axis.text.x = element_text(angle = 45, hjust=1),
          axis.text.y = element_text(angle = 0, hjust=1),
          axis.title.x = element_blank(),
          axis.title.y = element_blank()) +
    theme(axis.line.x = element_line(size = 0.25), axis.line.y = element_line(size = 0.25)) +
    theme(panel.grid.major = element_line(colour="grey90", size = (0.1)))
  gg <- gg + scale_y_discrete(limits = rev(levels(data3$CellGroup)))
  gg <- gg + scale_fill_manual(values = ggplot2::alpha(color.use, alpha = dot.alpha), drop = FALSE, na.value = "white")
  gg <- gg + scale_colour_manual(values = color.use, drop = FALSE, na.value = "white")
  gg <- gg + guides(colour=FALSE) + guides(fill=FALSE)
  gg <- gg + theme(legend.title = element_text(size = 10), legend.text = element_text(size = 8))
  gg
  return(gg)
}


#' 2D visualization of the learned manifold of signaling networks
#'
#' @param object CellChat object
#' @param slot.name the slot name of object that is used to compute centrality measures of signaling networks
#' @param type "functional","structural"
#' @param pathway.labeled a char vector giving the signaling names to show when labeling each point
#' @param top.label the fraction of signaling pathways to label
#' @param pathway.remove a character vector defining the signaling to remove
#' @param pathway.remove.show whether show the removed signaling names
#' @param color.use defining the color for each cell group
#' @param dot.size a range defining the size of the symbol
#' @param dot.alpha transparency
#' @param xlabel label of x-axis
#' @param ylabel label of y-axis
#' @param title main title of the plot
#' @param font.size font size of the text
#' @param font.size.title font size of the title
#' @param label.size font size of the text
#' @param do.label label the each point
#' @param show.legend whether show the legend
#' @param show.axes whether show the axes
#' @import ggplot2
#' @importFrom ggrepel geom_text_repel
#' @importFrom methods slot
#' @return
#' @export
#'
#' @examples
netVisual_embedding <- function(object, slot.name = "netP", type = c("functional","structural"), color.use = NULL, pathway.labeled = NULL, top.label = 1, pathway.remove = NULL, pathway.remove.show = TRUE, dot.size = c(2, 6), label.size = 2, dot.alpha = 0.5,
                                xlabel = "Dim 1", ylabel = "Dim 2", title = NULL,
                                font.size = 10, font.size.title = 12, do.label = T, show.legend = T, show.axes = T) {
  type <- match.arg(type)
  comparison <- "single"
  comparison.name <- paste(comparison, collapse = "-")

  Y <- methods::slot(object, slot.name)$similarity[[type]]$dr[[comparison.name]]
  Groups <- methods::slot(object, slot.name)$similarity[[type]]$group[[comparison.name]]
  prob <- methods::slot(object, slot.name)$prob
  if (is.null(pathway.remove)) {
    similarity <- methods::slot(object, slot.name)$similarity[[type]]$matrix[[comparison.name]]
    pathway.remove <- rownames(similarity)[which(colSums(similarity) == 1)]
  }

  if (length(pathway.remove) > 0) {
    pathway.remove.idx <- which(dimnames(prob)[[3]] %in% pathway.remove)
    prob <- prob[ , , -pathway.remove.idx]
  }

  prob_sum <- apply(prob, 3, sum)
  df <- data.frame(x = Y[,1], y = Y[, 2], Commun.Prob. = prob_sum/max(prob_sum), labels = as.character(unlist(dimnames(prob)[3])), Groups = as.factor(Groups))
  if (is.null(color.use)) {
    color.use <- ggPalette(length(unique(Groups)))
  }
  gg <- ggplot(data = df, aes(x, y)) +
    geom_point(aes(size = Commun.Prob.,fill = Groups, colour = Groups), shape = 21) +
    CellChat_theme_opts() +
    theme(text = element_text(size = font.size), legend.key.height = grid::unit(0.15, "in"))+
    guides(colour = guide_legend(override.aes = list(size = 3)))+
    labs(title = title, x = xlabel, y = ylabel) + theme(plot.title = element_text(size= font.size.title, face="plain"))+
    scale_size_continuous(limits = c(0,1), range = dot.size, breaks = c(0.1,0.5,0.9)) +
    theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank()) +
    theme(axis.line.x = element_line(size = 0.25), axis.line.y = element_line(size = 0.25))
  gg <- gg + scale_fill_manual(values = ggplot2::alpha(color.use, alpha = dot.alpha), drop = FALSE)
  gg <- gg + scale_colour_manual(values = color.use, drop = FALSE)
  if (do.label) {
    if (is.null(pathway.labeled)) {
      if (top.label < 1) {
        if (length(comparison) == 2) {
          g.t <- rankSimilarity(object, slot.name = slot.name, type = type, comparison1 = comparison)
          pathway.labeled <- as.character(g.t$data$name[(nrow(g.t$data)-ceiling(top.label * nrow(g.t$data))+1):nrow(g.t$data) ])
          data.label <- df[df$labels %in% pathway.labeled, , drop = FALSE]
        }
      } else {
        data.label <- df
      }

    } else {
      data.label <- df[df$labels %in% pathway.labeled, , drop = FALSE]
    }
    gg <- gg + ggrepel::geom_text_repel(data = data.label, mapping = aes(label = labels, colour = Groups), size = label.size, show.legend = F,segment.size = 0.2, segment.alpha = 0.5) + scale_alpha_discrete(range = c(1, 0.6))

    # gg <- gg + ggrepel::geom_text_repel(mapping = aes(label = labels, colour = Groups), size = label.size, show.legend = F,segment.size = 0.2, segment.alpha = 0.5)
  }

  if (length(pathway.remove) > 0 & pathway.remove.show) {
    gg <- gg + annotate(geom = 'text', label =  paste("Isolate pathways: ", paste(pathway.remove, collapse = ', ')), x = -Inf, y = Inf, hjust = 0, vjust = 1, size = label.size,fontface="italic")
  }
  if (!show.legend) {
    gg <- gg + theme(legend.position = "none")
  }

  if (!show.axes) {
    gg <- gg + theme_void()
  }
  gg
}


#' Zoom into the 2D visualization of the learned manifold learning of the signaling networks
#'
#' @param object CellChat object
#' @param slot.name the slot name of object that is used to compute centrality measures of signaling networks
#' @param type "functional","structural"
#' @param pathway.remove a character vector defining the signaling to remove
#' @param color.use defining the color for each cell group
#' @param nCol the number of columns of the plot
#' @param dot.size a range defining the size of the symbol
#' @param dot.alpha transparency
#' @param xlabel label of x-axis
#' @param ylabel label of y-axis
#' @param label.size font size of the text
#' @param do.label label the each point
#' @param show.legend whether show the legend
#' @param show.axes whether show the axes
#' @import ggplot2
#' @importFrom ggrepel geom_text_repel
#' @importFrom cowplot plot_grid
#' @importFrom methods slot
#' @return
#' @export
#'
#' @examples
netVisual_embeddingZoomIn <- function(object, slot.name = "netP", type = c("functional","structural"), color.use = NULL, pathway.remove = NULL,  nCol = 1, dot.size = c(2, 6), label.size = 2.8, dot.alpha = 0.5,
                                      xlabel = NULL, ylabel = NULL, do.label = T, show.legend = F, show.axes = T) {
  comparison <- "single"
  comparison.name <- paste(comparison, collapse = "-")
  Y <- methods::slot(object, slot.name)$similarity[[type]]$dr[[comparison.name]]
  clusters <- methods::slot(object, slot.name)$similarity[[type]]$group[[comparison.name]]
  prob <- methods::slot(object, slot.name)$prob
  if (is.null(pathway.remove)) {
    similarity <- methods::slot(object, slot.name)$similarity[[type]]$matrix[[comparison.name]]
    pathway.remove <- rownames(similarity)[which(colSums(similarity) == 1)]
  }

  if (length(pathway.remove) > 0) {
    pathway.remove.idx <- which(dimnames(prob)[[3]] %in% pathway.remove)
    prob <- prob[ , , -pathway.remove.idx]
  }

  prob_sum <- apply(prob, 3, sum)
  df <- data.frame(x = Y[,1], y = Y[, 2], Commun.Prob. = prob_sum/max(prob_sum), labels = as.character(unlist(dimnames(prob)[3])), clusters = as.factor(clusters))

  if (is.null(color.use)) {
    color.use <- ggPalette(length(unique(clusters)))
  }

  # zoom into each cluster and do labels
  ggAll <- vector("list", length(unique(clusters)))
  for (i in 1:length(unique(clusters))) {
    clusterID = i
    title <- paste0("Group ",  clusterID)
    df2 <- df[df$clusters %in% clusterID,]
    gg <- ggplot(data = df2, aes(x, y)) +
      geom_point(aes(size = Commun.Prob.), shape = 21, colour = alpha(color.use[clusterID], alpha = 1), fill = alpha(color.use[clusterID], alpha = dot.alpha)) +
      CellChat_theme_opts() +
      theme(text = element_text(size = 10), legend.key.height = grid::unit(0.15, "in"))+
      labs(title = title, x = xlabel, y = ylabel) + theme(plot.title = element_text(size=12))+
      scale_size_continuous(limits = c(0,1), range = dot.size, breaks = c(0.1,0.5,0.9)) +
      theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank()) +
      theme(axis.line.x = element_line(size = 0.25), axis.line.y = element_line(size = 0.25))
    if (do.label) {
      gg <- gg + ggrepel::geom_text_repel(mapping = aes(label = labels), colour = color.use[clusterID], size = label.size, segment.size = 0.2, segment.alpha = 0.5)
    }

    if (!show.legend) {
      gg <- gg + theme(legend.position = "none")
    }

    if (!show.axes) {
      gg <- gg + theme_void()
    }
    ggAll[[i]] <- gg
  }
  gg.combined <- cowplot::plot_grid(plotlist = ggAll, ncol = nCol)

  gg.combined

}



#' 2D visualization of the joint manifold learning of signaling networks from two datasets
#'
#' @param object CellChat object
#' @param slot.name the slot name of object that is used to compute centrality measures of signaling networks
#' @param type "functional","structural"
#' @param comparison a numerical vector giving the datasets for comparison. Default are all datasets when object is a merged object
#' @param pathway.labeled a char vector giving the signaling names to show when labeling each point
#' @param top.label the fraction of signaling pathways to label
#' @param pathway.remove a character vector defining the signaling to remove
#' @param pathway.remove.show whether show the removed signaling names
#' @param color.use defining the color for each cell group
#' @param point.shape a numeric vector giving the point shapes. By default point.shape <- c(21, 0, 24, 23, 25, 10, 12), see available shapes at http://www.sthda.com/english/wiki/r-plot-pch-symbols-the-different-point-shapes-available-in-r
#' @param dot.size a range defining the size of the symbol
#' @param dot.alpha transparency
#' @param xlabel label of x-axis
#' @param ylabel label of y-axis
#' @param title main title of the plot
#' @param label.size font size of the text
#' @param do.label label the each point
#' @param show.legend whether show the legend
#' @param show.axes whether show the axes
#' @import ggplot2
#' @importFrom ggrepel geom_text_repel
#' @importFrom methods slot
#' @return
#' @export
#'
#' @examples
netVisual_embeddingPairwise <- function(object, slot.name = "netP", type = c("functional","structural"), comparison = NULL, color.use = NULL, point.shape = NULL, pathway.labeled = NULL, top.label = 1, pathway.remove = NULL, pathway.remove.show = TRUE, dot.size = c(2, 6), label.size = 2.5, dot.alpha = 0.5,
                                        xlabel = "Dim 1", ylabel = "Dim 2", title = NULL,do.label = T, show.legend = T, show.axes = T) {
  type <- match.arg(type)
  if (is.null(comparison)) {
    comparison <- 1:length(unique(object@meta$datasets))
  }
  cat("2D visualization of signaling networks from datasets", as.character(comparison), '\n')
  comparison.name <- paste(comparison, collapse = "-")

  Y <- methods::slot(object, slot.name)$similarity[[type]]$dr[[comparison.name]]
  clusters <- methods::slot(object, slot.name)$similarity[[type]]$group[[comparison.name]]
  object.names <- setdiff(names(methods::slot(object, slot.name)), "similarity")[comparison]
  prob <- list()
  for (i in 1:length(comparison)) {
    object.net <- methods::slot(object, slot.name)[[comparison[i]]]
    prob[[i]] = object.net$prob
  }

  if (is.null(point.shape)) {
    point.shape <- c(21, 0, 24, 23, 25, 10, 12)
  }

  if (is.null(pathway.remove)) {
    similarity <- methods::slot(object, slot.name)$similarity[[type]]$matrix[[comparison.name]]
    pathway.remove <- rownames(similarity)[which(colSums(similarity) == 1)]
    # pathway.remove <- sub("--.*", "", pathway.remove)
  }

  if (length(pathway.remove) > 0) {
    for (i in 1:length(prob)) {
      probi <- prob[[i]]
      pathway.remove.idx <- which(paste0(dimnames(probi)[[3]],"--",object.names[i]) %in% pathway.remove)
      #  pathway.remove.idx <- which(dimnames(probi)[[3]] %in% pathway.remove)
      if (length(pathway.remove.idx) > 0) {
        probi <- probi[ , , -pathway.remove.idx]
      }
      prob[[i]] <- probi
    }
  }
  prob_sum.each <- list()
  signalingAll <- c()
  for (i in 1:length(prob)) {
    probi <- prob[[i]]
    prob_sum.each[[i]] <- apply(probi, 3, sum)
    signalingAll <- c(signalingAll, paste0(names(prob_sum.each[[i]]),"--",object.names[i]))
  }
  prob_sum <- unlist(prob_sum.each)
  names(prob_sum) <- signalingAll

  group <- sub(".*--", "", names(prob_sum))
  labels = sub("--.*", "", names(prob_sum))

  df <- data.frame(x = Y[,1], y = Y[, 2], Commun.Prob. = prob_sum/max(prob_sum),
                   labels = as.character(labels), clusters = as.factor(clusters), group = factor(group, levels = unique(group)))
  # color dots (light inside color and dark border) based on clustering and no labels
  if (is.null(color.use)) {
    color.use <- ggPalette(length(unique(clusters)))
  }
  gg <- ggplot(data = df, aes(x, y)) +
    geom_point(aes(size = Commun.Prob.,fill = clusters, colour = clusters, shape = group)) +
    CellChat_theme_opts() +
    theme(text = element_text(size = 10), legend.key.height = grid::unit(0.15, "in"))+
    guides(colour = guide_legend(override.aes = list(size = 3)))+
    labs(title = title, x = xlabel, y = ylabel) +
    scale_size_continuous(limits = c(0,1), range = dot.size, breaks = c(0.1,0.5,0.9)) +
    theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank()) +
    theme(axis.line.x = element_line(size = 0.25), axis.line.y = element_line(size = 0.25))
  gg <- gg + scale_fill_manual(values = ggplot2::alpha(color.use, alpha = dot.alpha), drop = FALSE) #+ scale_alpha(group, range = c(0.1, 1))
  gg <- gg + scale_colour_manual(values = color.use, drop = FALSE)
  gg <- gg + scale_shape_manual(values = point.shape[1:length(prob)])
  if (do.label) {
    gg <- gg + ggrepel::geom_text_repel(mapping = aes(label = labels, colour = clusters, alpha=group), size = label.size, show.legend = F,segment.size = 0.2, segment.alpha = 0.5) + scale_alpha_discrete(range = c(1, 0.6))
  }

  if (length(pathway.remove) > 0 & pathway.remove.show) {
    gg <- gg + annotate(geom = 'text', label =  paste("Isolate pathways: ", paste(pathway.remove, collapse = ', ')), x = -Inf, y = Inf, hjust = 0, vjust = 1, size = label.size,fontface="italic")
  }

  if (!show.legend) {
    gg <- gg + theme(legend.position = "none")
  }

  if (!show.axes) {
    gg <- gg + theme_void()
  }
  gg
}



#' Zoom into the 2D visualization of the joint manifold learning of signaling networks from two datasets
#'
#' @param object CellChat object
#' @param slot.name the slot name of object that is used to compute centrality measures of signaling networks
#' @param type "functional","structural"
#' @param comparison a numerical vector giving the datasets for comparison. Default are all datasets when object is a merged object
#' @param pathway.remove a character vector defining the signaling to remove
#' @param color.use defining the color for each cell group
#' @param nCol number of columns in the plot
#' @param point.shape a numeric vector giving the point shapes. By default point.shape <- c(21, 0, 24, 23, 25, 10, 12), see available shapes at http://www.sthda.com/english/wiki/r-plot-pch-symbols-the-different-point-shapes-available-in-r
#' @param dot.size a range defining the size of the symbol
#' @param dot.alpha transparency
#' @param xlabel label of x-axis
#' @param ylabel label of y-axis
#' @param label.size font size of the text
#' @param do.label label the each point
#' @param show.legend whether show the legend
#' @param show.axes whether show the axes
#' @import ggplot2
#' @importFrom ggrepel geom_text_repel
#' @importFrom methods slot
#' @return
#' @export
#'
#' @examples
netVisual_embeddingPairwiseZoomIn <- function(object, slot.name = "netP", type = c("functional","structural"), comparison = NULL, color.use = NULL, nCol = 1, point.shape = NULL, pathway.remove = NULL, dot.size = c(2, 6), label.size = 2.8, dot.alpha = 0.5,
                                              xlabel = NULL, ylabel = NULL, do.label = T, show.legend = F, show.axes = T) {

  type <- match.arg(type)
  if (is.null(comparison)) {
    comparison <- 1:length(unique(object@meta$datasets))
  }
  cat("2D visualization of signaling networks from datasets", as.character(comparison), '\n')
  comparison.name <- paste(comparison, collapse = "-")

  Y <- methods::slot(object, slot.name)$similarity[[type]]$dr[[comparison.name]]
  clusters <- methods::slot(object, slot.name)$similarity[[type]]$group[[comparison.name]]
  object.names <- setdiff(names(methods::slot(object, slot.name)), "similarity")[comparison]
  prob <- list()
  for (i in 1:length(comparison)) {
    object.net <- methods::slot(object, slot.name)[[comparison[i]]]
    prob[[i]] = object.net$prob
  }

  if (is.null(point.shape)) {
    point.shape <- c(21, 0, 24, 23, 25, 10, 12)
  }

  if (is.null(pathway.remove)) {
    similarity <- methods::slot(object, slot.name)$similarity[[type]]$matrix[[comparison.name]]
    pathway.remove <- rownames(similarity)[which(colSums(similarity) == 1)]
    # pathway.remove <- sub("--.*", "", pathway.remove)
  }

  if (length(pathway.remove) > 0) {
    for (i in 1:length(prob)) {
      probi <- prob[[i]]
      pathway.remove.idx <- which(paste0(dimnames(probi)[[3]],"--",object.names[i]) %in% pathway.remove)
      #  pathway.remove.idx <- which(dimnames(probi)[[3]] %in% pathway.remove)
      if (length(pathway.remove.idx) > 0) {
        probi <- probi[ , , -pathway.remove.idx]
      }
      prob[[i]] <- probi
    }
  }

  prob_sum.each <- list()
  signalingAll <- c()
  for (i in 1:length(prob)) {
    probi <- prob[[i]]
    prob_sum.each[[i]] <- apply(probi, 3, sum)
    signalingAll <- c(signalingAll, paste0(names(prob_sum.each[[i]]),"--",object.names[i]))
  }
  prob_sum <- unlist(prob_sum.each)
  names(prob_sum) <- signalingAll

  group <- sub(".*--", "", names(prob_sum))
  labels = sub("--.*", "", names(prob_sum))

  df <- data.frame(x = Y[,1], y = Y[, 2], Commun.Prob. = prob_sum/max(prob_sum),
                   labels = as.character(labels), clusters = as.factor(clusters), group = factor(group, levels = unique(group)))
  if (is.null(color.use)) {
    color.use <- ggPalette(length(unique(clusters)))
  }

  # zoom into each cluster and do labels
  ggAll <- vector("list", length(unique(clusters)))
  for (i in 1:length(unique(clusters))) {
    clusterID = i
    title <- paste0("Cluster ",  clusterID)
    df2 <- df[df$clusters %in% clusterID,]
    gg <- ggplot(data = df2, aes(x, y)) +
      geom_point(aes(size = Commun.Prob., shape = group),fill = alpha(color.use[clusterID], alpha = dot.alpha), colour = alpha(color.use[clusterID], alpha = 1)) +
      CellChat_theme_opts() +
      theme(text = element_text(size = 10), legend.key.height = grid::unit(0.15, "in"))+
      guides(colour = guide_legend(override.aes = list(size = 3)))+
      labs(title = title, x = xlabel, y = ylabel) +
      scale_size_continuous(limits = c(0,1), range = dot.size, breaks = c(0.1,0.5,0.9)) +
      theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank()) +
      theme(axis.line.x = element_line(size = 0.25), axis.line.y = element_line(size = 0.25))
    idx <- match(unique(df2$group), levels(df$group), nomatch = 0)
    gg <- gg + scale_shape_manual(values= point.shape[idx])
    if (do.label) {
      gg <- gg + ggrepel::geom_text_repel(mapping = aes(label = labels), colour = color.use[clusterID], size = label.size, show.legend = F,segment.size = 0.2, segment.alpha = 0.5) + scale_alpha_discrete(range = c(1, 0.6))
    }

    if (!show.legend) {
      gg <- gg + theme(legend.position = "none")
    }

    if (!show.axes) {
      gg <- gg + theme_void()
    }
    ggAll[[i]] <- gg
  }
  gg.combined <- cowplot::plot_grid(plotlist = ggAll, ncol = nCol)

  gg.combined

}


#' Stacked bar plot showing the proportion of cells across certrain cell groups
#'
#' @param object SpatialCellChat object or seurat object
#' @param n.colors Number of colors when setting colors.use to be a palette name from brewer.pal
#' @param title.name Name of the main title
#' @param legend.title Name of legend
#' @param xlabel Name of x label
#' @param ylabel Name of y label
#' @param y.group grouping variable used to define rows of the topic composition matrix
#' @param x.group grouping variable used to define columns of the topic composition matrix
#' @param x.levels a character vector specifying the order of column groups
#' @param y.levels a character vector specifying the order of row groups
#' @param scale.rows whether to normalize proportions across rows.
#'   When TRUE (default), values are scaled within each row to sum to 1.
#'   When FALSE, values are scaled within each column.
#' @param cutoff.prop numeric threshold for filtering low-proportion groups
#' @param trim the fraction (0 to 0.5) of observations to be trimmed from each end of x before the mean is computed
#' @param type visualization type. Either "heatmap" (default) or "dot"
#' @param color.use defining the color for each column group
#' @param color.heatmap a character string or vector indicating the colormap option to use. It can be the avaibale color palette in brewer.pal() or viridis_pal()
#' @param dot.size numeric vector of length two specifying the minimum and maximum point sizes when `type = "dot"`
#' @param angle.x angle for x-axis text rotation
#' @param font.size font size in heatmap
#' @param font.size.title font size of the title
#' @param cluster.rows whether cluster rows
#' @param cluster.cols whether cluster columns
#' @param clustering_distance_rows distance metric used for row clustering in the heatmap
#' @param row.show,col.show a vector giving the index or the name of row or columns to show in the heatmap
#' @param remove.isolate whether remove the isolate nodes in the communication network
#' @param x.lab.rot rotation angle of x-axis tick labels
#'
#' @return ggplot2 object
#' @export
#' @examples
#'
#' @import ggplot2
netVisual_TopicComposition <- function(
    object,
    y.group,
    x.group = NULL,
    x.levels = NULL,
    y.levels = NULL,
    scale.rows = TRUE,
    cutoff.prop = 0.1,
    trim = 0.1,
    type = "heatmap",
    color.use = NULL,
    color.heatmap = "RdPu",
    n.colors = 8,
    dot.size = c(1, 6),
    angle.x = 45,
    xlabel = NULL,
    ylabel = NULL,
    title.name = NULL,
    legend.title = NULL,
    font.size = 8,
    font.size.title = 10,
    cluster.rows = TRUE,
    cluster.cols = FALSE,
    clustering_distance_rows = "euclidean",
    x.lab.rot = 45,
    row.show = NULL,
    col.show = NULL,
    remove.isolate = TRUE
) {
  if (is.character(x.group) | is.null(x.group)) {
    if (is.null(x.group)) {
      if (inherits(object,what = c("SpatialCellChat"))) {
        x.cell.state <- object@idents
        x.cell.state.level <- levels(x.cell.state)
      } else if(inherits(object,what = c("Seurat"))) {
        x.cell.state <- Seurat::Idents(object)
        x.cell.state.level <- levels(x.cell.state)
      }
    } else if (x.group %in% colnames(object@meta.data) == TRUE) {
      x.cell.state <- object@meta.data[, x.group]
      x.cell.state.level <- levels(x.cell.state)
    } else if (x.group %in% colnames(object@meta.data) == FALSE) {
      stop("'x.group' is not the column of `object@meta.data`! \n")
    }

    if (is.character(y.group)) {
      if (y.group %in% c("topicOut", "topicIn","inmfNorm") == TRUE) {
        if (inherits(object,what = c("SpatialCellChat"))) {
          topic.key <- list("topicIn"="incoming","topicOut"="outgoing")
          pattern <- topic.key[[y.group]]
          pred <- t(methods::slot(object, "net")$topic[[pattern]]$cell)
        }else if(inherits(object,what = c("Seurat"))) {
          pred <- t(object@reductions[[y.group]]@cell.embeddings)
        }

        if (y.group == "topicOut") {
          if (is.null(title.name)) {
            title.name <- paste0("outgoing", " CCC topics")
          }

          #ylabel <- "Outgoing communication topics"
        } else if (y.group == "topicIn") {
          if (is.null(title.name)) {
            title.name <- paste0("incoming", " CCC topics")
          }
          #ylabel <- "Incoming communication topics"
        }
      } else if (y.group %in% colnames(object@meta.data) == TRUE) {
        y.cell.state <- object@meta.data[, y.group]
        y.cell.state.level <- levels(y.cell.state)
        pred <-
          matrix(0,
                 nrow = length(y.cell.state.level),
                 ncol = length(y.cell.state))
        for (i in 1:length(y.cell.state.level)) {
          pred[i, y.cell.state == y.cell.state.level[i]] <- 1
        }
        rownames(pred) <- y.cell.state.level
      } else {
        stop("'y.group' is not the column of `object@meta.data`! \n")
      }
    } else if (
      inherits(y.group,
               what = c("data.frame", "matrix", "Matrix", "dgCMatrix"))
    ) {
      pred <- y.group
      cat(
        "The rownames of the input `y.group` are ",
        toString(rownames(y.group)),
        " which will be used as group names. \n"
      )
    } else {
      stop("Please check your input `y.group`! \n")
    }

    labels <- x.cell.state
    labels.level <- x.cell.state.level
    prop <-
      matrix(0, nrow = length(labels.level), ncol = nrow(pred))
    for (i in 1:length(labels.level)) {
      cell.use <- which(labels == labels.level[i])
      if (length(cell.use) > 0) {
        prop[i,] <- apply(pred[, cell.use, drop = FALSE], 1, function(x) thresholdedMean(x, trim = trim, na.rm = TRUE))
      }
    }
    colnames(prop) <- rownames(pred)
    rownames(prop) <- labels.level
    prop <- t(prop)

  } else if (inherits(x.group, what = c("data.frame", "matrix", "Matrix", "dgCMatrix"))) {
    pred <- x.group # this is the prediction score of each cell type for each spot for 10X visium
    cat("The rownames of the input `x.group` are ", toString(rownames(x.group)), " which will be used as group names. \n")
    if (all(y.group %in% colnames(object@meta.data)) == TRUE) {
      y.cell.state <- object@meta.data[, y.group]
      y.cell.state.level <- levels(y.cell.state)

      labels <- y.cell.state
      labels.level <- y.cell.state.level
      prop <-
        matrix(0, nrow = length(labels.level), ncol = nrow(pred))
      for (i in 1:length(labels.level)) {
        cell.use <- which(labels == labels.level[i])
        prop[i,] <-
          apply(pred[, cell.use, drop = FALSE], 1, function(x)
            thresholdedMean(x, trim = trim, na.rm = TRUE))
      }
      colnames(prop) <- rownames(pred)
      rownames(prop) <- labels.level
    } else if (
      inherits(y.group,
               what = c("data.frame", "matrix", "Matrix", "dgCMatrix"))
    ) {
      pred <- y.group
      cat(
        "The rownames of the input `y.group` are ",
        toString(rownames(y.group)),
        ", which will be used as the names of another group. \n"
      )

      y.cell.state <- rownames(y.group)
      y.cell.state <- factor(y.cell.state,levels = unique(y.cell.state))
      y.cell.state.level <- levels(y.cell.state)

      labels <- y.cell.state
      labels.level <- y.cell.state.level
      prop <- Matrix::as.matrix(x.group) %*% Matrix::as.matrix(t(y.group))
      colnames(prop) <- rownames(pred)
      rownames(prop) <- rownames(x.group)
      prop <- t(prop)
    }
  } else {
    stop("Please check your input `x.group`! \n")
  }

  mat <- prop

  if (scale.rows) {
    mat <- sweep(mat, 1L, apply(mat, 1, sum), '/', check.margin = FALSE)
  } else {
    mat <- sweep(mat, 2L, apply(mat, 2, sum), '/', check.margin = FALSE)
  }
  mat[is.na(mat)] <- 0
  if (remove.isolate) {
    idx <- which(apply(mat, 1, max) < cutoff.prop)
    if (length(idx) > 0) {
      mat <- mat[,-idx]
    }
  }

  if (!is.null(x.levels)) {
    x.levels <- x.levels[x.levels %in% colnames(mat)]
    mat <- mat[, order(factor(colnames(mat), levels = x.levels)), drop = FALSE]
  }
  if (!is.null(y.levels)) {
    y.levels <- y.levels[y.levels %in% rownames(mat)]
    mat <-
      mat[order(factor(rownames(mat), levels = y.levels)), , drop = FALSE]
  }

  if (is.null(color.use)) {
    color.use <- scPalette(NCOL(mat))
  }
  names(color.use) <- colnames(mat)
  if (!is.null(row.show)) {
    mat <- mat[row.show,]
  }
  if (!is.null(col.show)) {
    mat <- mat[, col.show]
    color.use <- color.use[col.show]
  }

  if (type == "dot") {
    df <- as.data.frame(as.table(mat))
    colnames(df) <- c("x", "y", "proportion")
    gg <-
      ggplot(df, aes(x, y, size = proportion, color = proportion)) +
      geom_point() +
      scale_size_continuous(range = dot.size) +
      theme_linedraw() +
      labs(x = xlabel,
           y = ylabel,
           size = "Proportion") +
      scale_color_viridis_c(option = "D") +
      theme(legend.key.height = grid::unit(0.15, "in")) +
      guides(color = guide_colourbar(barwidth = 0.5, title = "Proportion"))
    # gg <- gg + guides(color = guide_colorbar(barwidth = legend.width, title = "Scaled expression"),size = guide_legend(title = 'Percent expressed'))
    gg <- gg + theme(
      text = element_text(size = 10),
      axis.text.x = element_text(angle = x.lab.rot, hjust =
                                   1),
      axis.text.y = element_text(angle = 0, hjust = 1),
      axis.title.x = element_blank(),
      axis.title.y = element_blank()
    ) +
      theme(
        axis.line.x = element_line(linewidth = 0.25),
        axis.line.y = element_line(linewidth = 0.25)
      ) +
      theme(panel.grid.major = element_line(colour = "grey90", linewidth = (0.1)))
    return(gg)

  } else if (type == "heatmap") {
    color.heatmap.use = grDevices::colorRampPalette((RColorBrewer::brewer.pal(n = 9, name = color.heatmap)))(100)

    df <-
      data.frame(group = colnames(mat))
    rownames(df) <- colnames(mat)
    col_annotation <-
      ComplexHeatmap::HeatmapAnnotation(
        df = df,
        col = list(group = color.use),
        which = "column",
        show_legend = FALSE,
        show_annotation_name = FALSE,
        simple_anno_size = grid::unit(0.2, "cm")
      )

    # ha1 = rowAnnotation(Strength = anno_barplot(rowSums(abs(mat)), border = FALSE,gp = grid::gpar(fill = color.use, col=color.use)), show_annotation_name = FALSE)
    # mat[mat == 0] <- NA
    color.heatmap.use = c("white", color.heatmap.use)
    ht1 = ComplexHeatmap::Heatmap(
      mat,
      col = color.heatmap.use,
      na_col = "white",
      name = "Proportion",
      bottom_annotation = col_annotation,
      cluster_rows = cluster.rows,
      cluster_columns = cluster.cols,
      clustering_distance_rows = clustering_distance_rows,
      row_names_side = "left",
      row_names_rot = 0,
      row_names_gp = grid::gpar(fontsize = font.size),
      column_names_gp = grid::gpar(fontsize = font.size),
      # width = grid::unit(width, "cm"), height = grid::unit(height, "cm"),
      row_title = ylabel,
      row_title_gp = grid::gpar(fontsize = font.size.title),
      column_title = paste0("Compositions of ", title.name),
      column_title_gp = grid::gpar(fontsize = font.size.title),
      column_names_rot = x.lab.rot,
      heatmap_legend_param = list(
        title = "Proportion",
        title_gp = grid::gpar(fontsize = 8, fontface = "plain"),
        title_position = "leftcenter-rot",
        border = NA,
        legend_height = grid::unit(20, "mm"),
        labels_gp = grid::gpar(fontsize = 8),
        grid_width = grid::unit(2, "mm")
      )
    )

    return(ht1)
  }
}


#' @title netVisual_TopicSignaling
#' @description
#' Ranking the signalings(LR pairs or pathways) and show the top marker signalings in each factor
#'
#' @param object Seurat object or SpatialCellChat object
#' @param topics NULL or Numeric. Determine which factors should be displayed.
#' By default, use `NULL` to display all factors' feature ranking plot.
#' @param slot.name "net" or "netP".
#' @param pattern "outgoing" or "incoming"
#' @param ncol number of columns in plot
#' @param topic.key character. By default, "Topic_"
#' @param nfeatures numeric. Determine how many top marker signalings in each factor to show
#' @param nTop numeric. Determine how many top marker signalings in each factor to show with a text annotation
#' @param feature_cat a character string specifying the feature category to summarize loadings, one of \code{"LR-pair"}, \code{"ligand"}, or \code{"signaling"}
#' @param col numeric of length 2. Colors of points to use, default is \code{c('neg'="#377EB8", 'pos'="#E41A1C")}
#' @param combine combine plots into a single patchwork ggplot object. If FALSE, return a list of ggplot objects
#' @param ylabel name of y label
#' @param balanced return an equal number of signalings with + and - scores. If FALSE (default), returns the top signalings
#'  ranked by the scores absolute values
#' @param text_nudge_x horizontal and vertical adjustments to nudge the starting position of each text label. The units for
#' nudge_x and nudge_y are the same as for the data units on the x-axis and y-axis
#' @param text_size numeric. Fontsize of text annotation
#' @param point_size numeric. Size of points
#' @param line_size numeric. Width of lines
#' @param species a character string specifying the species used to select the default CellChat database, either \code{"mouse"} or \code{"human"}
#' @param db an optional interaction database
#' @param simplify whether to combine results from all topics into a single data frame (TRUE) or return a list split by topic (FALSE)
#' @import ggplot2
#'
#' @return ggplot2 object
#' @export
netVisual_TopicSignaling <- function(
    object,
    pattern = c("incoming","outgoing"),
    slot.name = "net",
    topics = NULL,
    topic.key = "Topic_",
    nfeatures = 30,
    nTop = 10,
    feature_cat = c("LR-pair","ligand", "signaling"),
    col = c("black","black"),
    ncol = NULL,
    combine = T,
    ylabel = NULL,
    balanced = FALSE,
    text_nudge_x = -8,
    text_size = 3,
    point_size = 2,
    line_size = 1,
    species = NULL, # c('mouse','human')
    db = NULL,
    simplify=TRUE
){
  pattern <- match.arg(pattern)
  feature_cat <- match.arg(feature_cat)
  species <- match.arg(species)
  if (inherits(object,what = c("SpatialCellChat"))){
    pattern.res <- methods::slot(object, slot.name)$topic[[pattern]]
    loadings <- pattern.res$signaling  # signaling.name x topic
    # if(is.null(colnames(loadings)))
    colnames(loadings) <- stringr::str_c(topic.key,1:NCOL(loadings)) # signaling.name x topic
  } else if(inherits(object,what = c("Seurat"))){
    # topic.key <- list("incoming"="topicIn","outgoing"="topicOut")
    # topic.key.name <- topic.key[[pattern]]
    topic.key.name <- paste0(topic.key,".",pattern)
    loadings <- object@reductions[[topic.key.name]]@feature.loadings

    # loadings: nFeatures x nPCs
    loadings <- Loadings(object = object[[reduction]], projected = projected)

    # if(is.null(colnames(loadings)))
    colnames(loadings) <- stringr::str_c(topic.key,1:NCOL(loadings)) # signaling.name x topic
  } else if(inherits(object,what = c("matrix","array"))){
    loadings <- object
    # if(is.null(colnames(loadings)))
    colnames(loadings) <- stringr::str_c(topic.key,1:NCOL(loadings)) # signaling.name x topic
  }
  if (inherits(object,what = c("SpatialCellChat"))) {
    CellChatDB <- object@DB
  } else {
    if (is.null(db)) {
      if (species == "mouse") {
        CellChatDB <- CellChatDB.mouse
      } else if (species == 'human') {
        CellChatDB <- CellChatDB.human
      } else {
        stop("Only mouse and human are supported currently. Please provide a `db` instead! ")
      }
    } else if (!is.null(db)) {
      CellChatDB <- db
      if (all(c("ligand","receptor") %in% colnames(db)) == FALSE) {
        stop("The input `db` must contain at least two columns named as ligand,receptor")
      }
      if (all(c("interaction_name") %in% colnames(db)) == FALSE) {
        db$interaction_name <- paste0(toupper(db$ligand), "_", toupper(db$receptor))
      }
      if (all(c("interaction_name_2") %in% colnames(db)) == FALSE) {
        db$interaction_name_2 <- paste0(db$ligand, " - ", db$receptor)
      }
      if (feature_cat == "signaling") {
        if (all(c("pathway_name") %in% colnames(db)) == FALSE) {
          stop("The input `db` must contain a column named as pathway_name when feature_cat == signaling !")
        }
      }
    } else {
      stop("Please provide either `species` or db` when the input is not a SpatialCellChat object! ")
    }
  }


  if (is.null(topics)) {
    topics <- seq_len(NCOL(loadings))
  }

  # Column number of combined plots
  if (is.null(x = ncol)) {
    ncol <- 2
    if (length(x = topics) == 1) {
      ncol <- 1
    }
    if (length(x = topics) > 6) {
      ncol <- 3
    }
    if (length(x = topics) > 9) {
      ncol <- 4
    }
  }


  # define `Top`
  MyTop <- function (data, num, balanced)
  {
    nr <- nrow(x = data)
    if (num > nr) {
      # warning("Requested number is larger than the number of available items (",
      #         nr, "). Setting to ", nr, ".", call. = FALSE)
      num <- nr
    }
    if (num == 1) {
      balanced <- FALSE
    }
    top <- if (balanced) {
      num <- round(x = num/2)
      data <- data[order(data, decreasing = TRUE), , drop = FALSE]
      positive <- head(x = rownames(x = data), n = num)
      negative <- rev(x = tail(x = rownames(x = data), n = num))
      if (positive[num] == negative[num]) {
        negative <- negative[-num]
      }
      list(positive = positive, negative = negative)
    }
    else {
      data <- data[rev(x = order(abs(x = data))), , drop = FALSE]
      top <- head(x = rownames(x = data), n = num)
      top[order(data[top, ])]
    }
    return(top)
  }
  if (balanced) {
    if (length(col) == 1) {
      stop("Please input `col` as a vector with length 2 when `balanced` is TRUE! \n")
    }
  }

  # Get features of each plot
  features <- lapply(
    X = topics,
    FUN = function(topics){
      data.fa <- loadings[, topics, drop = FALSE]
      MyTop(data = data.fa, num = nfeatures, balanced = balanced)
    }
  )
  features <- lapply(X = features, FUN = unlist, use.names = FALSE)


  # Subset loadings: mFeatures x mPCs
  loadings <- loadings[unlist(x = features), topics, drop = FALSE]
  names(features) <- colnames(loadings) <- as.character(topics)

  plots <- lapply(
    X = as.character(topics),
    FUN = function(i) {

      # get loading of PC i
      df <- as.data.frame(loadings[features[[i]],i, drop = FALSE])
      colnames(df) <- "score"
      #rownames(df) <- gsub("-", "_", rownames(df))
      df$feature <- rownames(df)

      if (feature_cat %in% c("ligand","signaling") ) {
        db.use <- CellChatDB$interaction[df$feature, c("ligand",'receptor','pathway_name')]
        df <- cbind(df, db.use)
        # summarize the loading score
        if (feature_cat == "ligand") {
          data.plot <- df %>% group_by(ligand) %>% summarize(avg = mean(score))  # total = sum(score),
        } else {
          data.plot <- df %>% group_by(pathway_name) %>% summarize(avg = mean(score))  # avg = mean(score),
        }
        colnames(data.plot) <- c("feature", "score")
      } else {
        idx <- match(df$feature, CellChatDB$interaction$interaction_name)
        if (sum(is.na(idx)) > 0) {
          stop(paste0("Please check the feature names or the input database! \n",
                      "The feature names are ", df$feature, ", and the names in the database are ", CellChatDB$interaction$interaction_name[idx[!is.na(idx)]], "! \n"))
        }
        df$feature <- CellChatDB$interaction$interaction_name_2[idx]

        data.plot <- df
      }

      data.plot$positive <- if_else(data.plot[,"score",drop=T]>0,"pos","neg")
      data.plot$absolute.v <- abs(data.plot[,"score",drop=T])

      data.plot <- data.plot %>%
        arrange(desc(absolute.v)) %>%
        mutate(rank=(seq_len(n())))

      data.topFeature <- data.plot %>%
        head(nTop)

      plot <- ggplot(data.plot, aes(rank, absolute.v)) +
        geom_line(colour = "grey80",linewidth = line_size) +
        theme(text = element_text(size = 10), axis.text.x = element_blank(),axis.ticks.x = element_blank()) +
        xlab("Ranking") + #ylab("Weight")+
        ylab(ylabel)+
        geom_point(data = data.topFeature,mapping = aes(color=positive),size = point_size, shape = 1)+
        scale_color_manual(values = col)+
        guides(color="none")+

        ggrepel::geom_text_repel(aes(label = feature),
                                 data = data.topFeature,
                                 segment.color = "grey50",
                                 segment.alpha = 1,
                                 direction = "y", max.overlaps=nrow(data.topFeature),
                                 nudge_x = text_nudge_x,
                                 hjust = 1,
                                 size = text_size,
                                 segment.size = 0.3)+
        theme_classic()+
        theme(
          panel.grid = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()
        ) +
        theme(strip.background = element_rect(fill="grey80")) +
        ggtitle(paste0(topic.key,i))+theme(plot.title = element_text(hjust = 0.5, vjust = 0, size = 10))

      return(plot)
    })
  if (combine) {
    plots <- patchwork::wrap_plots(plots, ncol = ncol)
  }
  return(plots)
}


#' @title vizDimLoadings
#'
#' @param object Seurat object
#' @param dims an integer vector specifying which dimensions to analyze
#' @param nfeatures an integer specifying the number of top features to extract from each dimension
#' @param nTop an integer specifying the number of top-ranked features to show with a text annotation
#' @param feature_cat a character string specifying the feature category to summarize loadings, one of \code{"LR-pair"}, \code{"ligand"}, or \code{"signaling"}
#' @param col numeric of length 2. Colors of points to use, default is \code{c('neg'="#377EB8", 'pos'="#E41A1C")}
#' @param reduction a character string specifying the dimensional reduction method stored in the Seurat object
#' @param projected whether to use the projected feature loadings
#' @param balanced whether to balance positive and negative features
#' @param species a character string specifying the species used to select the default CellChat database, either \code{"mouse"} or \code{"human"}
#' @param db an optional interaction database
#' @param ncol number of columns to show in the plot
#' @param combine whether outputting an combined plot using patchwork
#' @param ylabel name of y label
#' @param text_nudge_x horizontal and vertical adjustments to nudge the starting position of each text label. The units for
#' nudge_x and nudge_y are the same as for the data units on the x-axis and y-axis
#' @param text_size fontsize of text annotation
#' @param point_size numeric. Size of points
#' @param line_size numeric. Width of lines
#'
#' @return ggplot2 object
#' @import Seurat
#' @import ggplot2
#' @export
vizDimLoadings <- function (
    object,
    dims = 1:5,
    nfeatures = 30,
    nTop = 10,
    feature_cat = c("LR-pair","ligand", "signaling"),
    col = c('neg'="#377EB8", 'pos'="#E41A1C"),
    reduction = "pca",
    projected = FALSE,
    balanced = FALSE,
    species = c('mouse','human'), db = NULL,
    ncol = NULL,
    combine = TRUE,
    ylabel = "Signaling score",
    text_nudge_x = -8,
    text_size = 3,
    point_size = 2,
    line_size = 1
){
  feature_cat <- match.arg(feature_cat)
  species <- match.arg(species)
  if (is.null(x = ncol)) {
    ncol <- 2
    if (length(x = dims) == 1) {
      ncol <- 1
    }
    if (length(x = dims) > 6) {
      ncol <- 3
    }
    if (length(x = dims) > 9) {
      ncol <- 4
    }
  }
  # loadings: nFeatures x nPCs
  loadings <- Loadings(object = object[[reduction]], projected = projected)
  # get the features of each PC
  features <- lapply(X = dims, FUN = TopFeatures, object = object[[reduction]],
                     nfeatures = nfeatures, projected = projected, balanced = balanced)
  features <- lapply(X = features, FUN = unlist, use.names = FALSE)

  # subset loadings: mFeatures x mPCs
  loadings <- loadings[unlist(x = features), dims, drop = FALSE]
  names(features) <- colnames(loadings) <- as.character(dims)
  plots <- lapply(X = as.character(dims), FUN = function(i) {
    # get loading of PC i
    df <- as.data.frame(loadings[features[[i]],i, drop = FALSE])
    #colnames(data.plot) <- paste0(Key(object = object[[reduction]]),i)
    colnames(df) <- "score"
    rownames(df) <- gsub("-", "_", rownames(df))
    df$feature <- rownames(df)

    if (feature_cat %in% c("ligand","signaling") ) {
      if (is.null(db)) {
        if (species == "mouse") {
          CellChatDB <- CellChatDB.mouse
        } else if (species == 'human') {
          CellChatDB <- CellChatDB.human
        } else {
          stop("Only mouse and human are supported currently. Please provide a `db` instead! ")
        }
      } else {
        CellChatDB <- db
      }
      db.use <- CellChatDB$interaction[df$feature, c("ligand",'receptor','pathway_name')]
      df <- cbind(df, db.use)
      # summarize the loading score
      if (feature_cat == "ligand") {
        data.plot <- df %>% group_by(ligand) %>% summarize(avg = mean(score))  # total = sum(score),
      } else {
        data.plot <- df %>% group_by(pathway_name) %>% summarize(avg = mean(score))  # avg = mean(score),
      }
      colnames(data.plot) <- c("feature", "score")
    } else {
      data.plot <- df
    }

    data.plot$positive <- if_else(data.plot[,"score",drop=T]>0,"pos","neg")
    data.plot$absolute.v <- abs(data.plot[,"score",drop=T])

    data.plot <- data.plot %>%
      arrange(desc(absolute.v)) %>%
      mutate(rank=(seq_len(n())))

    data.topFeature <- data.plot %>%
      head(nTop)

    plot <- ggplot(data.plot, aes(rank, absolute.v)) +
      geom_line(colour = "grey80",linewidth = line_size) +
      theme(text = element_text(size = 10), axis.text.x = element_blank(),axis.ticks.x = element_blank()) +
      xlab("Ranking") + #ylab("Weight")+
      ylab(ylabel)+
      geom_point(data = data.topFeature,mapping = aes(color=positive),size = point_size, shape = 1)+
      scale_color_manual(values = col)+
      guides(color="none")+

      ggrepel::geom_text_repel(aes(label = feature),
                               data = data.topFeature,
                               segment.color = "grey50",
                               segment.alpha = 1,
                               direction = "y", max.overlaps=nrow(data.topFeature),
                               nudge_x = text_nudge_x,
                               hjust = 1,
                               size = text_size,
                               segment.size = 0.3)+
      theme_classic()+
      theme(
        panel.grid = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
      ) +
      theme(strip.background = element_rect(fill="grey80")) +
      ggtitle(paste0(Key(object = object[[reduction]]),i))+theme(plot.title = element_text(hjust = 0.5, vjust = 0, size = 10))

    return(plot)
  })
  if (combine) {
    plots <- patchwork::wrap_plots(plots, ncol = ncol)
  }
  return(plots)
}


#' @title netVisual_CommunField
#' @param object SpatialCellChat object
#' @param signaling a signaling pathway or ligand-receptor pair to visualize
#' @param slot.name the slot name of object. Set is to be "netP" if input signaling is a pathway name; Set is to be "net" if input signaling is a ligand-receptor pair
#' @param pattern "outgoing" or "incoming"
#' @param zoom.to an integer specifying the zoom level for spatial visualization.
#'   If \code{0}, no zooming is applied and the full field is shown.
#'   If greater than 0, the spatial area is divided into \code{zoom.to + 1}
#'   equal intervals along each axis
#' @param zoom.mode a logical value indicating how zooming is handled.
#'   If \code{FALSE}, only the selected zoomed region is shown.
#'   If \code{TRUE}, the selected region is highlighted while the rest of the
#'   field remains visible with reduced transparency
#' @param anchor.x a numeric value between 0 and 1 specifying the relative horizontal position of the zoom anchor, used to select the zoomed region along the x-axis
#' @param anchor.y a numeric value between 0 and 1 specifying the relative vertical position of the zoom anchor, used to select the zoomed region along the y-axis
#' @param color.use defining the color for each cell group
#' @param point.size the size of spots
#' @param min.mag a numeric value specifying the minimum vector magnitude required for arrows to be displayed
#' @param alpha.image the transparency of individual spots
#' @param arrow.angle the angle of arrows
#' @param arrow.length the length of arrows
#' @param arrow.size the size of arrows
#' @param arrow.color the color of arrows
#' @param legend.size the size of legend
#' @param legend.text.size the text size on the legend
#' @param plot.title.name title name of the plot
#' @param plot.title.text.size the text size of the title
#' @param ... other parameters passing to metR::geom_arrow
#' @import ggplot2
#'
#' @return ggplot2 object
#' @export
netVisual_CommunField <- function(
    object,
    signaling,
    slot.name = "netP",
    pattern = c("incoming","outgoing"),
    zoom.to = 0L,
    zoom.mode = T,
    anchor.x = 0.5,
    anchor.y = 0.5,
    color.use = NULL,
    point.size = 2,
    min.mag=0,
    alpha.image = 0.2,
    arrow.angle = 15,
    arrow.length = 0.6,
    arrow.size = 1,
    arrow.color = NULL,
    ...,
    legend.size = 0.5,
    legend.text.size = 8,
    plot.title.name = NULL,
    plot.title.text.size = 8
){
  pattern <- match.arg(pattern)
  if (pattern == "outgoing") {
    arrow.ends = "last"
    legend.title <- "Sources"
    if (is.null(plot.title.name)) {
      plot.title.name <- paste0("Outgoing communication field of ", signaling)
    }
  } else {
    arrow.ends = "last"
    legend.title <- "Targets"
    if (is.null(plot.title.name)) {
      plot.title.name <- paste0("Incoming communication field of ", signaling)
    }
  }

  net0 <- methods::slot(object, slot.name)
  field <- my_as_sparse3Darray(net0$field[[pattern]])
  rm(net0)
  gc()

  siganling_name = dimnames(field)[[3]]
  if (!(signaling %in% siganling_name)) {
    stop(
      "Please check the input 'signaling' and make sure it has been computed via `computeCommunField`."
    )
  }

  # x->y;y->x
  field.use = data.frame(field[, , signaling])
  temp_field.use = field.use
  field.use[, 1] = temp_field.use[, 2]
  field.use[, 2] = temp_field.use[, 1]
  coordinates <- object@images$coordinates
  temp_coordinates = coordinates
  coordinates[, 1] = temp_coordinates[, 2]
  coordinates[, 2] = temp_coordinates[, 1]
  colnames(field.use) = c('dx', 'dy')
  field.mag = sqrt(field.use$dx ^ 2 + field.use$dy ^ 2)
  colnames(coordinates) <- c("x_cent", "y_cent")

  labels <- object@idents
  cells.level <- levels(labels)
  df.field <-
    data.frame(coordinates, field.use, mag = field.mag, labels)

  if (is.integer(zoom.to) & (zoom.to >= 0)) {
    if (zoom.to >= 1) {
      n.breaks <- zoom.to + 1
      x_range <- range(df.field[["x_cent"]])
      x_breaks <- seq(x_range[1], x_range[2], length.out = n.breaks + 1)
      y_range <- range(df.field[["y_cent"]])
      y_breaks <- seq(y_range[1], y_range[2], length.out = n.breaks + 1)
      df.field[["x_cent_cut"]] <- cut(df.field[["x_cent"]], breaks = x_breaks)
      levels_x_cent <- levels(df.field[["x_cent_cut"]])
      df.field[["y_cent_cut"]] <- cut(df.field[["y_cent"]], breaks = y_breaks)
      levels_y_cent <- levels(df.field[["y_cent_cut"]])

      anchor.Interval <- seq(0, 1, len = (n.breaks + 1))
      anchor.Interval <- anchor.Interval[-1]
      anchor.Interval <- anchor.Interval[-length(anchor.Interval)]

      anchor.x.Interval <- findInterval(anchor.x, anchor.Interval)
      anchor.y.Interval <- findInterval(anchor.y, anchor.Interval)

      if (!zoom.mode) {
        df.field <-
          dplyr::filter(df.field,
                        (x_cent_cut == levels_x_cent[[anchor.x.Interval + 1L]]) &
                          (y_cent_cut == levels_y_cent[[anchor.y.Interval + 1L]]))
      } else{
        # do highlight
        df.field <-
          dplyr::mutate(df.field,
                        alpha_image = dplyr::if_else(
                          (x_cent_cut == levels_x_cent[[anchor.x.Interval + 1L]]) &
                            (y_cent_cut == levels_y_cent[[anchor.y.Interval + 1L]]),
                          true = min((alpha.image + 0.1), 1),
                          false = alpha.image
                        ))
      }
    } else {
      # zoom to 0x, do nothing
    }
  } else {
    stop(cli.symbol("fail"), "Please check your `zoom.to`!")
  }

  idx.remove <- which(df.field$mag == 0)
  if(length(idx.remove)>0){
    df.field$dx[idx.remove] <- NA
    df.field$dy[idx.remove] <- NA
    df.field.use <- df.field[-idx.remove,]
  } else { #length(idx.remove)=0
    df.field.use <- df.field
  }

  if (is.null(color.use)) {
    color.use = scPalette(length(cells.level))
    names(color.use) <- cells.level
  }

  if (length(setdiff(levels(df.field.use$labels), unique(df.field.use$labels))) > 0) {
    color.use.field <-
      color.use[names(color.use) %in% unique(df.field.use$labels)]
    df.field.use$labels <-
      droplevels(df.field.use$labels, exclude = setdiff(
        levels(df.field.use$labels),
        unique(df.field.use$labels)
      ))
  } else {
    color.use.field <- color.use
  }

  # Add points based on zoom.to and zoom.mode
  if (zoom.to == 0L) {
    gg <- ggplot(data = df.field, aes(x_cent, y_cent)) +
      geom_point(
        mapping = aes(colour = labels),
        size = point.size,
        show.legend = TRUE,
        alpha = alpha.image
      )
  } else {
    if (zoom.mode) {
      gg <- ggplot(data = df.field, aes(x_cent, y_cent)) +
        geom_point(
          mapping = aes(colour = labels,alpha = alpha_image),
          size = point.size,
          show.legend = FALSE
        )
    } else {
      gg <- ggplot(data = df.field, aes(x_cent, y_cent)) +
        geom_point(
          mapping = aes(colour = labels),
          size = point.size,
          show.legend = TRUE,
          alpha = alpha.image
        )
    }
  }

  # Add arrows
  if (is.null(arrow.color)) {
    gg <-
      gg + metR::geom_arrow(
        data = df.field.use,
        mapping = aes(
          dx = dx,
          dy = -dy,
          colour = labels
          # alpha = alpha_image
        ),
        min.mag = min.mag,
        pivot = 0,
        ...,
        arrow = grid::arrow(
          arrow.angle,
          grid::unit(arrow.length, "lines"),
          ends = arrow.ends,
          type = "closed"
        )
      ) +  metR::scale_mag(guide = "none")
  } else {
    gg <-
      gg + metR::geom_arrow(
        data = df.field.use,
        mapping = aes(
          dx = dx,
          dy = -dy
          # alpha = alpha_image
        ),
        min.mag = min.mag,
        pivot = 0,
        colour = arrow.color[1],
        ...,
        arrow = grid::arrow(
          arrow.angle,
          grid::unit(arrow.length, "lines"),
          ends = arrow.ends,
          type = "closed"
        )
      ) + metR::scale_mag(guide = "none")
  }

  # Finalize the plot
  if(isTRUE(zoom.mode) & (zoom.to > 0)){
    gg <- gg + geom_vline(xintercept = x_breaks[2:n.breaks], color = "gray", linetype = "dashed")
    gg <- gg + geom_hline(yintercept = y_breaks[2:n.breaks], color = "gray", linetype = "dashed")

    x_breaks_prop <- (x_breaks - x_range[1]) / (x_range[2] - x_range[1])
    y_breaks_prop <- (y_breaks - y_range[1]) / (y_range[2] - y_range[1])

    df_vline_text <-
      data.frame(x = x_breaks[2:n.breaks],
                 y = y_range[1],
                 label = round(x_breaks_prop[2:n.breaks], 2))
    gg <- gg + geom_text(
      data = df_vline_text,
      aes(x = x, y = y, label = label),
      vjust = 1.5, hjust = 0.5, size = 3
    )

    df_vline_text <-
      data.frame(x = x_range[1],
                 y = y_breaks[2:n.breaks],
                 label = round(y_breaks_prop[2:n.breaks], 2))
    gg <- gg + geom_text(
      data = df_vline_text,
      aes(x = x, y = y, label = label),
      hjust = -0.1, vjust = 0.5, size = 3
    )
  }

  gg <- gg + scale_color_manual(values = color.use.field) +
    theme(
      legend.key = element_blank(),
      legend.text = element_text(size = legend.text.size),
      legend.title = element_text(size = legend.text.size)
    ) +
    theme(
      panel.background = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank()
    ) + xlab(NULL) + ylab(NULL)
  gg <- gg + coord_fixed() + scale_y_reverse()
  # gg <- gg + coord_fixed(ratio = 1*diff(range(coordinates$x_cent))/diff(range(coordinates$y_cent))) + scale_y_reverse()
  gg <-
    gg + guides(color = guide_legend(
      title = legend.title,
      override.aes = list(size = legend.size)
    ))
  gg <-
    gg + ggtitle(plot.title.name) +
    theme(plot.title = element_text(
      hjust = 0.5,
      vjust = 0,
      size = plot.title.text.size
    ))

  return(gg)
}


#' @title netVisual_CommunFieldGrid
#' @param object SpatialCellChat object
#' @param signaling a signaling pathway or ligand-receptor pair to visualize
#' @param slot.name the slot name of object. Set is to be "netP" if input signaling is a pathway name; Set is to be "net" if input signaling is a ligand-receptor pair
#' @param pattern "outgoing" or "incoming"
#' @param cellsize an integer specifying the spatial resolution (in pixels) used to construct the grid. If \code{NULL}, the minimum pairwise distance between spatial coordinates is used as the default cell size.
#' @param what a character string specifying the type of geometry returned by \code{sf::st_make_grid}. Default is "polygons"
#' @param square a logical value indicating whether the grid cells should be square
#' @param grid.resolution an integer specifying the scaling factor applied to \code{cellsize} when constructing the grid
#' @param grid.field.mean a character string specifying how vector fields within each grid cell are aggregated. Options include \code{"median"} or \code{"sum"}
#' @param color.use defining the color for each cell group
#' @param point.size the size of spots
#' @param image.alpha the transparency of individual spots
#' @param min.mag a numeric value specifying the minimum vector magnitude required for arrows to be displayed
#' @param arrow.line.color the color of arrows
#' @param arrow.line.alpha the transparency of arrows
#' @param arrow.angle the angle of arrows
#' @param arrow.length the length of arrows
#' @param legend.size the size of legend
#' @param legend.text.size the text size on the legend
#' @param title.name title name of the plot
#' @param ... other parameters passing to metR::geom_arrow
#' @import ggplot2
#'
#' @return ggplot2 object
#' @export
netVisual_CommunFieldGrid <- function(
    object,
    signaling,
    slot.name = "netP",
    pattern = c("incoming","outgoing"),
    cellsize=NULL,
    what = "polygons",
    square = T,
    grid.resolution=NULL,
    grid.field.mean = c("median","sum"),
    color.use = NULL,
    point.size = 2,
    image.alpha = 0.3,
    min.mag = 0,
    arrow.line.color = "black",
    # arrow.line.width = 0.7,
    arrow.line.alpha = 0.9,
    arrow.angle = 22.5,
    arrow.length = 0.5,
    ...,
    legend.size = 2,
    legend.text.size = 8,
    title.name = NULL
){
  pattern <- match.arg(pattern)
  grid.field.mean <- match.arg(grid.field.mean)

  if (pattern == "outgoing") {
    arrow.ends = "last"
    legend.title <- "Sources"
    if (is.null(title.name)) {
      title.name <- paste0("Outgoing communication flow of ", signaling)
    }
  } else if (pattern == "incoming") {
    arrow.ends = "last"
    legend.title <- "Targets"
    if (is.null(title.name)) {
      title.name <- paste0("Incoming communication flow of ", signaling)
    }
  }

  field <- my_as_sparse3Darray(methods::slot(object, slot.name)$field[[pattern]])

  siganling_name = dimnames(field)[[3]]
  if (!(signaling %in% siganling_name)) {
    stop("Please check the input 'signaling' and make sure it has been computed via `computeCommunField`.")
  }

  field.use = data.frame(field[, , signaling,drop=T]) # return a matrix
  temp_field.use = field.use
  field.use[, 1] = temp_field.use[, 2]
  field.use[, 2] = temp_field.use[, 1]

  coordinates <- object@images$coordinates
  temp_coordinates = coordinates
  coordinates[, 1] = temp_coordinates[, 2]
  coordinates[, 2] = temp_coordinates[, 1]
  colnames(field.use) = c('dx', 'dy')
  field.mag = sqrt(field.use$dx ^ 2 + field.use$dy ^ 2)
  colnames(coordinates) <- c("x_cent", "y_cent")

  labels <- object@idents
  cells.level <- levels(labels)
  df.field <- data.frame(coordinates, field.use, mag = field.mag, labels)

  # some hints about contact.range and spot size in the new grid SpatialCellChat
  spatial.factors <- object@images$spatial.factors

  if(( !is.null(cellsize) ) & is.integer(cellsize)){
    spot.size <- cellsize*spatial.factors[["ratio"]] %>% round(digits = 4)
    cat(cli.symbol(3),"The `cellsize` you input is ",cellsize," units(pixels). It is about ",spot.size," um in International System of Units.\n")
  }
  if(is.null(cellsize)){
    cell2cellDist <- Rfast::Dist(coordinates)
    diag(cell2cellDist) <- NA
    cellsize <- min(cell2cellDist,na.rm = T)
    spot.size <- cellsize*spatial.factors[["ratio"]] %>% round(digits = 4)
    cat(cli.symbol(3),"The default `cellsize` in the SpatialCellChat object is ",cellsize," units(pixels). It is about ",spot.size," um in International System of Units.\n")
  }

  if(is.null(grid.resolution)) {grid.resolution <- 2}
  newcellsize <- cellsize*grid.resolution
  spot.size <- newcellsize*spatial.factors[["ratio"]] %>% round(digits = 4)
  cat(cli.symbol(),"Do grid to show communication flow, the new `cellsize` will be ",newcellsize," units(pixels). It is about ",spot.size," um in International System of Units.\n")


  df.field <-
    sf::st_as_sf(df.field,
                 coords = c("x_cent", "y_cent"),
                 remove = FALSE)
  sf::st_crs(df.field) <- 3857

  square_grid <-
    sf::st_make_grid(df.field,
                     cellsize = newcellsize,
                     what = what,
                     square = square)

  ## convert `square_grid` into sf and add grid IDs
  square_grid_sf = sf::st_sf(square_grid) %>%
    # add grid ID
    dplyr::mutate(grid_id = seq_len(length(lengths(square_grid))))

  if(grid.field.mean=="median"){
    field.mean <- median
  } else if (grid.field.mean=="sum"){
    field.mean <- sum
  }

  field.gridded = square_grid_sf %>%
    mutate(id = 1:n(),
           contained = lapply(sf::st_intersects(square_grid,df.field),identity),
           obs = sapply(contained, length),
           u = sapply(contained, function(x) {field.mean( df.field[x,,drop=F]$dx, na.rm = TRUE )}),
           v = sapply(contained, function(x) {field.mean( df.field[x,,drop=F]$dy, na.rm = TRUE )}))

  # field.gridded = field.gridded %>% select(obs, u, v) %>% na.omit()

  ## obtain the centroid coordinates from the grid as table
  field.coordinates = field.gridded %>%
    sf::st_centroid() %>%
    sf::st_coordinates() %>%
    dplyr::as_tibble() %>%
    rename(x = X, y = Y)
  # sf::st_geometry(field.gridded) = NULL
  current.gridded = field.coordinates %>%
    bind_cols(field.gridded) #%>% na.omit()
  current.gridded.na <- current.gridded[current.gridded$obs==0,,drop=F]
  current.gridded <- current.gridded %>% na.omit()

  ## interpolate the U component
  u.se = oce::interpBarnes(x = current.gridded$x, y = current.gridded$y, z = current.gridded$u)

  ## obtain dimension that determine the width (ncol) and length (nrow) for tranforming wide into long format table
  dimension = data.frame(lon = u.se$xg, u.se$zg) %>% dim()

  ## make a U component data table from interpolated matrix
  u.tb = data.frame(lon = u.se$xg, u.se$zg) %>%
    tidyr::gather(key = "lata", value = "u", 2:dimension[2]) %>%
    dplyr::mutate(lat = rep(u.se$yg, each = dimension[1])) %>%
    dplyr::select(lon,lat, u) %>%
    dplyr::as_tibble()

  ## interpolate the V component
  v.se = oce::interpBarnes(x = current.gridded$x, y = current.gridded$y, z = current.gridded$v)

  ## make the V component data table from interpolated matrix
  v.tb = data.frame(lon = v.se$xg, v.se$zg) %>%
    tidyr::gather(key = "lata", value = "v", 2:dimension[2]) %>%
    dplyr::mutate(lat = rep(v.se$yg, each = dimension[1])) %>%
    dplyr::select(lon,lat, v) %>%
    dplyr::as_tibble()

  ## stitch now the V component intot the U data table and compute the velocity
  uv.se = u.tb %>%
    bind_cols(v.tb %>% select(v)) %>%
    mutate(vel = sqrt(u^2+v^2))

  # omit uvs which are not in non-zero grid
  uv.se <-
    sf::st_as_sf(uv.se,
                 coords = c("lon", "lat"),
                 remove = FALSE)
  sf::st_crs(uv.se) <- 3857

  uv.se.use <- sf::st_intersects(uv.se$geometry,current.gridded$square_grid,sparse = T)
  uv.se.use <- purrr::map_dbl(.x = uv.se.use,.f = function(x){return(length(x))})
  uv.se.omit <- which(uv.se.use==0)

  uv.se[uv.se.omit,c("u","v","vel")] <- 0

  if (is.null(color.use)) {
    color.use = scPalette(length(cells.level))
    names(color.use) <- cells.level
  }

  scChat_theme_opts <- theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.grid.major = element_line(colour = "transparent"),
    panel.grid.minor = element_line(colour = "transparent"),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank()
  ) + theme(legend.key = element_blank())

  gg <- ggplot() +
    # metR::geom_contour_fill(data = uv.se, aes(x = lon, y = lat, z = vel),
    # na.fill = TRUE, bins = 70) +
    # !y=-lat;dy=-v  =>  scale_y_reverse()!
    geom_point(
      data = df.field,
      mapping = aes(x = x_cent, y = -y_cent, color = labels),
      size = point.size,
      alpha = image.alpha
    ) +
    scale_color_manual(values = color.use)+
    metR::geom_arrow(
      mapping = aes(
        x = lon,
        y = -lat,
        dx = u,
        dy = -v
      ),
      data = uv.se,
      color = arrow.line.color,
      min.mag = min.mag,
      # linewidth = arrow.line.width,
      alpha = arrow.line.alpha,
      ... = ...,
      arrow.angle = arrow.angle,
      arrow.length = arrow.length,
      arrow.ends = arrow.ends
    )+
    metR::scale_mag(guide = "none")+
    # geom_sf(data =current.gridded.na,mapping = aes(geometry=square_grid),color="white",fill="white",linewidth=0)+
    # scale_fill_gradientn(name = "Current",colours = oceColorsVelocity(120),
    #                      limits = c(0,1.6), breaks = seq(0.1,1.6,.3))+
    scChat_theme_opts+#coord_fixed()+  #scale_y_reverse()+
    theme(legend.position = "right") + xlab(NULL) + ylab(NULL) +
    coord_fixed()
  # coord_fixed(ratio = 1*diff(range(coordinates$x_cent))/diff(range(coordinates$y_cent)))


  gg <-
    gg + guides(color = guide_legend(
      title = legend.title,
      override.aes = list(size = legend.size)
    ))
  gg <-
    gg + ggtitle(title.name) + theme(plot.title = element_text(
      hjust = 0.5,
      vjust = 0,
      size = 10
    ))

  return(gg)
}


#' @title netVisual_CommunFlow
#' @param object SpatialCellChat object
#' @param signaling a signaling pathway or ligand-receptor pair to visualize
#' @param slot.name the slot name of object. Set is to be "netP" if input signaling is a pathway name; Set is to be "net" if input signaling is a ligand-receptor pair
#' @param pattern "outgoing" or "incoming"
#' @param cellsize an integer specifying the spatial resolution (in pixels) used to construct the grid. If \code{NULL}, the minimum pairwise distance between spatial coordinates is used as the default cell size.
#' @param what a character string specifying the type of geometry returned by \code{sf::st_make_grid}. Default is "polygons"
#' @param square a logical value indicating whether the grid cells should be square
#' @param grid.resolution an integer specifying the scaling factor applied to \code{cellsize} when constructing the grid
#' @param grid.field.mean a character string specifying how vector fields within each grid cell are aggregated. Options include \code{"median"} or \code{"sum"}
#' @param color.use defining the color for each cell group
#' @param point.size the size of spots
#' @param image.alpha the transparency of individual spots
#' @param L typical length of a streamline
#' @param min.L minimum length of segments to show
#' @param res resolution parameter. Higher numbers increases the resolution used for streamline integration
#' @param n numeric indicating the number of points to draw
#' @param jitter amount of jitter of the starting points
#' @param arrow.line.color the color of arrows
#' @param arrow.line.width the width of arrows
#' @param arrow.line.alpha the transparency of arrows
#' @param arrow.angle the angle of arrows
#' @param arrow.length the length of arrows
#' @param legend.size the size of legend
#' @param legend.text.size the text size on the legend
#' @param title.name title name of the plot
#' @param ... other parameters passing to metR::geom_streamline
#' @import ggplot2
#'
#' @return ggplot2 object
#' @export
netVisual_CommunFlow <- function(
    object,
    signaling,
    slot.name = "netP",
    pattern = c("incoming","outgoing"),
    cellsize=NULL,
    what = "polygons",
    square = T,
    grid.resolution=NULL,
    grid.field.mean = c("median","sum"),
    color.use = NULL,
    point.size = 2,
    image.alpha = 0.9,
    L = 300,
    min.L = 0,
    res = 0.5,
    n = 30,
    jitter=3,
    arrow.line.color="grey40",
    arrow.line.width=0.7,
    arrow.line.alpha=0.6,
    arrow.angle = 12,
    arrow.length = 0.3,
    ...,
    legend.size = 2,
    legend.text.size = 8,
    title.name = NULL
){
  pattern <- match.arg(pattern)
  grid.field.mean <- match.arg(grid.field.mean)


  if (pattern == "outgoing") {
    arrow.ends = "last"
    legend.title <- "Sources"
    if (is.null(title.name)) {
      title.name <- paste0("Outgoing communication flow of ", signaling)
    }
  } else if (pattern == "incoming") {
    arrow.ends = "last"
    legend.title <- "Targets"
    if (is.null(title.name)) {
      title.name <- paste0("Incoming communication flow of ", signaling)
    }
  }

  field <- my_as_sparse3Darray(methods::slot(object, slot.name)$field[[pattern]])

  siganling_name = dimnames(field)[[3]]
  if (!(signaling %in% siganling_name)) {
    stop("Please check the input 'signaling' and make sure it has been computed via `computeCommunField`.")
  }

  field.use = data.frame(field[, , signaling,drop=T]) # return a matrix
  temp_field.use = field.use
  field.use[, 1] = temp_field.use[, 2]
  field.use[, 2] = temp_field.use[, 1]

  coordinates <- object@images$coordinates
  temp_coordinates = coordinates
  coordinates[, 1] = temp_coordinates[, 2]
  coordinates[, 2] = temp_coordinates[, 1]
  colnames(field.use) = c('dx', 'dy')
  field.mag = sqrt(field.use$dx ^ 2 + field.use$dy ^ 2)
  colnames(coordinates) <- c("x_cent", "y_cent")

  labels <- object@idents
  cells.level <- levels(labels)
  df.field <- data.frame(coordinates, field.use, mag = field.mag, labels)

  # some hints about contact.range and spot size in the new grid SpatialCellChat
  spatial.factors <- object@images$spatial.factors

  if(( !is.null(cellsize) ) & is.integer(cellsize)){
    spot.size <- cellsize*spatial.factors[["ratio"]] %>% round(digits = 4)
    cat(cli.symbol(3),"The `cellsize` you input is ",cellsize," units(pixels). It is about ",spot.size," um in International System of Units.\n")
  }
  if(is.null(cellsize)){
    cell2cellDist <- Rfast::Dist(coordinates)
    diag(cell2cellDist) <- NA
    cellsize <- min(cell2cellDist,na.rm = T)
    spot.size <- cellsize*spatial.factors[["ratio"]] %>% round(digits = 4)
    cat(cli.symbol(3),"The default `cellsize` in the SpatialCellChat object is ",cellsize," units(pixels). It is about ",spot.size," um in International System of Units.\n")
  }

  if(is.null(grid.resolution)) {grid.resolution <- 2}
  newcellsize <- cellsize*grid.resolution
  spot.size <- newcellsize*spatial.factors[["ratio"]] %>% round(digits = 4)
  cat(cli.symbol(),"Do grid to show communication flow, the new `cellsize` will be ",newcellsize," units(pixels). It is about ",spot.size," um in International System of Units.\n")


  df.field <-
    sf::st_as_sf(df.field,
                 coords = c("x_cent", "y_cent"),
                 remove = FALSE)
  sf::st_crs(df.field) <- 3857

  square_grid <-
    sf::st_make_grid(df.field,
                     cellsize = newcellsize,
                     what = what,
                     square = square)

  ## convert `square_grid` into sf and add grid IDs
  square_grid_sf = sf::st_sf(square_grid) %>%
    # add grid ID
    dplyr::mutate(grid_id = seq_len(length(lengths(square_grid))))

  if(grid.field.mean=="median"){
    field.mean <- median
  } else if (grid.field.mean=="sum"){
    field.mean <- sum
  }

  field.gridded = square_grid_sf %>%
    mutate(id = 1:n(),
           contained = lapply(sf::st_intersects(square_grid,df.field),identity),
           obs = sapply(contained, length),
           u = sapply(contained, function(x) {field.mean( df.field[x,,drop=F]$dx, na.rm = TRUE )}),
           v = sapply(contained, function(x) {field.mean( df.field[x,,drop=F]$dy, na.rm = TRUE )}))

  # field.gridded = field.gridded %>% select(obs, u, v) %>% na.omit()

  ## obtain the centroid coordinates from the grid as table
  field.coordinates = field.gridded %>%
    sf::st_centroid() %>%
    sf::st_coordinates() %>%
    dplyr::as_tibble() %>%
    rename(x = X, y = Y)
  # sf::st_geometry(field.gridded) = NULL
  current.gridded = field.coordinates %>%
    bind_cols(field.gridded) #%>% na.omit()
  current.gridded.na <- current.gridded[current.gridded$obs==0,,drop=F]
  current.gridded <- current.gridded %>% na.omit()

  ## interpolate the U component
  u.se = oce::interpBarnes(x = current.gridded$x, y = current.gridded$y, z = current.gridded$u)

  ## obtain dimension that determine the width (ncol) and length (nrow) for tranforming wide into long format table
  dimension = data.frame(lon = u.se$xg, u.se$zg) %>% dim()

  ## make a U component data table from interpolated matrix
  u.tb = data.frame(lon = u.se$xg, u.se$zg) %>%
    tidyr::gather(key = "lata", value = "u", 2:dimension[2]) %>%
    dplyr::mutate(lat = rep(u.se$yg, each = dimension[1])) %>%
    dplyr::select(lon,lat, u) %>%
    dplyr::as_tibble()

  ## interpolate the V component
  v.se = oce::interpBarnes(x = current.gridded$x, y = current.gridded$y, z = current.gridded$v)

  ## make the V component data table from interpolated matrix
  v.tb = data.frame(lon = v.se$xg, v.se$zg) %>%
    tidyr::gather(key = "lata", value = "v", 2:dimension[2]) %>%
    dplyr::mutate(lat = rep(v.se$yg, each = dimension[1])) %>%
    dplyr::select(lon,lat, v) %>%
    dplyr::as_tibble()

  ## stitch now the V component intot the U data table and compute the velocity
  uv.se = u.tb %>%
    bind_cols(v.tb %>% select(v)) %>%
    mutate(vel = sqrt(u^2+v^2))

  # omit uvs which are not in non-zero grid
  uv.se <-
    sf::st_as_sf(uv.se,
                 coords = c("lon", "lat"),
                 remove = FALSE)
  sf::st_crs(uv.se) <- 3857

  uv.se.use <- sf::st_intersects(uv.se$geometry,current.gridded$square_grid,sparse = T)
  uv.se.use <- purrr::map_dbl(.x = uv.se.use,.f = function(x){return(length(x))})
  uv.se.omit <- which(uv.se.use==0)

  uv.se[uv.se.omit,c("u","v","vel")] <- 0

  if (is.null(color.use)) {
    color.use = scPalette(length(cells.level))
    names(color.use) <- cells.level
  }

  scChat_theme_opts <- theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.grid.major = element_line(colour = "transparent"),
    panel.grid.minor = element_line(colour = "transparent"),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank()
  ) + theme(legend.key = element_blank())

  gg <- ggplot() +
    # metR::geom_contour_fill(data = uv.se, aes(x = lon, y = lat, z = vel),
    # na.fill = TRUE, bins = 70) +
    # !y=-lat;dy=-v  =>  scale_y_reverse()!
    geom_point(data = df.field,mapping = aes(x=x_cent,y=-y_cent,color=labels),size=point.size,alpha=image.alpha)+
    scale_color_manual(values = color.use)+
    # ggnewscale::new_scale_color()+
    metR::geom_streamline(
      mapping = aes(
        x = lon,
        y = -lat,
        dx = u,
        dy = -v
      ),
      data = uv.se,
      color = arrow.line.color,
      linewidth = arrow.line.width,
      alpha = arrow.line.alpha,
      ... = ...,
      L = L,
      min.L = min.L,
      res = res,
      n = n,
      jitter = jitter,
      arrow.angle = arrow.angle,
      arrow.length = arrow.length,
      arrow.ends = arrow.ends
    ) +
    # geom_sf(data =current.gridded.na,mapping = aes(geometry=square_grid),color="white",fill="white",linewidth=0)+
    # scale_fill_gradientn(name = "Current",colours = oceColorsVelocity(120),
    #                      limits = c(0,1.6), breaks = seq(0.1,1.6,.3))+
    scChat_theme_opts+#coord_fixed()+  #scale_y_reverse()+
    theme(legend.position = "right") + xlab(NULL) + ylab(NULL) +
    coord_fixed()
  # coord_fixed(ratio = 1*diff(range(coordinates$x_cent))/diff(range(coordinates$y_cent))) #+ scale_y_reverse()


  gg <-
    gg + guides(color = guide_legend(
      title = legend.title,
      override.aes = list(size = legend.size)
    ))
  gg <-
    gg + ggtitle(title.name) + theme(plot.title = element_text(
      hjust = 0.5,
      vjust = 0,
      size = 10
    ))

  return(gg)
}



#' @title netVisual_CommunFieldZoomIn
#' @param object SpatialCellChat object
#' @param zoom.pt.size a numeric value specifying the point size used in the zoomed-in view
#' @inheritParams netVisual_CommunField
#' @import ggplot2
#'
#' @return ggplot2 object
#' @export
netVisual_CommunFieldZoomIn <- function(
    object,
    signaling,
    slot.name = "netP",
    pattern = c("incoming", "outgoing"),
    color.use = NULL,
    point.size = 2,
    zoom.to = 1L,
    anchor.x = 0.5,
    anchor.y = 0.5,
    zoom.pt.size = point.size+zoom.to*0.5,
    alpha.image = 0.2,
    arrow.angle = 15,
    arrow.length = 0.6,
    arrow.size = 1,
    arrow.color = NULL,
    ...,
    legend.size = 0.5,
    legend.text.size = 8,
    plot.title.name = NULL,
    plot.title.text.size = 8
) {
  gg1 <-
    netVisual_CommunField(
      object = object,
      signaling = signaling,
      slot.name = slot.name,
      pattern = pattern,
      zoom.to = zoom.to,
      zoom.mode = T,
      anchor.x = anchor.x,
      anchor.y = anchor.y,
      color.use = color.use,
      point.size = point.size,
      alpha.image = alpha.image,
      arrow.angle = arrow.angle,
      arrow.length = arrow.length,
      arrow.size = arrow.size,
      arrow.color = arrow.color,
      ...,
      legend.size = legend.size,
      legend.text.size = legend.text.size,
      plot.title.name = plot.title.name,
      plot.title.text.size = plot.title.text.size
    )
  title.gg2 <- paste0(plot.title.name,
                      "\n",
                      stringr::str_remove_all(
                        string = stringr::str_c("zoom.to=", zoom.to, "x,anchor.x=", anchor.x, ",anchor.y=", anchor.y),
                        pattern = " "
                      ))

  gg2 <-
    netVisual_CommunField(
      object = object,
      signaling = signaling,
      slot.name = slot.name,
      pattern = pattern,
      zoom.to = zoom.to,
      zoom.mode = F,
      anchor.x = anchor.x,
      anchor.y = anchor.y,
      color.use = color.use,
      point.size = zoom.pt.size,
      alpha.image = alpha.image,
      arrow.angle = arrow.angle,
      arrow.length = arrow.length,
      arrow.size = arrow.size,
      arrow.color = arrow.color,
      ...,
      legend.size = legend.size,
      legend.text.size = legend.text.size,
      plot.title.name = title.gg2,
      plot.title.text.size = plot.title.text.size
    )
  gg <- cowplot::plot_grid(gg1, gg2, align = "h", nrow = 1)
  return(gg)
}


#' Show the description of CellChatDB databse
#'
#' @param CellChatDB CellChatDB databse
#' @param nrow the number of rows in the plot
#' @importFrom dplyr group_by summarise n %>%
#'
#' @return
#' @export
#'
showDatabaseCategory <- function(CellChatDB, nrow = 1) {
  interaction_input <- CellChatDB$interaction
  geneIfo <- CellChatDB$geneInfo
  df <- interaction_input %>% group_by(annotation) %>% summarise(value=n())
  df$group <- factor(df$annotation, levels = c("Secreted Signaling","ECM-Receptor","Cell-Cell Contact"))
  gg1 <- pieChart(df)
  binary <- (interaction_input$ligand %in% geneIfo$Symbol) & (interaction_input$receptor %in% geneIfo$Symbol)
  df <- data.frame(group = rep("Heterodimers", dim(interaction_input)[1]),stringsAsFactors = FALSE)
  df$group[binary] <- rep("Others",sum(binary),1)
  df <- df %>% group_by(group) %>% summarise(value=n())
  df$group <- factor(df$group, levels = c("Heterodimers","Others"))
  gg2 <- pieChart(df)

  kegg <- grepl("KEGG", interaction_input$evidence)
  df <- data.frame(group = rep("Literature", dim(interaction_input)[1]),stringsAsFactors = FALSE)
  df$group[kegg] <- rep("KEGG",sum(kegg),1)
  df <- df %>% group_by(group) %>% summarise(value=n())
  df$group <- factor(df$group, levels = c("KEGG","Literature"))
  gg3 <- pieChart(df)

  gg <- cowplot::plot_grid(gg1, gg2, gg3, nrow = nrow, align = "h", rel_widths = c(1, 1,1))
  return(gg)
}


#' Plot pie chart
#'
#' @param df a dataframe
#' @param label.size a character
#' @param color.use the name of the variable in CellChatDB interaction_input
#' @param title the title of plot
#' @import ggplot2
#' @importFrom scales percent
#' @importFrom dplyr arrange desc mutate
#' @importFrom ggrepel geom_text_repel
#' @return
#' @export
#'
pieChart <- function(df, label.size = 2.5, color.use = NULL, title = "") {
  df %>% arrange(dplyr::desc(value)) %>%
    mutate(prop = scales::percent(value/sum(value))) -> df

  gg <- ggplot(df, aes(x="", y=value, fill=forcats::fct_inorder(group))) +
    geom_bar(stat="identity", width=1) +
    coord_polar("y", start=0)+theme_void() +
    ggrepel::geom_text_repel(aes(label = prop), size= label.size, show.legend = F, position = position_stack(vjust=0.5))
  #  ggrepel::geom_text_repel(aes(label = prop), size= label.size, show.legend = F, nudge_x = 0)
  gg <- gg + theme(legend.position="bottom", legend.direction = "vertical")

  if(!is.null(color.use)) {
    gg <- gg + scale_fill_manual(values=color.use)
    # gg <- gg + scale_color_manual(color.use)
  }

  if (!is.null(title)) {
    gg <- gg + guides(fill = guide_legend(title = title))
  }
  gg
}




#' A Seurat wrapper function for plotting gene expression using violin plot, dot plot or bar plot
#'
#' This function create a Seurat object from an input CellChat object, and then plot gene expression distribution using a modified violin plot or dot plot based on Seurat's function or a bar plot.
#' Please check \code{\link{StackedVlnPlot}},\code{\link{dotPlot}} and \code{\link{barPlot}}for detailed description of the arguments.
#'
#' USER can extract the signaling genes related to the inferred L-R pairs or signaling pathway using \code{\link{extractEnrichedLR}}, and then plot gene expression using Seurat package.
#'
#' @param object CellChat object
#' @param features Features to plot gene expression
#' @param signaling a char vector containing signaling pathway names for searching
#' @param enriched.only whether only return the identified enriched signaling genes in the database. Default = TRUE, returning the significantly enriched signaling interactions
#' @param type violin plot or dot plot
#' @param color.use defining the color for each cell group
#' @param group.by Name of one metadata columns to group (color) cells. Default is the defined cell groups in CellChat object
#' @param ... other arguments passing to either VlnPlot or DotPlot from Seurat package
#' @return
#' @export
#'
#' @examples

plotGeneExpression <- function(object, features = NULL, signaling = NULL, enriched.only = TRUE, type = c("violin", "dot","bar"), color.use = NULL, group.by = NULL, ...) {
  type <- match.arg(type)
  meta <- object@meta
  if (is.list(object@idents)) {
    meta$group.cellchat <- object@idents$joint
  } else {
    meta$group.cellchat <- object@idents
  }
  if (!identical(rownames(meta), colnames(object@data.signaling))) {
    cat("The cell barcodes in 'meta' is ", head(rownames(meta)),'\n')
    warning("The cell barcodes in 'meta' is different from those in the used data matrix.
              We now simply assign the colnames in the data matrix to the rownames of 'mata'!")
    rownames(meta) <- colnames(object@data.signaling)
  }

  w10x <- Seurat::CreateSeuratObject(counts = object@data.signaling, meta.data = meta)
  if (is.null(group.by)) {
    group.by <- "group.cellchat"
  }
  Seurat::Idents(w10x) <- group.by
  if (!is.null(features) & !is.null(signaling)) {
    warning("`features` will be used when inputing both `features` and `signaling`!")
  }
  if (!is.null(features)) {
    feature.use <- features
  } else if (!is.null(signaling)) {
    res <- extractEnrichedLR(object, signaling = signaling, geneLR.return = TRUE, enriched.only = enriched.only)
    feature.use <- res$geneLR
  }
  if (type == "violin") {
    gg <- StackedVlnPlot(w10x, features = feature.use, color.use = color.use, ...)
  } else if (type == "dot") {
    gg <- dotPlot(w10x, features = feature.use, color.use = color.use, ...)
  } else if (type == "bar") {
    gg <- barPlot(w10x, features = feature.use, color.use = color.use, ...)
  }
  return(gg)
}


#' Dot plot
#'
#'The size of the dot encodes the percentage of cells within a class, while the color encodes the AverageExpression level across all cells within a class
#'
#' @param object seurat object
#' @param features Features to plot (gene expression, metrics)
#' @param rotation whether rotate the plot
#' @param colormap RColorbrewer palette to use (check available palette using RColorBrewer::display.brewer.all()). default will use customed color palette
#' @param color.direction Sets the order of colours in the scale. If 1, the default, colours are as output by RColorBrewer::brewer.pal(). If -1, the order of colours is reversed.
#' @param color.use defining the color for each condition/dataset
#' @param idents Which classes to include in the plot (default is all)
#' @param group.by Name of one or more metadata columns to group (color) cells by
#' (for example, orig.ident); pass 'ident' to group by identity class
#' @param split.by Name of a metadata column to split plot by;
#' @param legend.width legend width
#' @param scale whther show x-axis text
#' @param col.min Minimum scaled average expression threshold (everything smaller will be set to this)
#' @param col.max Maximum scaled average expression threshold (everything larger will be set to this)
#' @param dot.scale Scale the size of the points, similar to cex
#' @param assay Name of assay to use, defaults to the active assay
#' @param angle.x angle for x-axis text rotation
#' @param hjust.x adjust x axis text
#' @param angle.y angle for y-axis text rotation
#' @param hjust.y adjust y axis text
#' @param show.legend whether show the legend
#' @param ... Extra parameters passed to DotPlot from Seurat package
#' @return ggplot2 object
#' @export
#'
#' @examples
#' @import ggplot2
dotPlot <- function(object, features, rotation = TRUE, colormap = "OrRd", color.direction = 1,  color.use = c("#F8766D","#00BFC4"), scale = TRUE, col.min = -2.5, col.max = 2.5, dot.scale = 6, assay = "RNA",
                    idents = NULL, group.by = NULL, split.by = NULL, legend.width = 0.5,
                    angle.x = 45, hjust.x = 1, angle.y = 0, hjust.y = 0.5, show.legend = TRUE, ...) {

  gg <- Seurat::DotPlot(object, features = features, assay = assay, cols = color.use,
                        scale = scale, col.min = col.min, col.max = col.max, dot.scale = dot.scale,
                        idents = idents, group.by = group.by, split.by = split.by,...)
  gg <- gg + theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
    theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10), axis.line = element_line(colour = 'black')) +
    theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5))+
    theme(axis.text.x = element_text(angle = angle.x, hjust = hjust.x), axis.text.y = element_text(angle = angle.y, hjust = hjust.y))

  gg <- gg + theme(legend.title = element_text(size = 10), legend.text = element_text(size = 8))
  if (is.null(split.by)) {
    gg <- gg + guides(color = guide_colorbar(barwidth = legend.width, title = "Scaled expression"),size = guide_legend(title = 'Percent expressed'))
  }

  if (rotation) {
    gg <- gg + coord_flip()
  }
  if (!is.null(colormap)) {
    if (is.null(split.by)) {
      gg <- gg + scale_color_distiller(palette = colormap, direction = color.direction, guide = guide_colorbar(title = "Scaled Expression", ticks = T, label = T, barwidth = legend.width), na.value = "lightgrey")
    }
  }
  if (!show.legend) {
    gg <- gg + theme(legend.position = "none")
  }
  return(gg)
}



#' Stacked Violin plot
#'
#' @param object seurat object
#' @param features Features to plot (gene expression, metrics)
#' @param color.use defining the color for each cell group
#' @param colors.ggplot whether use ggplot color scheme; default: colors.ggplot = FALSE
#' @param split.by Name of a metadata column to split plot by;
#' @param idents Which classes to include in the plot (default is all)
#' @param show.median whether show the median value
#' @param median.size the shape size of the median
#' @param show.text.y whther show y-axis text
#' @param line.size line width in the violin plot
#' @param pt.size size of the dots
#' @param plot.margin adjust the white space between each plot
#' @param angle.x angle for x-axis text rotation
#' @param vjust.x adjust x axis text
#' @param hjust.x adjust x axis text
#' @param ... Extra parameters passed to VlnPlot from Seurat package
#' @return ggplot2 object
#' @export
#'
#' @examples
#' @import ggplot2
#' @importFrom  patchwork wrap_plots
StackedVlnPlot<- function(object, features, idents = NULL, split.by = NULL,
                          color.use = NULL, colors.ggplot = FALSE,show.median = FALSE, median.size = 1,
                          angle.x = 90, vjust.x = NULL, hjust.x = NULL, show.text.y = TRUE, line.size = NULL,
                          pt.size = 0,
                          plot.margin = margin(0, 0, 0, 0, "cm"),
                          ...) {
  options(warn=-1)
  if (is.null(color.use)) {
    numCluster <- length(levels(Seurat::Idents(object)))
    if (colors.ggplot) {
      color.use <- NULL
    } else {
      color.use <- scPalette(numCluster)
    }
  }
  if (is.null(vjust.x) | is.null(hjust.x)) {
    angle=c(0, 45, 90)
    hjust=c(0, 1, 1)
    vjust=c(0, 1, 0.5)
    vjust.x = vjust[angle == angle.x]
    hjust.x = hjust[angle == angle.x]
  }

  plot_list<- purrr::map(features, function(x) modify_vlnplot(object = object, features = x, idents = idents, split.by = split.by, cols = color.use, show.median = show.median, median.size = median.size, pt.size = pt.size,
                                                              show.text.y = show.text.y, line.size = line.size, ...))

  # Add back x-axis title to bottom plot. patchwork is going to support this?
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(), axis.ticks.x = element_line()) +
    theme(axis.text.x = element_text(angle = angle.x, hjust = hjust.x, vjust = vjust.x)) +
    theme(axis.text.x = element_text(size = 10))

  # change the y-axis tick to only max value
  ymaxs<- purrr::map_dbl(plot_list, extract_max)
  plot_list<- purrr::map2(plot_list, ymaxs, function(x,y) x +
                            scale_y_continuous(breaks = c(y)) +
                            expand_limits(y = y))

  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1) + patchwork::plot_layout(guides = "collect")
  return(p)
}


#' modified vlnplot
#' @param object Seurat object
#' @param features Features to plot (gene expression, metrics)
#' @param split.by Name of a metadata column to split plot by;
#' @param idents Which classes to include in the plot (default is all)
#' @param cols defining the color for each cell group
#' @param show.median whether show the median value
#' @param median.size the shape size of the median
#' @param show.text.y whther show y-axis text
#' @param line.size line width in the violin plot
#' @param pt.size size of the dots
#' @param plot.margin adjust the white space between each plot
#' @param ... pass any arguments to VlnPlot in Seurat
#' @import ggplot2
#'
modify_vlnplot<- function(object,
                          features,
                          idents = NULL,
                          split.by = NULL,
                          cols = NULL,
                          show.median = FALSE,
                          median.size = 1,
                          show.text.y = TRUE,
                          line.size = NULL,
                          pt.size = 0,
                          plot.margin = margin(0, 0, 0, 0, "cm"),
                          ...) {
  options(warn=-1)
  p<- Seurat::VlnPlot(object, features = features, cols = cols, pt.size = pt.size, idents = idents, split.by = split.by,  ... )  +
    xlab("") + ylab(features) + ggtitle("")
  if (show.median) {
    p <- p + stat_summary(fun.y=median, geom="point", shape=3, size=median.size)
  }
  p <- p + theme(text = element_text(size = 10)) + theme(axis.line = element_line(size=line.size)) +
    theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 8), axis.line.x = element_line(colour = 'black', size=line.size),axis.line.y = element_line(colour = 'black', size= line.size))
  # theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5))
  p <- p + theme(plot.title= element_blank(), # legend.position = "none",
                 axis.title.x = element_blank(),
                 axis.text.x = element_blank(),
                 axis.ticks.x = element_blank(),
                 axis.title.y = element_text(size = rel(1), angle = 0),
                 axis.text.y = element_text(size = rel(1)),
                 plot.margin = plot.margin ) +
    theme(axis.text.y = element_text(size = 8))
  p <- p + theme(element_line(size=line.size))

  if (!show.text.y) {
    p <- p + theme(axis.ticks.y=element_blank(), axis.text.y=element_blank())
  }
  return(p)
}

#' extract the max value of the y axis
#' @param p ggplot object
#' @importFrom  ggplot2 ggplot_build
extract_max<- function(p){
  ymax<- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
}


#' Bar plot for average gene expression
#'
#' Please check \code{\link{barplot_internal}}for detailed description of the arguments.
#'
#' @param object seurat object
#' @param features Features to plot (gene expression, metrics)
#' @param color.use defining the color for each condition/dataset
#' @param group.by Name of one or more metadata columns to group (color) cells by
#' (for example, orig.ident); pass 'ident' to group by identity class
#' @param method methods for computing the average gene expression per cell group. By default = "truncatedMean", where a value should be assigned to 'trim;
#' @param trim the fraction (0 to 0.5) of observations to be trimmed from each end of x before the mean is computed
#' @param split.by Name of a metadata column to split plot by;
#' @param assay Name of assay to use, defaults to the active assay
#' @param x.lab.rot whether do rotation for the x.tick.label
#' @param ncol number of columns to show in the plot
#' @param ... Extra parameters passed to barplot_internal
#' @return ggplot2 object
#' @export
#'
#' @examples
#' @import ggplot2
barPlot <- function(object, features, group.by = NULL, split.by = NULL, color.use = NULL, method = c("truncatedMean", "triMean","median"),trim = 0.1, assay = "RNA",
                    x.lab.rot = FALSE, ncol = 1, ...) {
  method <- match.arg(method)
  if (is.null(group.by)) {
    labels = Idents(object)
  } else {
    labels = object@meta.data[,group.by]
  }
  FunMean <- switch(method,
                    truncatedMean = function(x) mean(x, trim = trim, na.rm = TRUE),
                    triMean = triMean,
                    median = function(x) median(x, na.rm = TRUE))

  if (!is.null(split.by)) {
    group = object@meta.data[,split.by]
    group.levels <- levels(group)
    df <- data.frame()
    for (i in 1:length(group.levels)) {
      data = GetAssayData(object, slot = "data", assay = assay)[, group == group.levels[i]]
      labels.use <- labels[group == group.levels[i]]
      dataavg <- aggregate(t(data[features, ]), list(labels.use) , FUN = FunMean)
      dataavg <- t(dataavg[,-1])
      colnames(dataavg) <- levels(labels.use)
      dataavg <- as.data.frame(dataavg)
      dataavg$gene = rownames(dataavg)
      df1 = reshape2::melt(dataavg, id.vars = c("gene"))
      colnames(df1) <- c("gene","labels","value")
      df1$condition = group.levels[i]
      df = rbind(df, df1)
    }
    df$labels <- factor(df$labels, levels = levels(labels))
    df$condition <- factor(df$condition, levels = group.levels)

  } else {
    data = GetAssayData(object, slot = "data", assay = assay)
    dataavg <- aggregate(t(data[features, ]), list(labels) , FUN = FunMean)
    dataavg <- t(dataavg[,-1])
    colnames(dataavg) <- levels(labels)
    dataavg$gene = rownames(dataavg)
    df1 = reshape2::melt(dataavg, id.vars = c("gene"))
    colnames(df1) <- c("gene","labels","value")
    df1$condition = df1[,"labels"]
    df = df1
  }
  gg <- list()
  for (i in 1:length(features)) {
    if (i < length(features)) {
      df.use = subset(df, gene == features[i])
      gg[[i]] <- barplot_internal(df.use, x = "labels", y = "value", fill = "condition",color.use = color.use,ylabel = features[i],remove.xtick = TRUE,x.lab.rot = x.lab.rot,...)
    }else {
      gg[[i]] <- barplot_internal(df.use, x = "labels", y = "value", fill = "condition",color.use = color.use,ylabel = features[i],remove.xtick = FALSE,x.lab.rot = x.lab.rot,...)
    }
  }

  p<- patchwork::wrap_plots(plotlist = gg, ncol = ncol)+ patchwork::plot_layout(guides = "collect")
  return(p)

}

#' Bar plot for dataframe
#'
#' @param df a dataframe
#' @param x Name of one column to show on the x-axis
#' @param y Name of one column to show on the y-axis
#' @param fill Name of one column to compare the values
#' @param color.use defining the color of bar plot;
#' @param percent.y whether showing y-values as percentage
#' @param width bar width
#' @param legend.title Name of legend
#' @param xlabel Name of x label
#' @param ylabel Name of y label
#' @param remove.xtick whether remove x tick
#' @param title.name Name of the main title
#' @param stat.add whether adding statistical test
#' @param stat.method,label.x parameters for ggpubr::stat_compare_means
#' @param show.legend Whether show the legend
#' @param x.lab.rot Whether rorate the xtick labels
#' @param size.text font size

#' @import ggplot2
#' @importFrom ggpubr stat_compare_means
#'
#' @return ggplot2 object
#' @export
barplot_internal <- function(df, x = "cellType", y = "value", fill = "condition", legend.title = NULL, width=0.6, title.name = NULL,
                             xlabel = NULL, ylabel = NULL, color.use = NULL,remove.xtick = FALSE,
                             stat.add = FALSE, stat.method = "wilcox.test", percent.y = FALSE, label.x = 1.5,
                             show.legend = TRUE, x.lab.rot = FALSE, size.text = 10) {

  gg <- ggplot(df, aes_string(x=x, y=y, fill = fill, color = fill)) + geom_bar(stat="identity", width=width, position=position_dodge()) +
    theme_classic() + scale_x_discrete(limits = (levels(df$x))) + theme(axis.text.x = element_text(angle = 45, hjust = 1,size=10))

  gg <- gg + ylab(ylabel) + xlab(xlabel) + theme_classic() +
    labs(title = title.name) +  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) +
    theme(text = element_text(size = size.text), axis.text = element_text(colour="black"))
  if (!is.null(color.use)) {
    gg <- gg + scale_fill_manual(values = alpha(color.use, alpha = 1), drop = FALSE)
    gg <- gg + scale_color_manual(values = alpha(color.use, alpha = 1), drop = FALSE) + guides(colour = FALSE)
  }
  if (stat.add) {
    gg <- gg + ggpubr::stat_compare_means(mapping = aes_string(group = fill), method = stat.method, label.x = label.x,
                                          label = "p.format", size = 3)
  }
  # if (show.mean) {
  #   gg <- gg + stat_summary(fun.y=mean, geom="point", shape=20, size=10, color="red", fill="red")
  # }
  if (remove.xtick) {
    gg <- gg + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.x=element_blank())
  }
  if (percent.y) {
    gg <- gg + scale_y_continuous(labels = scales::percent_format(accuracy = 1))
  }
  if (is.null(legend.title)) {
    gg <- gg + theme(legend.title = element_blank())
  } else {
    gg <- gg + guides(fill=guide_legend(legend.title))
  }
  if (!show.legend) {
    gg <- gg + theme(legend.position = "none")
  }
  if (x.lab.rot) {
    gg <- gg + theme(axis.text.x = element_text(angle = 45, hjust = 1, size=size.text))
  }
  gg
  return(gg)
}


#' Visualize the spatial distribution of cell type proportion in a geom scatterpie plot
#'
#' @param proportion Data frame, cell type proportion estimated by CARD in either original resolution or enhanced resolution.
#' @param spatial_location Data frame, spatial location information.
#' @param colors Vector of color names that you want to use, if NULL, we will use the color palette "Spectral" from RColorBrewer package.
#' @param radius Numeric value about the radius of each pie chart, if NULL, we will calculate it inside the function.
#' @param seed Seed number about generating the colors if users do not provide the colors, if NULL, we will generate it inside the function
#'
#' @import ggplot2
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette
#' @return Returns a ggplot2 figure.
#'
#' @export
#'
# This code is adapted from CARD package
visualizePie <- function(proportion,spatial_location,colors = NULL,radius = NULL,seed = NULL){
  res = as.data.frame(proportion)
  res = res[,gtools::mixedsort(colnames(res))]
  location = as.data.frame(spatial_location)
  if(sum(rownames(res)==rownames(location))!= nrow(res)){
    stop("The rownames of proportion data does not match with the rownames of spatial location data")
  }
  colorCandidate = c("#1e77b4","#ff7d0b","#ceaaa3","#2c9f2c","#babc22","#d52828","#9267bc",
                     "#8b544c","#e277c1","#d42728","#adc6e8","#97df89","#fe9795","#4381bd","#f2941f","#5aa43a","#cc4d2e","#9f83c8","#91675a",
                     "#da8ec8","#929292","#c3c237","#b4e0ea","#bacceb","#f7c685",
                     "#dcf0d0","#f4a99f","#c8bad8",
                     "#F56867", "#FEB915", "#C798EE", "#59BE86", "#7495D3",
                     "#D1D1D1", "#6D1A9C", "#15821E", "#3A84E6", "#997273",
                     "#787878", "#DB4C6C", "#9E7A7A", "#554236", "#AF5F3C",
                     "#93796C", "#F9BD3F", "#DAB370", "#877F6C", "#268785",
                     "#f4f1de","#e07a5f","#3d405b","#81b29a","#f2cc8f","#a8dadc","#f1faee","#f08080")
  if(is.null(colors)){
    #colors = brewer.pal(11, "Spectral")
    if(ncol(res) > length(colorCandidate)){
      colors = colorRampPalette(colorCandidate)(ncol(res))
    }else{

      if(is.null(seed)){
        iseed = 12345
      }else{
        iseed = seed
      }
      set.seed(iseed)
      colors = colorCandidate[sample(1:length(colorCandidate),ncol(res))]
    }
  }else{
    colors = colors
  }
  data = cbind(res,location)
  ct.select = colnames(res)
  if(is.null(radius)){
    radius = (max(data$x) - min(data$x)) * (max(data$y) - min(data$y))
    radius = radius / nrow(data)
    radius = radius / pi
    radius = sqrt(radius) * 0.85
  }else{
    #### avoid the situation when the radius does not generate the correct figure
    radius = radius
  }
  p = suppressMessages(ggplot() + scatterpie::geom_scatterpie(aes(x=x_cent, y=y_cent,r = radius),data=data,
                                                  cols=ct.select,color=NA) + coord_fixed()+ #coord_fixed(ratio = 1*max(data$x)/max(data$y)) +
                         scale_fill_manual(values =  colors)+
                         theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
                               panel.background = element_blank(),
                               plot.background = element_blank(),
                               panel.border = element_rect(colour = "grey89", fill=NA, size=0.5),
                               axis.text =element_blank(),
                               axis.ticks =element_blank(),
                               axis.title =element_blank(),
                               legend.title=element_text(size = 16,face="bold"),
                               legend.text=element_text(size = 15),
                               legend.key = element_rect(colour = "transparent", fill = "white"),
                               legend.key.size = unit(0.45, 'cm'),
                               strip.text = element_text(size = 16,face="bold"),
                               legend.position="bottom")+
                         guides(fill=guide_legend(title="Cell Type")))
  return(p)
}






#' @title SpatialTopicPlot
#'
#' @param object SpatialCellChat object
#' @param topics a vector of numeric indices or character names specifying which spatial topics to visualize. If \code{NULL}, all topics available are plotted
#' @param slot.name the slot name of object that is used for topic analysis
#' @param pattern "outgoing" or "incoming"
#' @param topics.key a character string used as a prefix to generate topic names
#' @param color.heatmap a character string or vector indicating the colormap option to use. It can be the avaibale color palette in brewer.pal() or viridis_pal()
#' @param n.colors number of basic colors to generate from color palette
#' @param direction sets the order of colors in the scale. If 1, the default colors are used. If -1, the order of colors is reversed
#' @param cutoff a numeric value used to filter low loading values
#' @param color.use defining the color for each cell group
#' @param alpha the transparency of individual spot
#' @param point.size the size of spots
#' @param legend.size the size of legend
#' @param legend.text.size the text size on the legend
#' @param shape.by the shape of individual spot
#' @param ncol number of columns if plotting multiple plots
#' @param combine whether outputting an combined plot using patchwork
#' @param show.legend whether show each figure legend
#' @param show.legend.combined whether show the figure legend for the last plot
#' @import ggplot2
#'
#' @return a ggplot2 object
#' @export
spatialTopicPlot <- function(
    object,
    topics = NULL,
    slot.name = "net",
    pattern = c("incoming","outgoing"),
    topics.key = "Topic_",
    color.heatmap = "Reds",
    n.colors = 8,
    direction = 1,
    cutoff = NULL,
    color.use = NULL,
    alpha = 1,
    point.size = 1,
    legend.size = 3,
    legend.text.size = 8,
    shape.by = 16,
    ncol = NULL,
    combine = T,
    show.legend = TRUE,
    show.legend.combined = FALSE
){
  # Column number of combined plots
  if (is.null(x = ncol)) {
    ncol <- 2
    if (length(x = topics) == 1) {
      ncol <- 1
    }
    if (length(x = topics) > 6) {
      ncol <- 3
    }
    if (length(x = topics) > 9) {
      ncol <- 4
    }
  }


  coords <- object@images$coordinates
  if (ncol(coords) == 2) {
    colnames(coords) <- c("x_cent", "y_cent")
    temp_coord = coords
    coords[, 1] = temp_coord[, 2]
    coords[, 2] = temp_coord[, 1]
  } else {
    stop("Please check the input 'coordinates' and make sure it is a two column matrix.")
  }

  pattern <- match.arg(pattern)
  pattern.res <- methods::slot(object, slot.name)$topic[[pattern]]
  cell.embeddings <- pattern.res$cell %>% Matrix::t() # topic x cell
  rownames(cell.embeddings) <- paste0(topics.key,seq_len(NROW(cell.embeddings)))
  meta <- object@meta

  if (length(color.heatmap) == 1) {
    colormap <- tryCatch({
      RColorBrewer::brewer.pal(n = n.colors, name = color.heatmap)
    }, error = function(e) {
      scales::viridis_pal(option = color.heatmap, direction = -1)(n.colors)
    })
    if (direction == -1) {
      colormap <- rev(colormap)
    }
    colormap <- colorRampPalette(colormap)(99)
    colormap[1] <- "#E5E5E5"
  } else {
    colormap <- color.heatmap
  }

  if (is.null(topics)) {
    feature.use <- rownames(cell.embeddings)
  } else {
    feature.use <- topics
  }


  df <- data.frame(x = coords[, 1], y = coords[, 2])
  data.use <- cell.embeddings[feature.use, , drop = FALSE]


  if (!is.null(cutoff)) {
    cat("Applying a cutoff of ", cutoff, "to the values...", '\n')
    data.use[data.use <= cutoff] <- NA
  }

  numFeature = length(feature.use)
  gg <- vector("list", numFeature)
  for (i in seq_len(numFeature)) {
    feature.name <- rownames(data.use)[i]

    df$feature.data <- data.use[i, , drop = T]
    g <- ggplot(data = df, aes(x, y)) +
      geom_point(
        aes(colour = feature.data),
        alpha = alpha,
        size = point.size,
        shape = shape.by
      ) +
      scale_colour_gradientn(
        colours = colormap,
        guide = guide_colorbar(
          title = NULL,
          ticks = T,
          label = T,
          barwidth = 0.5
        ),
        na.value = "grey90"
      ) +
      theme(legend.position = "right") +
      theme(
        legend.title = element_blank(),
        legend.text = element_text(size = legend.text.size),
        legend.key.size = grid::unit(0.15, "inches")
      )  + # , legend.key.size = grid::unit(0.4, "inches")
      ggtitle(feature.name) + theme(plot.title = element_text(
        hjust = 0.5,
        vjust = 0,
        size = 10
      )) +
      theme(
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank()
      ) + xlab(NULL) + ylab(NULL) +
      theme(legend.key = element_blank())
    g <- g + coord_fixed() + scale_y_reverse()
    #g <- g + coord_fixed(ratio = 1*diff(range(coords$x_cent))/diff(range(coords$y_cent))) + scale_y_reverse()
    if (!show.legend) {
      g <- g + theme(legend.position = "none")
    }
    if (show.legend.combined & i == numFeature) {
      g <-
        g + theme(
          legend.position = "right",
          legend.key.height = grid::unit(0.15, "in"),
          legend.key.width = grid::unit(0.5, "in"),
          legend.title = element_blank(),
          legend.key = element_blank()
        )
    }
    gg[[i]] <- g
  }
  if (combine) {
    gg <- patchwork::wrap_plots(gg, ncol = ncol)
  }
  return(gg)
}


#' @title plotly_spatialLRpairPlot
#'
#' @param object SpatialCellChat object
#' @param features a char vector containing features to visualize. `features` can be genes or column names of `object@meta`.
#' @param pairLR.use a char vector or a data frame consisting of one column named "interaction_name", defining the L-R pairs of interest
#' @param cells.highlight cells/spots to be highlighted. It can be cell IDs (character vector) or indices (numeric vector)
#' @param enriched.only whether only return the identified enriched signaling genes in the database. Default = TRUE, returning the significantly enriched signaling interactions
#' @param thresh threshold of the p-value for determining significant interaction
#' @param do.group set `do.group = TRUE` when computing scores of each cell group; set `do.group = FALSE` when computing scores of each individual cell
#' @param colormap RColorbrewer palette to use (check available palette using RColorBrewer::display.brewer.all()). default will use customed color palette
#' @param reversescale whether to reverse the direction of the color scale
#' @param point.size the size of spots
#' @param nohighlight.alpha the opacity of non-highlighted cells/spots
#' @param do.plot whether to generate and return plots
#'
#' @return a list of ggplot2/plotly objects (`do.plot = TRUE`) or the number of features (`do.plot = FALSE`)
#' @export
#'
#' @examples
plotly_spatialLRpairPlot <- function (
    object,
    features = NULL,
    pairLR.use = NULL,
    cells.highlight = NULL,
    enriched.only = F,
    thresh = 0.05,
    do.group = F,
    colormap = "Spectral",
    # n.colors = 8,
    reversescale = T,
    point.size = 4,
    nohighlight.alpha = 0.3,
    do.plot=T
){
  coords <- object@images$coordinates
  if (ncol(coords) == 2) {
    colnames(coords) <- c("x_cent", "y_cent")
    temp_coord = coords
    coords[, 1] = temp_coord[, 2]
    coords[, 2] = temp_coord[, 1]
  }
  else {
    stop("Please check the input 'coordinates' and make sure it is a two column matrix.")
  }

  # add idents info
  spot_labels <- object@idents

  data <- as.matrix(object@data)
  meta <- object@meta

  if (!is.null(features) & !is.null(pairLR.use)) {
    stop("Please don't input features or pairLR.use simultaneously.")
  }
  if (!is.null(features) & length(features) > 2) {
    stop("Please don't input more than 2 features.")
  }
  if (!is.null(pairLR.use) & length(pairLR.use) > 1) {
    stop("Please don't input more than 1 LR pair.")
  }

  df <- coords


  if (!is.null(pairLR.use)) {
    if (is.character(pairLR.use)) {
      pairLR.use <- data.frame(interaction_name = pairLR.use)
    }

    if (enriched.only) {
      if (do.group) {
        object@net$prob[object@net$pval > thresh] <- 0
        pairLR.use.name <- pairLR.use$interaction_name[pairLR.use$interaction_name %in% dimnames(object@net$prob)[[3]]]
        prob <- object@net$prob[, , pairLR.use.name,drop = F] # return an array
        prob.sum <- apply(prob > 0, 3L, sum)
        names(prob.sum) <- pairLR.use.name
        signaling.includes <- names(prob.sum)[prob.sum > 0]
        pairLR.use <- pairLR.use[pairLR.use$interaction_name %in% signaling.includes, , drop = F]
      } else {
        pairLR.use.name <-
          pairLR.use$interaction_name[pairLR.use$interaction_name %in% dimnames(object@net$prob.cell)[[3]]]
        prob.cell <-
          object@net$prob.cell[, , pairLR.use.name,drop = FALSE]
        prob.sum <- apply(prob.cell > 0, 3L, sum)
        names(prob.sum) <- pairLR.use.name
        signaling.includes <- names(prob.sum)[prob.sum > 0]
        pairLR.use <-
          pairLR.use[pairLR.use$interaction_name %in% signaling.includes, , drop = FALSE]
      }
      if (length(pairLR.use$interaction_name) == 0) {
        stop(
          paste0(
            "There is no significant communication related with the input `pairLR.use`. Set `enriched.only = FALSE` to show non-significant signaling."
          )
        )
      }
    }
    LR.pair <- object@LR$LRsig[pairLR.use$interaction_name, c("ligand", "receptor")]
    geneL <- unique(LR.pair$ligand)
    geneR <- unique(LR.pair$receptor)
    # geneL <- extractGeneSubset(geneL, object@DB$complex,
    #                            object@DB$geneInfo)
    # geneR <- extractGeneSubset(geneR, object@DB$complex,
    #                            object@DB$geneInfo)
    feature.use <- c(geneL, geneR)
  } else {
    feature.use <- features
  }

  if (length(intersect(feature.use, rownames(data))) > 0) {
    feature.use <- feature.use[feature.use %in% rownames(data)]
    data.use <- data[feature.use, , drop = FALSE]
  } else if (length(intersect(feature.use, colnames(meta))) > 0) {
    feature.use <- feature.use[feature.use %in% colnames(meta)]
    data.use <- t(meta[, feature.use, drop = FALSE])
  } else {
    stop("Please check your input! ")
  }

  numFeature = length(feature.use)
  list_pl <- vector("list", numFeature)
  max.feature.data <- max(data.use)

  if(do.plot){
    # Custom Hover Text: https://plotly.com/r/text-and-annotations/
    for (i in seq_len(numFeature)) {
      df$feature.data <- data.use[i,,drop=T]
      list_pl[[i]] <-
        plotFeatures(
          df_score = df,
          cells.highlight = cells.highlight,
          method = c("2d"),
          plot.title = feature.use[i],
          point.size = point.size,
          nohighlight.alpha = nohighlight.alpha,
          show.colorscale = T,
          color.heatmap=colormap,
          color.0 = "default",
          color.NA = "#E5E5E5",
          n.colors = 8,
          direction=if_else(reversescale,-1,1),
          legend.group = "score",
          legend.group.title = "score's value",
          scene.name = "scene1",
          colorbar.y = 0.5
        )
    }
    return(list_pl)
    # do.plot
  } else {
    return(numFeature)
    # not do.plot
  }
}


#' @title plotFeatures
#'
#' @param df_score a data frame with columns ordered as "x", "y", and "score".
#' @param group.by name of one or more metadata columns to group (color) cells by
#' @param group.highlight cell groups to be highlighted
#' @param show.group whether to display group annotations
#' @param cells.highlight cells/spots to be highlighted
#' @param method plotting mode, one of `"normal"` (ggplot2), `"2d`, or `"3d"` (Plotly)
#' @param plot.title title of the plot
#' @param plot.ncol the number of columns in the grid when `show.group = TRUE`
#' @param point.size the size of spots
#' @param highlight.alpha the opacity of highlighted cells/spots
#' @param nohighlight.alpha the opacity of non-highlighted cells/spots
#' @param show.colorscale whether to display the colorbar
#' @param color.heatmap A character string or vector indicating the colormap option to use. It can be the avaibale color palette in brewer.pal() or viridis_pal()
#' @param color.0 color used for 0 in the colormap
#' @param color.NA color used to represent missing values (NA) in the plot
#' @param n.colors number of basic colors to generate from color palette
#' @param direction sets the order of colors in the scale. If 1, the default colors are used. If -1, the order of colors is reversed.
#' @param color.use defining the color for each cell group
#' @param legend.group the parameter `legendgroup` in `plotly::plot_ly` and `plotly::add_trace`. Name of the legend group for traces, used to combine multiple traces under the same legend.
#' @param legend.group.title the parameter `legendgrouptitle` in `plotly::plot_ly` and `plotly::add_trace`. The title text displayed for the legend group.
#' @param scene.name name of the Plotly 3D scene
#' @param colorbar.y the relative position of the color bar along the y-axis
#' @param three.dim.z the fixed z-axis coordinate in 3D mode
#' @import ggplot2
#'
#' @return ggplot2 object or plotly object
#' @export
#'
#' @examples
plotFeatures <- function(
    df_score,
    group.by=NULL,
    group.highlight=NULL,
    show.group=F,
    cells.highlight=NULL,
    method = c("normal","3d","2d"),
    plot.title = NULL,
    plot.ncol=2,
    point.size = 3,
    # point.stroke=2,
    highlight.alpha = 1,
    nohighlight.alpha = 0.65,
    show.colorscale = T,
    color.heatmap="Blues",
    color.0 = NULL,
    color.NA = "#E5E5E5",
    n.colors = 8,
    direction=1,
    color.use = NULL,

    legend.group = "score",
    legend.group.title = "score's value",
    scene.name = "scene1",
    colorbar.y = 0.5,
    three.dim.z = 0
){
  method <- match.arg(method)

  # rename the cols
  colnames(df_score) <- c("x", "y", "score")

  # set the title of plot
  if (is.null(plot.title)) {
    plot.title <- plot.title #NULL
  } else{
    plot.title <- labs(title = plot.title)
  }
  plot.title.theme.opts <- theme(plot.title = element_text(hjust = 0.5,vjust = 0, size = 10,face="bold"))

  if(!is.null(cells.highlight)){
    if(is.character(cells.highlight) & (!is.null(rownames(df_score))) ){
      if(method == "normal"){
        df_score$spot_highlight <- dplyr::if_else(rownames(df_score)%in%cells.highlight,highlight.alpha,nohighlight.alpha)
      } else if (method == "3d"){
        df_score$spot_highlight <- rownames(df_score)%in%cells.highlight
      } else if (method == "2d"){
        df_score$spot_highlight <- rownames(df_score)%in%cells.highlight
      }
    }
    if(is.numeric(cells.highlight)){
      if(method == "normal"){
        df_score$spot_highlight <- dplyr::if_else(seq_len(NROW(df_score)) %in% cells.highlight,highlight.alpha,nohighlight.alpha)
      } else if (method == "3d"){
        df_score$spot_highlight <- seq_len(NROW(df_score)) %in% cells.highlight
      } else if (method == "2d"){
        df_score$spot_highlight <- seq_len(NROW(df_score)) %in% cells.highlight
      }
    }
  }
  if((!is.null(group.highlight)) & (!is.null(group.by))){
    if(method == "normal"){
      cells.highlight <- which(group.by==group.highlight)
      df_score$spot_highlight <- dplyr::if_else(seq_len(NROW(df_score)) %in% cells.highlight,highlight.alpha,nohighlight.alpha)
    } else if (method == "3d"){
      cells.highlight <- which(group.by==group.highlight)
      df_score$spot_highlight <- seq_len(NROW(df_score)) %in% cells.highlight
    } else if (method == "2d"){
      cells.highlight <- which(group.by==group.highlight)
      df_score$spot_highlight <- seq_len(NROW(df_score)) %in% cells.highlight
    }
  }

  if(length(color.heatmap) == 1){

    if(color.heatmap=="RedsBlues"){

      if(direction==-1) {switch.RedsBlues=T} else {switch.RedsBlues=F}
      if(!is.null(color.0)){
        if(color.0=="default"){
          colormap <- generate_RedsBlues_colormap(df_score$score,mid.color = "#E5E5E5",switch.RedsBlues = switch.RedsBlues)
        } else if(color.0!="default" & is.character(color.0)){
          colormap <- generate_RedsBlues_colormap(df_score$score,mid.color = color.0,switch.RedsBlues = switch.RedsBlues)
        }
      } else { # color.0=NULL
        colormap <- generate_RedsBlues_colormap(df_score$score,mid.color = "#E5E5E5",switch.RedsBlues = switch.RedsBlues)
      }

    } else {

      colormap <- tryCatch({
        RColorBrewer::brewer.pal(n = n.colors, name = color.heatmap)
      }, error = function(e) {
        (scales::viridis_pal(option = color.heatmap, direction = -1))(n.colors)
      })
      if (direction == -1) {
        colormap <- rev(colormap)
      }
      colormap <- colorRampPalette(colormap)(99)
      if(!is.null(color.0)){
        if(color.0=="default"){
          colormap[1] <- "#E5E5E5"
        } else if(color.0!="default" & is.character(color.0)){
          colormap[1] <- color.0
        }
      }

    }

  } else {
    colormap <- color.heatmap
    if (direction == -1) {
      colormap <- rev(colormap)
    }
  }

  my.guide = guide_colorbar(title = NULL, ticks = T, label = T, barwidth = 0.5)
  df_score$cell_id <- rownames(df_score)
  if (method == "normal") {

    scChat_theme_opts <- theme(
      panel.background = element_blank(),
      panel.grid = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank()
    ) + theme(aspect.ratio = 1, legend.key = element_blank())

    if(!show.group){

      if(is.null(df_score$spot_highlight)){
        g <- ggplot(data = df_score) +
          geom_point(aes(
            x = x,
            y = y,
            color = score #,shape = point_shape
          ), size = point.size)
      } else {
        g <- ggplot(data = df_score) +
          geom_point(aes(
            x = x,
            y = y,
            alpha = spot_highlight,
            color = score #,shape = point_shape
          ), size = point.size)
      }

      g <- g + scale_colour_gradientn(colours = colormap, guide = my.guide,na.value = color.NA) +
        # scale_color_distiller(palette = "RdYlBu",direction = -1)+
        guides(alpha = guide_none()) + plot.title + plot.title.theme.opts +
        coord_fixed() + scale_y_reverse() + scChat_theme_opts

    }else{
      if(!is.factor(group.by)) group.by <- factor(group.by)
      df_score$group.label <- group.by
      if (is.null(color.use)) {
        color.use <- scPalette(nlevels(df_score$group.label))
        names(color.use) <- levels(df_score$group.label)
      }

      if(is.null(df_score$spot_highlight)){
        g.dim <- ggplot(data = df_score) +
          geom_point(mapping = aes(x = x, y = y, colour = group.label), size = point.size) +
          scale_color_manual(values = color.use, na.value = color.NA,guide=guide_legend(title = NULL))

        g.feature <- ggplot(data = df_score) +
          geom_point(aes(
            x = x,
            y = y,
            color = score #,shape = point_shape
          ), size = point.size)
      } else {
        g.dim <- ggplot(data = df_score) +
          geom_point(
            mapping = aes(
              x = x,
              y = y,
              alpha = spot_highlight,
              colour = group.label
            ),
            size = point.size
          ) +
          scale_color_manual(values = color.use,
                             na.value = color.NA,
                             guide = guide_legend(title = NULL))

        g.feature <- ggplot(data = df_score) +
          geom_point(aes(
            x = x,
            y = y,
            alpha = spot_highlight,
            color = score #,shape = point_shape
          ), size = point.size)
      }

      #   # theme(legend.position = legend.position)+
      #   # theme(legend.title = element_blank(), legend.text = element_text(size = legend.text.size)) +
      #   # guides(color = guide_legend(override.aes = list(size = legend.size),ncol = ncol, byrow = byrow)) +
      #   # theme(panel.background = element_blank(),axis.ticks = element_blank(), axis.text = element_blank()) +
      #   # xlab(NULL) + ylab(NULL) + coord_fixed() + theme(aspect.ratio = 1) +
      #   # theme(legend.key = element_blank()) + scale_y_reverse()
      g.dim <- g.dim +guides(alpha = guide_none())+ ggtitle("Cell annotation") +
        plot.title.theme.opts +
        coord_fixed() + scale_y_reverse() + scChat_theme_opts

      g.feature <- g.feature +
        scale_color_gradientn(colours = colormap, guide = my.guide ,na.value = color.NA) +
        # scale_fill_manual(values = color.use, na.value = color.NA)+
        guides(alpha = guide_none()) + plot.title + plot.title.theme.opts +
        coord_fixed() + scale_y_reverse() + scChat_theme_opts

      g <- patchwork::wrap_plots(g.dim,g.feature,ncol = plot.ncol)
    }
    return(g)

  } else if (method == "3d") {
    df_score$z <- three.dim.z
    if(is.null(df_score$spot_highlight)){
      df_score$spot_highlight <- T
    }

    pl <- plotly::plot_ly(
      data = df_score[df_score$spot_highlight, ],
      x = ~ x,
      y = ~ y,
      z = ~ z,
      scene = scene.name,
      color =  ~ score,
      colors = colormap,
      opacity = 1,
      symbol = I(19),
      # span = I(2),
      legendgroup = legend.group,
      legendgrouptitle = list(text = legend.group.title),
      type = "scatter3d",
      mode = "markers",
      hoverinfo = 'text',
      customdata = ~ cell_id,
      text = ~ paste(
        "<i><b>",
        cell_id,
        "</b></i>",
        '<br> x:',
        x,
        ", y:",
        y,
        ", z:",
        z,
        # '<br>cell type: <b>',
        # spot_labels,
        '<br>score:',
        score,
        '<br>'
      ),
      marker = list(size = point.size)
    ) %>% plotly::add_trace(
      data = df_score[!df_score$spot_highlight, ],
      x = ~ x,
      y = ~ y,
      z = ~ z,
      scene = scene.name,
      color =  ~ score,
      colors = colormap,
      opacity = nohighlight.alpha,
      symbol = I(19),
      # span = I(2),
      legendgroup = legend.group,
      legendgrouptitle = list(text = legend.group.title),
      type = "scatter3d",
      mode = "markers",
      hoverinfo = 'text',
      customdata = ~ cell_id,
      text = ~ paste(
        "<i><b>",
        cell_id,
        "</b></i>",
        '<br> x:',
        x,
        ", y:",
        y,
        ", z:",
        z,
        # '<br>cell type: <b>',
        # spot_labels,
        '<br>score:',
        score,
        '<br>'
      ),
      marker = list(size = point.size)
    ) %>% plotly::layout(
      title = plot.title[["title"]],
      grid = generate_grid_nrows(1),
      scene = generate_custom_scenes3d(
        row = 0,
        projection = "orthographic",
        zoom = 2.5
      ),
      legend = custom_legend
    ) %>% plotly::plotly_build()

    # change the style of colorbar
    pl[["x"]][["data"]][[1]][["marker"]][["colorbar"]] <- list(
      orientation = "v",
      len = 1 / 4,
      y = colorbar.y,
      lenmode = "fraction",
      thickness = 8,
      # default's 1/2
      title = list(text = "")
    )
    # show the colorbar
    pl[["x"]][["data"]][[1]][["marker"]][["showscale"]] <-
      show.colorscale
    # not show other elements' scales
    pl[["x"]][["data"]][[2]][["marker"]][["showscale"]] <- F
    pl[["x"]][["data"]][[3]][["marker"]][["showscale"]] <- F

    # not show "trace" legend
    pl[["x"]][["layout"]][["showlegend"]] <- F
    return(pl)


  } else if (method == "2d"){

    if(is.null(df_score$spot_highlight)){
      df_score$spot_highlight <- T
    }

    pl <- plotly::plot_ly(
      data = df_score[df_score$spot_highlight, ],
      x = ~ x,
      y = ~ y,

      color =  ~ score,
      colors = colormap,
      opacity = 1,
      symbol = I(19),
      # span = I(2),
      legendgroup = legend.group,
      legendgrouptitle = list(text = legend.group.title),
      type = "scatter",
      mode = "markers",
      hoverinfo = 'text',
      customdata = ~ cell_id,
      text = ~ paste(
        "<i><b>",
        plot.title[["title"]],
        "</b></i>",
        '<br> x:',
        x,
        ", y:",
        y,
        '<br> cell_id:',
        cell_id,
        # '<br>cell type: <b>',
        # spot_labels,
        '<br>score:',
        score,
        '<br>'
      ),
      marker = list(size = point.size)
    ) %>% plotly::add_trace(
      data = df_score[!df_score$spot_highlight, ],
      x = ~ x,
      y = ~ y,

      color =  ~ score,
      colors = colormap,
      opacity = nohighlight.alpha,
      symbol = I(19),
      # span = I(2),
      legendgroup = legend.group,
      legendgrouptitle = list(text = legend.group.title),
      type = "scatter",
      mode = "markers",
      hoverinfo = 'text',
      customdata = ~ cell_id,
      text = ~ paste(
        "<i><b>",
        plot.title[["title"]],
        "</b></i>",
        '<br> x:',
        x,
        ", y:",
        y,
        '<br> cell_id:',
        cell_id,
        # '<br>cell type: <b>',
        # spot_labels,
        '<br>score:',
        score,
        '<br>'
      ),
      marker = list(size = point.size)
    ) %>% plotly::layout(
      title = plot.title[["title"]],
      yaxis = list(
        autorange = "reversed",
        title = "",
        showgrid = FALSE,
        ticks = "",
        # ticktext = "",
        tickvals = "",
        scaleanchor = "x",
        tickformat = "%0f",
        showline = FALSE,
        visible = FALSE
      ),
      xaxis = list(
        title = "",
        showgrid = FALSE,
        ticks = "",
        # ticktext = "",
        tickvals = "",
        tickformat = "%0f",
        showline = FALSE,
        visible = FALSE
      ),
      legend = custom_legend
    ) %>% plotly::plotly_build()

    # change the style of colorbar
    pl[["x"]][["data"]][[1]][["marker"]][["colorbar"]] <- list(
      orientation = "v",
      len = 1 / 4,
      y = colorbar.y,
      lenmode = "fraction",
      thickness = 8,
      # default's 1/2
      title = list(text = "")
    )
    # show the colorbar
    pl[["x"]][["data"]][[1]][["marker"]][["showscale"]] <-
      show.colorscale
    # not show other elements' scales
    pl[["x"]][["data"]][[2]][["marker"]][["showscale"]] <- F
    pl[["x"]][["data"]][[3]][["marker"]][["showscale"]] <- F

    # not show "trace" legend
    pl[["x"]][["layout"]][["showlegend"]] <- F
    return(pl)

  }
}



#' @title plotly_spatialLRpairPlot_shiny
#'
#' @param object SpatialCellChat object
#' @param features a char vector containing features to visualize. `features` can be genes or column names of `object@meta`.
#' @param pairLR.use a char vector or a data frame consisting of one column named "interaction_name", defining the L-R pairs of interest
#' @param point.size.dimplot point size used in the spatialDimPlot
#' @param point.size.LRplot point size used in the spatialLRpairPlot
#' @param nohighlight.alpha the opacity of non-highlighted cells/spots
#' @param port port number used to launch the Shiny application
#'
#' @return A Shiny application object
#' @export
#'
#' @examples
plotly_spatialLRpairPlot_shiny <- function(
    object,
    features = NULL,
    pairLR.use = NULL,
    point.size.dimplot = 2,
    point.size.LRplot = 6,
    nohighlight.alpha = 0.45,
    port = 10333
){

  # object <- cortexChat
  # features = "Tgfb1"
  # pairLR.use = NULL
  # point.size.dimplot = 6
  # point.size.LRplot = 4
  # nohighlight.alpha = 0.45
  # port = 10333

  n.pls <- plotly_spatialLRpairPlot(
    object = object,
    features = features,
    pairLR.use = pairLR.use,
    enriched.only = F,
    do.group = F,
    do.plot=F
  )

  ui <- shiny::fluidPage(
    shiny::tags$head(
      shiny::tags$style(shiny::HTML("
      /* Your CSS code here */
      .plot-column {
        padding: 5px;
      }

      .plotly {
        width: 100% !important;
        height: 100% !important;
      }

      @media (max-width: 768px) {
        .plot-column {
          width: 100% !important;
        }
      }
    "))
    ),
    shiny::fluidRow(
      shiny::column(width = 6, class = "plot-column",
                    shiny::div(style = "height: 250px;"),
                    plotly::plotlyOutput(outputId = "Celltype", width = "100%"),
                    # rev button
                    shiny::actionButton("rev", "Reverse"),
                    # update button
                    shiny::actionButton("update", "Update"),
                    # clear button
                    shiny::actionButton("reset", "Clear"),
                    # done button
                    shiny::actionButton("done", "Done")),
      shiny::column(width = 6, class = "plot-column",
                    shiny::uiOutput("dynamic_outputs"))
    ),
    shiny::br(),
    shiny::p("Selected Cells by Tools Bar:"),
    shiny::verbatimTextOutput("selecting"),
    shiny::p("Selected Cells by Highlight:"),
    shiny::verbatimTextOutput("selecting2")
  )

  server <- function(input, output, session) {
    vals <- shiny::reactiveValues(
      keepcells = rownames(object@images$coordinates),
      Ncells = NROW(object@images$coordinates),
      rev = F
    )

    output_list <- shiny::reactive({
      # generate multiple Plotly plots by `lapply` function
      lapply(seq_len(n.pls), function(i) {
        plotly::plotlyOutput(outputId = paste0('LRExpr', i), width = "100%")
      })
    })

    # render dynamic output
    output$dynamic_outputs <- shiny::renderUI({
      # convert a list of Plotly plots into a list of tags
      shiny::tagList(output_list())
    })

    observe({
      p <- plotly_spatialDimPlot(
        object,
        point.size = point.size.dimplot,
        cells.highlight = vals$keepcells,
        nohighlight.alpha = nohighlight.alpha
      ) %>%
        plotly::layout(dragmode = "select") %>%
        plotly::event_register("plotly_selecting")
      new_env1 <- new.env()
      new_env1$p <- p
      output$Celltype <- plotly::renderPlotly(expr = {p},env = new_env1)

      pls <- plotly_spatialLRpairPlot(
        object = object,
        features = features,
        pairLR.use = pairLR.use,
        cells.highlight = vals$keepcells,
        nohighlight.alpha = nohighlight.alpha,
        point.size = point.size.LRplot,
        reversescale = vals$rev,
        enriched.only = F,
        do.group = F,
        do.plot=T
      )
      for (i in seq_len(n.pls)) {
        new_env2 <- new.env()
        new_env2$idx <- i
        new_env2$pls <- pls
        output[[stringr::str_c("LRExpr", i)]] <- plotly::renderPlotly(
          expr = {
            pls[[idx]]%>%
              plotly::layout(dragmode = "select") %>%
              plotly::event_register("plotly_selecting")
          },
          env = new_env2
        )
      }
    })

    output$selecting <- shiny::renderPrint({
      d <- plotly::event_data("plotly_selecting")
      if (is.null(d))
        "Please select cells in Spatial Dim Plot."
      else {
        cells.highlight <- d$customdata
        c_ <- cells.highlight %>% paste0("\'",.,"\',")
        for (i in 1:length(c_)) {
          if (i%%6==0) {
            cat(c_[[i]],sep = "\n")
          }else{
            cat(c_[[i]],sep = "")
          }
        }
      }

    })
    output$selecting2 <- shiny::renderPrint({
      c_ <- vals$keepcells %>% paste0("\'",.,"\',")
      for (i in 1:length(c_)) {
        if (i%%6==0) {
          cat(c_[[i]],sep = "\n")
        }else{
          cat(c_[[i]],sep = "")
        }
      }
    })

    # Toggle points that are clicked/selected
    shiny::observeEvent(plotly::event_data("plotly_click"), {
      d <- plotly::event_data("plotly_click")
      if (!is.null(d)){
        if(length(vals$keepcells)>=vals$Ncells){
          vals$keepcells <- d$customdata
        } else {
          vals$keepcells <- c(vals$keepcells,d$customdata)
        }
      }
    })
    shiny::observeEvent(plotly::event_data("plotly_selecting"), {
      d <- plotly::event_data("plotly_selecting")
      if (!is.null(d)){
        shiny::observeEvent(input$update, {
          vals$keepcells <- d$customdata
        })
      }
    })

    # Reset all points
    shiny::observeEvent(input$reset, {
      vals$keepcells = rownames(object@images$coordinates)
    })

    # reverse colorscale
    shiny::observeEvent(input$rev, {
      vals$rev = !vals$rev
    })
    # close
    shiny::observe({
      if(input$done > 0){
        shiny::stopApp(0)
      }
    })

  }

  shiny::shinyApp(ui, server, options = list(port = port))
}


#' @title interactive spatialDimPlot
#'
#' @param object a SpatialCellChat object
#' @param color.use defining the color for each cell group
#' @param group.by name of one or more metadata columns to group (color) cells by
#' @param group.highlight cell groups to be highlighted
#' @param cells.highlight cells/spots to be highlighted
#' @param sources.use a vector giving the index or the name of source cell groups of interest
#' @param targets.use a vector giving the index or the name of target cell groups of interest
#' @param idents.use a vector giving the index or the name of cell groups of interest
#' @param nohighlight.alpha the opacity of non-highlighted cells/spots
#' @param title.name set the title of this plot
#' @param scene.name name of the Plotly 3D scene
#' @param point.size the size of spots
#' @param show.legend whether show each figure legend
#' @param legend.group the parameter `legendgroup` in `plotly::plot_ly` and `plotly::add_trace`. Name of the legend group for traces, used to combine multiple traces under the same legend.
#' @param legend.group.title the parameter `legendgrouptitle` in `plotly::plot_ly` and `plotly::add_trace`. The title text displayed for the legend group.
#' @param three.dim.z the fixed z-axis coordinate in 3D mode
#' @param method plotting mode, one of `"2d`, or `"3d"` (Plotly)
#'
#' @return a plotly plot
#' @export
plotly_spatialDimPlot <- function (
    object,
    group.by = NULL,
    group.highlight = NULL,
    cells.highlight = NULL,
    sources.use = NULL,
    targets.use = NULL,
    idents.use = NULL,
    method = c("2d", "3d"),
    color.use = NULL,
    nohighlight.alpha = 0.65,
    title.name = NULL,
    scene.name = "scene",
    point.size = 2,
    show.legend = T,
    legend.group.title = "labels",
    three.dim.z = 0
){
  method = match.arg(method)
  coordinates <- object@images$coordinates
  if (ncol(coordinates) == 2) {
    colnames(coordinates) <- c("x_cent", "y_cent")
    temp_coordinates = coordinates
    coordinates[, 1] = temp_coordinates[, 2]
    coordinates[, 2] = temp_coordinates[, 1]
  }
  else {
    stop("Please check the input 'coordinates' and make sure it is a two column matrix.")
  }

  if (is.null(group.by)) {
    labels <- object@idents
  }
  else {
    labels <-  object@meta[, group.by]
    labels <- factor(labels)
  }
  cells.level <- levels(labels)

  coordinates$spot_labels <- labels
  coordinates$cell_id <- BiocGenerics::rownames(coordinates)

  if (!is.null(idents.use)) {
    if (is.numeric(idents.use)) {
      idents.use <- cells.level[idents.use]
    }
    group <- rep("Others", length(labels))
    group[(labels %in% idents.use)] <- idents.use
    group = factor(group, levels = c(idents.use, "Others"))

    if (is.null(color.use)) {
      color.use.all <- scPalette(nlevels(labels))

      # get the first 2 colors, then convert the 3th color to "grey90"
      color.use <- color.use.all[match(c(idents.use), levels(labels))]
      color.use[nlevels(group)] <- "grey90"
    }
    labels <- group
  }

  if (is.null(sources.use) & is.null(targets.use)) {
    if (is.null(color.use)) {
      color.use <- scPalette(nlevels(labels))
    }
  } else {
    if (is.numeric(sources.use)) {
      sources.use <- cells.level[sources.use]
    }
    if (is.numeric(targets.use)) {
      targets.use <- cells.level[targets.use]
    }
    group <- rep("Others", length(labels))
    group[(labels %in% sources.use)] <- sources.use
    group[(labels %in% targets.use)] <- targets.use
    group = factor(group, levels = c(sources.use, targets.use,
                                     "Others"))
    if (is.null(color.use)) {
      color.use.all <- scPalette(nlevels(labels))
      # get the first 2 colors, then convert the 3th color to "grey
      color.use <- color.use.all[match(c(sources.use,targets.use), levels(labels))]
      color.use[nlevels(group)] <- "grey90"
    }
    labels <- group
  }
  coordinates$color_labels <- labels


  if (!is.null(group.highlight) & !is.null(cells.highlight)) {
    stop("Please don't input `group.highlight` or `cells.highlight` simultaneously.")
  }

  if (!is.null(group.highlight)) {
    if (all(group.highlight %in% levels(coordinates$spot_labels))) {
      coordinates <- dplyr::mutate(
        coordinates,
        spot_highlight = if_else(.data[["spot_labels"]] %in% group.highlight, F, T)
      )
    } else{
      stop("`group.highlight` not in cell groups you used! Please check your input.")
    }

  }
  if (!is.null(cells.highlight)) {
    if(is.numeric(cells.highlight)){
      coordinates$spot_highlight <- dplyr::if_else(seq_len(NROW(coordinates)) %in% cells.highlight,F,T)
    } else{
      if (all(cells.highlight %in% coordinates$cell_id)) {
        coordinates <- dplyr::mutate(
          coordinates,
          spot_highlight = if_else(.data[["cell_id"]] %in% cells.highlight, F, T)
        )
      } else{
        stop("`cells.highlight` has some wrong cell_ids! Please check your input.")
      }
    }
  }
  if(is.null(group.highlight) & is.null(cells.highlight)){
    coordinates <- dplyr::mutate(coordinates, spot_highlight = F)
  }

  if (method == "3d") {
    pl <- plotly::plot_ly(
      data = coordinates[!coordinates$spot_highlight, ],
      x = ~ x_cent,
      y = ~ y_cent,
      z = three.dim.z,
      scene = scene.name,
      color =  ~ color_labels,
      colors = color.use,
      opacity = 1,
      legendgroup = legend.group.title,
      legendgrouptitle = list(text = legend.group.title),
      type = "scatter3d",
      mode = "markers",
      hoverinfo = 'text',
      text = ~ paste(
        "<i><b>cell type: ",
        spot_labels,
        "</b></i>",
        '<br> x:',
        x_cent,
        ", y:",
        y_cent,
        '<br> cell id:',
        cell_id
      ),
      marker = list(size = point.size,
                    symbol = "0")
    ) %>% plotly::add_trace(
      data = coordinates[coordinates$spot_highlight, ],
      x = ~ x_cent,
      y = ~ y_cent,
      z = three.dim.z,
      scene = scene.name,
      color =  ~ color_labels,
      colors = color.use,
      opacity = nohighlight.alpha,
      legendgroup = legend.group.title,
      legendgrouptitle = list(text = legend.group.title),
      type = "scatter3d",
      mode = "markers",
      hoverinfo = 'text',
      customdata = ~ cell_id,
      text = ~ paste(
        "<i><b>cell type: ",
        spot_labels,
        "</b></i>",
        '<br> x:',
        x_cent,
        ", y:",
        y_cent,
        '<br> cell id:',
        cell_id
      ),
      marker = list(size = point.size,
                    symbol = "0")
    ) %>%
      plotly::layout(
        title = title.name,
        showlegend = show.legend,
        grid = generate_grid_nrows(1),
        scene = generate_custom_scenes3d(
          row = 0,
          projection = "orthographic",
          zoom = 2.5
        ),
        legend = custom_legend
      )

  } else if (method == "2d") {
    # pl.key <- plotly::highlight_key(coordinates,  ~ cell_id)
    pl <- plotly::plot_ly(
      data = coordinates[!coordinates$spot_highlight, ],
      x = ~ x_cent,
      y = ~ y_cent,
      color =  ~ color_labels,
      colors = color.use,
      opacity = 1,
      legendgroup = legend.group.title,
      legendgrouptitle = list(text = legend.group.title),
      hoverinfo = 'text',
      text = ~ paste(
        "<i><b>cell type: ",
        spot_labels,
        "</b></i>",
        '<br> x:',
        x_cent,
        ", y:",
        y_cent,
        '<br>cell id:',
        cell_id
      ),
      customdata = ~ cell_id,
      marker = list(size = point.size+4,
                    symbol = "0"),
      type = "scatter",
      mode = 'markers'
    ) %>% plotly::add_trace(
      data = coordinates[coordinates$spot_highlight, ],
      x = ~ x_cent,
      y = ~ y_cent,
      color =  ~ color_labels,
      colors = color.use,
      opacity = nohighlight.alpha,
      legendgroup = legend.group.title,
      legendgrouptitle = list(text = legend.group.title),
      hoverinfo = 'text',
      text = ~ paste(
        "<i><b>cell type: ",
        spot_labels,
        "</b></i>",
        '<br> x:',
        x_cent,
        ", y:",
        y_cent,
        '<br>cell id:',
        cell_id
      ),
      customdata = ~ cell_id,
      marker = list(size = point.size+4,
                    symbol = "0"),
      type = "scatter",
      mode = 'markers'
    ) %>%
      plotly::layout(
        title = title.name,
        yaxis = list(
          autorange = "reversed",
          title = "",
          showgrid = FALSE,
          ticks = "",
          # ticktext = "",
          tickvals = "",
          scaleanchor = "x",
          tickformat = "%0f",
          showline = FALSE,
          visible = FALSE
        ),
        xaxis = list(
          title = "",
          showgrid = FALSE,
          ticks = "",
          # ticktext = "",
          tickvals = "",
          tickformat = "%0f",
          showline = FALSE,
          visible = FALSE
        ),
        legend = custom_legend
      )
    # %>%plotly::highlight(
    #   # on = "plotly_click",
    #   # off = "plotly_doubleclick",
    #   selectize = TRUE,
    #   dynamic = TRUE,
    #   color = RColorBrewer::brewer.pal(length(cells.level), "Dark2"),
    #   persistent = TRUE,
    #   selected = plotly::attrs_selected(opacity = nohighlight.alpha)
    # )


  } else {
    stop(
      "The parameter \"method\" in this function only supports for \"normal\",\"3d\"!\n Please check your input."
    )
  }
  return(pl)
}


#' @title themeFeaturePlot
#' @description Defines a reusable ggplot2 theme for feature plots
#' @param legend.text.size the text size on the legend
#' @import ggplot2
#' @return a ggplot2 theme object
#' @export
themeFeaturePlot <- function(legend.text.size = 8) {
  theme(legend.position = "right") +
    theme(
      legend.title = element_blank(),
      legend.text = element_text(size = legend.text.size),
      legend.key.size = grid::unit(0.15, "inches")
    )  +
    theme(
      panel.background = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank()
    ) + xlab(NULL) + ylab(NULL) +
    theme(legend.key = element_blank())
}


#' @title spatialLeePlot
#' @param object SpatialCellChat object
#' @param group.by1 a character string specifying a column name in `object@meta` used to define the grouping variable for y.
#' If NULL, `object@idents` is used by default. Alternatively, a vector of labels can be provided directly
#' @param group.by2 a character string specifying a column name in `object@meta` used to define the grouping variable for x.
#' Alternatively, a vector of labels can be provided directly
#' @param lw.type select the type of distance weight matrix based on the communication mode, one of "Itself","contact.range", or "interaction.range"
#' @param lw spatial distance weight matrix. If `lw = NULL`, it will be automatically computed based on `lw.type`. Default is `NULL`
#' @param interaction.range numeric. By default 250um. The maximum interaction/diffusion length of ligands (Unit: microns)
#' This hard threshold is used to filter out the connections between spatially distant individual cells
#' @param contact.range numeric. By default 10um. The interaction range (Unit: microns) to restrict the contact-dependent signaling.
#' For spatial transcriptomics in a single-cell resolution, `contact.range` is approximately equal to the estimated cell diameter (i.e., the cell center-to-center distance), which means that contact-dependent and juxtacrine signaling can only happens when the two cells are contact to each other
#' @inheritParams plotStatistics_Lee
#'
#' @import ComplexHeatmap
#' @import grid
#'
#' @return a ggplot2 object
#' @export
spatialLeePlot <- function(
    object,
    group.by2,
    group.by1 = NULL,
    x.levels = NULL,
    y.levels = NULL,
    lw.type = c("Itself","contact.range","interaction.range"),
    lw = NULL,
    interaction.range = 250,
    contact.range = 10,
    cutoff = TRUE,
    cutoff.Lee = 0,
    type = "heatmap",
    color.use = NULL,
    color.heatmap = "Reds",
    reverse.colorscale = F,
    n.colors = 6,
    dot.size = c(1, 6),
    angle.x = 45,
    xlabel = NULL,
    ylabel = NULL,
    title.name = NULL,
    legend.title = NULL,
    font.size = 10,
    font.size.title = 10,
    x.lab.rot = 45,
    row.show = NULL,
    col.show = NULL
){
  # get coordinates
  coords <- object@images$coordinates
  # check coordinates
  if (ncol(coords) == 2) {
    colnames(coords) <- c("x_cent", "y_cent")
    temp_coord = coords
    coords[, 1] = temp_coord[, 2]
    coords[, 2] = temp_coord[, 1]
  } else {
    stop("Please check the input 'coordinates' and make sure it is a two column matrix.")
  }

  if (is.null(group.by1)) {
    y.labels <- object@idents
  } else {
    if(is.character(group.by1)){
      y.labels = object@meta[,group.by1]
      # avoid loss of factor levels
      if(!is.factor(y.labels)) y.labels <- factor(y.labels)
    } else {
      y.labels = group.by1
    }
  }

  if(is.character(group.by2)){
    x.labels = object@meta[,group.by2]

  } else {
    x.labels = group.by2
  }

  if(is.null(lw)) {
    lw.type <- match.arg(lw.type)
    if(lw.type == "Itself"){
      print("lw.type = Itself")
      plotStatistics_Lee(
        x.labels = x.labels,
        y.labels = y.labels,
        coords = coords,
        lw = NULL,
        x.levels = x.levels,
        y.levels = y.levels,
        cutoff = cutoff,
        cutoff.Lee = cutoff.Lee,
        type = type,
        color.use = color.use,
        color.heatmap = color.heatmap,
        reverse.colorscale = reverse.colorscale,
        n.colors = n.colors,
        dot.size = dot.size,
        angle.x = angle.x,
        xlabel = xlabel,
        ylabel = ylabel,
        title.name = title.name,
        legend.title = legend.title,
        font.size = font.size,
        font.size.title = font.size.title,
        x.lab.rot = x.lab.rot,
        row.show = row.show,
        col.show = col.show
      )
    } else {
      if(!is.null(object@images[["result.computeCellDistance"]])){
        resCellDistance <- object@images[["result.computeCellDistance"]]
      } else {
        resCellDistance <- SpatialCellChat::computeCellDistance(
          coordinates = coords,
          interaction.range = interaction.range,
          contact.range = contact.range,
          ratio = object@images[["spatial.factors"]]$ratio,
          tol = object@images[["spatial.factors"]]$tol
        )
      }

      if(lw.type == "contact.range"){
        print("lw.type = contact.range")
        weight.Mat <- resCellDistance[["adj.contact"]]
        lw <- spdep::mat2listw(weight.Mat,style = "B") # all weights are 1.
      } else if (lw.type == "interaction.range"){
        print("lw.type = interaction.range")
        d.spatial <- resCellDistance[["d.spatial"]]
        P.spatial <- SpatialCellChat::createPspatialFrom_dspatial(d.spatial,distance.use = T)
        Matrix::diag(P.spatial) <- 1
        lw <- spdep::mat2listw(P.spatial,style = "W") # all weights.
      }
      # View(lw)
      plotStatistics_Lee(
        x.labels = x.labels,
        y.labels = y.labels,
        coords = NULL,
        lw = lw,
        x.levels = x.levels,
        y.levels = y.levels,
        cutoff = cutoff,
        cutoff.Lee = cutoff.Lee,
        type = type,
        color.use = color.use,
        color.heatmap = color.heatmap,
        reverse.colorscale = reverse.colorscale,
        n.colors = n.colors,
        dot.size = dot.size,
        angle.x = angle.x,
        xlabel = xlabel,
        ylabel = ylabel,
        title.name = title.name,
        legend.title = legend.title,
        font.size = font.size,
        font.size.title = font.size.title,
        x.lab.rot = x.lab.rot,
        row.show = row.show,
        col.show = col.show
      )
    }
  } else {
    print("lw is given")
    plotStatistics_Lee(
      x.labels = x.labels,
      y.labels = y.labels,
      coords = NULL,
      lw = lw,
      x.levels = x.levels,
      y.levels = y.levels,
      cutoff = cutoff,
      cutoff.Lee = cutoff.Lee,
      type = type,
      color.use = color.use,
      color.heatmap = color.heatmap,
      reverse.colorscale = reverse.colorscale,
      n.colors = n.colors,
      dot.size = dot.size,
      angle.x = angle.x,
      xlabel = xlabel,
      ylabel = ylabel,
      title.name = title.name,
      legend.title = legend.title,
      font.size = font.size,
      font.size.title = font.size.title,
      x.lab.rot = x.lab.rot,
      row.show = row.show,
      col.show = col.show
    )
  }
}


#' @title plotStatistics_Lee
#' @param x.labels labels defining the categories for the x variables. Can be a factor/vector, a numeric vector, or a matrix
#' @param y.labels labels defining the categories for the y variables. Can be a factor/vector, a numeric vector, or a matrix
#' @param coords a two-column matrix or data frame specifying spatial coordinates
#' @param lw spatial distance weight matrix
#' @param x.levels a character vector specifying the order of x categories
#' @param y.levels a character vector specifying the order of y categories
#' @param cutoff whether to apply a threshold to Lee's statistics
#' @param cutoff.Lee numeric threshold for Lee's statistic. Values less than or equal to this cutoff are removed when `cutoff = TRUE`
#' @param type the type of the output, one of "data", "dot" or "heatmap"
#' @param color.use defining the color for each group
#' @param color.heatmap A character string or vector indicating the colormap option to use. It can be the avaibale color palette in brewer.pal() or viridis_pal() (e.g., "Spectral","viridis")
#' @param reverse.colorscale  whether to reverse the colorbar
#' @param n.colors number of basic colors to generate from color palette
#' @param dot.size a range defining the size of the symbol
#' @param angle.x angle for x-axis text rotation
#' @param xlabel label of x-axis
#' @param ylabel label of y-axis
#' @param title.name title name
#' @param legend.title name of legend
#' @param font.size font size in the plot
#' @param font.size.title font size of the title
#' @param x.lab.rot rotation angle of x-axis tick labels
#' @param row.show,col.show a vector giving the index or the name of row or columns to show in the heatmap
#'
#' @import ComplexHeatmap
#' @import grid
#' @return a data frame, a ggplot2 object or a Heatmap-class object
#' @export
plotStatistics_Lee <- function(
    x.labels,
    y.labels,
    coords = NULL,
    lw = NULL,
    x.levels = NULL,
    y.levels = NULL,
    cutoff = TRUE,
    cutoff.Lee = 0,
    type = "heatmap",
    color.use = NULL,
    color.heatmap = "Reds",
    reverse.colorscale = F,
    n.colors = 6,
    dot.size = c(1, 6),
    angle.x = 45,
    xlabel = NULL,
    ylabel = NULL,
    title.name = NULL,
    legend.title = NULL,
    font.size = 10,
    font.size.title = 10,
    x.lab.rot = 45,
    row.show = NULL,
    col.show = NULL
) {
  if(inherits(x = x.labels,what = c("matrix", "Matrix", "dgCMatrix"))){
    labels1.Mat <- as.matrix(x.labels)
    x.labels <- factor(colnames(labels1.Mat),levels = colnames(labels1.Mat))
  } else {
    if(is.numeric(x.labels)){
      labels1.Mat <- matrix(x.labels,ncol = 1)
      colnames(labels1.Mat) <- xlabel%||%"variable"
      x.labels <- factor(colnames(labels1.Mat),levels = colnames(labels1.Mat))
      type = "dot"
    } else {
      if(!is.factor(x.labels)){
        x.labels <- factor(x.labels)
      }
      labels1.Mat <- model.matrix(~x.labels-1)
    }
  }


  if(inherits(x = y.labels,what = c("matrix", "Matrix", "dgCMatrix"))){
    labels2.Mat <- as.matrix(y.labels)
    y.labels <- factor(colnames(labels2.Mat),levels = colnames(labels2.Mat))
  } else {
    if(is.numeric(y.labels)){
      labels2.Mat <- matrix(y.labels,ncol = 1)
      colnames(labels2.Mat) <- ylabel%||%"variable"
      y.labels <- factor(colnames(labels2.Mat),levels = colnames(labels2.Mat))
      type = "dot"
    } else {
      if(!is.factor(y.labels)){
        y.labels <- factor(y.labels)
      }
      labels2.Mat <- model.matrix(~y.labels-1)
    }
  }

  if(is.null(coords)){
    if(is.null(lw)) stop("Please provide either coords or lw for Lee statistics.")
  } else {
    if(!is.null(lw)) stop("Please don't provide both coords and lw at the same time.")
    if(NCOL(coords)!=2) stop("Please check your coordinates. It must be 2 columns.")
    colnames(coords) <- c("x", "y")

    weight.Mat <- Matrix::sparseMatrix(
      i = seq_len(NROW(coords)),
      j = seq_len(NROW(coords)),
      x = rep.int(1,times = NROW(coords)),
      dims = c(NROW(coords),NROW(coords)),
      index1 = T,
      dimnames = NULL
    )
    lw <- spdep::mat2listw(weight.Mat,style = "B")
  }


  dfLee <- purrr::map_dfc(
    .x = seq_len(NCOL(labels1.Mat)),
    .f = function(L1) {
      label1 <- labels1.Mat[, L1, drop = T]
      cooccurrenceL1 <- purrr::map_dbl(
        .x = seq_len(NCOL(labels2.Mat)),
        .f = function(L2) {
          label2 <- labels2.Mat[, L2, drop = T]
          Lee.S <- spdep::lee(
            label1,
            label2,
            listw = lw ,
            n = length(lw$neighbours),
            zero.policy = TRUE
          )
          return(Lee.S$L)
        }
      )
      cooccurrenceL1 <- data.frame(cooccurrenceL1)
      colnames(cooccurrenceL1) <- colnames(labels1.Mat)[[L1]]
      return(cooccurrenceL1)
    }
  )

  colnames(dfLee) <- levels(x.labels)
  rownames(dfLee) <- levels(y.labels)
  # View(dfLee);cat(x.labels);cat("\n");cat(y.labels)

  if (cutoff) {
    dfLee[dfLee<=cutoff.Lee] <- NA
  }

  mat <- Matrix::t(Matrix::as.matrix(dfLee))

  if (!is.null(x.levels)) {
    x.levels <- x.levels[x.levels %in% colnames(mat)]
    mat <-
      mat[, order(factor(colnames(mat), levels = x.levels)), drop = FALSE]
  }
  if (!is.null(y.levels)) {
    y.levels <- y.levels[y.levels %in% rownames(mat)]
    mat <-
      mat[order(factor(rownames(mat), levels = y.levels)), , drop = FALSE]
  }

  if (is.null(color.use)) {
    color.use <- scPalette(NCOL(mat))
  }
  names(color.use) <- colnames(mat)


  if (type == "data") {
    df <- as.data.frame(Matrix::t(mat))
    return(df)
  } else if (type == "dot") {
    df <- as.data.frame(as.table(Matrix::t(mat[nrow(mat):1,])))
    colnames(df) <- c("x", "y", "proportion")
    color.heatmap.use = grDevices::colorRampPalette((RColorBrewer::brewer.pal(n = n.colors, name = color.heatmap)))(100)
    if(reverse.colorscale){
      color.heatmap.use <- rev(color.heatmap.use)
    }
    gg <-
      ggplot(df, aes(x, y, size = proportion, color = proportion)) +
      geom_point() +
      scale_size_continuous(range = dot.size) +
      theme_linedraw() +
      labs(size = "cooccurrence") +
      scale_color_gradientn(colors = color.heatmap.use) +
      theme(legend.key.height = grid::unit(0.15, "in")) +
      guides(color = guide_colourbar(barwidth = 0.5, title = "Cooccurrence"))
    # gg <- gg + guides(color = guide_colorbar(barwidth = legend.width, title = "Scaled expression"),size = guide_legend(title = 'Percent expressed'))
    gg <- gg + theme(
      text = element_text(size = 10),
      axis.text.x = element_text(angle = x.lab.rot, hjust = 1,size=font.size-2),
      axis.text.y = element_text(angle = 0, hjust = 1,size=font.size),
      axis.title.x = element_blank(),
      axis.title.y = element_blank()
    ) +
      theme(
        axis.line.x = element_line(linewidth = 0.25),
        axis.line.y = element_line(linewidth = 0.25)
      ) +
      theme(panel.grid.major = element_line(colour = "grey90", linewidth = (0.1)))
    return(gg)

  } else if (type == "heatmap") {

    if (!is.null(row.show)) {
      mat <- mat[row.show,,drop=F]
    }
    if (!is.null(col.show)) {
      mat <- mat[, col.show,drop=F]
      color.use <- color.use[col.show]
    }
    color.heatmap.use = grDevices::colorRampPalette((RColorBrewer::brewer.pal(n = n.colors, name = color.heatmap)))(100)

    if(reverse.colorscale){
      color.heatmap.use <- rev(color.heatmap.use)
    }

    df <-
      data.frame(group = colnames(mat))
    rownames(df) <- colnames(mat)
    col_annotation <-
      ComplexHeatmap::HeatmapAnnotation(
        df = df,
        col = list(group = color.use),
        which = "column",
        show_legend = FALSE,
        show_annotation_name = FALSE,
        simple_anno_size = grid::unit(0.2, "cm")
      )

    # ha1 = rowAnnotation(Strength = anno_barplot(rowSums(abs(mat)), border = FALSE,gp = grid::gpar(fill = color.use, col=color.use)), show_annotation_name = FALSE)
    # mat[mat == 0] <- NA
    # color.heatmap.use = c("white", color.heatmap.use)

    # View(mat)
    ht1 = ComplexHeatmap::Heatmap(
      mat,
      col = color.heatmap.use,
      na_col = "white",
      name = "Proportion",
      bottom_annotation = col_annotation,
      cluster_rows = F,
      cluster_columns = F,
      # clustering_distance_rows = clustering_distance_rows,
      column_order = colnames(mat),
      row_names_side = "left",
      row_names_rot = 0,
      row_names_gp = grid::gpar(fontsize = font.size),
      column_names_gp = grid::gpar(fontsize = font.size),
      # width = grid::unit(width, "cm"), height = grid::unit(height, "cm"),
      row_title = xlabel,
      row_title_gp = grid::gpar(fontsize = font.size.title),
      column_title = paste0(title.name,"Spatial co-occurrence ", "(Lee statistic)"),
      column_title_gp = grid::gpar(fontsize = font.size.title),
      column_names_rot = x.lab.rot,
      heatmap_legend_param = list(
        title = "Cooccurrence",
        title_gp = grid::gpar(fontsize = 10, fontface = "plain"),
        title_position = "leftcenter-rot",
        border = NA,
        legend_height = grid::unit(15, "mm"),
        labels_gp = grid::gpar(fontsize = 10),
        grid_width = grid::unit(2, "mm")
      ),
      cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
        if (!is.na(mat[i, j])) {
          grid::grid.text(round(mat[i, j],digits = 2), x, y, gp = grid::gpar(col = "black", fontsize = 7))
        }
      }
    )

    return(ht1)
  }
}


#' @title spatialGiPlot
#' @param object SpatialCellChat object
#' @param slot.name the slot name of object that is used for analysis
#' @param signaling.name alternative signaling pathway name to show on the plot
#' @param do.binary whether to plot the scores in a binary mode
#' @param group.highlight cell groups to be highlighted
#' @param measure "indeg" or "outdeg". "indeg": in-degree of the communication network. "outdeg": out-degree of the communication network
#' @param measure.name the measure names to show, which should correspond to the defined `measure`
#' @param plot.title title of the plot
#' @param n the number of nearest neighbours for each spot, used to compute the Gi statistic
#' @param method plotting mode. one of "scoring", "normal", "3d", "density" or "contour"
#' @param gstar logical; If true, compute the Gi* statistic, else compute the Gi statistic
#' @param p.value significance threshold for p-values
#' @param three.dim.z the fixed z-axis coordinate in 3D mode
#' @param point.size the size of spots
#' @param nohighlight.alpha the opacity of non-highlighted cells/spots
#' @param show.colorscale whether to display the colorbar
#' @param return.object return either SpatialCellChat object or ggplot2 object
#' @param show.plot if `return.object = TRUE`, whether to display the plot
#'
#' @return SpatialCellChat object, a ggplot2 object or a plotly object
#' @export
spatialGiPlot <- function(
    object,
    slot.name = "netP",
    signaling.name=NULL,
    do.binary = T,
    group.highlight = NULL,
    measure = c("indeg", "outdeg"),
    measure.name = c("incoming", "outgoing"),
    plot.title = NULL,
    n = 12,
    method = "normal",
    gstar = FALSE,
    p.value = NULL,
    three.dim.z = 0,
    point.size = 2,
    nohighlight.alpha = 0.65,
    show.colorscale = T,
    return.object = T,
    show.plot = T
) {
  for (i in 1:length(measure)) {
    if (measure[i] == "outdeg") {
      measure.name[i] = "outgoing"
    } else if (measure[i] == "indeg") {
      measure.name[i] = "incoming"
    } else{
      measure.name[i] = measure[i]
    }
  }

  # get coordinates
  coords <- object@images$coordinates
  # check coordinates
  if (ncol(coords) == 2) {
    colnames(coords) <- c("x_cent", "y_cent")
    temp_coord = coords
    coords[, 1] = temp_coord[, 2]
    coords[, 2] = temp_coord[, 1]
  } else {
    stop("Please check the input 'coordinates' and make sure it is a two column matrix.")
  }


  # check the slot's existence
  if ( is.null(methods::slot(object, slot.name)$centr.cell) ) {
    stop(cli.symbol(2),"Please run `netAnalysis_computeCentrality` to compute the network centrality scores! ")
  }

  numFeature = length(measure)

  # get the corresponding scores of each spot according to 'measure'
  # centr's dim: centr.measure x nCell x nLR/nPathway
  centr <- slot(object, slot.name)$centr.cell[measure, , ,drop=F]
  if (!is.null(signaling.name)) {
    signaling.use <- dimnames(centr)[[3]]
    if (is.numeric(signaling.name)) signaling.name <- signaling.use[signaling.name]
    centr <- centr[ , ,signaling.name,drop=F]
  }

  if (do.binary) {
    centr[centr>0] <- 1
  } else {
    NULL
  }

  # simplify = F => ordered list; simplify = T => array
  centr.sum <- apply(centr,c(1,2),sum,simplify = T) %>% Matrix::t()

  # combine the coordinates and communication scores by columns
  df_score <- cbind(coords, centr.sum)
  colnames(df_score) <- c("x_cent", "y_cent", measure)

  # add idents info
  spot_labels <- object@idents

  # highlight the selected groups
  if (is.null(group.highlight)) {
    spot_highlight <- rep(F, length(spot_labels))
  } else{
    if (all(group.highlight %in% levels(spot_labels))) {
      spot_highlight <- sapply(
        spot_labels,
        FUN = function(x) {
          if (x %in% group.highlight) {
            return(F)
          } else{
            return(T)
          }
        },
        simplify = T
      )

    } else{
      stop("`group.highlight` not %in% object@idents!")
    }
  }


  # no more than 2
  plot_ncol <- min(numFeature, 2)

  # df <- data.frame(x = coords[, 1], y = coords[, 2])
  gg <- vector("list", numFeature)

  if (numFeature == 1) {
    colorbar.y <- c(0.5)
  } else{
    colorbar.y <- c(0.75, 0.25)
  }

  # plot...
  for (i in seq_along(measure)) {
    if(is.null(plot.title)){
      feature.name <-
        paste0("Hotspots of ",signaling.name,"\'s ",measure.name[i], " scoring")
    } else {
      feature.name <- plot.title
    }

    # only 3 cols, because knn2nb can only accept 3 cols
    df_sub <- df_score[c("x_cent", "y_cent", measure[i])]

    g <- df_sub %>% plotStatistics_Gi(
      df_score = .,
      group.by = spot_labels,
      group.highlight = spot_highlight,
      n = n,
      method = method,
      gstar = gstar,
      p.value = p.value,
      plot.title = feature.name,
      point.size = point.size,
      nohighlight.alpha = nohighlight.alpha,
      show.colorscale = show.colorscale,
      legend.group = paste0("Gi", i),
      legend.group.title = "Getis-Ord Gi",
      scene.name = paste0("scene", ifelse(i == 1, "", i)),
      colorbar.y = colorbar.y[[i]],
      three.dim.z = three.dim.z
    )

    g.info <- df_sub %>% plotStatistics_Gi(
      df_score = .,
      group.by = spot_labels,
      group.highlight = spot_highlight,
      n = n,
      method = "computation",
      gstar = gstar,
      p.value = p.value,
      plot.title = feature.name,
      point.size = point.size,
      nohighlight.alpha = nohighlight.alpha,
      show.colorscale = show.colorscale,
      legend.group = paste0("Gi", i),
      legend.group.title = "Getis-Ord Gi",
      scene.name = paste0("scene", ifelse(i == 1, "", i)),
      colorbar.y = colorbar.y[[i]],
      three.dim.z = three.dim.z
    )
    Gi.feature.name <- paste0("Getis-Ord Gi ",feature.name) %>% make.names()
    object@meta[[Gi.feature.name]] <- g.info[["statistic_G"]]
    Gip.feature.name <- paste0("Getis-Ord Gi P ",feature.name) %>% make.names()
    object@meta[[Gip.feature.name]] <- g.info[["Pr(z != E(Gi))"]]

    gg[[i]] <- g
  }

  if(return.object){
    if(show.plot){
      if (class(gg[[1]])[[1]] %in% c('ggplot2::ggplot', 'ggplot', 'ggplot2::gg', 'gg')) {
        gg <- patchwork::wrap_plots(gg, ncol = plot_ncol)
        # do plot
        show(gg)
      } else{
        # combine the plotly plots
        if (length(gg) > 1) {
          gg <- plotly::subplot(gg[[1]], gg[[2]], margin = 0) %>%
            plotly::layout(
              grid = generate_grid_nrows(2),
              scene = generate_custom_scenes3d(row = 0),
              scene2 = generate_custom_scenes3d(row = 1),
              showlegend = F
            )
        } else{
          gg <- gg[[1]] %>%
            plotly::layout(
              showlegend = F,
              grid = generate_grid_nrows(1),
              scene = generate_custom_scenes3d(row = 0,
                                               z.space = 1)
            )
        }
        # do plot
        show(gg)
      }
    }
    return(object)
  } else { # return plot
    if (class(gg[[1]])[[1]] %in% c('ggplot2::ggplot', 'ggplot', 'ggplot2::gg', 'gg')) {
      gg <- patchwork::wrap_plots(gg, ncol = plot_ncol)
      return(gg)
    } else{
      # combine the plotly plots
      if (length(gg) > 1) {
        gg <- plotly::subplot(gg[[1]], gg[[2]], margin = 0) %>%
          plotly::layout(
            grid = generate_grid_nrows(2),
            scene = generate_custom_scenes3d(row = 0),
            scene2 = generate_custom_scenes3d(row = 1),
            showlegend = F
          )
      } else{
        gg <- gg[[1]] %>%
          plotly::layout(
            showlegend = F,
            grid = generate_grid_nrows(1),
            scene = generate_custom_scenes3d(row = 0,
                                             z.space = 1)
          )
      }
      return(gg)
    }
  }
}


#' @title plotStatistics_Gi
#' @param df_score a data frame with columns ordered as "x", "y", and "score"
#' @param group.by name of one or more metadata columns to group (color) cells by
#' @param group.highlight cell groups to be highlighted
#' @param n the number of nearest neighbours for each spot, used to compute the Gi statistic
#' @param method plotting mode. one of "scoring", "normal", "3d", "density" or "contour"
#' @param gstar logical; If true, compute the Gi* statistic, else compute the Gi statistic
#' @param p.value significance threshold for p-values
#' @param plot.title the title of the plot
#' @param point.size the size of spots
#' @param nohighlight.alpha the opacity of non-highlighted cells/spots
#' @param show.colorscale whether to display the colorbar
#' @param legend.group the parameter `legendgroup` in `plotly::plot_ly` and `plotly::add_trace`. Name of the legend group for traces, used to combine multiple traces under the same legend
#' @param legend.group.title the parameter `legendgrouptitle` in `plotly::plot_ly` and `plotly::add_trace`. The title text displayed for the legend group
#' @param scene.name name of the Plotly 3D scene
#' @param colorbar.y the relative position of the color bar along the y-axis
#' @param three.dim.z the fixed z-axis coordinate in 3D mode
#'
#' @return a ggplot2 object, a plotly object or a data frame
#' @export
plotStatistics_Gi <- function(df_score,
                              group.by,
                              group.highlight,
                              n,
                              method = "normal",
                              gstar = FALSE,
                              p.value = NULL,
                              plot.title = NULL,
                              point.size = 2,
                              nohighlight.alpha = 0.65,
                              show.colorscale = T,
                              legend.group = "Gi",
                              legend.group.title = "Gi\'s value",
                              scene.name = "scene1",
                              colorbar.y = 0.5,
                              three.dim.z = 0) {
  # rename the cols
  colnames(df_score) <- c("x", "y", "score")

  # set the title of plot
  if (is.null(plot.title)) {
    plot.title <- plot.title #NULL
  } else{
    plot.title <- labs(title = plot.title)
  }

  # calculate the statistic-Gi(Gi*)
  if (gstar) {
    knn_nb <-
      spdep::knearneigh(df_score[c("x", "y")],
                        k = n,
                        longlat = FALSE,
                        use_kd_tree = T) %>%
      spdep::knn2nb() %>%
      spdep::include.self()
    p.value.label <- "Pr(z != E(G*i))"
    gstar.label <- "Gi*"
  } else{
    knn_nb <-
      spdep::knearneigh(df_score[c("x", "y")],
                        k = n,
                        longlat = FALSE,
                        use_kd_tree = T) %>%
      spdep::knn2nb()
    p.value.label <- "Pr(z != E(Gi))"
    gstar.label <- "Gi"
  }
  G_knn <-
    spdep::localG(df_score$score, spdep::nb2listw(knn_nb, style = "B"))

  # update the score-data.frame, with the values of Gi, Gi's p.value and so on
  df_score$statistic_G <- c(G_knn)
  df_score <- cbind(df_score, attr(G_knn, "internals"))

  df_score$spot_labels <- group.by
  df_score$spot_highlight <- group.highlight
  my.guide = guide_colorbar(title = NULL, ticks = T, label = T, barwidth = 0.5)

  if (method == "scoring") {
    colors.use <- generate_custom_colorscale("Spectral",n.colors = 8,reversescale = T,plotly.format = F)

    g <- df_score %>% ggplot() +
      geom_point(aes(
        x = x,
        y = y,
        color = score
      ), size = point.size) +
      scale_colour_gradientn(colours = colors.use,guide = my.guide) +
      # scale_color_distiller(palette = "RdYlBu",direction = -1)+
      guides(alpha = guide_none()) +
      labs(title = NULL,colour = "score") + plot.title +
      theme(
        plot.title = element_text(hjust = 0.5),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank()
      ) +
      coord_fixed() + theme(aspect.ratio = 1, legend.key = element_blank()) + scale_y_reverse()
    return(g)

  } else if (method == "normal") {
    if (is.null(p.value)) {
      df_score <- df_score %>% mutate(point_shape = I(19))
    } else{
      df_score <-
        df_score %>% mutate(point_shape = if_else(.data[[p.value.label]] < p.value, I(19), I(1)))
    }

    colors.use <-
      generate_RedsBlues_colormap(df_score$statistic_G)

    g <- df_score %>% ggplot() +
      geom_point(aes(
        x = x,
        y = y,
        color = statistic_G,
        shape = point_shape
      ), size = point.size) +
      scale_colour_gradientn(colours = colors.use,guide = my.guide) +
      # scale_color_distiller(palette = "RdYlBu",direction = -1)+
      guides(alpha = guide_none()) +
      labs(title = paste0(
        "Spatial statistics with ",
        gstar.label,
        ifelse(
          is.null(p.value),
          "",
          paste0("(Normal, p.value=", p.value, ")")
        )
      ),colour = "Gi/Gi*") + plot.title +
      theme(
        plot.title = element_text(hjust = 0.5),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank()
      ) +
      coord_fixed() + theme(aspect.ratio = 1, legend.key = element_blank()) + scale_y_reverse()
    return(g)

  } else if (method == "3d") {
    if (is.null(p.value)) {
      df_score <- df_score %>% mutate(point_shape = I(19),
                                      p_value = "")
    } else{
      # shape: https://zhuanlan.zhihu.com/p/104143758
      df_score <-
        df_score %>% mutate(
          point_shape = if_else(.data[[p.value.label]] < p.value, I(19), I(1)),
          p_value = if_else(
            .data[[p.value.label]] < p.value,
            "Hypothesis testing: <b>Accept</b>",
            "Hypothesis testing: <b>Reject</b>"
          )
        )
    }

    colors.use <-
      generate_RedsBlues_colormap(df_score$statistic_G)

    df_score$z <- three.dim.z

    pl <- plotly::plot_ly(
      data = df_score[!df_score$spot_highlight, ],
      x = ~ x,
      y = ~ y,
      z = three.dim.z,
      scene = scene.name,
      color =  ~ statistic_G,
      colors = colors.use,
      opacity = 1,
      symbol = ~ point_shape,
      span = I(2),
      legendgroup = legend.group,
      legendgrouptitle = list(text = legend.group.title),
      type = "scatter3d",
      mode = "markers",
      hoverinfo = 'text',
      text = ~ paste(
        "<i><b>",
        plot.title,
        "</b></i>",
        '<br> x:',
        x,
        ", y:",
        y,
        ", z:",
        z,
        '<br>cell type: <b>',
        spot_labels,
        '</b><br>Gi/Gi* statistic: ',
        statistic_G,
        '<br>',
        p_value
      ),
      marker = list(size = point.size)
    ) %>% plotly::add_trace(
      data = df_score[df_score$spot_highlight, ],
      x = ~ x,
      y = ~ y,
      z = three.dim.z,
      scene = scene.name,
      color =  ~ statistic_G,
      colors = colors.use,
      symbol = ~ point_shape,
      span = I(2),
      opacity = nohighlight.alpha,
      legendgroup = legend.group,
      legendgrouptitle = list(text = legend.group.title),
      type = "scatter3d",
      mode = "markers",
      hoverinfo = 'text',
      text = ~ paste(
        "<i><b>",
        plot.title,
        "</b></i>",
        '<br> x:',
        x,
        ", y:",
        y,
        ", z:",
        z,
        '<br>cell type: <b>',
        spot_labels,
        '</b><br>Gi/Gi* statistic: ',
        statistic_G,
        '<br>',
        p_value
      ),
      marker = list(size = point.size)
    ) %>% plotly::plotly_build()

    # change the style of colorbar
    pl[["x"]][["data"]][[1]][["marker"]][["colorbar"]] <- list(
      orientation = "v",
      len = 1 / 4,
      y = colorbar.y,
      lenmode = "fraction",
      thickness = 8,
      # default's 1/2
      title = list(text = "Gi/Gi*")
    )
    # show the colorbar
    pl[["x"]][["data"]][[1]][["marker"]][["showscale"]] <-
      show.colorscale
    # not show other elements' scales
    pl[["x"]][["data"]][[2]][["marker"]][["showscale"]] <- F
    pl[["x"]][["data"]][[3]][["marker"]][["showscale"]] <- F

    return(pl)

  } else if (method == "density") {
    if (is.null(p.value)) {
      df_score <- df_score %>% mutate(point_shape = I(19))
    } else{
      df_score <-
        df_score %>% mutate(point_shape = if_else(.data[[p.value.label]] < p.value, I(19), I(1)))
    }

    # use ks method to calculate density!
    G_dens <-
      calculate_density(df_score[["statistic_G"]], df_score[c("x", "y")], method = "wkde")
    df_score$statistic_G.d <- G_dens

    colors.use <-
      generate_RedsBlues_colormap(df_score$statistic_G)

    g <- df_score %>% ggplot() +
      geom_point(aes(
        x = x,
        y = y,
        color = statistic_G.d,
        shape = point_shape
      ),
      size = point.size) +
      scale_colour_gradientn(colours = colors.use,guide = my.guide) +
      # scale_color_distiller(palette = "RdYlBu",direction = -1)+
      guides(alpha = guide_none()) +
      labs(title = paste0(
        "Spatial statistics with ",
        gstar.label,
        ifelse(
          is.null(p.value),
          "",
          paste0("(Density, p.value=", p.value, ")")
        )
      ),colour = "") + plot.title +
      theme(
        plot.title = element_text(hjust = 0.5),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank()
      ) +
      coord_fixed() + theme(aspect.ratio = 1, legend.key = element_blank()) + scale_y_reverse()
    return(g)

  } else if (method == "contour") {
    cat("It takes some time to compute and plot...please hang on a moment.")
    if (is.null(p.value)) {
      df_score <- df_score %>% mutate(point_shape = I(19))
    } else{
      df_score <-
        df_score %>% mutate(point_shape = if_else(.data[[p.value.label]] < p.value, I(19), I(1)))
    }

    colors.use <-
      generate_RedsBlues_colormap(df_score$statistic_G)

    g <- df_score %>% ggplot() +
      metR::geom_contour_fill(aes(x, y, z = statistic_G), kriging = T) +
      scale_fill_gradientn(colors = colors.use,guide = my.guide) +
      # scale_fill_distiller(palette = "RdYlBu",direction = -1)+
      geom_point(aes(x = x, y = y, shape = point_shape),
                 color = "grey10",
                 size = point.size) +
      guides(shape = guide_none()) +
      labs(title = paste0(
        "Spatial statistics with ",
        gstar.label,
        ifelse(
          is.null(p.value),
          "",
          paste0("(Contour, p.value=", p.value, ")")
        )
      ),colour = "") + plot.title +
      theme(
        plot.title = element_text(hjust = 0.5),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank()
      ) +
      coord_fixed() + theme(aspect.ratio = 1, legend.key = element_blank()) + scale_y_reverse()

    return(g)

  } else if (method == "computation") {
    return(df_score)
  } else{
    stop(
      "The parameter \"method\" in plotStatistics_Gi-function only supports for \"normal\",\"scoring\",\"3d\",\"density\",\"contour\"!\n Please check your input."
    )
  }
}


#### Heatmap ####
#' Feature expression heatmap
#'
#' Draws a heatmap of single cell feature expression.
#'
#' @param object Seurat object
#' @param features A vector of features to plot, defaults to \code{VariableFeatures(object = object)}
#' @param cells A vector of cells to plot
#' @param group.bar Add a color bar showing group status for cells
#' @param group.by A vector of variables to group cells by; pass 'ident' to group by cell identity classes
#' @param additional.group.by A vector of variables to group cells by.
#' @param additional.group.sort.by which variable from additional.group.by to use for sorting cells
#' @param colors.use A list of colors to use for the color bar. e.g., colors.use = list(ident = c("red","blue"))
#' @param colors.ggplot whether use ggplot color scheme; default: colors.ggplot = FALSE
#' @param color.heatmap A character string or vector indicating the colormap option to use. It can be the avaibale color palette in viridis_pal() or brewer.pal()
#' @param direction Sets the order of colors in the scale. If 1, the default colors are used. If -1, the order of colors is reversed.
#' @param n.colors number of basic colors to generate from color palette
#' @param slot Data slot to use, choose from 'raw.data', 'data', or 'scale.data'
#' @param disp.min Minimum display value (all values below are clipped)
#' @param disp.max Maximum display value (all values above are clipped); defaults to 2.5
#' if \code{slot} is 'scale.data', 6 otherwise
#' @param assay Assay to pull from
# @param check.plot Check that plotting will finish in a reasonable amount of time
#' @param text.size text size of gene names
#' @param label Label the cell identies above the color bar
#' @param size Size of text above color bar
#' @param hjust Horizontal justification of text above color bar
#' @param angle Angle of text above color bar
#' @param raster If true, plot with geom_raster, else use geom_tile. geom_raster may look blurry on
#' some viewing applications such as Preview due to how the raster is interpolated. Set this to FALSE
#' if you are encountering that issue (note that plots may take longer to produce/render).
#' @param draw.lines Include white lines to separate the groups
#' @param lines.width Integer number to adjust the width of the separating white lines.
#' Corresponds to the number of "cells" between each group.
#' @param group.bar.height Scale the height of the color bar
#' @param legend.remove whether remove legend
#' @param legend.title legend title
#' @param combine Combine plots into a single \code{\link[patchwork]{patchwork}ed}
#' ggplot object. If \code{FALSE}, return a list of ggplot objects
#'
#' @return A \code{\link[patchwork]{patchwork}ed} ggplot object if
#' \code{combine = TRUE}; otherwise, a list of ggplot objects
#'
#' @importFrom stats median
#' @importFrom scales hue_pal
#' @importFrom ggplot2 annotation_raster coord_cartesian scale_color_manual
#' ggplot_build aes_string
#' @importFrom patchwork wrap_plots
#' @importFrom grid textGrob gpar
#' @importFrom RColorBrewer brewer.pal
#' @import rlang
#' @import Seurat
#' @export
doHeatmap <- function (object,
                       features = NULL,
                       cells = NULL,
                       group.by = "ident",
                       colors.use = NULL,
                       colors.ggplot = FALSE,
                       group.bar = TRUE,
                       additional.group.by = NULL,
                       additional.group.sort.by = NULL,
                       color.heatmap = "RdBu",
                       n.colors = 11,
                       direction = -1,
                       disp.min = -2.1,
                       disp.max = 2.1,
                       slot = "scale.data",
                       assay = NULL,
                       label = TRUE,
                       text.size = 8,
                       size = 4,
                       hjust = 0,
                       angle = 45,
                       raster = TRUE,
                       draw.lines = TRUE,
                       lines.width = NULL,
                       group.bar.height = 0.02,
                       legend.remove = FALSE,
                       legend.title = NULL,
                       combine = TRUE)
{
  # if (is.null(color.heatmap)) {
  #   color.heatmap.use <- PurpleAndYellow()
  # } else

  group.colors <- colors.use
  if (!is.null(additional.group.by)) {
    if (is.null(additional.group.sort.by)) {
      additional.group.sort.by <- additional.group.by[1]
    }
  }

  cells <- cells %||% colnames(x = object)
  if (is.numeric(x = cells)) {
    cells <- colnames(x = object)[cells]
  }
  assay <- assay %||% DefaultAssay(object = object)
  DefaultAssay(object = object) <- assay
  features <- features %||% VariableFeatures(object = object)
  ## Why reverse???
  features <- rev(x = unique(x = features))
  disp.max <- disp.max %||% ifelse(test = slot == "scale.data",
                                   yes = 2.5,
                                   no = 6)
  possible.features <- rownames(x = GetAssayData(object = object,
                                                 slot = slot))
  if (any(!features %in% possible.features)) {
    bad.features <- features[!features %in% possible.features]
    features <- features[features %in% possible.features]
    if (length(x = features) == 0) {
      stop("No requested features found in the ",
           slot,
           " slot for the ",
           assay,
           " assay.")
    }
    warning(
      "The following features were omitted as they were not found in the ",
      slot,
      " slot for the ",
      assay,
      " assay: ",
      paste(bad.features,
            collapse = ", ")
    )
  }

  if (!is.null(additional.group.sort.by)) {
    if (any(!additional.group.sort.by %in% additional.group.by)) {
      bad.sorts <-
        additional.group.sort.by[!additional.group.sort.by %in% additional.group.by]
      additional.group.sort.by <-
        additional.group.sort.by[additional.group.sort.by %in% additional.group.by]
      if (length(x = bad.sorts) > 0) {
        warning(
          "The following additional sorts were omitted as they were not a subset of additional.group.by : ",
          paste(bad.sorts, collapse = ", ")
        )
      }
    }
  }

  data <-
    as.data.frame(x = as.matrix(x = t(x = GetAssayData(
      object = object,
      slot = slot
    )[features, cells, drop = FALSE])))

  object <- suppressMessages(expr = StashIdent(object = object,
                                               save.name = "ident"))
  group.by <- group.by %||% "ident"

  groups.use <-
    object[[c(group.by, additional.group.by[!additional.group.by %in% group.by])]][cells, , drop = FALSE]
  plots <- list()
  for (i in group.by) {
    data.group <- data
    if (!purrr::is_null(additional.group.by)) {
      additional.group.use <- additional.group.by[additional.group.by != i]
      if (!purrr::is_null(additional.group.sort.by)) {
        additional.sort.use = additional.group.sort.by[additional.group.sort.by != i]
      } else {
        additional.sort.use = NULL
      }
    } else {
      additional.group.use = NULL
      additional.sort.use = NULL
    }

    group.use <-
      groups.use[, c(i, additional.group.use), drop = FALSE]

    for (colname in colnames(group.use)) {
      if (!is.factor(x = group.use[[colname]])) {
        group.use[[colname]] <- factor(x = group.use[[colname]])
      }
    }

    if (draw.lines) {
      lines.width <- lines.width %||% ceiling(x = nrow(x = data.group) *
                                                0.0025)
      placeholder.cells <-
        sapply(
          X = 1:(length(x = levels(x = group.use[[i]])) *
                   lines.width),
          FUN = function(x) {
            return(RandomName(length = 20))
          }
        )
      placeholder.groups <-
        data.frame(rep(x = levels(x = group.use[[i]]), times = lines.width))
      group.levels <- list()
      group.levels[[i]] = levels(x = group.use[[i]])
      for (j in additional.group.use) {
        group.levels[[j]] <- levels(x = group.use[[j]])
        placeholder.groups[[j]] = NA
      }

      colnames(placeholder.groups) <- colnames(group.use)
      rownames(placeholder.groups) <- placeholder.cells

      group.use <- sapply(group.use, as.vector)
      rownames(x = group.use) <- cells

      group.use <- rbind(group.use, placeholder.groups)

      for (j in names(group.levels)) {
        group.use[[j]] <-
          factor(x = group.use[[j]], levels = group.levels[[j]])
      }

      na.data.group <-
        matrix(
          data = NA,
          nrow = length(x = placeholder.cells),
          ncol = ncol(x = data.group),
          dimnames = list(placeholder.cells,
                          colnames(x = data.group))
        )
      data.group <- rbind(data.group, na.data.group)
    }

    order_expr <-
      paste0('order(', paste(c(i, additional.sort.use), collapse = ','), ')')
    group.use = with(group.use, group.use[eval(parse(text = order_expr)), , drop =
                                            F])

    plot <- SingleRasterMap(
      data = data.group,
      raster = raster,
      disp.min = disp.min,
      disp.max = disp.max,
      feature.order = features,
      cell.order = rownames(x = group.use),
      group.by = group.use[[i]]
    )

    if (group.bar) {
      pbuild <- ggplot_build(plot = plot)
      group.use2 <- group.use
      cols <- list()
      na.group <- RandomName(length = 20)
      for (colname in rev(x = colnames(group.use2))) {
        if (colname == i) {
          colid = paste0('Identity (', colname, ')')
        } else {
          colid = colname
        }

        # Default
        if (colors.ggplot) {
          cols[[colname]] <-
            c(scales::hue_pal()(length(x = levels(x = group.use[[colname]]))))
        } else {
          cols[[colname]] <-
            scPalette(length(x = levels(x = group.use[[colname]])))
          if (colname %in% additional.group.by) {
            cols[[colname]] <-
              ggPalette(length(x = levels(x = group.use[[colname]])))
          }
        }


        #Overwrite if better value is provided
        if (!purrr::is_null(group.colors[[colname]])) {
          req_length = length(x = levels(group.use))
          if (length(group.colors[[colname]]) < req_length) {
            warning(
              "Cannot use provided colors for ",
              colname,
              " since there aren't enough colors."
            )
          } else {
            if (!purrr::is_null(names(group.colors[[colname]]))) {
              if (all(levels(group.use[[colname]]) %in% names(group.colors[[colname]]))) {
                cols[[colname]] <-
                  as.vector(group.colors[[colname]][levels(group.use[[colname]])])
              } else {
                warning(
                  "Cannot use provided colors for ",
                  colname,
                  " since all levels (",
                  paste(levels(group.use[[colname]]), collapse = ","),
                  ") are not represented."
                )
              }
            } else {
              cols[[colname]] <-
                as.vector(group.colors[[colname]])[c(1:length(x = levels(x = group.use[[colname]])))]
            }
          }
        }

        # Add white if there's lines
        if (draw.lines) {
          levels(x = group.use2[[colname]]) <-
            c(levels(x = group.use2[[colname]]), na.group)
          group.use2[placeholder.cells, colname] <- na.group
          cols[[colname]] <- c(cols[[colname]], "#FFFFFF")
        }
        names(x = cols[[colname]]) <-
          levels(x = group.use2[[colname]])

        y.range <- diff(x = pbuild$layout$panel_params[[1]]$y.range)
        y.pos <-
          max(pbuild$layout$panel_params[[1]]$y.range) + y.range * 0.015
        y.max <- y.pos + group.bar.height * y.range
        pbuild$layout$panel_params[[1]]$y.range <-
          c(pbuild$layout$panel_params[[1]]$y.range[1], y.max)

        plot <- suppressMessages(
          plot +
            annotation_raster(
              raster = t(x = cols[[colname]][group.use2[[colname]]]),
              xmin = -Inf,
              xmax = Inf,
              ymin = y.pos,
              ymax = y.max
            ) +
            annotation_custom(
              grob = grid::textGrob(
                label = colid,
                hjust = 0,
                gp = grid::gpar(cex = 0.75)
              ),
              ymin = mean(c(y.pos, y.max)),
              ymax = mean(c(y.pos, y.max)),
              xmin = Inf,
              xmax = Inf
            ) +
            coord_cartesian(ylim = c(0, y.max), clip = "off")
        )

        #  same time run  or not   ggplot version    ?
        if ((colname == i) && label) {
          x.max <- max(pbuild$layout$panel_params[[1]]$x.range)
          # x.divs <- pbuild$layout$panel_params[[1]]$x.major
          x.divs <-
            pbuild$layout$panel_params[[1]]$x.major %||% attr(x = pbuild$layout$panel_params[[1]]$x$get_breaks(),
                                                              which = "pos")
          # group.use$x <- x.divs
          # label.x.pos <- tapply(X = group.use$x, INDEX = group.use[[colname]],
          #                       FUN = median) * x.max
          # label.x.pos <- data.frame(group = names(x = label.x.pos),
          #                           label.x.pos)

          x <-
            data.frame(group = sort(x = group.use[[colname]]), x = x.divs)
          label.x.pos <-
            tapply(
              X = x$x,
              INDEX = x$group,
              FUN = function(y) {
                if (isTRUE(x = draw.lines)) {
                  mean(x = y[-length(x = y)])
                } else {
                  mean(x = y)
                }
              }
            )
          label.x.pos <-
            data.frame(group = names(x = label.x.pos), label.x.pos)

          plot <- plot + geom_text(
            stat = "identity",
            data = label.x.pos,
            aes_string(label = "group",
                       x = "label.x.pos"),
            y = y.max + y.max *
              0.03 * 0.5,
            angle = angle,
            hjust = hjust,
            size = size
          )
          plot <-
            suppressMessages(plot + coord_cartesian(
              ylim = c(0,
                       y.max + y.max * 0.002 * max(nchar(
                         x = levels(x = group.use[[colname]])
                       )) *
                         size),
              clip = "off"
            ))
        }
      }
    }
    plot <- plot + theme(line = element_blank()) +
      theme(axis.text.y = element_text(size = text.size))
    if (!is.null(color.heatmap)) {
      if(color.heatmap=="MyRdBu"){
        data_ <- MinMax(data = data, min = disp.min, max = disp.max)
        color.heatmap.use <- generate_RedsBlues_colormap(vec = as.matrix(data_),mid.color = "white")
      } else{
        if (length(color.heatmap) == 1) {
          color.heatmap.use <- tryCatch({
            RColorBrewer::brewer.pal(n = n.colors, name = color.heatmap)
          }, error = function(e) {
            scales::viridis_pal(option = color.heatmap, direction = -1)(n.colors)
          })
        } else if (length(color.heatmap) > 1) {
          color.heatmap.use <- color.heatmap
        }
        if (direction == -1) {
          color.heatmap.use <- rev(color.heatmap.use)
        }
        color.heatmap.use <- colorRampPalette(color.heatmap.use)(99)

      }
      plot <-
        plot + scale_fill_gradientn(colors = color.heatmap.use)
    }

    if (legend.remove) {
      plot <-  plot + NoLegend()
    }
    plots[[i]] <- plot
  }
  if (combine) {
    plots <- wrap_plots(plots)
  }
  return(plots)
}


#' generate_custom_colorscale
#'
#' @param color.heatmap a character string or vector indicating the colormap option to use. It can be the avaibale color palette in brewer.pal() or viridis_pal() (e.g., "Spectral","viridis")
#' @param n.colors number of basic colors to generate from color palette
#' @param reversescale whether to reverse the direction of the color scale
#' @param plotly.format whether to convert the colormap into a plotly-compatible format
#'
#' @return List. See the parameter `scene` in \code{plotly::layout}
#' @export
generate_custom_colorscale <- function(
    color.heatmap,
    n.colors = 8,
    reversescale = F,
    plotly.format = T
){
  if (length(color.heatmap) == 1) {
    colormap <- tryCatch({
      RColorBrewer::brewer.pal(n = n.colors, name = color.heatmap)
    }, error = function(e) {
      scales::viridis_pal(option = color.heatmap, direction = -1)(n.colors)
    })
    if (reversescale) {
      colormap <- rev(colormap)
    }
    colormap <- colorRampPalette(colormap)(99)
    colormap[[1]] <- "#E5E5E5"
  } else {
    colormap <- color.heatmap
  }
  if (n.colors < 2) {
    stop("`n.colors` must be greater than 1.")
  }
  if (plotly.format) {
    len.colormap <- length(colormap)
    color.index <- seq(0, 1, length.out = len.colormap)
    colormap <- BiocGenerics::sapply(
      X = 1:len.colormap,
      FUN = function(x) {
        return(list(color.index[[x]], colormap[[x]]))
      },
      simplify = F
    )

    return(colormap)
  } else {
    return(colormap)
  }
}


#' generate RedsBlues colormap for a numeric vector
#' @param vec A numeric vector. Positive values will be drawn red,
#' negative values will be drawn blue.
#' @param mid.color A string representing the color mapping to the `0` in the numeric vector.
#' @param switch.RedsBlues Whether to reverse the colormap
#'
#' @return Function. Generate RedsBlues colors.
#' @export
generate_RedsBlues_colormap <- function(
    vec,
    mid.color = "#E5E5E5",
    # mid.color = "#fffbff",
    switch.RedsBlues = F
){
  min_vec <- min(vec,na.rm = T)
  max_vec <- max(vec,na.rm = T)
  Colors.length <- 6
  Reds <-  c('#fee5d9','#fcbba1','#fc9272','#fb6a4a','#de2d26','#a50f15')
  Blues <- c('#eff3ff','#c6dbef','#9ecae1','#6baed6','#3182bd','#08519c')
  if (switch.RedsBlues) {
    Blues <- c('#fee5d9','#fcbba1','#fc9272','#fb6a4a','#de2d26','#a50f15')
    Reds <-  c('#eff3ff','#c6dbef','#9ecae1','#6baed6','#3182bd','#08519c')
  }
  if (0 <= min_vec) {
    # positive <- Reds
    col_fun <-
      colorRamp3(breaks = seq(0, max_vec ,length.out = (Colors.length+1) ),
                 colors = c(mid.color, Reds))
  } else if (0 >= max_vec) {
    # negative <- Blues
    col_fun <-
      colorRamp3(breaks = seq(min_vec, 0 ,length.out = (Colors.length+1) ),
                 colors = c(rev(Blues), mid.color))
  } else{
    # both <- both Blues and Reds
    Blues.breaks <- seq(min_vec, 0 ,length.out = (Colors.length+1) )
    Reds.breaks <- seq(0, max_vec ,length.out = (Colors.length+1) )

    col_fun <-
      colorRamp3(breaks = c(Blues.breaks,Reds.breaks[-1]), # onl one zero
                 colors = c(rev(Blues), mid.color, Reds))
  }
  return( col_fun(seq(from = min_vec,to = max_vec,length.out = 99)) )
}
