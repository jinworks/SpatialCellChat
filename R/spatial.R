#' Plot the cell-cell communication distance distribution
#'
#' @param object CellChat object
#' @param signaling.type the type of signaling
#' @param enriched.only whether to only show the communication distance distribution of the identified significant communication
#' @param density.alpha the transparence of the density plot
#'
#' @return
#' @export
#'
communicationDistPlot2 <-  function(
    object,
    signaling.type = c("All","Secreted","Contact"),
    enriched.only = TRUE,
    density.alpha=0.5
) {
  signaling.type <- match.arg(signaling.type)
  if (object@options$parameter[["all.contact.dependent"]] == TRUE) {
    if (signaling.type != "Contact") {
      message("All the signaling are the contact-dependent, and please set `signaling.type` as `Contact`! \n")
    }
  }
  interaction.range <- object@options$parameter$interaction.range
  res <- object@images$result.computeCellDistance

  # long-range distance
  d.spatial <- res$d.spatial
  # short-range distance adjacent matrix for contact-dependent and juxtacrine signaling
  adj.contact <- d.spatial*res$adj.contact

  if (enriched.only) {
    prob.cell <- object@net$prob.cell
    #LRsig.use.idx <- object@net$tmp$LRsig.use.idx
    d.spatial.all <- my_future_lapply(
      X = seq_len(length(prob.cell$i)),
      FUN = function(x){
        c(d.spatial[prob.cell$i[x], prob.cell$j[x]], adj.contact[prob.cell$i[x], prob.cell$j[x]])
      },
      simplify = F
    )
    d.spatial.all <- do.call(rbind, d.spatial.all)
    d.spatial.all[d.spatial.all == 0] <- NA
    df.Secreted <- data.frame(x = d.spatial.all[,1])
    #df.Secreted <- df.Secreted[df.Secreted$x > 0, , drop = FALSE]
    df.Contact <- data.frame(x = d.spatial.all[,2])
    #df.Contact <- df.Contact[df.Contact$x > 0, , drop = FALSE]
  } else {
    df.Secreted <- data.frame(
      x = d.spatial@x
    )
    df.Contact <- data.frame(
      x = adj.contact@x
    )
  }

  my.theme <- theme_bw() +theme(
    # axis.ticks = element_blank(),
    # axis.text = element_blank(),
    # plot.background = element_rect(linetype = "transparent")
  ) +
    theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) +
    theme(text = element_text(size = 10))

  p.Secreted <- ggplot(df.Secreted) +
    geom_density( aes(x = x, y = ..density..), fill="#69b3a2",alpha= density.alpha, na.rm = TRUE)+
    scale_y_continuous(expand = expansion(c(0, 0)), breaks = NULL, labels = NULL) +
    scale_x_continuous(limits = c(0, (interaction.range+10)),breaks = seq(0, (interaction.range+10), by = 20))+
    labs(title = "Secreting-dependent signaling", x="Distance between cell pairs (um)", y = "Density") +my.theme

  p.Contact <- ggplot(df.Contact) +
    geom_density( aes(x = x, y = ..density..),fill= "#404080",alpha=density.alpha, na.rm = TRUE)+
    scale_y_continuous(expand = expansion(c(0, 0)), breaks = NULL, labels = NULL) +
    scale_x_continuous(limits = c(0, (interaction.range+10)),breaks = seq(0, (interaction.range+10), by = 20))+
    labs(title = "Contact-dependent signaling", x="Distance between cell pairs (um)", y = "Density")+my.theme

  if(signaling.type == "All"){
    p <- patchwork::wrap_plots(p.Secreted,p.Contact,nrow = 2,ncol = 1)
    return(p)
  } else if(signaling.type == "Secreted"){
    return(p.Secreted)
  } else if(signaling.type == "Contact"){
    return(p.Contact)
  }
}

#' Plot the spatial distance distribution of inferred cell-cell communication
#'
#' @param object CellChat object
#' @param enriched.only whether to only show the communication distance distribution of the identified significant communication
#' @param density.alpha the transparence of the density plot
#'
#' @return
#' @export
#'
communicationDistPlot <-  function(
    object,
    enriched.only = TRUE,
    density.alpha=0.5
) {
  interaction.range <- object@options$parameter$interaction.range
  res <- object@images$result.computeCellDistance

  # long-range distance
  d.spatial <- res$d.spatial
  # short-range distance adjacent matrix for contact-dependent and juxtacrine signaling
  adj.contact <- d.spatial*res$adj.contact

  df.spatial <- vector("list", 2)
  if (enriched.only) {
    prob.cell <- object@net$prob.cell
    #LRsig.use.idx <- object@net$tmp$LRsig.use.idx
    d.spatial.all <- my_future_lapply(
      X = seq_len(length(prob.cell$i)),
      FUN = function(x){
        d.spatial[prob.cell$i[x], prob.cell$j[x]]
      },
      simplify = T
    )

    # d.spatial.all <- sapply(
    #   X = 1:length(prob.cell$i),
    #   FUN = function(x){
    #     d.spatial[prob.cell$i[x], prob.cell$j[x]]
    #   }
    # )
    d.spatial.all[d.spatial.all == 0] <- NA

    if (object@options$parameter[["all.diffusible"]] == TRUE) {
      cat(cli.symbol(),"All the signaling are diffusible! \n")
      df.spatial[[1]] <- data.frame(x = d.spatial.all)
    } else if (object@options$parameter[["all.contact.dependent"]] == TRUE) {
      cat(cli.symbol(),"All the signaling are contact-dependent! \n")
      df.spatial[[2]] <- data.frame(x = d.spatial.all)
    } else {
      cat(cli.symbol(),"The enriched signaling include both diffusible signaling and contact-dependent signaling!  \n")
      pairLRsig <- object@LR$LRsig
      nLR1 <- max(which(pairLRsig$annotation %in% c("Secreted Signaling", "ECM-Receptor", "Non-protein Signaling")))
      df.spatial[[1]] <- data.frame(x = d.spatial.all[prob.cell$k %in% seq_len(nLR1)])
      df.spatial[[2]] <- data.frame(x = d.spatial.all[prob.cell$k %in% seq(nLR1+1, nrow(pairLRsig))])
    }

  } else {
    df.spatial[[1]] <- data.frame(
      x = d.spatial@x
    )
    df.spatial[[2]] <- data.frame(
      x = adj.contact@x
    )
  }

  my.theme <- theme_bw() +theme(
    # axis.ticks = element_blank(),
    # axis.text = element_blank(),
    # plot.background = element_rect(linetype = "transparent")
  ) +
    theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) +
    theme(text = element_text(size = 10))

  gg <- list(NA, NA)
  if (!is.null(df.spatial[[1]])) {
    gg[[1]] <- ggplot(df.spatial[[1]]) +
      geom_density( aes(x = x, y = ..density..), fill="#69b3a2",alpha= density.alpha, na.rm = TRUE)+
      scale_y_continuous(expand = expansion(c(0, 0)), breaks = NULL, labels = NULL) +
      scale_x_continuous(limits = c(0, (interaction.range+10)),breaks = seq(0, (interaction.range+10), by = 20))+
      labs(title = "Diffusible signaling", x="Distance between cell pairs (um)", y = "Density") +my.theme
  }
  if (!is.null(df.spatial[[2]]) ) {
    gg[[2]] <- ggplot(df.spatial[[2]]) +
      geom_density( aes(x = x, y = ..density..),fill= "#404080",alpha=density.alpha, na.rm = TRUE)+
      scale_y_continuous(expand = expansion(c(0, 0)), breaks = NULL, labels = NULL) +
      scale_x_continuous(limits = c(0, (interaction.range+10)),breaks = seq(0, (interaction.range+10), by = 20))+
      labs(title = "Contact-dependent signaling", x="Distance between cell pairs (um)", y = "Density")+my.theme
  }
  gg <- gg[!is.na(gg)]
  p <- patchwork::wrap_plots(gg, nrow = length(gg),ncol = 1)

  return(p)
}


#' @title computeGridSize
#' @description
#' Use this function to view SpatialCellChat object after the "grid" operation and visualize it.
#'
#' @param object CellChat object
#' @param grid.resolution Numeric. By default, it is set to be 2.
#' @param do.plot Boolean. If FALSE, only print hints in terminal.
#' @param cellsize NULL or Numeric. If NULL, the function will tell you the default `cellsize` in the SpatialCellChat object.
#' If not NULL and be Numeric, the new cellsize will be equal to cellsize*grid.resolution.
#' Note: We use [sf::st_make_grid()] to make grid data, so the param is a numeric vector
#' of length 1 or 2 with target cellsize, please see details in [sf::st_make_grid()]:
#' for square or rectangular cells the width and height, for hexagonal cells the distance between
#' opposite edges (edge length is cellsize/sqrt(3)). A length units object can be passed,
#' or an area unit object with area size of the square or hexagonal cell.
#' @param what Character. One of: "polygons", "corners", or "centers". See details in [sf::st_make_grid()]
#' @param square Boolean. If FALSE, create hexagonal grid. See details in [sf::st_make_grid()]
#'
#' @return ggplot
#' @export
computeGridSize <- function(
    object,
    cellsize=NULL,
    grid.resolution=NULL,
    do.plot=T,
    what = "polygons",
    square = T
){

  # object <- cortexChat
  # cellsize=NULL;grid.resolution=NULL
  # what = "polygons";
  # square = T
  coordinates <- object@images$coordinates

  # if (ncol(coordinates) == 2) {
  #   colnames(coordinates) <- c("x_cent","y_cent")
  #   temp_coord = coordinates
  #   coordinates[,1] = temp_coord[,2]
  #   coordinates[,2] = temp_coord[,1]
  # } else {
  #   stop("Please check the input 'coordinates' and make sure it is a two column matrix.")
  # }

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
  cat(cli.symbol(),"If you do grid, the new `cellsize` will be ",newcellsize," units(pixels). It is about ",spot.size," um in International System of Units.\n")

  if(do.plot){
    df.sf <-
      sf::st_as_sf(coordinates,
                   coords = c("x_cent", "y_cent"),
                   remove = FALSE)
    sf::st_crs(df.sf) <- 3857

    square_grid <-
      sf::st_make_grid(df.sf,
                       cellsize = newcellsize,
                       what = what,
                       square = square)

    # convert `square_grid` into sf and add grid IDs
    square_grid_sf = sf::st_sf(square_grid) %>%
      # add grid ID
      dplyr::mutate(grid_id = seq_len(length(lengths(square_grid))))

    WithinMat <- sf::st_intersects(df.sf, square_grid_sf, sparse = F)

    within_nGrid <- rowSums(WithinMat)
    grid_nspots <- colSums(WithinMat)
    cat(cli.symbol(),"The \"Grid\" operation will generate about ",sum(grid_nspots!=0)," new spots.\n")

    if(all(within_nGrid>0)){
      cat(cli.symbol(symbol = "success"),"All spots have been assigned to grids respectively.\n")
    } else {
      cat(cli.symbol(symbol = "fail"),"Some spots have not been assigned to grids respectively.\n")
    }

    if(is.null(object@idents)){
      df.sf$cell_type <- factor("SpatialCellChatObj",levels = c("SpatialCellChatObj"))
    } else {
      df.sf$cell_type <- object@idents
    }

    color.use <- scPalette(nlevels(df.sf$cell_type))
    names(color.use) <- levels(df.sf$cell_type)

    p <- ggplot() +
      ggplot2::geom_sf(data = square_grid)+
      ggplot2::geom_sf(data = df.sf,mapping = aes(color=cell_type))+
      scale_color_manual(values = color.use, na.value = "grey40")+
      # scale_y_reverse()+
      theme_minimal()+theme(
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        # axis.title.y = element_blank(),
        # plot.background = element_blank(),
        # panel.grid.major = element_blank(),
        # panel.grid.minor = element_blank(),
      )+ theme(legend.key = element_blank())+
      labs(color="Cell Type")
    return(p)
  }
}


#' @title makeGridSpatialCellChat
#'
#' @param object CellChat object
#' @param data.slot A CellChat object may have 5 slots to store data matrix:
#' data.raw,data,data.signaling,data.scale,data.project, choose one data slot
#' to make a new grid CellChat object. By default, use "data".
# #' @param re.normalize Boolean. Wether to normalize the data in the data.slot you choose
# #' with the function `normalizeData`. See details in [SpatialCellChat::normalizeData()]
#' @param cellsize Numeric vector of length 1 or 2. We use [sf::st_make_grid()] to make grid data. It is a numeric vector
#' of length 1 or 2 with target cellsize, please see details in [sf::st_make_grid()]:
#' for square or rectangular cells the width and height, for hexagonal cells the distance between
#' opposite edges (edge length is cellsize/sqrt(3)). A length units object can be passed,
#' or an area unit object with area size of the square or hexagonal cell.
#'
#' Tips: Please run `computeGridSize` to get a proper `cellsize`.the unit of the `cellsize` you input is "pixel" or a unit used in old CellChat object's coordinates
#' @param what Character. One of: "polygons", "corners", or "centers". See details in [sf::st_make_grid()]
#' @param square Boolean. If FALSE, create hexagonal grid. See details in [sf::st_make_grid()]
#' @param idents.ties.method To determine the ident/cell group of each grid, we use [base::max.col()]
#' to choose a ident containing the most cells per grid. See `ties.method` in [base::max.col()].
#'
#' @return CellChat object
#' @export
#'
#' @examples
makeGridSpatialCellChat <- function(object,
                                data.slot=c("data","data.raw","data.signaling","data.scale","data.project"),
                                # re.normalize=F,
                                cellsize = c(5, 5),
                                what = "polygons",
                                square = T,
                                idents.ties.method = "first")
{
  data.slot <- match.arg(data.slot)
  coordinates <- object@images$coordinates
  colnames(coordinates) <- c("x_cent", "y_cent")
  coordinates$cell_type <- object@idents

  # some hints about contact.range and spot size in the new grid SpatialCellChat
  spatial.factors <- object@images$spatial.factors
  spot.size <- cellsize[[1]]*spatial.factors[["ratio"]]
  cat(cli.symbol(3),"The `cellsize` you input is ",cellsize[[1]]," units(pixels). It is about ",spot.size," um in International System of Units.\n")

  # `tol` can be the the half value of the minimum center-to-center distance(cellsize)
  spatial.factors[["tol"]] <- spot.size/2
  cat(cli.symbol(3),"A grid in the new object will have a spot diameter equal to ",spot.size," um. So the new `tol` in `spatial.factors` will be ",spatial.factors[["tol"]]," um.\n")


  df.sf <-
    sf::st_as_sf(coordinates,
                 coords = c("x_cent", "y_cent"),
                 remove = FALSE)
  sf::st_crs(df.sf) <- 3857

  square_grid <-
    sf::st_make_grid(df.sf,
                     cellsize = cellsize,
                     what = what,
                     square = square)

  # convert `square_grid` into sf and add grid IDs
  square_grid_sf = sf::st_sf(square_grid) %>%
    # add grid ID
    dplyr::mutate(grid_id = 1:length(lengths(square_grid)))

  # nCell x nGrid
  WithinMat <- sf::st_intersects(df.sf, square_grid_sf, sparse = F)

  within_nGrid <- Matrix::rowSums(WithinMat)

  if(all(within_nGrid>0)){
    cat(cli.symbol(symbol = "success"),"All spots have been assigned to grids respectively.\n")
  } else {
    cat(cli.symbol(symbol = "fail"),"Some spots have not been assigned to grids respectively.\n")
  }
  # Sum.Weight <- 1/within_nGrid

  SpotsCounts_perGrid <- Matrix::colSums(WithinMat)

  # create new idents
  cat(cli.symbol(),"create new idents...\n")
  LevelsIdents <- levels(object@idents)
  PreviousIdents <- object@idents
  within_Idents <- purrr::map(
    .x = 1:NCOL(WithinMat),
    .f = function(col) {
      index.within <- WithinMat[, col, drop = T] # return a Boolean vec
      new_col <- dplyr::if_else(index.within, PreviousIdents, NA)
      return(new_col)
    },
    .progress = T
  )

  # Convert WithinMat into list object
  within_Spots <- purrr::map(
    .x = 1:NCOL(WithinMat),
    .f = function(col) {
      index.within <- WithinMat[, col, drop = T] # return a vec
      new_col <- which(index.within == T)
      return(new_col)
    },
    .progress = T
  )

  # create new meta data
  cat(cli.symbol(),"create new meta data...\n")

  NewMeta <- pbapply::pbsapply(
    X = seq_along(within_Idents),
    FUN = function(i) {
      table_ <- within_Idents[[i]] %>% table() %>% as.vector()
      return(table_)
    }
  ) %>% t()
  colnames(NewMeta) <- LevelsIdents

  PresentIdents <-
    LevelsIdents[max.col(NewMeta, ties.method = idents.ties.method)]

  NewMeta <- as.data.frame(NewMeta)
  NewMeta[["cell.type"]] <- PresentIdents
  NewMeta[["spots.counts"]] <- SpotsCounts_perGrid
  NewMeta <- NewMeta[SpotsCounts_perGrid > 0, ]

  # NewCoordinates <- square_grid_sf[SpotsCounts_perGrid > 0, ]

  # new cell names
  NewSpotsNames <-
    stringr::str_c("Grid", square_grid_sf$grid_id[SpotsCounts_perGrid > 0])

  # NewCoordinates[["cell.type"]] <- NewMeta$cell.type

  # which spots within a specific grid
  within_Spots <- within_Spots[SpotsCounts_perGrid > 0]
  cat(cli.symbol(),"create new coordinates...\n")
  NewCoordinates_ <- pbapply::pbsapply(
    X = seq_along(within_Spots),
    FUN = function(grid_) {
      idx <- within_Spots[[grid_]] # idx in previous coordinate
      coor_ <- coordinates[idx, c("x_cent", "y_cent"),drop=F] %>% colMeans()
      return(coor_)
    }
  ) %>% t()
  # NewCoordinates <- cbind(NewCoordinates, NewCoordinates_)

  PreviousDataInput <- methods::slot(object, data.slot)
  cat(cli.symbol(),"create new expression data input...\n")
  NewDataInput <- pbapply::pbsapply(
    X = seq_along(within_Spots),
    FUN = function(grid_) {
      idx <- within_Spots[[grid_]] # idx in previous coordinate
      # idx.weight <- Sum.Weight[idx]

      # gene x cell
      # maybe need weight-sum!!!
      # now use `rowMeans`, calculate mean for each gene's expression.
      # expr_: nGene x 1
      expr_ <- PreviousDataInput[, idx, drop = F] %>% Matrix::rowMeans()
      # expr_ <- sapply(
      #   X = 1:NCOL(expr_),
      #   FUN = function(col){
      #     expr_[,col] <- expr_[,col]*idx.weight[[col]]
      #   },
      #   simplify = T
      # ) %>% Matrix::rowSums()
      return(expr_)
    }
  )

  NewDataInput <- as(NewDataInput,Class = "CsparseMatrix")

  # NewCoordinates_dist <- Rfast::Dist(NewCoordinates_)
  # diag(NewCoordinates_dist) <- NA
  # min(NewCoordinates_dist, na.rm = T)

  # set names of all objects
  rownames(NewDataInput) <- rownames(PreviousDataInput)
  colnames(NewDataInput) <- NewSpotsNames
  rownames(NewMeta) <- NewSpotsNames
  rownames(NewCoordinates_) <- NewSpotsNames
  colnames(NewCoordinates_) <- c("x_cent", "y_cent")
  NewMeta$cell.type <- factor(NewMeta$cell.type, levels = LevelsIdents)

  # create a SpatialCellChat Object
  data.input <- NewDataInput
  # if(re.normalize){
  #   # seems unnecessary
  #   data.input <- normalizeData(NewDataInput)
  # }else{
  #   data.input <- NewDataInput
  # }


  meta <- NewMeta
  coordinates <- NewCoordinates_ %>% as.data.frame()


  object_grid <-
    createSpatialCellChat(
      object = data.input,
      meta = meta,
      group.by = "cell.type",
      datatype = "spatial",
      coordinates = coordinates,
      spatial.factors = spatial.factors
    )
  object_grid@images[["within_nGrid"]] <- within_nGrid
  object_grid@images[["recommended.contact.range"]] <- spot.size

  cat(cli.symbol(1),"Making grid is done. A recommended contact.range is stored in object@images.\n")
  return(object_grid)
}


#' compute cell-cell distance and cell-cell contact adjacency matrix
#'
#' @description
#' Compute cell-cell distance based on the spatial coordinates and
#' generate cell-cell contact adjacency matrix with a contact.range restriction
#'
#' @param coordinates a data matrix in which each row gives the spatial locations/coordinates of each cell/spot
#' @param interaction.range Numeric(must positive). The maximum interaction/diffusion range of ligands. This hard threshold is used to filter out the connections between spatially distant cells
#' @param contact.range Numeric. The interaction range (Unit: microns) to restrict the contact-dependent signaling.
#' For spatial transcriptomics in a single-cell resolution, `contact.range` is approximately equal to the estimated cell diameter (i.e., the cell center-to-center distance), which means that contact-dependent and juxtacrine signaling can only happens when the two cells are contact to each other.
#' Typically, `contact.range = 10`, which is a typical human cell size. However, for low-resolution spatial data such as 10X visium, it should be the cell center-to-center distance (i.e., `contact.range = 100` for visium data).
#' Users can run the function `computeCellDistance` to get the center-to-center distance in the result's "d.spatial" key, which will help decide the value of `contact.range`.
#' @param tol Numeric. set a distance tolerance when computing cell-cell distances and contact adjacent matrix. Typically, `tol` should equal to the half value of cell/spot size in the unit of um.
#' For example, for 10X visium, `tol` can be set as `65/2`; for slide-seq, `tol` can be set as `10/2`.
#' If the cell/spot size is not known, we provide a function `computeCellDistance` to compute the center-to-center distance. `tol` can be the the half value of the minimum center-to-center distance.
#' By default `tol = 10/2`.
#'
#' @param ratio NULL or Numeric. ratio The conversion factor when converting spatial coordinates from Pixels or other units to Micrometers (i.e.,Microns).
#'
#' For example, setting `ratio = 0.18` indicates that 1 pixel equals 0.18um in the coordinates.
#' For 10X visium, it is the ratio of the theoretical spot size (i.e., 65um) over the number of pixels that span the diameter of a theoretical spot size in the full-resolution image (i.e., 'spot.size.fullres' in the 'scalefactors_json.json' file).
#'
#' @return List. A list has two keys: "d.spatial" and "adj.contact", storing cell-cell distance matrix and cell-cell contact adjacency matrix
#' @export
#'
#' @examples
computeCellDistance <- function (
    coordinates,
    interaction.range = 250,
    contact.range = 10,
    ratio = NULL,
    tol = 10/2
)
{
  if (ncol(coordinates) == 2) {
    colnames(coordinates) <- c("x_cent", "y_cent")
  }
  else {
    stop(cli.symbol(2),"Please check the input 'coordinates' and make sure it is a two column matrix.")
  }
  # calculate the distances
  d.spatial <- Rfast::Dist(coordinates)
  NC <- NCOL(d.spatial);gc()

  cat(cli.symbol(),paste0("Apply a predefined spatial distance threshold based on the interaction length(=",interaction.range,"um)...\n"))
  # ((interaction.range+tol)/ratio) will convert the interaction range's unit
  # into `pixel/specific unit` used in coordinates from `um`
  # d.spatial[d.spatial > ((interaction.range+tol)/ratio)] <- 0
  interaction.range.threshold <- ((interaction.range+tol)/ratio)

  SparseMatSlots <- purrr::map_dfr(
    .x = 1:NC,
    .f = function(j){
      i=which(d.spatial[,j]<=interaction.range.threshold & d.spatial[,j]>0)
      x=d.spatial[i,j,drop=T]
      j=rep.int(j,length(i))
      return(data.frame("i"=i,"j"=j,"x"=x))
    },
    .progress=T
  )

  # long-range distance
  d.spatial <- Matrix::sparseMatrix(
    i=SparseMatSlots$i,
    j=SparseMatSlots$j,
    x = SparseMatSlots$x,
    repr = "C",
    symmetric = F,
    index1 = T  # i and j are interpreted as 1-based indices
  )
  gc()

  # scale the distances
  if (!is.null(ratio)) {
    d.spatial@x <- d.spatial@x * ratio
  }

  # short-range distance based on contact.range
  adj.contact <- createCellCellContactMatrixFrom_dspatial(d.spatial = d.spatial,
                                                          tol = tol,
                                                          contact.threshold = contact.range)

  res <- list("d.spatial" = d.spatial, "adj.contact" = adj.contact)
  return(res)
}

#' create a cell-cell contact matrix from `d.spatial` object
#' @description
#' A function defined for `computeCellDistance`. We use the function `createSparseMatrixFrom_distObj` to convert the dist object into a CsparseMatrix object when running `computeCellDistance`.
#' Here we design the function to generate cell-cell contact adjacency matrix with a `contact.range` restriction and a `tol` tolerance passed from `computeCellDistance`'s parameters respectively
#'
#' @param d.spatial a CsparseMatrix (dsC Matrix) object, see in `computeCellDistance`:
#' d.spatial <- stats::dist(coordinates);
#' d.spatial <- createSparseMatrixFrom_distObj(d.spatial)
#' @param tol Numeric. set a distance tolerance when computing cell-cell distances and contact adjacent matrix,
#' passed from `computeCellDistance`'s same parameter.
#' By default `tol = NULL` means `tol` equals to the half value of cell/spot size in the unit of um.
#' @param contact.threshold Numeric. The interaction range threshold(Unit: microns) to restrict the contact-dependent signaling,
#' passed from `computeCellDistance`'s parameter `contact.range`.
#' If the distance of a pair of cells is larger than contact.threshold, regard the two cells as not in contact. Otherwise, regard the two cells as in contact.
#'
#' @return A `CsparseMatrix` object in `Matrix` Package
#' @export
#'
#' @examples
createCellCellContactMatrixFrom_dspatial <- function(d.spatial,contact.threshold,tol){

  adj.contact <- as(d.spatial,Class = "TsparseMatrix")
  # `min(adj.contact@x)` is the smallest value in the d.spatial matrix.
  # Take the larger value in `contact.threshold` and `min(adj.contact@x)` as the
  # interaction range threshold (Unit: microns) to restrict the contact-dependent signaling.
  # contact.range <- max(contact.threshold,min(adj.contact@x))
  contact.range <- contact.threshold

  cat(paste0(cli.symbol(),"Contact range is set as ",round(contact.range,3),"um, to restrict the contact-dependent signaling.\n"))
  cellcell.contact.index <- which(adj.contact@x <= (contact.range+tol))

  # update
  adj.contact@i <- adj.contact@i[cellcell.contact.index]
  adj.contact@j <- adj.contact@j[cellcell.contact.index]
  adj.contact@x <- rep.int(1,times = length(cellcell.contact.index))

  adj.contact <- as(adj.contact,Class = "CsparseMatrix")

  # fill up the diagonal, because d.spatial's diagonal is all-zero
  Matrix::diag(adj.contact) <- 1
  return(adj.contact)

}



#' Assessment of the colocalization between any cell groups
#'
#' @param coordinates a data matrix in which each row gives the spatial locations/coordinates of each cell/spot
#' @param group a factor vector defining the labels of each cell/spot
#' @param idents.use the identity to use
#' @param symmetric whether producing a symmetric colocalization matrix
#' @param nboot number of permutation when assessing the colocalization
#' @param seed.use set a random seed. By default, set the seed to 1.
#'
#' @return A square matrix giving the computed p-values after permuting cell labels
#'
#' @export

computeColocalization <- function(coordinates, group = NULL, idents.use = NULL, symmetric = TRUE, nboot = 100, seed.use = 1L) {
  if (ncol(coordinates) == 2) {
    colnames(coordinates) <- c("x_cent","y_cent")
  } else {
    stop("Please check the input 'coordinates' and make sure it is a two column matrix.")
  }
  if (!is.factor(group)) {
    stop("Please input the `group` as a factor!")
  }
  numCluster <- nlevels(group)
  level.use <- levels(group)
  level.use0 <- level.use
  if (is.null(idents.use)) {
    level.use <- level.use[level.use %in% unique(group)]
  } else {
    numCluster <- length(idents.use)
    level.use <- level.use[level.use %in% idents.use]
    level.use0 <- level.use
  }

  fdr = matrix(NaN,numCluster,numCluster)
  for (i in c(1:numCluster)){
    for (j in c(1:numCluster)){
      label_1 = level.use[i]
      label_2 = level.use[j]
      data_label_1 = coordinates[group==label_1,]
      data_label_2 = coordinates[group==label_2,]

      shuffle_sel = !(group %in% c(label_1,label_2))
      data_nont = coordinates[shuffle_sel,]

      coord2_dist = matrix(0,nrow=nrow(data_label_1),ncol=nrow(data_label_2))
      coord1 = cbind(data_label_1$x_cent,data_label_1$y_cent) # every row is a cell in label1
      coord2_shuf_dist_array = array(dim=c(nrow(data_label_1),nrow(data_label_2), nboot))

      for (idx2 in c(1:nrow(data_label_2))){
        coord2 = c(data_label_2$x_cent[idx2], data_label_2$y_cent[idx2])
        coord2_dist[,idx2] = raster::pointDistance(coord1,coord2,lonlat=F)
      }

      set.seed(seed.use)
      permutation <- replicate(nboot, sample.int(nrow(data_nont), size = nrow(data_label_1)))
      for (idx3 in c(1:nboot)){
        #data_label_1_shuf = data_nont[sample(c(1:nrow(data_nont)),nrow(data_label_1)),]
        data_label_1_shuf = data_nont[permutation[, idx3], ]
        coord1_shuf = cbind(data_label_1_shuf$x_cent,data_label_1_shuf$y_cent)
        for (idx2 in c(1:nrow(data_label_2))){
          coord2 = c(data_label_2$x_cent[idx2], data_label_2$y_cent[idx2])
          coord2_shuf_dist_array[,idx2,idx3] = raster::pointDistance(coord1_shuf,coord2,lonlat=F)
        }
      }


      label2_mean_dists = rowMeans(coord2_dist)
      label2_mean_dists_median = median(label2_mean_dists) # observed median of mean dists
      label2_shuf_medians = rep(0,nboot)

      fdr_count= c(0,0)
      for (idx4 in c(1: nboot)){
        label2_mean_dists_shuf = rowMeans(coord2_shuf_dist_array[,,idx4])
        label2_shuf_medians[idx4] = median(label2_mean_dists_shuf)
        if (label2_shuf_medians[idx4]>label2_mean_dists_median){
          fdr_count[1] = fdr_count[1]+1
        }
        else{
          fdr_count[2] = fdr_count[2]+1
        }
      }
      fdr[i,j] = fdr_count[2]/nboot
    }
  }
  if (symmetric) {
    for (i in 1:(numCluster-1)){
      for (j in (i+1):numCluster){
        if (fdr[i,j] > fdr[j,i]) {
          fdr[i,j] <- fdr[j,i]
        }
        fdr[j,i] <- fdr[i,j]
      }
    }

  }
  rownames(fdr) <- level.use0; colnames(fdr) <- level.use0
  return(fdr)
}




