
#' The SpatialCellChat Class
#'
#' The SpatialCellChat object is created from a single-cell transcriptomic data matrix, Seurat V3 or SingleCellExperiment object.
#' When inputting an data matrix, it takes a digital data matrices as input. Genes should be in rows and cells in columns. rownames and colnames should be included.
#' The class provides functions for data preprocessing, intercellular communication network inference, communication network analysis, and visualization.
#'
#'
#'# Class definitions
#' @importFrom methods setClassUnion
#' @importClassesFrom Matrix dgCMatrix
setClassUnion(name = 'AnyMatrix', members = c("matrix", "dgCMatrix"))
setClassUnion(name = 'AnyFactor', members = c("factor", "list"))

#' The key slots used in the SpatialCellChat object are described below.
#'
#' @slot data.raw raw count data matrix
#' @slot data normalized data matrix for SpatialCellChat analysis (Genes should be in rows and cells in columns)
#' @slot data.signaling a subset of normalized matrix only containing signaling genes
#' @slot data.scale scaled data matrix
#' @slot data.project projected data
#' @slot images a list of spatial image objects
#' @slot net a three-dimensional array P (K×K×N), where K is the number of cell groups and N is the number of ligand-receptor pairs. Each row of P indicates the communication probability originating from the sender cell group to other cell groups.
#' @slot netP a three-dimensional array representing cel-cell communication networks on a signaling pathway level
#' @slot DB ligand-receptor interaction database used in the analysis (a subset of CellChatDB)
#' @slot LR a list of information related with ligand-receptor pairs
#' @slot meta data frame storing the information associated with each cell
#' @slot idents a factor defining the cell identity used for all analysis. It becomes a list for a merged SpatialCellChat object
#' @slot var.features A list: one element is a vector consisting of the identified over-expressed signaling genes; one element is a data frame returned from the differential expression analysis
#' @slot dr List of the reduced 2D coordinates, one per method, e.g., umap/tsne/dm
#' @slot options List of miscellaneous data, such as parameters used throughout analysis, and a indicator whether the SpatialCellChat object is a single or merged
#'
#' @exportClass SpatialCellChat
#' @importFrom Rcpp evalCpp
#' @importFrom methods setClass
# #' @useDynLib SpatialCellChat
SpatialCellChat <- methods::setClass("SpatialCellChat",
                              slots = c(data.raw = 'AnyMatrix',
                                        data = 'AnyMatrix',
                                        data.signaling = "AnyMatrix",
                                        data.scale = "matrix",
                                        data.project = "AnyMatrix",
                                        images = "list",
                                        net = "list",
                                        netP = "list",
                                        meta = "data.frame",
                                        idents = "AnyFactor",
                                        DB = "list",
                                        LR = "list",
                                        var.features = "list",
                                        dr = "list",
                                        options = "list")
)
#' show method for SpatialCellChat
#'
#' @param SpatialCellChat object
#' @param show show the object
#' @param object object
#' @docType methods
#'
setMethod(f = "show", signature = "SpatialCellChat", definition = function(object) {
  if (object@options$mode == "single") {
    cat("An object of class", class(object), "created from a single dataset", "\n", nrow(object@data), "genes.\n",  ncol(object@data), "cells. \n")
  } else if (object@options$mode == "merged") {
    cat("An object of class", class(object), "created from a merged object with multiple datasets", "\n", nrow(object@data.signaling), "signaling genes.\n",  ncol(object@data.signaling), "cells. \n")
  }
  if (object@options$datatype == "RNA") {
    cat("SpatialCellChat analysis of single cell RNA-seq data! \n")
  } else {
    cat("SpatialCellChat analysis of", object@options$datatype, "data! The input spatial locations are \n")
    print(head(object@images$coordinates))
  }


  invisible(x = NULL)
})



#' Create a new SpatialCellChat object from a data matrix, Seurat or SingleCellExperiment object
#'
#' @param object a normalized (NOT count) data matrix (genes by cells), Seurat or SingleCellExperiment object
#' @param meta a data frame (rows are cells with rownames) consisting of cell information, which will be used for defining cell groups.
#' If input is a Seurat or SingleCellExperiment object, the meta data in the object will be used
#' @param group.by a char name of the variable in meta data, defining cell groups.
#' If input is a data matrix and group.by is NULL, the input `meta` should contain a column named 'labels',
#' If input is a Seurat or SingleCellExperiment object, USER must provide `group.by` to define the cell groups. e.g, group.by = "ident" for Seurat object
#' @param datatype By default datatype = "RNA"; when running SpatialCellChat on spatial imaging/transcriptomics data, set type = "spatial" and input `spatial.factors`
#'
#' @param coordinates a data matrix in which each row gives the spatial locations/coordinates of each cell/spot
#' @param spatial.factors a list containing two distance factors `ratio` and `tol`, which is dependent on spatial transcriptomics technologies (and specific datasets).
#'
#' USER must input this list when datatype = "spatial". spatial.factors must contain an element named `ratio`, which is the conversion factor when converting spatial coordinates from Pixels or other units to Micrometers (i.e.,Microns). For example, setting `ratio = 0.18` indicates that 1 pixel equals 0.18um in the coordinates,
#'
#' and another element named `tol`, which is the tolerance factor to increase the robustness when comparing the center-to-center distance against the `interaction.range`. This can be the half value of cell/spot size in the unit of um. If the cell/spot size is not known, we provide a function `computeCellDistance` to compute the cell center-to-center distance.
#' `tol` can be the the half value of the minimum center-to-center distance. Of note, CellChat does not need an accurate tolerance factor, which is used for determining whether considering the cell-pair as spatially proximal if their distance is greater than `interaction.range` but smaller than "`interaction.range` + `tol`".
#'
#' @param assay Assay to use when the input is a Seurat object. NB: The data in the `integrated` assay is not suitable for SpatialCellChat analysis because it contains negative values.
#' @param do.sparse whether use sparse format
#'
#' @return
#' @export
#' @importFrom methods as new
#' @examples
#' \dontrun{
#' Create a SpatialCellChat object from single-cell transcriptomics data
#' # Input is a data matrix
#' ## create a dataframe consisting of the cell labels
#' meta = data.frame(labels = cell.labels, row.names = names(cell.labels))
#' SpatialCellChat <- createSpatialCellChat(object = data.input, meta = meta, group.by = "labels")
#'
#' # input is a Seurat object
#' ## use the default cell identities of Seurat object
#' SpatialCellChat <- createSpatialCellChat(object = seurat.obj, group.by = "ident", assay = "RNA")
#' ## use other meta information as cell groups
#' SpatialCellChat <- createSpatialCellChat(object = seurat.obj, group.by = "seurat.clusters")
#'
#' # input is a SingleCellExperiment object
#' SpatialCellChat <- createSpatialCellChat(object = sce.obj, group.by = "sce.clusters")
#'
#' Create a SpatialCellChat object from spatial imaging data
#' # Input is a data matrix
#' SpatialCellChat <- createSpatialCellChat(object = data.input, meta = meta, group.by = "labels",
#'                            datatype = "spatial", coordinates = coordinates, spatial.factors = spatial.factors)
#'
#' # input is a Seurat object
#' SpatialCellChat <- createSpatialCellChat(object = seurat.obj, group.by = "ident", assay = "SCT",
#'                            datatype = "spatial", spatial.factors = spatial.factors)
#'
#' }
createSpatialCellChat <- function(object, meta = NULL, group.by = NULL,
                           datatype = c("RNA", "spatial"), coordinates = NULL, spatial.factors = NULL,
                           assay = NULL, do.sparse = T) {
  datatype <- match.arg(datatype)
  # data matrix as input
  if (inherits(x = object, what = c("matrix", "Matrix", "dgCMatrix"))) {
    cat(cli.symbol(3),"Create a SpatialCellChat object from a data matrix")
    data <- object
    if (is.null(group.by)) {
      group.by <- "labels"
    }
  }

  # Seurat object as input
  .error_if_no_Seurat <- function() {
    if (!requireNamespace("Seurat", quietly = TRUE)) {
      stop("Seurat installation required for working with Seurat objects")
    }
  }
  if (is(object,"Seurat")) {
    .error_if_no_Seurat()
    cat(cli.symbol(3),"Create a SpatialCellChat object from a Seurat object.\n")
    if (is.null(assay)) {
      assay = Seurat::DefaultAssay(object)
      if (assay == "integrated") {
        warning("The data in the `integrated` assay is not suitable for SpatialCellChat analysis! Please use the `RNA` or `SCT` assay! ")
      }
      cat(paste0("The `data` slot in the default assay is used. The default assay is ", assay),'\n')
    }

    data <- Seurat::GetAssayData(object, assay = assay, slot = "data") # normalized data matrix
    if (min(data) < 0) {
      stop("The data matrix contains negative values. Please ensure the normalized data matrix is used.")
    }
    if (is.null(meta)) {
      cat("The `meta.data` slot in the Seurat object is used as cell meta information",'\n')
      meta <- object@meta.data
      meta$ident <- Seurat::Idents(object)
    }
    if (is.null(group.by)) {
      group.by <- "ident"
    }
    if (is.null(coordinates)) {
      coordinates <- Seurat::GetTissueCoordinates(object, scale = NULL, cols = c("imagerow", "imagecol"))
      # spatial.factors <- object@images[["slice1"]]@spatial.factors
    }

  }
  # SingleCellExperiment object as input
  if (is(object,"SingleCellExperiment")) {
    cat(cli.symbol(3),"Create a SpatialCellChat object from a SingleCellExperiment object")
    if ("logcounts" %in% SummarizedExperiment::assayNames(object)) {
      cat("The `logcounts` assay is used",'\n')
      data <- SingleCellExperiment::logcounts(object)
    } else {
      stop("SingleCellExperiment object must contain an assay named `logcounts`")
    }
    if (is.null(meta)) {
      cat("The `colData` assay in the SingleCellExperiment object is used as cell meta information",'\n')
      meta <- as.data.frame(SingleCellExperiment::colData(object))
    }
    if (is.null(group.by)) {
      stop("`group.by` should be defined!")
    }
  }

  if (!inherits(x = data, what = c("dgCMatrix")) & do.sparse) {
    data <- as(data, "dgCMatrix")
  }

  if (!is.null(meta)) {
    if (inherits(x = meta, what = c("matrix", "Matrix"))) {
      meta <- as.data.frame(x = meta)
    }
    if (!is.data.frame(meta)) {
      stop("The input `meta` should be a data frame")
    }
    if (!identical(rownames(meta), colnames(data))) {
      cat("The cell barcodes in 'meta' is ", head(rownames(meta)),'\n')
      warning("The cell barcodes in 'meta' is different from those in the used data matrix.
              We now simply assign the colnames in the data matrix to the rownames of 'meta'!")
      rownames(meta) <- colnames(data)
    }
  } else {
    meta <- data.frame()
  }
  if (datatype %in% c("spatial")) {
    if (ncol(coordinates) == 2) {
      colnames(coordinates) <- c("x_cent","y_cent")
    } else {
      stop("Please check the input 'coordinates' and make sure it is a two column matrix.")
    }
    if (is.null(spatial.factors) | !("ratio" %in% names(spatial.factors)) | !("tol" %in% names(spatial.factors))) {
      stop("spatial.factors with elements named `ratio` and `tol` should be provided!")
    } else {
      images = list("coordinates" = coordinates,
                    "spatial.factors" = spatial.factors)
    }
    cat("Create a SpatialCellChat object from spatial imaging data...",'\n')
  } else {
    images <- list()
  }

  object <- methods::new(Class = "SpatialCellChat",
                         data = data,
                         images = images,
                         meta = meta)

  if (!is.null(meta) & nrow(meta) > 0) {
    cat("Set cell identities for the new SpatialCellChat object", '\n')
    if (!(group.by %in% colnames(meta))) {
      stop("The 'group.by' is not a column name in the `meta`, which will be used for cell grouping.")
    }
    object <- setIdent(object, ident.use = group.by) # set "labels" as default cell identity
    cat("The cell groups used for SpatialCellChat analysis are ",cli::col_red( levels(object@idents) ), '\n')
  }
  object@options$mode <- "single"
  object@options$datatype <- datatype
  return(object)
}




#' Merge SpatialCellChat objects
#'
#' @param object.list  A list of multiple SpatialCellChat objects.
#' @param add.names A vector containing the name of each dataset.
#' @param merge.data Whether to merge the `data` slot for ALL genes. Default only merges the data of signaling genes.
#' @param cell.prefix Whether to prefix cell names.
#' @param show.plot Whether to show merged plots.
#'
#' @importFrom methods slot new
#'
#' @return
#' @export
#'
#' @examples
mergeSpatialCellChat <- function(object.list, add.names = NULL, merge.data = FALSE, cell.prefix = FALSE,show.plot=T) {
  if(length(object.list)<2){
    stop(cli.symbol("fail")," Make sure the `object.list` has at least 2 objects.")
  }
  if (is.null(add.names)) {
    add.names <- paste("Dataset",1:length(object.list),sep = "_")
  }
  # initialization
  slot.name <- c("net", "netP", "idents" ,"LR", "var.features", "images")
  slot.combined <- vector("list", length(slot.name))
  names(slot.combined) <- slot.name
  for ( i in seq_len(length(slot.name)) ) {
    object.slot <- vector("list", length(object.list))
    for (j in seq_len(length(object.list))) {
      object.slot[[j]] <- methods::slot(object.list[[j]], slot.name[i])
    }
    slot.combined[[i]] <- object.slot
    names(slot.combined[[i]]) <- add.names
  }

  if (cell.prefix) {
    warning("Prefix cell names!")
    for (i in seq_len(length(object.list))) {
      colnames(object.list[[i]]@data) <- paste(add.names[i], sep = "_",colnames(object.list[[i]]@data))
    }
  } else {
    cell.names <- lapply(
      X = seq_len(length(object.list)),
      FUN = function(i){
        rownames(object.list[[i]]@meta)
      }
    )
    cell.names <- do.call(c,args = cell.names)

    if (sum(duplicated(cell.names))>0) {
      stop(cli.symbol("fail")," Duplicated cell names were detected across datasets!! Please set cell.prefix = TRUE")
    }
  }

  meta.use <- colnames(object.list[[1]]@meta)
  for (i in 2:length(object.list)) {
    meta.use <- meta.use[meta.use %in% colnames(object.list[[i]]@meta)]
  }

  dataset.name <- c()
  cell.names <- c()
  meta.joint <- data.frame()
  for (i in 1:length(object.list)) {
    dataset.name <- c(dataset.name, rep(add.names[i], length(colnames(object.list[[i]]@data))))
    cell.names <- c(cell.names, colnames(object.list[[i]]@data))
    meta.joint <- rbind(meta.joint, object.list[[i]]@meta[ , meta.use, drop = FALSE])
  }
  if (!identical(rownames(meta.joint), cell.names)) {
    cat("The cell barcodes in merged 'meta' is ", cli::col_blue(head(rownames(meta.joint))),'\n')
    cat(cli::col_red("The cell barcodes in merged 'meta' is different from those in the used data matrix.\nWe now simply assign the colnames in the data matrix to the rownames of merged 'meta'!\n"))
    rownames(meta.joint) <- cell.names
    cat("Now the cell barcodes in merged 'meta' is ", cli::col_blue(head(rownames(meta.joint))),'\n')
  }

  # dataset.name <- data.frame(dataset.name = dataset.name, row.names = cell.names)
  meta.joint$datasets <- factor(dataset.name, levels = add.names)

  genes.use <- rownames(object.list[[1]]@data)
  for (i in 2:length(object.list)) {
    genes.use <- genes.use[genes.use %in% rownames(object.list[[i]]@data)]
  }
  data.joint <- c()
  for (i in 1:length(object.list)) {
    data.joint <- cbind(data.joint, object.list[[i]]@data[genes.use, ])
  }

  gene.signaling.joint = unique( do.call(what = c,args = lapply(object.list, function(x) rownames(x@data.signaling))) )
  data.signaling.joint <- data.joint[rownames(data.joint) %in% gene.signaling.joint, ]

  idents.joint <- c()
  idents.levels <- c()
  for (i in 1:length(object.list)) {
    idents.joint <- c(idents.joint, as.character(object.list[[i]]@idents))
    idents.levels <- union(idents.levels, levels(object.list[[i]]@idents))
  }
  names(idents.joint) <- cell.names
  idents.joint <- factor(idents.joint, levels = idents.levels)
  slot.combined$idents$joint <- idents.joint

  if (merge.data) {
    cat("Merge the following slots:",cli::col_blue(" 'data','data.signaling','images','net', 'netP','meta', 'idents', 'var.features' , 'DB', and 'LR'.\n"))
    merged.object <- methods::new(
      Class = "SpatialCellChat",
      data = data.joint,
      data.signaling = data.signaling.joint,
      images = slot.combined$images,
      net = slot.combined$net,
      netP = slot.combined$netP,
      meta = meta.joint,
      idents = slot.combined$idents,
      var.features = slot.combined$var.features,
      LR = slot.combined$LR,
      DB = object.list[[1]]@DB)
  } else {
    cat("Merge the following slots:",cli::col_blue(" 'data.signaling','images','net', 'netP','meta', 'idents', 'var.features' , 'DB', and 'LR'.\n"))
    merged.object <- methods::new(
      Class = "SpatialCellChat",
      data.signaling = data.signaling.joint,
      images = slot.combined$images,
      net = slot.combined$net,
      netP = slot.combined$netP,
      meta = meta.joint,
      idents = slot.combined$idents,
      var.features = slot.combined$var.features,
      LR = slot.combined$LR,
      DB = object.list[[1]]@DB)
  }
  merged.object@options$mode <- "merged"
  datatype.joint <- c()
  for (j in 1:length(object.list)) {
    datatype.joint <- union(datatype.joint, slot(object.list[[j]], "options")$datatype)
  }
  if (length(datatype.joint) == 1){
    merged.object@options$datatype <- datatype.joint
  } else {
    message("The data types in these objects are  ", datatype.joint,'\n')
    stop("Comparison analysis is not suggested for different types of data.\n")
  }

  gg1 <- spatialDimPlot(object.list[[1]],title.name = add.names[[1]],legend.text.size = 10)
  gg2 <- spatialDimPlot(object.list[[2]],title.name = add.names[[2]],legend.text.size = 10)

  df.labels <- data.frame(
    "group.cellchat"=merged.object@idents$joint,
    "dataset"=merged.object@meta$datasets
  )
  df.stat <- df.labels %>%
    group_by(group.cellchat,dataset) %>%
    summarise(cell.num=n()) %>%
    group_by(dataset) %>%
    mutate(cell.num2=sum(cell.num)) %>%
    mutate(cell.fraction=cell.num/cell.num2)

  dataset.levels <- levels(df.labels$dataset)
  df.stat <- df.stat %>%
    mutate(flag=if_else(dataset==dataset.levels[[1]],-1,1)) %>%
    mutate(cell.num3=cell.num*flag)

  max.cell.num <- max(abs(df.stat$cell.num))
  brks <- round(quantile(0:max.cell.num,probs = c(0,0.5,1))[-1])
  brks.lab <- as.character(c(rev(brks),0,brks))
  brks <- c(-rev(brks),0,brks)

  gg3 <- ggplot(df.stat, aes(x = group.cellchat, y = cell.num3, fill = dataset)) +   # Fill column
    geom_bar(stat = "identity", width = .6) +   # draw the bars
    geom_text(aes(label = cell.num,hjust = -flag),
              # hjust=0,
              # vjust = -0.5,
              color = "black")+
    scale_y_continuous(expand=expansion(mult = c(0.1, 0.1)),
                       limits = c(-max.cell.num, max.cell.num),
                       breaks = brks,labels = brks.lab)+
    scale_fill_manual(values = c("#d6604d","#4393c3"),guide=F) +
    coord_flip() +  # Flip axes
    labs(title=NULL,x=NULL,y="Cell Group Size") +
    theme_bw() +
    theme(
      # plot.title = element_text(hjust = .5, size = 14),
      panel.grid.major.y = element_blank(),
      # axis.ticks = element_blank(),
      legend.title = element_blank(),
      axis.title.x = element_text(size = 10),
      axis.text.y = element_text(
        colour = scPalette(nlevels(df.labels$group.cellchat)),
        size = 12,
        face = "bold"
      )
    )
  # df.stat.order <- df.stat %>%
  #   group_by(dataset) %>%
  #   mutate(group.cellchat=forcats::fct_reorder(group.cellchat,cell.num))
  # df.stat.order$group.cellchat <- as.character(df.stat.order$group.cellchat)
  gg4 <- ggplot(df.stat, aes(x = cell.fraction, y = group.cellchat)) +
    geom_segment(aes(yend = group.cellchat), xend = 0, colour = "grey50") +
    geom_point(size = 3, aes(colour = dataset)) +
    scale_colour_manual(values = c("#d6604d","#4393c3")) +
    theme_bw() +
    theme(panel.grid.major.y = element_blank()) +
    labs(title = NULL,x="Cell Group Fraction",y=NULL)+
    # facet_grid(dataset ~ ., scales = "free_y", space = "free_y")+
    theme(
      legend.title = element_blank(),
      axis.title.x = element_text(size = 10),
      axis.text.y = element_text(colour=scPalette(nlevels(df.labels$group.cellchat)),
                                 size = 12,face="bold"))

  pp <- cowplot::plot_grid(gg1,gg2,gg3,gg4,ncol=2,nrow=2,
                           rel_heights = c(1.5,1),
                           rel_widths = c(1.5,1.5))
  misc(merged.object,key="merge.stat") <- pp
  if(show.plot){
    show(pp)
  }

  return(merged.object)
}



#' Update a single SpatialCellChat object
#'
#' Update a single previously calculated SpatialCellChat object (version < 1.6.0)
#'
#' version < 0.5.0: `object@var.features` is now `object@var.features$features`; `object@net$sum` is now `object@net$weight` if `aggregateNet` has been run.
#'
#' version 1.6.0: a `object@images` slot is added and `datatype` is added in `object@options$datatype`
#'
#' version 2.0.0: `images$scale.factors` is changed to `images$spatial.factors` for spatial transcriptomics data analysis.
#'
#' @param object SpatialCellChat object
#'
#' @return a updated SpatialCellChat object
#' @export
#'
updateSpatialCellChat <- function(object) {

  DB <- object@DB
  # interaction_input <- DB$interaction
  # if (("category" %in% colnames(interaction_input) == FALSE) & ("annotation" %in% colnames(interaction_input) == TRUE)) {
  #   message("Change the column name `annotation` in object@DB$interaction to `category` since CellChat v2")
  #   colnames(interaction_input) <- plyr::mapvalues(colnames(interaction_input),from = c("annotation"), to = c("category"), warn_missing = TRUE)
  #   DB$interaction <- interaction_input
  # }
  if (is.character(object@var.features)) {
    message("Update slot 'var.features' from a vector to a list")
    var.features.new <- list(features = object@var.features)
  } else {
    var.features.new <- object@var.features
  }
  if ("sum" %in% names(object@net)) {
    net <- object@net
    net$weight <- net$sum
  } else {
    net <- object@net
  }
  if (!("mode" %in% names(object@options))) {
    object@options$mode <- "single"
  }
  if (!("datatype" %in% names(object@options))) {
    object@options$datatype <- "RNA"
    images = list()
  } else {
    images = object@images
  }
  meta = object@meta
  if ("slices" %in% colnames(meta)) {
    meta$samples <-  meta$slices
    meta$slices = NULL
  }
  if (!("samples" %in% colnames(meta))) {
    warning("The 'meta' data does not have a column named `samples`. We now add this column and all cells are assumed to belong to `sample1`!")
    meta$samples <- "sample1"
    meta$samples <- factor(meta$samples)
  } else if (is.factor(meta$samples) == FALSE) {
    warning("The 'meta$samples' is not a factor. We now force it as a factor!")
    meta$samples <- factor(meta$samples)
  }
  if (object@options$datatype %in% c("spatial")) {
    if ("scale.factors" %in% names(object@images)) {
      images$spatial.factors <- as.data.frame(images$scale.factors)
      images$scale.factors <- NULL
      cat(cli.symbol(3),"Please check whether `ratio` and `tol` are in the `spatial.factors`
          when using spatial transcriptomics data.")
    }
  }
  object.new <- methods::new(
    Class = "SpatialCellChat",
    data.raw = object@data.raw,
    data = object@data,
    data.signaling = object@data.signaling,
    data.scale = object@data.scale,
    data.project = object@data.project,
    images = images,
    net = net,
    netP = object@netP,
    meta = meta,
    idents = object@idents,
    DB = DB,
    LR = object@LR,
    var.features = var.features.new,
    dr = object@dr,
    options = object@options
  )
  return(object.new)
}

#' Update a SpatialCellChat object by lifting up the cell groups to the same cell labels across all datasets
#'
#' This function is useful when comparing inferred communications across different datasets with different cellular compositions
#'
#' @param object A single or merged SpatialCellChat object
#' @param group.new A char vector giving the cell labels to lift up. The order of cell labels in the vector will be used for setting the new cell identity.
#'
#'  If the input is a merged SpatialCellChat object and group.new = NULL, it will use the cell labels from one dataset with the maximum number of cell groups
#'
#'  If the input is a single SpatialCellChat object, `group.new` must be defined.
#'
#' @return a updated SpatialCellChat object
#'
#' @export
#'
liftSpatialCellChat <- function(object, group.new = NULL) {
  if (object@options$mode == "merged") {
    idents <- object@idents[1:(length(object@idents)-1)]
    if (is.null(group.new)) {
      group.max.all <- unique(unlist(sapply(idents, levels)))
      group.num <- sapply(idents, nlevels)
      group.num.max <- max(group.num)
      group.max <- levels(idents[[which(group.num == group.num.max)]])
      if (length(group.max) != length(group.max.all)) {
        stop("SpatialCellChat object cannot lift up due to the missing cell groups in any dataset. Please define the parameter `group.new`!")
      }
    } else {
      group.max <- group.new
      group.num.max <- length(group.new)
    }
    message(paste0("The SpatialCellChat object will be lifted up using the cell labels ", paste(group.max, collapse=", ")))
    for (i in 1:length(idents)) {
      cat("Update slots object@net, object@netP, object@idents in dataset ", names(object@idents)[i],'\n')
      # cat("Update slot object@net...", '\n')
      net <- object@net[[i]]
      group.i <- levels(idents[[i]])
      # group.existing <- group.max[group.max %in% group.i]
      group.existing <- group.i[group.i %in% group.max]
      group.existing.index <- which(group.max %in% group.existing)
      for (net.j in names(net)) {
        values <- net[[net.j]]
        if (net.j %in% c("prob","pval")) {
          values.new <- array(data = 0, dim = c(group.num.max, group.num.max, dim(values)[3]),
                              dimnames = list(group.max, group.max, dimnames(values)[[3]]))
          values.new[group.existing.index, group.existing.index, ] <- values
          net[[net.j]] <- values.new
        }
        if (net.j %in% c("count","sum","weight")) {
          values.new <- array(data = 0, dim = c(group.num.max, group.num.max),
                              dimnames = list(group.max, group.max))
          values.new[group.existing.index, group.existing.index] <- values
          net[[net.j]] <- values.new
        }
        if (net.j %in% c("pairwiseRank")) {
          for (k in 1:length(values)) {
            values.new1 <- vector("list", group.num.max)
            values.new1[group.existing.index] <- values[[k]]
            temp <- values[[k]][[1]]
            temp$prob  <- 0; temp$pval <- 1
            for (kk in setdiff(1:group.num.max, group.existing.index)) {
              values.new1[[kk]] <- temp
            }
            names(values.new1) <- group.max
            values[[k]] <- values.new1
          }
          values.new <- vector("list", group.num.max)
          values.new[group.existing.index] <- values
          temp <- lapply(values.new1, function(x) {
            x$prob  <- 0; x$pval <- 1
            return(x)
          })
          for (kk in setdiff(1:group.num.max, group.existing.index)) {
            values.new[[kk]] <- temp
          }
          names(values.new) <- group.max
        }
        net[[net.j]] <- values.new
      }
      object@net[[i]] <- net

      # cat("Update slot object@netP...", '\n')
      netP <- object@netP[[i]]
      for (netP.j in names(netP)) {
        values <- netP[[netP.j]]
        if (netP.j %in% c("pathways")) {
          values.new <- values
          netP[[netP.j]] <- values.new
        }
        if (netP.j %in% c("prob")) {
          values.new <- array(data = 0, dim = c(group.num.max, group.num.max, dim(values)[3]),
                              dimnames = list(group.max, group.max, dimnames(values)[[3]]))
          values.new[group.existing.index, group.existing.index, ] <- values
          netP[[netP.j]] <- values.new
        }

        if (netP.j %in% c("centr")) {
          values.new <- array(data = 0, dim = c(dim(values)[1], group.num.max, dim(values)[3]),
                              dimnames = list(dimnames(values)[[1]], group.max, dimnames(values)[[3]]))
          values.new[ , group.existing.index, ] <- values
          netP[[netP.j]] <- values.new
        }
      }
      object@netP[[i]] <- netP
      # cat("Update slot object@idents...", '\n')
      # idents[[i]] <- factor(group.max, levels = group.max)
      idents[[i]] <- factor(idents[[i]], levels = group.max)
    }
    object@idents[1:(length(object@idents)-1)] <- idents
  } else {
    if (is.null(group.new)) {
      stop("Please define the parameter `group.new`!")
    } else {
      group.max <- as.character(group.new)
      group.num.max <- length(group.new)
      message(paste0("The SpatialCellChat object will be lifted up using the cell labels ", paste(group.max, collapse=", ")))
    }
    cat("Update slots object@net, object@netP, object@idents in a single dataset...", '\n')
    # cat("Update slot object@net...", '\n')
    net <- object@net
    idents <- object@idents
    group.i <- levels(idents)
    # group.existing <- group.max[group.max %in% group.i]
    group.existing <- group.i[group.i %in% group.max]
    group.existing.index <- which(group.max %in% group.existing)
    for (net.j in names(net)) {
      values <- net[[net.j]]
      if (net.j %in% c("prob","pval")) {
        values.new <- array(data = 0, dim = c(group.num.max, group.num.max, dim(values)[3]),
                            dimnames = list(group.max, group.max, dimnames(values)[[3]]))
        values.new[group.existing.index, group.existing.index, ] <- values
        net[[net.j]] <- values.new
      }
      if (net.j %in% c("count","sum","weight")) {
        values.new <- array(data = 0, dim = c(group.num.max, group.num.max),
                            dimnames = list(group.max, group.max))
        values.new[group.existing.index, group.existing.index] <- values
        net[[net.j]] <- values.new
      }
      if (net.j %in% c("pairwiseRank")) {
        for (k in 1:length(values)) {
          values.new1 <- vector("list", group.num.max)
          values.new1[group.existing.index] <- values[[k]]
          temp <- values[[k]][[1]]
          temp$prob  <- 0; temp$pval <- 1
          for (kk in setdiff(1:group.num.max, group.existing.index)) {
            values.new1[[kk]] <- temp
          }
          names(values.new1) <- group.max
          values[[k]] <- values.new1
        }
        values.new <- vector("list", group.num.max)
        values.new[group.existing.index] <- values
        temp <- lapply(values.new1, function(x) {
          x$prob  <- 0; x$pval <- 1
          return(x)
        })
        for (kk in setdiff(1:group.num.max, group.existing.index)) {
          values.new[[kk]] <- temp
        }
        names(values.new) <- group.max
      }
      net[[net.j]] <- values.new
    }
    object@net <- net


    # cat("Update slot object@netP...", '\n')
    netP <- object@netP
    for (netP.j in names(netP)) {
      values <- netP[[netP.j]]
      if (netP.j %in% c("pathways")) {
        values.new <- values
        netP[[netP.j]] <- values.new
      }
      if (netP.j %in% c("prob")) {
        values.new <- array(data = 0, dim = c(group.num.max, group.num.max, dim(values)[3]),
                            dimnames = list(group.max, group.max, dimnames(values)[[3]]))
        values.new[group.existing.index, group.existing.index, ] <- values
        netP[[netP.j]] <- values.new
      }
      if (netP.j %in% c("centr")) {
        values.new <- array(data = 0, dim = c(dim(values)[1], group.num.max, dim(values)[3]),
                            dimnames = list(dimnames(values)[[1]], group.max, dimnames(values)[[3]]))
        values.new[ , group.existing.index, ] <- values
        netP[[netP.j]] <- values.new
      }
    }
    object@netP <- netP

    # cat("Update slot object@idents...", '\n')
    idents <- factor(idents, levels = group.max)
    object@idents <- idents
  }

  return(object)
}


#' Subset SpatialCellChat object using a portion of cells
#'
#' @param object  A SpatialCellChat object (either an object from a single dataset or a merged objects from multiple datasets)
#' @param cells.use a char vector giving the cell barcodes to subset. If cells.use = NULL, USER must define `idents.use`
#' @param idents.use a subset of cell groups used for analysis
#' @param group.by cell group information; default is `object@idents`; otherwise it should be one of the column names of the meta slot
#' @param invert whether invert the idents.use
#' @param thresh threshold of the p-value for determining significant interaction. A parameter as an input of the function \link{computeCommunProbPathway}
#'
#' @importFrom methods slot new
#'
#' @return
#' @export
#'

subsetSpatialCellChat <- function(object, cells.use = NULL, idents.use = NULL, group.by = NULL, invert = FALSE, thresh = 0.05) {
  # code for test
  # object <- Chat
  # cells.use = NULL
  # idents.use = c("Basal","Spinous","Supraspinous")
  # group.by=NULL
  # invert=F
  # thresh = 0.05

  if (!is.null(idents.use)) {
    if (is.null(group.by)) {
      labels <- object@idents
      if (object@options$mode == "merged") {
        message("Use the joint cell labels from the merged SpatialCellChat object")
        labels <- object@idents$joint
      }
    } else {
      labels <- object@meta[[group.by]]
    }
    if (!is.factor(labels)) {
      labels <- factor(labels)
    }

    # View(rownames(object@meta))
    names(labels) <- dimnames(object@meta)[[1]] # key
    level.use0 <- levels(labels)
    level.use <- levels(labels)[levels(labels) %in% unique(labels)]

    if (invert) {
      level.use <- level.use[!(level.use %in% idents.use)]
    } else {
      level.use <- level.use[level.use %in% idents.use]
    }
    cells.use.index <- which(as.character(labels) %in% level.use)
    cells.use <- names(labels)[cells.use.index]
  } else if (!is.null(cells.use)) {
    labels <- object@idents
    if (object@options$mode == "merged") {
      message("Use the joint cell labels from the merged SpatialCellChat object")
      labels <- object@idents$joint
    }
    names(labels) <- rownames(object@meta) # key
    level.use0 <- levels(labels)
    level.use <- levels(labels)[levels(labels) %in% unique(as.character(labels[cells.use]))]
    cells.use.index <- which(names(labels) %in% cells.use)
  } else {
    stop("USER should define either `cells.use` or `idents.use`!")
  }
  cat("The subset of cell groups used for SpatialCellChat analysis are ", cli::col_red(level.use), '\n')

  if (nrow(object@data) > 0) {
    data.subset <- object@data[, cells.use.index]
  } else {
    data.subset <- matrix(0, nrow = 0, ncol = 0)
  }
  if (nrow(object@data.project) > 0) {
    data.project.subset <- object@data.project[, cells.use.index]
  } else {
    data.project.subset <- matrix(0, nrow = 0, ncol = 0)
  }
  data.signaling.subset <- object@data.signaling[, cells.use.index]

  meta.subset <- object@meta[cells.use.index, , drop = FALSE]


  if (object@options$mode == "merged") {
    idents <- object@idents[1:(length(object@idents)-1)]
    group.existing <- level.use0[level.use0 %in% level.use]
    group.existing.index <- which(level.use0 %in% level.use)
    net.subset <- vector("list", length = length(object@net))
    netP.subset <- vector("list", length = length(object@netP))
    idents.subset <- vector("list", length = length(idents))
    names(net.subset) <- names(object@net)
    names(netP.subset) <- names(object@netP)
    names(idents.subset) <- names(object@idents[1:(length(object@idents)-1)])
    images.subset <- vector("list", length = length(idents))
    names(images.subset) <- names(object@idents[1:(length(object@idents)-1)])

    for (i in seq_len(length(idents))) {
      cat(cli.symbol(),"Update slots object@images, object@net, object@netP, object@idents in dataset ", cli::col_red(names(object@idents)[i]),'\n')

      images <- object@images[[i]]
      for (images.j in names(images)) {
        values <- images[[images.j]]
        if (images.j %in% c("coordinates")) {
          # values.new <- values[cells.use.index, ,drop=F]
          images[[images.j]] <- NULL
        }
        if (images.j %in% c("distance")) {
          values.new <- values[group.existing.index, group.existing.index, drop = FALSE]
          images[[images.j]] <- values.new
        }
      }
      images.subset[[i]] <- images

      # cat("Update slot object@net...", '\n')

      net <- methods::slot(object,"net")[[i]]
      for (net.j in names(net)) {
        values <- net[[net.j]]
        if (net.j %in% c("prob","pval")) {
          values.new <- values[group.existing.index, group.existing.index, ,drop = FALSE]
          net[[net.j]] <- values.new
        }
        if (net.j %in% c("count","sum","weight")) {
          values.new <- values[group.existing.index, group.existing.index,drop = FALSE]
          net[[net.j]] <- values.new
        }
        # net[[net.j]] <- values.new
        if (net.j %in% c("centr")) {
          # Nvar x Ngroup x NLR
          values.new <- values[ , group.existing.index, ,drop = FALSE]
          net[[net.j]] <- values.new
        }
        if (net.j %in% c("count.cell","sum.cell","weight.cell")) {
          net[[net.j]] <- NULL
        }
        if (net.j %in% c("centr.cell")) {
          net[[net.j]] <- NULL
        }

        if (net.j %in% c("tmp")) {
          net[[net.j]] <- NULL
          net[["prob.cell"]] <- NULL
        }
      }
      net.subset[[i]] <- net

      # cat("Update slot object@netP...", '\n')
      # netP <- object@netP[[i]]
      # for (netP.j in names(netP)) {
      #   values <- netP[[netP.j]]
      #   if (netP.j %in% c("pathways")) {
      #     values.new <- values
      #     netP[[netP.j]] <- values.new
      #   }
      #   if (netP.j %in% c("prob")) {
      #     values.new <- values[group.existing.index, group.existing.index, ]
      #     netP[[netP.j]] <- values.new
      #   }
      #   if (netP.j %in% c("centr")) {
      #     for (k in 1:length(values)) {
      #       values.new <- lapply(values, function(x) {
      #         values.new2 <- lapply(x, function(x) {
      #           values.new1 <- x[group.existing.index]
      #           names(values.new1) <- group.existing
      #           return(values.new1)
      #         })
      #         names(values.new2) <- names(x)
      #         return(values.new2)
      #       })
      #       names(values.new) <- names(values)
      #     }
      #   }
      #   netP[[netP.j]] <- values.new
      # }
      netP = computeCommunProbPathway(net = net.subset[[i]], pairLR.use = object@LR[[i]]$LRsig, thresh = thresh)
      netP[["centr"]] = netAnalysis_computeCentrality(net =  netP[["prob"]])
      netP.subset[[i]] <- netP
      idents.subset[[i]] <- idents[[i]][names(idents[[i]]) %in% cells.use]
      idents.subset[[i]] <- factor(idents.subset[[i]], levels = levels(idents[[i]])[levels(idents[[i]]) %in% level.use])
    }
    idents.subset$joint <- factor(object@idents$joint[cells.use.index], levels = level.use)

  } else {
    cat(cli.symbol(),"Update slots object@images, object@net, object@netP in a single dataset...", '\n')

    group.existing <- level.use0[level.use0 %in% level.use]
    group.existing.index <- which(level.use0 %in% level.use)

    images <- object@images
    for (images.j in names(images)) {
      values <- images[[images.j]]
      if (images.j %in% c("coordinates")) {
        values.new <- values[cells.use.index, ,drop = F]
        images[[images.j]] <- values.new
      }
      if (images.j %in% c("distance")) {
        values.new <- values[group.existing.index, group.existing.index, drop = F]
        images[[images.j]] <- values.new
      }
      if (images.j %in% c("result.computeCellDistance")) {
        values[["d.spatial"]] <- values[["d.spatial"]][cells.use.index,cells.use.index,drop = F]
        values[["adj.contact"]] <- values[["adj.contact"]][cells.use.index,cells.use.index,drop = F]
        images[[images.j]] <- values
      }

    }
    images.subset <- images


    # cat("Update slot object@net...", '\n')
    net <- object@net
    for (net.j in names(net)) {

      # code for test
      # net.j <- "centr.cell"
      values <- net[[net.j]]
      if (net.j %in% c("prob","pval")) {
        values.new <- values[group.existing.index, group.existing.index, ,drop = FALSE]
        net[[net.j]] <- values.new
      }
      if (net.j %in% c("count","sum","weight")) {
        values.new <- values[group.existing.index, group.existing.index, drop = FALSE]
        net[[net.j]] <- values.new
      }
      if (net.j %in% c("centr")) {
        # Nvar x Ngroup x NLR
        values.new <- values[ , group.existing.index, ,drop = FALSE]
        net[[net.j]] <- values.new
      }
      if (net.j %in% c("count.cell","sum.cell","weight.cell")) {
        values.new <- values[cells.use.index, cells.use.index, drop = FALSE]
        net[[net.j]] <- values.new
      }
      if (net.j %in% c("centr.cell")) {
        # Nvar x Ncell x NLR
        values.new <- values[ , cells.use.index, ,drop = FALSE]
        net[[net.j]] <- values.new
      }

      if (net.j %in% c("tmp")) {
        if("cell.type.decomposition" %in% names(values)){
          new.CTD <- values[["cell.type.decomposition"]][cells.use.index,,drop=F]
          # update
          net[[net.j]][["cell.type.decomposition"]] <- new.CTD
        }
        if("prob.cell" %in% names(values)){
          new.prob.cell_ <- values[["prob.cell"]]
          new.prob.cell_ <- my_future_lapply(
            X = names(new.prob.cell_),
            FUN = function(i){
              return(new.prob.cell_[[i]][cells.use.index,cells.use.index,drop=F])
            },
            hint.message = "updating prob.cell..."
          )
          # update
          names(new.prob.cell_) <- names(values[["prob.cell"]])
          new.prob.cell <- my_as_sparse3Darray(new.prob.cell_)
          cat(cli.symbol(),"The number of cells and L-R pairs in Dim(Prob.cell):",dim(new.prob.cell),"\n")

          # set `Prob.cell`'s names
          dimnames(new.prob.cell) <- list(cells.use, cells.use,
                                      names(new.prob.cell_))
          net[[net.j]][["prob.cell"]] <- new.prob.cell_
          net[["prob.cell"]] <- new.prob.cell
        }
      }
    } # end for
    net.subset <- net

    # cat("Update slot object@netP...", '\n')
    # netP <- object@netP
    # for (netP.j in names(netP)) {
    #   values <- netP[[netP.j]]
    #   if (netP.j %in% c("pathways")) {
    #     values.new <- values
    #     netP[[netP.j]] <- values.new
    #   }
    #   if (netP.j %in% c("prob")) {
    #     values.new <- values[group.existing.index, group.existing.index, ]
    #     netP[[netP.j]] <- values.new
    #   }
    #   if (netP.j %in% c("centr")) {
    #     for (k in 1:length(values)) {
    #       values.new <- lapply(values, function(x) {
    #         values.new2 <- lapply(x, function(x) {
    #           values.new1 <- x[group.existing.index]
    #           names(values.new1) <- group.existing
    #           return(values.new1)
    #         })
    #         names(values.new2) <- names(x)
    #         return(values.new2)
    #       })
    #       names(values.new) <- names(values)
    #     }
    #   }
    #   netP[[netP.j]] <- values.new
    # }
    netP = computeCommunProbPathway(net = net.subset, pairLR.use = object@LR$LRsig, thresh = thresh)
    netP[["centr"]] = netAnalysis_computeCentrality(net = netP[["prob"]])
    netP.subset <- netP
    idents.subset <- object@idents[cells.use.index]
    idents.subset <- factor(idents.subset, levels = level.use)
  }

  new.options <- object@options
  new.options[["parameter"]][["subset"]] <- T
  new.options[["parameter"]][["original.idents"]] <- levels(object@idents)

  object.subset <- methods::new(
    Class = "SpatialCellChat",
    data = data.subset,
    data.signaling = data.signaling.subset,
    data.project = data.project.subset,
    images = images.subset,
    net = net.subset,
    netP = netP.subset,
    meta = meta.subset,
    idents = idents.subset,
    var.features = object@var.features,
    LR = object@LR,
    DB = object@DB,
    options = new.options
  )
  return(object.subset)
}


#' @title getIdent
#'
#' @param object A SpatialCellChat object
#'
#' @return Character string. SpatialCellChat's idents saved in the `idents` slot.
#' @export
level.SpatialCellChat <- function(object){

  if (!is.null(object@idents)) {
    res <- levels(object@idents)
  } else {
    stop(cli.symbol("fail"),"Please check your idents in your SpatialCellChat object!")
  }

  cat(cli.symbol(),"Idents in your SpatialCellChat object are ",cli::col_red(paste0(levels(object@idents),collapse = ", ")),"\n")
  return(res)
}

#' Get and set miscellaneous data
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @return miscellaneous data
#'
#' @rdname misc
#' @export misc
#'
#' @concept data-access
#'
misc <- function(object, ...) {
  UseMethod(generic = 'misc', object = object)
}

#' @param value Data to add
#'
#' @return An object with miscellaneous data added
#'
#' @rdname misc
#' @export misc<-
#'
"misc<-" <- function(object, ..., value) {
  UseMethod(generic = 'misc<-', object = object)
}


#' @param key Name of specific bit of meta data to pull
#'
#' @rdname misc
#' @export
#' @method misc SpatialCellChat
#'
misc.SpatialCellChat <- function(object, key=NULL, ...) {
  if(!is.null(key)){
    object@options[["misc"]][[key]]
  } else {
    object@options[["misc"]]
  }
}

#' @param key Name of specific bit of meta data to pull
#'
#' @rdname misc
#' @export
#' @method misc<- SpatialCellChat
#'
"misc<-.SpatialCellChat" <- function(object, key=NULL, ..., value) {
  if(is.null(key)){
    if(!is.list(value)) stop("Please use a named `list` to store key information.")
    else object@options[["misc"]] <- value
  } else object@options[["misc"]][[key]] <- value
  return(object)
}
