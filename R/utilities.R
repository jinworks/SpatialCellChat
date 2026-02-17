#' @title savePlotWithCow
#' @description
#' Save a plot directly from `Plots` window in Rstudio with the help of \code{cowplot}.
#' See details in \code{\link[cowplot]{save_plot}}.
#'
#' @inheritParams cowplot::save_plot
#'
#' @return NULL
#' @export
savePlotWithCow <- function(
    filename,
    plot=NULL,
    ncol = 1,
    nrow = 1,
    base_height = 3.71,
    base_asp = 1.618,
    base_width = NULL,
    ...,
    cols,
    rows,
    base_aspect_ratio,
    width,
    height
){
  if(is.null(plot)){
    plot <- grDevices::recordPlot()
    plot <- cowplot::plot_grid(plot)
  }
  cowplot::save_plot(
    filename = filename,
    plot = plot,
    ncol = ncol,
    nrow = nrow,
    base_height = base_height,
    base_asp = base_asp,
    base_width = base_width,
    ...=...,
    cols=cols,
    rows=rows,
    base_aspect_ratio=base_aspect_ratio,
    width=width,
    height=height
  )
}

#' Merge SpatialCellChat objects
#'
#' @param object.list  A list of multiple SpatialCellChat objects.
#' @param slot.name Which slot in SpatialCellChat used to merge.
#' @param pattern SpatialCellChat pattern used to merge. In detail, it comes from
#' `centr.cell` in net or netP slot.
#' If the `object.list` is a SpatialCellChat object, the functions will return a seurat
#' object combining both "indeg" and "outdeg" in `centr.cell`.
#' If the `object.list` is a SpatialCellChat object list, the functions will return a seurat
#' object combining "indeg" or "outdeg" in `centr.cell` from objects respectively.
#' @param normalize.chat Whether to normalize SpatialCellChat pattern with \code{\link[Seurat]{NormalizeData.Seurat}}.
#' @param normalization.method Inherit from \code{\link[Seurat]{NormalizeData.Seurat}}.
#' @param merge.chat Use it when the `object.list` is a SpatialCellChat object, whether to
#' merge both "indeg" and "outdeg" in `centr.cell` into only one assay data in seurat.
#'
#' @importFrom methods slot new
#'
#' @return SeuratObject.
#' @export
mergeChatIntoSeurat <- function(
    object.list,
    slot.name = "net",
    pattern = c("incoming","outgoing"),
    normalize.chat = T,
    normalization.method = "LogNormalize",
    merge.chat = F
){
  if(inherits(x = object.list, what = c("SpatialCellChat"))){
    object <- object.list
    if(slot.name%in%c("net","netP")){
      if (length(methods::slot(object, slot.name)[["centr.cell"]]) == 0) {
        stop("Please run `netAnalysis_computeCentrality` with `do.group=F` to compute the network centrality scores! \n")
      }
      centr.cell <- methods::slot(object, slot.name)[["centr.cell"]]

      seu.list <- lapply(
        X = pattern,
        FUN = function(i){
          if (i == "outgoing") {
            # data_sender is matrix of features-by-samples: nCells x nPathways(nLRs), same as previous nGroups x nPathways(nLRs)
            data_sender <- centr.cell[c("outdeg"),,,drop=T]
            data_sr <- Matrix::t(data_sender)
            # use Seurat
            seu_<- CellChat2Seurat(object,counts = data_sr,assay = "Spatial",check.object = F)
            Seurat::DefaultAssay(seu_) <- "Spatial"
            if(normalize.chat){
              seu_ <- Seurat::NormalizeData(object = seu_,normalization.method = normalization.method, verbose = F)
            }
          } else if (i == "incoming") {
            # data_receiver is matrix of features-by-samples: nCells x nPathways(nLRs), same as previous nGroups x nPathways(nLRs)
            data_receiver <- centr.cell[c("indeg"),,,drop=T]
            data_sr <- Matrix::t(data_receiver)
            # use Seurat
            seu_<- CellChat2Seurat(object,counts = data_sr,assay = "Spatial",check.object = F)
            Seurat::DefaultAssay(seu_) <- "Spatial"
            if(normalize.chat){
              seu_ <- Seurat::NormalizeData(object = seu_,normalization.method = normalization.method, verbose = F)
            }
          }
          return(seu_)
        }
      )
      names(seu.list) <- pattern

      if(merge.chat){
        direction.keys <- c("outgoing"="-O","incoming"="-I")

        Chat.data1 <- seu.list[[1]]@assays[["Spatial"]]$data
        Chat.meta.features1 <- seu.list[[1]]@assays[["Spatial"]]@misc[["features.info"]]
        Chat.meta.features1[["orig.features"]] <- dimnames(centr.cell)[[3]]
        Chat.meta.features1[["srt.features"]] <- rownames(Chat.data1)
        rownames(Chat.data1) <- paste0(rownames(Chat.data1),direction.keys[[pattern[[1]]]])

        Chat.data2 <- seu.list[[2]]@assays[["Spatial"]]$data
        Chat.meta.features2 <- seu.list[[2]]@assays[["Spatial"]]@misc[["features.info"]]
        Chat.meta.features2[["orig.features"]] <- dimnames(centr.cell)[[3]]
        Chat.meta.features2[["srt.features"]] <- rownames(Chat.data2)
        rownames(Chat.data2) <- paste0(rownames(Chat.data2),direction.keys[[pattern[[2]]]])

        Chat.data <- rbind(Chat.data1,Chat.data2)
        Chat.meta.features <- rbind(Chat.meta.features1,Chat.meta.features2)

        seu <- CellChat2Seurat(object,counts = Chat.data,assay = "Spatial",check.object = F)

        seu@assays[["incoming"]] <- seu.list[["incoming"]]@assays[["Spatial"]]
        seu@assays[["outgoing"]] <- seu.list[["outgoing"]]@assays[["Spatial"]]
        Seurat::DefaultAssay(seu) <- "Spatial"
        rownames(Chat.meta.features) <- rownames(seu)
        Chat.meta.features[["feature.names"]] <- rownames(seu)
        seu@assays[["Spatial"]]@misc[["features.info"]] <- Chat.meta.features

        Seurat::Misc(seu,slot = "Chat.info") <-
          paste0(
            "slot.name=",
            slot.name,
            " ",
            "pattern=",
            stringr::str_c(pattern, collapse = ","),
            " ",
            "merge.chat=T",
            " ",
            "normalize.chat=",
            as.character(normalize.chat)
          )
        # condition.keys <- c("outgoing"="Sender","incoming"="Receiver")
        # seu.list[[1]]$Condition <- condition.keys[[pattern[[1]]]]
        # seu.list[[2]]$Condition <- condition.keys[[pattern[[2]]]]
        # seu <- merge(seu.list[[1]],seu.list[[2]], add.cell.ids = condition.keys)
        # Seurat::DefaultAssay(seu) <- "Spatial"
        # seu$Condition <- factor(seu$Condition, levels = unique(seu$Condition))
        # names(seu@images) <- pattern
      } else{
        seu <- seu.list[["incoming"]]
        seu@assays[["outgoing"]] <- seu.list[["outgoing"]]@assays[["Spatial"]]
        seu@assays[["incoming"]] <- seu@assays[["Spatial"]]
        Seurat::DefaultAssay(seu) <- "Spatial"
        Seurat::Misc(seu,slot = "Chat.info") <-
          paste0(
            "slot.name=",
            slot.name,
            " ",
            "pattern=",
            stringr::str_c(pattern, collapse = ","),
            " ",
            "merge.chat=F",
            " ",
            "normalize.chat=",
            as.character(normalize.chat)
          )
      }
      tmp.meta.features <- seu@assays[[Seurat::DefaultAssay(seu)]]@misc[["features.info"]]
      tmp.meta.features[["interaction_name"]] <- tmp.meta.features[["orig.features"]]
      tmp.meta.features <-
        dplyr::left_join(
          tmp.meta.features,
          object@DB$interaction[, c("interaction_name", "pathway_name", "ligand", "receptor")],
          by = "interaction_name"
        )
      seu@assays[[Seurat::DefaultAssay(seu)]]@misc[["features.info"]] <- tmp.meta.features
      return(seu)
    } else{
      stop("Wrong slot.name!")
    }
  } else if(inherits(x = object.list, what = c("list"))){
    if(is.null(names(object.list))){
      names(object.list) <- paste0("COND",seq_along(object.list))
    }

    if(slot.name%in%c("net","netP")){

      pattern <- match.arg(pattern)

      seu.list <- lapply(
        X = seq_along(object.list),
        FUN = function(i){
          object <- object.list[[i]]
          if (length(methods::slot(object, slot.name)[["centr.cell"]]) == 0) {
            stop("Please run `netAnalysis_computeCentrality` with `do.group=F` to compute the network centrality scores! \n")
          }
          centr.cell <- methods::slot(object, slot.name)[["centr.cell"]]

          if (pattern == "outgoing") {
            # data_sender is matrix of features-by-samples: nCells x nPathways(nLRs), same as previous nGroups x nPathways(nLRs)
            data_sender <- centr.cell[c("outdeg"),,,drop=T]
            data_sr <- Matrix::t(data_sender)

          } else if (pattern == "incoming") {
            # data_receiver is matrix of features-by-samples: nCells x nPathways(nLRs), same as previous nGroups x nPathways(nLRs)
            data_receiver <- centr.cell[c("indeg"),,,drop=T]
            data_sr <- Matrix::t(data_receiver)

          }

          # use Seurat
          seu_<- CellChat2Seurat(object,counts = data_sr,assay = "Spatial",check.object = F)
          Seurat::DefaultAssay(seu_) <- "Spatial"
          seu_$Condition <- names(object.list)[[i]]

          return(seu_)
        }
      )
      seu <- merge(seu.list[[1]],seu.list[[2]], add.cell.ids = names(object.list))
      Seurat::DefaultAssay(seu) <- "Spatial"
      seu$Condition <- factor(seu$Condition, levels = unique(seu$Condition))
      names(seu@images) <- names(object.list)
      Seurat::Misc(seu,slot = "Chat.info") <- paste0("slot.name=",slot.name," ","pattern=",pattern)

    } else if (slot.name%in%c("data","data.signaling","data.scale","data.project","data.raw")){
      seu.list <- lapply(
        X = seq_along(object.list),
        FUN = function(i){
          object <- object.list[[i]]
          if ( 0 %in% dim(methods::slot(object, slot.name)) ) {
            stop("Please check your CellChat/SpatialCellChat objects. The slot input is an empty matrix. \n")
          }
          data.use <- methods::slot(object, slot.name)

          # use Seurat
          seu_<- CellChat2Seurat(object,counts = data.use,assay = "Spatial",check.object = F)
          Seurat::DefaultAssay(seu_) <- "Spatial"
          seu_$Condition <- names(object.list)[[i]]

          return(seu_)
        }
      )
      seu <- merge(seu.list[[1]],seu.list[[2]], add.cell.ids = names(object.list))
      Seurat::DefaultAssay(seu) <- "Spatial"
      seu$Condition <- factor(seu$Condition, levels = unique(seu$Condition))
      names(seu@images) <- names(object.list)
      Seurat::Misc(seu,slot = "Chat.info") <- paste0("slot.name=",slot.name)
    }
    return(seu)
  }
}


#' @title RunSeurat
#'
#' @param Obejct A seurat object
#' @param normalization.method See details of `normalization.method` in \code{\link[Seurat]{NormalizeData.Seurat}}.
#' @param scale.factor See details of `scale.factor` in \code{\link[Seurat]{NormalizeData.Seurat}}.
#' @param selection.nfeatures See details of `nfeatures` in \code{\link[Seurat]{FindVariableFeatures}}.
#' @param selection.method See details of `selection.method` in \code{\link[Seurat]{FindVariableFeatures}}.
#' @param do.center See details of `do.center` in \code{\link[Seurat]{ScaleData}}.
#' @param do.scale See details of `do.scale` in \code{\link[Seurat]{ScaleData}}.
#' @param PCA.npcs See details of `npcs`  in \code{\link[Seurat]{RunPCA}}.
#' @param k_Neighbors See details of `k.param` in \code{\link[Seurat]{FindNeighbors}}.
#' @param dims_Neighbors See details of `dims` in \code{\link[Seurat]{FindNeighbors}}.
#' @param dims_UMAP See details of `dims` in \code{\link[Seurat]{RunUMAP}}.
#' @param seed.use Seed used for random number generation in PCA & UMAP.
#' @param verbose Controls verbosity.
#'
#' @return A seurat object
#' @import Seurat
#' @export
RunSeurat <- function(Obejct,
                      normalization.method = "LogNormalize",
                      scale.factor = 10000,
                      selection.nfeatures = 2000,
                      selection.method = "vst",
                      do.center=T,
                      do.scale=T,
                      PCA.npcs = 50,
                      k_Neighbors = 20,
                      dims_Neighbors = 1:40,
                      dims_UMAP = 1:10,
                      seed.use=666L,
                      verbose = TRUE){
  if(!is.null(normalization.method)) {Obejct <- Seurat::NormalizeData(object = Obejct, normalization.method = normalization.method, scale.factor = scale.factor, verbose = verbose)}
  Obejct <- Seurat::FindVariableFeatures(object = Obejct, nfeatures = selection.nfeatures,selection.method = selection.method, verbose = verbose)
  Obejct <- Seurat::ScaleData(object = Obejct, verbose = verbose,do.scale = do.scale,do.center = do.center)
  Obejct <- Seurat::RunPCA(object = Obejct, features = VariableFeatures(Obejct), npcs = PCA.npcs,verbose = verbose,seed.use = seed.use)
  Obejct <- Seurat::FindNeighbors(object = Obejct, k.param = k_Neighbors,dims = dims_Neighbors, verbose = verbose)
  Obejct <- Seurat::RunUMAP(object = Obejct, dims = dims_UMAP, verbose = verbose,seed.use = seed.use)
}


#' @title Fast `==` operation for `sparse3Darray` Class
#'
#' @param arr1 A sparse3Darray
#' @param arr2 A sparse3Darray
#'
#' @return Boolean.
#' @export
`==.sparse3Darray` <- function(arr1, arr2) {
  if (!inherits(arr1, "sparse3Darray") || !inherits(arr2, "sparse3Darray")) {
    stop("Both arguments must be of class 'sparse3Darray'")
  }
  i.equal = all(arr1$i==arr2$i)
  j.equal = all(arr1$j==arr2$j)
  k.equal = all(arr1$k==arr2$k)
  x.equal = all(arr1$x==arr2$x)
  dim.equal = all(arr1$dim==arr2$dim)

  return(all(c(i.equal,j.equal,k.equal,x.equal,dim.equal)))
}

#' @title my_obj_size
#' @description
#' works similarly to object.size, but counts more accurately and includes the size of environments.
#'
#' @param ... Any R object. See details in \code{\link[pryr]{object_size}}.
#'
#' @return String. Object size easy to read.
#' @export
#'
#' @examples
#' \dontrun{
#' my_obj_size(list(1,2,3))
#' }
my_obj_size <- function(...){
  size_ <- capture.output(print(pryr::object_size(... = ...)))
  return(size_)
}

#' Set working directory to the path of currently opened file.
#' @description
#' Set working directory to the path of currently opened file (usually an R script). You can use this function in both .R/.Rmd files and R Console. RStudio (version >= 1.2) is required for running this function.
#'
#' @return
#' @export
#'
setpwd <- function(){
  base::setwd(base::dirname(rstudioapi::getActiveDocumentContext()$path))
}

#' Generate 4 kinds of cli symbols widely used in CellChat
#' @description
#' Generate 4 kinds of cli symbols: 'success','fail','info','output'
#'
#' @param symbol A string in `c('success','fail','info','output')` or a integer from `1:4` respectively
#' @importFrom cli col_green col_red col_grey col_br_red symbol
#'
#' @return
#' @export
#'
#' @examples
#' cli.symbol();cli.symbol(1);cli.symbol("info")
cli.symbol <- function(symbol=NULL){
  list.symbol <- list(success = cli::col_green(cli::symbol$tick),
                      fail = cli::col_red(cli::symbol$cross),
                      info = cli::col_grey(cli::symbol$info),
                      output = cli::col_br_red(">>> "))
  names_ <- names(list.symbol)
  if(is.null(symbol)){
    return(list.symbol[["output"]])
  } else{
    if(symbol%in%names_ & length(symbol)==1){
      return(list.symbol[[symbol]])
    } else if(symbol%in%(1:length(names_)) & length(symbol)==1) {
      return(list.symbol[[symbol]])
    } else {
      stop("Please check your input. Parameter `symbol` should be one string in c('success','fail','info','output') or a integer from 1:4 !")
    }
  }
}

#' Set the computation environment
#'
#' @param workers A integer greater than 0. See details in \code{\link[future]{multisession}}.
#' When workers = 1, CellChat will run in sequential mode;
#' Otherwise, CellChat will run in parallel mode with a respective number of future parallel clusters.
#' @param gc Boolean. Whether the garbage collector run or not. It's an argument passed to Future().
#' See details in \code{\link[future]{Future}}
#' @param conda_env We use some functions coded with Python in CellChat. Users
#' may need to set Python Lib Path here or otherwhere explicitly. For example, setEnvironment(workers = 2L,conda_env = "path/to/.conda/envs/env_name/")
#' See details in \code{\link[reticulate]{use_condaenv}} or use `reticulate`'s other relevant functions to get a Python support.
#' @param future.globals.maxSize Run \code{?future::future.globals.maxSize} to see details.
#'
#' @export
#'
setEnvironment <- function(workers,
                            future.globals.maxSize = 10000*1024^2,
                            gc = T,
                            conda_env = NULL) {
  if (is.null(conda_env)) {
    NULL
  } else {
    reticulate::use_condaenv(conda_env)
  }

  if (is.numeric(workers) & (workers > 0)) {
    if (workers == 1) {
      future::plan("sequential", gc = gc)
    } else {
      future::plan("multisession", workers = workers, gc = gc)
    }
    options(future.globals.maxSize = future.globals.maxSize)
    future.strategy <- ifelse(
      test = future::nbrOfWorkers() == 1,
      yes = "sequential",
      no = paste0("parallel with ", future::nbrOfWorkers(), " workers")
    )
    cat(paste0(
      cli.symbol("info"),
      " Future strategy in use: `",
      future.strategy,
      "`\n"
    ))

  } else {
    cat(cli.symbol("fail"), "Please check your `workers`'s number.")
  }
}

#' my_as_sparse3Darray
#' @description
#' A faster realization in comparison to \code{\link[spatstat.sparse]{as.sparse3Darray}}
#'
#' @param x Data in another format (see Details in \code{\link[spatstat.sparse]{as.sparse3Darray}}.
#' @param ... Ignored.
#' @param strict Boolean. Logical value specifying whether to enforce the rule that each entry in i,j,k,x refers to a different cell. If strict=TRUE, entries which refer to the same cell in the array will be reduced to a single entry by summing the x values. Default is strict=FALSE.
#' @param nonzero Boolean. Logical value specifying whether to remove any entries of x which equal zero.
#' @importFrom spatstat.sparse sparse3Darray as.sparse3Darray
#' @return Sparse three-dimensional array (object of class "sparse3Darray").
#' @export
#'
my_as_sparse3Darray <- function (x, ... ,strict = FALSE, nonzero = FALSE)
{
  if (inherits(x, "sparse3Darray")) {
    y <- x
  }
  else if (inherits(x, c("matrix", "sparseMatrix"))) {
    z <- as(x, Class = "TsparseMatrix")
    dn <- dimnames(x)
    dn <- if (is.null(dn))
      NULL
    else c(dn, list(NULL))
    one <- if (length(z@i) > 0)
      1L
    else integer(0)
    y <- sparse3Darray(i = z@i + 1L, j = z@j + 1L, k = one,
                       x = z@x, dims = c(dim(x), 1L), dimnames = dn,strict = strict, nonzero = nonzero)
  }
  else if (is.array(x)) {
    stopifnot(length(dim(x)) == 3)
    dimx <- dim(x)
    if (prod(dimx) == 0) {
      y <- sparse3Darray(, dims = dimx, dimnames = dimnames(x))
    }
    else {
      ijk <- which(x != RelevantZero(x), arr.ind = TRUE)
      ijk <- cbind(as.data.frame(ijk), x[ijk])
      y <- sparse3Darray(i = ijk[, 1L], j = ijk[, 2L],
                         k = ijk[, 3L], x = ijk[, 4L], dims = dimx, dimnames = dimnames(x),strict = strict, nonzero = nonzero)
    }
  }
  else if (inherits(x, "sparseVector")) {
    one <- if (length(x@i) > 0)
      1L
    else integer(0)
    y <- sparse3Darray(i = x@i, j = one, k = one, x = x@x,
                       dims = c(x@length, 1L, 1L),strict = strict, nonzero = nonzero)
  }
  else if (is.null(dim(x)) && is.atomic(x)) {
    n <- length(x)
    dn <- names(x)
    if (!is.null(dn))
      dn <- list(dn, NULL, NULL)
    one <- if (n > 0)
      1L
    else integer(0)
    y <- sparse3Darray(i = seq_len(n), j = one, k = one,
                       x = x, dims = c(n, 1L, 1L), dimnames = dn,strict = strict, nonzero = nonzero)
  }
  else if (is.list(x) && length(x) > 0) {
    n <- length(x)
    if (all(sapply(x, is.matrix))) {
      z <- Reduce(abind, x)
      y <- as.sparse3Darray(z)
    }
    else if (all(sapply(x, inherits, what = "sparseMatrix"))) {
      dimlist <- unique(lapply(x, dim))
      if (length(dimlist) > 1)
        stop("Dimensions of matrices do not match")
      dimx <- c(dimlist[[1L]], n)
      dnlist <- lapply(x, dimnames)
      isnul <- sapply(dnlist, is.null)
      dnlist <- unique(dnlist[!isnul])
      if (length(dnlist) > 1)
        stop("Dimnames of matrices do not match")
      dn <- if (length(dnlist) == 0)
        NULL
      else c(dnlist[[1L]], list(NULL))

      # progressr::handlers(global = T)
      progressr::handlers(progressr::handler_progress(format = ":percent [:bar] :eta :message"))
      X <- seq_len(n)
      p <- progressr::progressor(along = X)
      df <- purrr::map_dfr(
        .x = X,
        .f = function(k) {
          mk <- as(x[[k]], "TsparseMatrix")
          kvalue <- if (length(mk@i) > 0) k else integer(0)
          dfk <- data.frame(
            i = mk@i + 1L,
            j = mk@j + 1L,
            k = kvalue,
            x = mk@x
          )
          p("as sparse3Darray...")
          return(dfk)
        }
      )

      y <- sparse3Darray(i = df$i, j = df$j, k = df$k,
                         x = df$x, dims = dimx, dimnames = dn,strict = strict, nonzero = nonzero)
    }
    else {
      warning("I don't know how to convert a list to a sparse array")
      return(NULL)
    }
  }
  else {
    warning("I don't know how to convert x to a sparse array")
    return(NULL)
  }
  return(y)
}



#' my_future_sapply
#'
#' @description
#' future_sapply with progress bar and stable random number generation
#'
#' @param X A vector-like object to iterate over.
#' @param FUN A function taking at least one argument.
#' @param ... (optional) Additional arguments passed to `FUN()`.
#' For `future_*apply()` functions and `replicate()`, any `future.*` arguments
#' part of `\ldots` are passed on to `future_lapply()` used internally.
#' @param simplify See \code{\link[base]{sapply}}, respectively.
#' @param USE.NAMES See \code{\link[base]{sapply}}.
#' @param future.envir An [environment] passed as argument `envir` to
#'        \code{\link[future]{future}} as-is.
#' @param future.seed A logical or an integer (of length one or seven),
#'        or a list of `length(X)` with pre-generated random seeds.
#'        For details, see \code{\link[future.apply]{future_lapply}}.
#' @param future.label If a character string, then each future is assigned
#'        a label `sprintf(future.label, chunk_idx)`.  If TRUE, then the
#'        same as `future.label = "future_sapply-%d"`.  If FALSE, no labels
#'        are assigned.
#' @param hint.message A hint message shown after the progress bar. By default,
#' the message is "Computing..."
#' @param strategy.message Whether to show the future stretegy information
#'
#' @return
#' For `future_sapply()`, a vector with same length and names as \code{X}.
#' See \code{\link[base]{sapply}} for reference.
#' @export
#'
#' @examples
#' my_future_sapply(X=1:10,FUN = function(x){Sys.sleep(11-x);sqrt(x)})
my_future_sapply <- function(X, FUN, ..., simplify = TRUE, USE.NAMES = TRUE, future.envir = parent.frame(),
                             future.seed = T, future.label = "future_sapply-%d", hint.message = "Computing...", strategy.message = TRUE) {
  future.strategy <- ifelse(test = future::nbrOfWorkers() == 1,
                            yes = "sequential",
                            no = paste0("parallel with ",future::nbrOfWorkers()," workers"))
  if (strategy.message) {
    cat(paste0(cli.symbol("info")," Future strategy in use: `",future.strategy,"`\n"))
  }

  # progressr::handlers(global = T)
  progressr::handlers(progressr::handler_progress(format = ":percent [:bar] :eta :message"))
  p <- progressr::progressor(along = X)

  FUN <- match.fun(FUN)
  answer <- future.apply::future_lapply(
    X = X,
    FUN = function(x){
      res <- FUN(x)
      p(hint.message)
      return(res)
    },
    ...,
    future.envir = future.envir,
    future.seed = future.seed,
    future.label = future.label
  )

  if (USE.NAMES && is.character(X) && is.null(names(answer)))
    names(answer) <- X
  if (!isFALSE(simplify) && length(answer))
    simplify2array(answer, higher = (simplify == "array"))
  else answer
}

#' my_future_lapply
#' @description
#' future_lapply with progress bar and stable random number generation
#' @param X A vector-like object to iterate over.
#' @param FUN A function taking at least one argument.
#' @param ... (optional) Additional arguments passed to `FUN()`.
#' For `future_*apply()` functions and `replicate()`, any `future.*` arguments
#' part of `\ldots` are passed on to `future_lapply()` used internally.
#' @param future.seed A logical or an integer (of length one or seven),
#'        or a list of `length(X)` with pre-generated random seeds.
#'        For details, see \code{\link[future.apply]{future_lapply}}.
#' @param simplify different from the `simplify` parameter in `my_future_sapply`.
#'        If simplify = TRUE, we run \code{\link[base]{unlist}} upon the result.
#' @param hint.message A hint message shown after the progress bar. By default,
#' the message is "Computing..."
#' @param strategy.message Whether to show the future stretegy information
#' @param future.label If a character string, then each future is assigned
#'        a label `sprintf(future.label, chunk_idx)`.  If TRUE, then the
#'        same as `future.label = "future_lapply-%d"`.  If FALSE, no labels
#'        are assigned.
#'
#' @return
#' For `future_lapply()`, a list with same length and names as `X`.
#' See \code{\link[base]{lapply}} for reference.
#' @export
#'
#' @examples
#' my_future_lapply(X=1:10,FUN = function(x){Sys.sleep(11-x);sqrt(x)})
my_future_lapply <- function(X, FUN, ..., future.seed = T, simplify = F, hint.message = "Computing...",strategy.message = TRUE, future.label = "future_lapply-%d") {
  future.strategy <- ifelse(test = future::nbrOfWorkers() == 1,
                            yes = "sequential",
                            no = paste0("parallel with ",future::nbrOfWorkers()," workers"))
  if (strategy.message) {
    cat(paste0(cli.symbol("info")," Future strategy in use: `",future.strategy,"`\n"))
  }

  # progressr::handlers(global = T)
  progressr::handlers(progressr::handler_progress(format = ":percent [:bar] :eta :message"))
  p <- progressr::progressor(along = X)
  res_ <- future.apply::future_lapply(
    X = X,
    FUN = function(x){
      res <- FUN(x)
      p(hint.message)
      return(res)
    },
    ...,
    future.seed = future.seed,
    future.label = future.label
  )
  if(simplify){
    res_ <- unlist(res_)
  } else{
    NULL
  }
  return(res_)
}

#' Convert a CellChat object into a Seurat object
#'
#' @param object CellChat object
#' @param counts Alternatively, use user-input `counts` for creating Seurat object instead of `object@data`
#' @param group.by Name of one metadata columns to group (color) cells. Default is the defined cell groups in CellChat object
#' @param assay Name of the initial assay. Default is "Spatial" when working on spatial imaging data
#' @param check.object whether checking the converted Seurat object. Set check.object = FALSE when not visualizing spatial clustering
#'
#' @return a Seurat object
#' @export
#'
CellChat2Seurat <- function(object,counts = NULL, group.by = NULL, assay = "Spatial", check.object = TRUE) {
  meta <- object@meta
  meta$group.cellchat <- object@idents
  if (!identical(rownames(meta), colnames(object@data))) {
    cat("The cell barcodes in 'meta' is ", head(rownames(meta)),'\n')
    warning("The cell barcodes in 'meta' is different from those in the used data matrix.
              We now simply assign the colnames in the data matrix to the rownames of 'meta'!")
    rownames(meta) <- colnames(object@data)
  }

  if(is.null(counts)){
    data.use <- object@data
  } else {
    if(is.null(colnames(counts))|is.null(rownames(counts))){
      stop("Please check the dimnames of `counts`!")
    }
    data.use <- counts
  }

  obj.new = Seurat::CreateSeuratObject(counts = data.use, meta.data = meta, assay = assay)
  if (assay == "Spatial") {
    spatial_locs <- object@images$coordinates
    obj.new@images$image =  new(
      Class = 'SlideSeq',
      assay = "Spatial",
      key = "image_",
      coordinates = spatial_locs
    )
  }
  if (is.null(group.by)) {
    group.by <- "group.cellchat"
  }
  Seurat::Idents(obj.new) <- group.by
  obj.new@assays[[SeuratObject::DefaultAssay(obj.new)]]@misc[["features.info"]] <-
    data.frame(
      "orig.features"=rownames(data.use)
    )

  if (check.object == TRUE) {
    cat("Check the converted Seurat object by visualizing spatial clustering...")
    pt.size.factor <- 3.5
    color.use <- scPalette(nlevels(obj.new))
    names(color.use) <- levels(obj.new)
    if (assay == "Spatial") {
      print(Seurat::SpatialDimPlot(obj.new, pt.size.factor = pt.size.factor, cols = color.use, label = FALSE) +
              theme(legend.title = element_blank(), legend.text = element_text(size = 10))+
              theme(legend.key = element_blank()))
    }
  }
  return(obj.new)
}

#' @title preProcessing
#' @description
#' Use this function to pre-process the data matrix input before running \link{createSpatialCellChat},
#' or pre-process the SpatialCellChat object following \link{subsetData} function.
#' @param object Matrix or SpatialCellChat object.
#' @param slot.name Character. By default "data.signaling". Use "data.signaling" or "data" when `object` is a SpatialCellChat object.
#'
#' @return Matrix or SpatialCellChat object.
#' @export
preProcessing <- function(object,slot.name=c("data.signaling","data")){
  if (inherits(x = object, what = c("matrix", "Matrix", "dgCMatrix"))) {
    cat(cli.symbol("info"),"Pre-processing from a data matrix.\n")
    data <- object %>% Matrix::t()
  } else if (is(object,"SpatialCellChat")) {
    cat(cli.symbol("info"),"Pre-processing from a SpatialCellChat object.\n")
    slot.name <- match.arg(slot.name)
    data <- methods::slot(object,slot.name) %>% Matrix::t()
    if(sum(dim(data))==0){
      stop(cli.symbol("fail"),"Please check the object's ",slot.name," slot. No data is in it.")
    }
  }

  k.boot <- sapply(X = 1:10,FUN = function(x){x*1;ALRA::choose_k(data)$k})
  k <- median(k.boot,na.rm = T)

  # a matrix where the cells are rows and genes are columns.
  capture.output({data.alra <- ALRA::alra(as.array(data),k = k,quantile.prob = 1e-5)})
  data.alra <- data.alra[[3]] %>% as(.,Class="CsparseMatrix")
  colnames(data.alra) <- colnames(data)
  rownames(data.alra) <- rownames(data)

  if (inherits(x = object, what = c("matrix", "Matrix", "dgCMatrix"))) {
    cat(cli.symbol("success"),"Pre-processing is done.\n")
    return(Matrix::t(data.alra))
  } else if (is(object,"SpatialCellChat")) {
    methods::slot(object,slot.name) <- Matrix::t(data.alra)
    object@options[["alra.k"]] <- k
    cat(cli.symbol("success"),"Pre-processing is done.\n")
    return(object)
  }
}


#' Normalize data using a scaling factor
#'
#' @param data.raw input raw data
#' @param scale.factor the scaling factor used for each cell
#' @param do.log whether do log transformation with pseudocount 1
#' @export
#'
normalizeData <- function(data.raw, scale.factor = 10000, do.log = TRUE) {
  # Scale counts within a sample
  library.size <- Matrix::colSums(data.raw)
  #scale.factor <- median(library.size)
  expr <- Matrix::t(Matrix::t(data.raw) / library.size) * scale.factor
  if (do.log) {
    data.norm <-log1p(expr)
  }
  return(data.norm)
}


#' Scale the data
#'
#' @param data.use input data
#' @param do.center whether center the values
#' @export
#'
scaleData <- function(data.use, do.center = T) {
  data.use <- Matrix::t(scale(Matrix::t(data.use), center = do.center, scale = TRUE))
  return(data.use)
}


#' Scale a data matrix
#'
#' @param x data matrix
#' @param scale the method to scale the data
#' @param na.rm whether remove na
#' @importFrom Matrix rowMeans colMeans rowSums colSums
#' @return
#' @export
#'
#' @examples
scaleMat <- function(x, scale, na.rm=TRUE){

  av <- c("none", "row", "column", 'r1', 'c1')
  i <- pmatch(scale, av)
  if(is.na(i) )
    stop("scale argument shoud take values: 'none', 'row' or 'column'")
  scale <- av[i]

  switch(
    scale,
    none = x,
    row = {
      x <-
        sweep(x, 1L, rowMeans(x, na.rm = na.rm), '-', check.margin = FALSE)
      sx <- apply(x, 1L, sd, na.rm = na.rm)
      sweep(x, 1L, sx, "/", check.margin = FALSE)
    },
    column = {
      x <-
        sweep(x, 2L, colMeans(x, na.rm = na.rm), '-', check.margin = FALSE)
      sx <- apply(x, 2L, sd, na.rm = na.rm)
      sweep(x, 2L, sx, "/", check.margin = FALSE)
    },
    r1 = sweep(x, 1L, rowSums(x, na.rm = na.rm), '/', check.margin = FALSE),
    c1 = sweep(x, 2L, colSums(x, na.rm = na.rm), '/', check.margin = FALSE)

  )
}

#' Downsampling single cell data using geometric sketching algorithm
#'
#' USERs need to install the python package `pip install geosketch` (https://github.com/brianhie/geosketch)
#'
#' @param object A data matrix (should have row names; samples in rows, features in columns) or a Seurat object.
#'
#' When object is a PCA or UMAP space, please set `do.PCA = FALSE`
#'
#' When object is a data matrix (cells in rows and genes in columns), it is better to use the highly variable genes. PCA will be done on this input data matrix.
#' @param percent the percent of data to sketch
#' @param idents A vector of identity classes to keep for sketching
#' @param do.PCA whether doing PCA on the input data
#' @param dimPC the number of components to use
#' @importFrom reticulate import
#' @return A vector of cell names to use for downsampling
#' @export
#'
sketchData <-
  function(object,
           percent,
           idents = NULL,
           do.PCA = TRUE,
           dimPC = 30) {
    # pip install geosketch
    geosketch <- reticulate::import('geosketch')
    if (is(object, "Seurat")) {
      sketch.size <- as.integer(percent * ncol(object))
      if (!is.null(idents)) {
        object <- subset(object, idents = idents)
      }
      object <- object %>% #Seurat::NormalizeData(verbose = FALSE) %>%
        FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
        RunPCA(pc.genes = object@var.genes,
               npcs = dimPC,
               verbose = FALSE)

      X.pcs <- object@reductions$pca@cell.embeddings
      cells.all <- Cells(object)

    } else {
      # Get top PCs
      if (do.PCA) {
        X.pcs <- runPCA(object, dimPC = dimPC)
      } else {
        X.pcs <- object
      }

      # Sketch percent of data.
      sketch.size <- as.integer(percent * nrow(X))
      cells.all <- rownames(object)
    }
    sketch.index <- geosketch$gs(X.pcs, sketch.size)
    sketch.index <- unlist(sketch.index) + 1
    sketch.cells <- cells.all[sketch.index]
    return(sketch.cells)
  }


#' Add the cell information into meta slot
#'
#' @param object CellChat object
#' @param meta cell information to be added
#' @param meta.name the name of column to be assigned
#'
#' @return
#' @export
#'
#' @examples
addMeta <- function(object, meta, meta.name = NULL) {
  if (is.null(x = meta.name) && is.atomic(x = meta)) {
    stop("'meta.name' must be provided for atomic meta types (eg. vectors)")
  }
  if (inherits(x = meta, what = c("matrix", "Matrix"))) {
    meta <- as.data.frame(x = meta)
  }

  if (is.null(x = meta.name)) {
    meta.name <- names(meta)
  } else {
    names(meta) <- meta.name
  }
  object@meta <- meta
  return(object)
}


#' Set the default identity of cells
#' @param object CellChat object
#' @param ident.use the name of the variable in object.meta;
#' @param levels set the levels of factor
#' @param display.warning whether display the warning message
#' @return
#' @export
#'
#' @examples
setIdent <- function(object, ident.use = NULL, levels = NULL, display.warning = TRUE){
  if (!is.null(ident.use)) {
    object@idents <- as.factor(object@meta[[ident.use]])
  }

  if (!is.null(levels)) {
    object@idents <- factor(object@idents, levels = levels)
  }
  if ("0" %in% as.character(object@idents)) {
    stop("Cell labels cannot contain `0`! ")
  }
  if (length(object@net) > 0) {
    if (all(dimnames(object@net$prob)[[1]] %in% levels(object@idents) )) {
      message("Reorder cell groups! ")
      cat("The cell group order before reordering is ", dimnames(object@net$prob)[[1]],'\n')
      # idx <- match(dimnames(object@net$prob)[[1]], levels(object@idents))
      idx <- match(levels(object@idents), dimnames(object@net$prob)[[1]])
      object@net$prob <- object@net$prob[idx, , ]
      object@net$prob <- object@net$prob[, idx, ]
      object@net$pval <- object@net$pval[idx, , ]
      object@net$pval <- object@net$pval[, idx, ]
      cat("The cell group order after reordering is ", dimnames(object@net$prob)[[1]],'\n')
    } else {
      message("Rename cell groups but do not change the order! ")
      cat("The cell group order before renaming is ", dimnames(object@net$prob)[[1]],'\n')
      dimnames(object@net$prob) <- list(levels(object@idents), levels(object@idents), dimnames(object@net$prob)[[3]])
      dimnames(object@net$pval) <- dimnames(object@net$prob)
      cat("The cell group order after renaming is ", dimnames(object@net$prob)[[1]],'\n')
    }
    if (display.warning) {
      warning("All the calculations after `computeCommunProb` should be re-run!!
    These include but not limited to `computeCommunProbPathway`,`aggregateNet`, and `netAnalysis_computeCentrality`.")
    }


  }
  return(object)
}


#' Update and re-order the cell group names after running `computeCommunProb`
#'
#' @param object CellChat object
#' @param old.cluster.name A vector defining old cell group labels in `object@idents`; Default = NULL, which will use `levels(object@idents)`
#' @param new.cluster.name A vector defining new cell group labels to rename
#' @param new.order reset order of cell group labels
#' @param new.cluster.metaname assign a name of the new labels, which will be the column name of new labels in `object@meta`
#' @return An updated CellChat object
#' @export
#'
updateClusterLabels <- function(object, old.cluster.name = NULL, new.cluster.name = NULL, new.order = NULL, new.cluster.metaname = "new.labels") {
  if (is.null(old.cluster.name)) {
    old.cluster.name <- levels(object@idents)
  }
  if (new.cluster.metaname %in% colnames(object@meta)) {
    stop("Please define another `new.cluster.metaname` as it exists in `colnames(object@meta)`!")
  }
  if (!is.null(new.cluster.name)) {
    labels.new <- plyr::mapvalues(object@idents, from = old.cluster.name, to = new.cluster.name)
    object@meta[[new.cluster.metaname]] <- labels.new
    object <- setIdent(object, ident.use = new.cluster.metaname, display.warning = FALSE)
  } else {
    new.cluster.metaname <- NULL
    cat("Only reorder cell groups but do not rename cell groups!")
  }

  if (!is.null(new.order)) {
    object <- setIdent(object, ident.use = new.cluster.metaname, levels = new.order, display.warning = FALSE)
  }
  message("We now re-run computeCommunProbPathway`,`aggregateNet`, and `netAnalysis_computeCentrality`...")
  object <- computeCommunProbPathway(object)
  ## calculate the aggregated network by counting the number of links or summarizing the communication probability
  object <- aggregateNet(object)
  # network importance analysis
  object <-netAnalysis_computeCentrality(object, slot.name = "netP")
  return(object)
}





#' Subset the expression data of signaling genes for saving computation cost
#'
#' @param object CellChat object
#' @param features default = NULL: subset the expression data of signaling genes in CellChatDB.use
#'
#' @return An updated CellChat object by assigning a subset of the data into the slot `data.signaling`
#' @export
#'
subsetData <- function(object, features = NULL) {
  interaction_input <- object@DB$interaction
  if (object@options$datatype != "RNA") {
    if ("annotation" %in% colnames(interaction_input) == FALSE) {
      warning("A column named `annotation` is required in `object@DB$interaction` when running CellChat on spatial transcriptomics! The `annotation` column is now automatically added and all L-R pairs are assigned as `Secreted Signaling`, which means that these L-R pairs are assumed to mediate diffusion-based cellular communication.")
      interaction_input$annotation <- "Secreted Signaling"
    }
  }
  if ("annotation" %in% colnames(interaction_input) == TRUE) {
    if (length(unique(interaction_input$annotation)) > 1) {
      interaction_input$annotation <- factor(interaction_input$annotation, levels = c("Secreted Signaling", "ECM-Receptor", "Non-protein Signaling", "Cell-Cell Contact"))
      interaction_input <- interaction_input[order(interaction_input$annotation), , drop = FALSE]
      interaction_input$annotation <- as.character(interaction_input$annotation)
    }
    object@DB$interaction <- interaction_input
  }

  if (is.null(features)) {
    DB <- object@DB
    gene.use_input <- extractGene(DB)
    gene.use <- intersect(gene.use_input, rownames(object@data))
  } else {
    gene.use <- intersect(features, rownames(object@data))
  }
  object@data.signaling <- object@data[rownames(object@data) %in% gene.use, ]
  return(object)
}


#' Identify over-expressed signaling genes associated with each cell group or spatially variable features independent of cell groups
#'
#' USERS can use customized gene set as over-expressed signaling genes by setting `object@var.features[[features.name]] <- features.sig`
#' The Bonferroni corrected/adjusted p value can be obtained via `object@var.features[[paste0(features.name, ".info")]]`. Note that by default `features.name = "features"`
#'
#' @param object CellChat object
#' @param do.grid Boolean. Whether to do "grid" operation to speed up computation.
#' @param data.use a customed data matrix. Default: data.use = NULL and the expression matrix in the slot 'data.signaling' is used
#' @param selection.method Method for selecting (spatially) variable features.
#' @param group.by cell group information; default is `object@idents`; otherwise it should be one of the column names of the meta slot
#' @param idents.use a subset of cell groups used for analysis
#' @param invert whether invert the idents.use
#' @param group.dataset dataset origin information in a merged CellChat object; set it as one of the column names of meta slot when identifying the highly enriched genes in one dataset for each cell group
#' @param pos.dataset the dataset name used for identifying highly enriched genes in this dataset for each cell group
#' @param features.name a char name used for storing the over-expressed signaling genes in `object@var.features[[features.name]]`
#' @param only.pos Only return positive markers
#' @param features features used for identifying Over Expressed genes. default use all features
#' @param return.object whether return the object; otherwise return a data frame consisting of over-expressed signaling genes associated with each cell group
#' @param thresh.pc Threshold of the percent of cells expressed in one cluster
#' @param thresh.fc Threshold of Log Fold Change
#' @param thresh.p Threshold of p-values
#' @inheritParams makeGridSpatialCellChat
#'
#' @importFrom stats sd wilcox.test
#' @importFrom stats p.adjust
#'
#' @return A CellChat object or a data frame. If returning a CellChat object, two new elements named 'features.name' and paste0(features.name, ".info") will be added into the list `object@var.features`
#' `object@var.features[[features.name]]` is a vector consisting of the identified over-expressed signaling genes;
#' `object@var.features[[paste0(features.name, ".info")]]` is a data frame returned from the differential expression analysis / the simple filtering
#' @export
#'
identifyOverExpressedGenes <- function(
    object,
    do.grid = F,
    cellsize = c(5, 5),
    what = "polygons",
    square = T,
    data.use = NULL,
    selection.method = c("wilcox", "moransi", "meringue"),
    features.name = "features",
    group.by = NULL,
    idents.use = NULL,
    invert = FALSE,
    group.dataset = NULL,
    pos.dataset = NULL,
    only.pos = TRUE,
    features = NULL,
    return.object = TRUE,
    thresh.pc = 0,
    thresh.fc = 0,
    thresh.p = 0.05

){
  selection.method <- match.arg(selection.method)

  # For code test
  # object = mousebrainGrid
  # data.use=NULL
  # # selection.method = c("wilcox", "moransi", "nnSVG", "simple")
  # selection.method = c("moransi")
  # features.name = "features"
  # group.by = NULL
  # idents.use = NULL
  # invert = FALSE
  # group.dataset = NULL
  # pos.dataset = NULL
  # only.pos = TRUE
  # features = NULL
  # return.object = TRUE
  # thresh.pc = 0
  # thresh.fc = 0
  # thresh.p = 0.05

  if(do.grid){
    # make grid SpatialCellChat
    grid.object <- makeGridSpatialCellChat(
      object = object,
      cellsize = cellsize,
      what = what,
      square = square
    )
    grid.object@DB <- object@DB
    grid.object <- subsetData(grid.object)


    # recursive definition
    grid.object <- identifyOverExpressedGenes(
      object = grid.object,
      do.grid=F,
      data.use = NULL,
      selection.method = selection.method,
      features.name = features.name,
      group.by = group.by,
      idents.use = idents.use,
      invert = invert,
      group.dataset = group.dataset,
      pos.dataset = pos.dataset,
      only.pos = only.pos,
      features = features,
      return.object = TRUE,
      thresh.pc = thresh.pc,
      thresh.fc = thresh.fc,
      thresh.p = thresh.p
    )

    if (return.object) {
      object@var.features <- grid.object@var.features
      misc(object,"grid.object") <- grid.object
      return(object)
    } else {
      markers.all = grid.object@var.features$features.info
      return(markers.all)
    }

  } else {
    if (!is.list(object@var.features)) {
      stop("Please update your CellChat object via `updateCellChat()`")
    }
    if (is.null(data.use)) {
      X <- object@data.signaling
      if (nrow(X) < 3) {stop("Please check `object@data.signaling` and ensure that you have run `subsetData` and that the data matrix `object@data.signaling` looks OK.")}
    } else {
      X <- data.use
    }

    if (is.null(features)) {
      features.use <- row.names(X)
    } else {
      features.use <- intersect(features, row.names(X))
    }
    data.use <- X[features.use,]
    data.use <- as.matrix(data.use)


    if (selection.method == "wilcox") {
      cat(cli.symbol(),"Choose variable features via Wilcox test using available cell group information...","\n")
      ### perform differential expression analysis using available cell group information
      if (is.null(group.by)) {
        labels <- object@idents
        if (!is.factor(labels)) {
          message("Use the joint cell labels from the merged CellChat object")
          labels <- object@idents$joint
        }
      } else {
        labels <- object@meta[[group.by]]
      }
      if (!is.factor(labels)) {
        labels <- factor(labels)
      }
      level.use <- levels(labels)[levels(labels) %in% unique(labels)]
      if (!is.null(idents.use)) {
        if (invert) {
          level.use <- level.use[!(level.use %in% idents.use)]
        } else {
          level.use <- level.use[level.use %in% idents.use]
        }
      }
      numCluster <- length(level.use)

      if (!is.null(group.dataset)) {
        labels.dataset <- as.character(object@meta[[group.dataset]])
        if (!(pos.dataset %in% unique(labels.dataset))) {
          cat(cli.symbol(),"Please set pos.dataset to be one of the following dataset names: ", unique(as.character(labels.dataset)))
          stop()
        }
      }

      mean.fxn <- function(x) {
        return(log(x = mean(x = expm1(x = x)) + 1))
      }
      labels <- as.character(labels)
      genes.de <- vector("list", length = numCluster)
      for (i in 1:numCluster) {
        features <- features.use
        if (is.null(group.dataset)) {
          cell.use1 <- which(labels == level.use[i])
          cell.use2 <- base::setdiff(1:length(labels), cell.use1)
        } else {
          cell.use1 <- which((labels == level.use[i]) & (labels.dataset == pos.dataset))
          cell.use2 <- which((labels == level.use[i]) & (labels.dataset != pos.dataset))
        }

        # feature selection (based on percentages)
        thresh.min <- 0
        pct.1 <- round(
          x = rowSums(data.use[features, cell.use1, drop = FALSE] > thresh.min) /
            length(x = cell.use1),
          digits = 3
        )
        pct.2 <- round(
          x = rowSums(data.use[features, cell.use2, drop = FALSE] > thresh.min) /
            length(x = cell.use2),
          digits = 3
        )
        data.alpha <- cbind(pct.1, pct.2)
        colnames(x = data.alpha) <- c("pct.1", "pct.2")
        alpha.min <- apply(X = data.alpha, MARGIN = 1, FUN = max)
        names(x = alpha.min) <- rownames(x = data.alpha)
        features <- names(x = which(x = alpha.min > thresh.pc))
        if (length(x = features) == 0) {
          #stop("No features pass thresh.pc threshold")
          next
        }

        # feature selection (based on average difference)
        data.1 <- apply(X = data.use[features, cell.use1, drop = FALSE],MARGIN = 1,FUN = mean.fxn)
        data.2 <- apply(X = data.use[features, cell.use2, drop = FALSE],MARGIN = 1,FUN = mean.fxn)
        FC <- (data.1 - data.2)
        if (only.pos) {
          features.diff <- names(which(FC > thresh.fc))
        } else {
          features.diff <- names(which(abs(FC) > thresh.fc))
        }

        features <- intersect(x = features, y = features.diff)
        if (length(x = features) == 0) {
          #  stop("No features pass thresh.fc threshold")
          next
        }

        data1 <- data.use[features, cell.use1, drop = FALSE]
        data2 <- data.use[features, cell.use2, drop = FALSE]

        pvalues <- unlist(
          x = my_future_lapply(
            X = 1:nrow(x = data1),
            FUN = function(x) {
              # return(wilcox.test(data1[x, ], data2[x, ], alternative = "greater")$p.value)
              return(wilcox.test(data1[x, ], data2[x, ])$p.value)
            }
          )
        )

        pval.adj = stats::p.adjust(
          p = pvalues,
          method = "bonferroni",
          n = nrow(X)
        )
        genes.de[[i]] <- data.frame(clusters = level.use[i], features = as.character(rownames(data1)), pvalues = pvalues, logFC = FC[features], data.alpha[features,, drop = F],pvalues.adj = pval.adj, stringsAsFactors = FALSE)
      }

      markers.all <- data.frame()
      for (i in 1:numCluster) {
        gde <- genes.de[[i]]
        if (!is.null(gde)) {
          gde <- gde[order(gde$pvalues, -gde$logFC), ]
          gde <- subset(gde, subset = pvalues < thresh.p)
          if (nrow(gde) > 0) {
            markers.all <- rbind(markers.all, gde)
          }
        }
      }
      if (only.pos & nrow(markers.all) > 0) {
        markers.all <- subset(markers.all, subset = logFC > 0)
      }
      if (!is.null(group.dataset)) {
        markers.all$datasets[markers.all$logFC > 0] <- pos.dataset
        markers.all$datasets[markers.all$logFC < 0] <- setdiff(unique(labels.dataset), pos.dataset)
        markers.all$datasets <- factor(markers.all$datasets, levels = levels(factor(object@meta[[group.dataset]])))
        markers.all <- markers.all[order(markers.all$datasets, markers.all$pvalues, -markers.all$logFC), ]
      }
      markers.all$features <- as.character(markers.all$features)
      features.sig <- markers.all$features

    } else if (selection.method == "moransi") {
      cat(cli.symbol(),"Choose spatially variable features using Moran's I value...","\n")
      coord <- object@images$coordinates

      markers.all <- Seurat::RunMoransI(data = data.use, pos = coord, verbose = T)
      # markers.all <- RunMoransI(data = data.use, pos = coord, verbose = FALSE)
      features.sig <- markers.all %>%
        dplyr::filter(p.value<thresh.p) %>%
        dplyr::filter(!is.nan(observed)) %>%
        rownames()
      # markers.all <- Seurat::RunMarkVario(spatial.location = coord, data = X)

    } else if (selection.method == "meringue") {
      cat(cli.symbol(),"Choose spatially variable features using 'MERINGUE' method...","\n")
      # refer to: https://github.com/PYangLab/SVGbench
      callSVG.MERINGUE = function(normMat, spatial_locs) {
        w <- MERINGUE::getSpatialNeighbors(spatial_locs, filterDist = NA)
        I <- MERINGUE::getSpatialPatterns(normMat, w)

        return(I)
      }

      coord <- object@images$coordinates

      markers.all <- callSVG.MERINGUE(normMat = data.use,spatial_locs = coord)
      features.sig <- markers.all %>%
        dplyr::filter(p.adj<thresh.p) %>%
        dplyr::filter(!is.nan(observed)) %>%
        rownames()

    }

    # else if (selection.method == "nnSVG") {
    #   cat(cli.symbol(),"Choose spatially variable features using 'nnSVG' method...","\n")
    #   set.seed(seed.use)
    #   coord <- object@images$coordinates
    #   n.cores <- min(future::availableCores(), n.cores)
    #   svg <-
    #     nnSVG::nnSVG(
    #       data.use,
    #       spatial_coords = as.matrix(coord),
    #       n_threads = n.cores,
    #       BPPARAM = BPPARAM,
    #       verbose = T
    #     )
    #   markers.all <- svg[, c("LR_stat","rank","pval","padj")]
    #   features.sig <- rownames(svg)[svg$pval < thresh.p]
    #
    # } else if (selection.method == "simple") {
    #   cat(cli.symbol(),"Choose spatially variable features using 'simple' method...","\n")
    #   nCell.perFeature <- Matrix::rowSums(data.use>0)# a named vector
    #
    #   markers.all <- as.data.frame(nCell.perFeature)
    #   colnames(markers.all) <- c("nNonZeroCell")
    #
    #   features.sig <- names(nCell.perFeature)[nCell.perFeature>thresh.nCell.perFeature]
    # }

    object@var.features[[features.name]] <- features.sig
    features.name <- paste0(features.name, ".info")
    object@var.features[[features.name]] <- markers.all
    cat(cli.symbol(),"The number of highly variable features is", length(features.sig), '\n')

    if (return.object) {
      return(object)
    } else {
      return(markers.all)
    }

  } # whether to do grid
}


#' Identify over-expressed ligands and (complex) receptors associated with each cell group
#'
#' This function identifies the over-expressed ligands and (complex) receptors based on the identified signaling genes from 'identifyOverExpressedGenes'.
#'
#' @param object CellChat object
#' @param features.name a char name used for storing the over-expressed ligands and receptors in `object@var.features[[paste0(features.name, ".LR")]]`
#' @param features a vector of features to use. default use all over-expressed genes in `object@var.features[[features.name]]`
#' @param return.object whether returning a CellChat object. If FALSE, it will return a data frame containing over-expressed ligands and (complex) receptors associated with each cell group
#' @importFrom future nbrOfWorkers
#' @importFrom future.apply future_sapply
#' @importFrom pbapply pbsapply
#' @importFrom dplyr select
#'
#' @return A CellChat object or a data frame. If returning a CellChat object, a new element named paste0(features.name, ".LR") will be added into the list `object@var.features`
#' @export
#'
identifyOverExpressedLigandReceptor <- function(object, features.name = "features", features = NULL, return.object = TRUE) {

  features.name.LR <- paste0(features.name, ".LR")
  features.name <- paste0(features.name, ".info")
  DB <- object@DB
  interaction_input <- DB$interaction
  complex_input <- DB$complex
  pairLR <- select(interaction_input, ligand, receptor)
  LR.use <- unique(c(pairLR$ligand, pairLR$receptor))
  if (is.null(features)) {
    if (is.list(object@var.features)) {
      markers.all <- object@var.features[[features.name]] # use the updated CellChat object 12/2020
    } else {
      stop("Please update your CellChat object via `updateCellChat()`")
    }

  } else {
    features.use <- features
    rm(features)
    markers.all <- subset(markers.all, subset = features %in% features.use)
  }


  complexSubunits <- complex_input[, grepl("subunit" , colnames(complex_input))]

  markers.all.new <- data.frame()
  for (i in 1:nrow(markers.all)) {
    if (markers.all$features[i] %in% LR.use) {
      markers.all.new <- rbind(markers.all.new, markers.all[i, , drop = FALSE])
    } else {
      index.sig <- unlist(
        x = my_future_lapply(
          X = 1:nrow(complexSubunits),
          FUN = function(x) {
            complexsubunitsV <- unlist(complexSubunits[x,], use.names = F)
            complexsubunitsV <- complexsubunitsV[complexsubunitsV != ""]
            if (markers.all$features[i] %in% complexsubunitsV) {
              return(x)
            }
          }
        )
      )
      complexSubunits.sig <- rownames(complexSubunits[index.sig,])
      markers.all.complex <- data.frame()
      for (j in 1:length(complexSubunits.sig)) {
        markers.all.complex <- rbind(markers.all.complex, markers.all[i, , drop = FALSE])
      }
      markers.all.complex$features <- complexSubunits.sig
      markers.all.new <- rbind(markers.all.new, markers.all.complex)
    }
  }

  object@var.features[[features.name.LR]] <- markers.all.new

  if (return.object) {
    return(object)
  } else {
    return(markers.all.new)
  }
}


#' Identify over-expressed ligand-receptor interactions (pairs) within the used CellChatDB
#'
#' @param object CellChat object
#' @param features.name a char name used for assess the results in `object@var.features[[features.name]]`
#' @param features a vector of features to use. default use all over-expressed genes in `object@var.features[[features.name]]`
#' @param variable.both variable.both = TRUE will require that both ligand and receptor from one pair are over-expressed;
#'
#' variable.both = FALSE will only require that either ligand or receptor from one pair is over-expressed, leading to more over-expressed ligand-receptor interactions (pairs) for further analysis.
#' @param return.object whether returning a CellChat object. If FALSE, it will return a data frame containing the over-expressed ligand-receptor pairs
#' @importFrom future nbrOfWorkers
#' @importFrom future.apply future_sapply
#' @importFrom pbapply pbsapply
#' @importFrom dplyr select
#'
#' @return A CellChat object or a data frame. If returning a CellChat object, a new element named 'LRsig' will be added into the list `object@LR`
#' @export
#'
identifyOverExpressedInteractions <- function(object, features.name = "features", variable.both = TRUE, features = NULL, return.object = TRUE) {
  gene.use <- row.names(object@data.signaling)
  DB <- object@DB
  if (is.null(features)) {
    if (is.list(object@var.features)) {
      features.sig <- object@var.features[[features.name]] # use the updated CellChat object 12/2020
    } else {
      stop("Please update your CellChat object via `updateCellChat()`")
    }

  } else {
    features.sig <- features
  }

  interaction_input <- DB$interaction
  complex_input <- DB$complex
  # my.sapply <- ifelse(
  #   test = future::nbrOfWorkers() == 1,
  #   yes = pbapply::pbsapply,
  #   no = future.apply::future_sapply
  # )
  complexSubunits <- complex_input[, grepl("subunit" , colnames(complex_input))]
  index.sig <- unlist(
    x = my_future_lapply(
      X = 1:nrow(complexSubunits),
      FUN = function(x) {
        complexsubunitsV <- unlist(complexSubunits[x,], use.names = F)
        complexsubunitsV <- complexsubunitsV[complexsubunitsV != ""]
        if (length(intersect(complexsubunitsV, features.sig)) > 0 & all(complexsubunitsV %in% gene.use)) {
          return(x)
        }
      }
    )
  )
  complexSubunits.sig <- complexSubunits[index.sig,]

  index.use <- unlist(
    x = my_future_lapply(
      X = 1:nrow(complexSubunits),
      FUN = function(x) {
        complexsubunitsV <- unlist(complexSubunits[x,], use.names = F)
        complexsubunitsV <- complexsubunitsV[complexsubunitsV != ""]
        if (all(complexsubunitsV %in% gene.use)) {
          return(x)
        }
      },
      strategy.message = FALSE
    )
  )
  complexSubunits.use <- complexSubunits[index.use,]

  pairLR <- select(interaction_input, ligand, receptor)

  if (variable.both) {
    index.sig <- unlist(
      x = my_future_lapply(
        X = 1:nrow(pairLR),
        FUN = function(x) {
          if (all(unlist(pairLR[x,], use.names = F) %in% c(features.sig, rownames(complexSubunits.sig)))) {
            return(x)
          }
        }
      )
    )
  } else {
    index.sig <- unlist(
      x = my_future_lapply(
        X = 1:nrow(pairLR),
        FUN = function(x) {
          # if (all(unlist(pairLR[x,], use.names = F) %in% c(features.sig, rownames(complexSubunits.sig)))) {
          if (all(unlist(pairLR[x,], use.names = F) %in% c(gene.use, rownames(complexSubunits.use))) & (length(intersect(unlist(pairLR[x,], use.names = F), c(features.sig, rownames(complexSubunits.sig)))) > 0)) {
            return(x)
          }
        },
        strategy.message = FALSE
      )
    )
  }

  pairLRsig <- interaction_input[index.sig, ]
  object@LR$LRsig <- pairLRsig
  cat(cli.symbol(),"The number of highly variable ligand-receptor pairs used for signaling inference is", nrow(pairLRsig), '\n')
  if (return.object) {
    return(object)
  } else {
    return(pairLRsig)
  }
}



#' Project gene expression data onto a protein-protein interaction network
#'
#' A diffusion process is used to smooth genes’ expression values based on their neighbors’ defined in a high-confidence experimentally validated protein-protein network.
#'
#' This function is useful when analyzing single-cell data with shallow sequencing depth because the projection reduces the dropout effects of signaling genes, in particular for possible zero expression of subunits of ligands/receptors
#'
#' @param object  CellChat object
#' @param adjMatrix adjacency matrix of protein-protein interaction network to use
#' @param alpha numeric in [0,1] alpha = 0: no smoothing; a larger value alpha results in increasing levels of smoothing.
#' @param normalizeAdjMatrix    how to normalize the adjacency matrix
#'                              possible values are 'rows' (in-degree)
#'                              and 'columns' (out-degree)
#' @return a projected gene expression matrix
#' @export
#'
# This function is adapted from https://github.com/BIMSBbioinfo/netSmooth
projectData <- function(object, adjMatrix, alpha=0.5, normalizeAdjMatrix=c('rows','columns')){
  data <- as.matrix(object@data.signaling)
  normalizeAdjMatrix <- match.arg(normalizeAdjMatrix)
  stopifnot(is(adjMatrix, 'matrix') || is(adjMatrix, 'sparseMatrix'))
  stopifnot((is.numeric(alpha) && (alpha > 0 && alpha < 1)))
  if(sum(Matrix::rowSums(adjMatrix)==0)>0) stop("PPI cannot have zero rows/columns")
  if(sum(Matrix::colSums(adjMatrix)==0)>0) stop("PPI cannot have zero rows/columns")
  if(is.numeric(alpha)) {
    if(alpha<0 | alpha > 1) {
      stop('alpha must be between 0 and 1')
    }
    data.projected <- projectAndRecombine(data, adjMatrix, alpha,normalizeAdjMatrix=normalizeAdjMatrix)
  } else stop("unsupported alpha value: ", class(alpha))
  object@data.project <- data.projected
  return(object)
}

#' Perform network projecting on network when the network genes and the
#' experiment genes aren't exactly the same.
#'
#' The gene network might be defined only on a subset of genes that are
#' measured in any experiment. Further, an experiment might not measure all
#' genes that are present in the network. This function projects the experiment
#' data onto the gene space defined by the network prior to projecting. Then,
#' it projects the projected data back into the original dimansions.
#'
#' @param gene_expression  gene expession data to be projected
#'                         [N_genes x M_samples]
#' @param adj_matrix  adjacenty matrix of network to perform projecting over.
#'                    Will be column-normalized.
#'                    Rownames and colnames should be genes.
#' @param alpha  network projecting parameter (1 - restart probability in random
#'                walk model.
#' @param projecting.function  must be a function that takes in data, adjacency
#'                            matrix, and alpha. Will be used to perform the
#'                            actual projecting.
#' @param normalizeAdjMatrix    which dimension (rows or columns) should the
#'                              adjacency matrix be normalized by. rows
#'                              corresponds to in-degree, columns to
#'                              out-degree.
#' @return  matrix with network-projected gene expression data. Genes that are
#'          not present in projecting network will retain original values.
#' @keywords internal
#'
projectAndRecombine <- function(gene_expression, adj_matrix, alpha,
                                projecting.function=randomWalkBySolve,
                                normalizeAdjMatrix=c('rows','columns')) {
  normalizeAdjMatrix <- match.arg(normalizeAdjMatrix)
  gene_expression_in_A_space <- projectOnNetwork(gene_expression,rownames(adj_matrix))
  gene_expression_in_A_space_project <- projecting.function(gene_expression_in_A_space, adj_matrix, alpha, normalizeAdjMatrix)
  gene_expression_project <- projectFromNetworkRecombine(gene_expression, gene_expression_in_A_space_project)
  return(gene_expression_project)
}


#' Project the gene expression matrix onto a lower space
#' of the genes defined in the projecting network
#' @param gene_expression    gene expression matrix
#' @param new_features       the genes in the network, on which to project
#'                           the gene expression matrix
#' @param missing.value      value to assign to genes that are in network,
#'                           but missing from gene expression matrix
#' @return the gene expression matrix projected onto the gene space defined by new_features
#' @keywords internal
projectOnNetwork <- function(gene_expression, new_features, missing.value=0) {
  # data_in_new_space = matrix(rep(0, length(new_features)*dim(gene_expression)[2]),nrow=length(new_features))
  data_in_new_space = matrix(0, ncol=dim(gene_expression)[2], nrow=length(new_features))
  rownames(data_in_new_space) <- new_features
  colnames(data_in_new_space) <- colnames(gene_expression)
  genes_in_both <- intersect(rownames(data_in_new_space),rownames(gene_expression))
  data_in_new_space[genes_in_both,] <- gene_expression[genes_in_both,]
  genes_only_in_network <- setdiff(new_features, rownames(gene_expression))
  data_in_new_space[genes_only_in_network,] <- missing.value
  return(data_in_new_space)
}

#' project data on graph by solving the linear equation (I - alpha*A) * E_sm = E * (1-alpha)

#' @param E      initial data matrix [NxM]
#' @param A      adjacency matrix of graph to network project on will be column-normalized.
#' @param alpha  projecting coefficient (1 - restart probability of random walk)
#' @return network-projected gene expression
#' @keywords internal
randomWalkBySolve <- function(E, A, alpha, normalizeAjdMatrix=c('rows','columns')) {
  normalizeAjdMatrix <- match.arg(normalizeAjdMatrix)
  if (normalizeAjdMatrix=='rows') {
    Anorm <- l1NormalizeRows(A)
  } else if (normalizeAjdMatrix=='columns') {
    Anorm <- l1NormalizeColumns(A)
  }
  eye <- diag(dim(A)[1])
  AA <- eye - alpha*Anorm
  BB <- (1-alpha) * E
  return(solve(AA, BB))
}

#' Column-normalize a sparse, symmetric matrix (using the l1 norm) so that each
#' column sums to 1.
#'
#' @param A matrix
#' @usage l1NormalizeColumns(A)
#' @return column-normalized sparse matrix object
#' @keywords internal
l1NormalizeColumns <- function(A) {
  return(Matrix::t(Matrix::t(A)/Matrix::colSums(A)))
}

#' Row-normalize a sparse, symmetric matrix (using the l1 norm) so that each
#' row sums to 1.
#'
#' @param A matrix
#' @usage l1NormalizeRows(A)
#' @return row-normalized sparse matrix object
#' @keywords internal
l1NormalizeRows <- function(A) {
  return(A/Matrix::rowSums(A))
}

#' Combine gene expression from projected space (that of the network) with the
#' expression of genes that were not projected (not present in network)
#' @keywords internal
#' @param original_expression    the non-projected expression
#' @param projected_expression    the projected gene expression, in the space
#'                               of the genes defined by the network
#' @return a matrix in the dimensions of original_expression, where values that
#'         are present in projected_expression are copied from there.
projectFromNetworkRecombine <- function(original_expression, projected_expression) {
  data_in_original_space <- original_expression
  genes_in_both <- intersect(rownames(original_expression),rownames(projected_expression))
  data_in_original_space[genes_in_both,] <- as.matrix(projected_expression[genes_in_both,])
  return(data_in_original_space)
}


#' Dimension reduction using PCA
#'
#' @param data.use input data (samples in rows, features in columns)
#' @param do.fast whether do fast PCA
#' @param dimPC the number of components to keep
#' @param seed.use set a seed
#' @param weight.by.var whether use weighted pc.scores
#' @importFrom stats prcomp
#' @importFrom irlba irlba
#' @return
#' @export
#'
#' @examples
runPCA <- function(data.use, do.fast = T, dimPC = 50, seed.use = 42, weight.by.var = T) {
  set.seed(seed = seed.use)
  if (do.fast) {
    dimPC <- min(dimPC, ncol(data.use) - 1)
    pca.res <- irlba::irlba(data.use, nv = dimPC)
    sdev <- pca.res$d/sqrt(max(1, nrow(data.use) - 1))
    if (weight.by.var){
      pc.scores <- pca.res$u %*% diag(pca.res$d)
    } else {
      pc.scores <- pca.res$u
    }
  } else {
    dimPC <- min(dimPC, ncol(data.use) - 1)
    pca.res <- stats::prcomp(x = data.use, rank. = dimPC)
    sdev <- pca.res$sdev
    if (weight.by.var) {
      pc.scores <- pca.res$x %*% diag(pca.res$sdev[1:dimPC]^2)
    } else {
      pc.scores <- pca.res$x
    }
  }
  rownames(pc.scores) <- rownames(data.use)
  colnames(pc.scores) <- paste0('PC', 1:ncol(pc.scores))
  return(pc.scores)
}


#' Run UMAP
#' @param data.use input data matrix
#' @param n_neighbors This determines the number of neighboring points used in
#' local approximations of manifold structure. Larger values will result in more
#' global structure being preserved at the loss of detailed local structure. In general this parameter should often be in the range 5 to 50.
#' @param n_components The dimension of the space to embed into.
#' @param metric This determines the choice of metric used to measure distance in the input space.
#' @param n_epochs the number of training epochs to be used in optimizing the low dimensional embedding. Larger values result in more accurate embeddings. If NULL is specified, a value will be selected based on the size of the input dataset (200 for large datasets, 500 for small).
#' @param learning_rate The initial learning rate for the embedding optimization.
#' @param min_dist This controls how tightly the embedding is allowed compress points together.
#' Larger values ensure embedded points are moreevenly distributed, while smaller values allow the
#' algorithm to optimise more accurately with regard to local structure. Sensible values are in the range 0.001 to 0.5.
#' @param spread he effective scale of embedded points. In combination with min.dist this determines how clustered/clumped the embedded points are.
#' @param set_op_mix_ratio Interpolate between (fuzzy) union and intersection as the set operation used to combine local fuzzy simplicial sets to obtain a global fuzzy simplicial sets.
#' @param local_connectivity The local connectivity required - i.e. the number of nearest neighbors
#' that should be assumed to be connected at a local level. The higher this value the more connected
#' the manifold becomes locally. In practice this should be not more than the local intrinsic dimension of the manifold.
#' @param repulsion_strength Weighting applied to negative samples in low dimensional embedding
#' optimization. Values higher than one will result in greater weight being given to negative samples.
#' @param negative_sample_rate The number of negative samples to select per positive sample in the
#' optimization process. Increasing this value will result in greater repulsive force being applied, greater optimization cost, but slightly more accuracy.
#' @param a More specific parameters controlling the embedding. If NULL, these values are set automatically as determined by min. dist and spread.
#' @param b More specific parameters controlling the embedding. If NULL, these values are set automatically as determined by min. dist and spread.
#' @param seed.use Set a random seed. By default, sets the seed to 42.
#' @param metric_kwds,angular_rp_forest,verbose other parameters used in UMAP
#' @import reticulate
#' @export
#'
runUMAP <- function(
  data.use,
  n_neighbors = 30L,
  n_components = 2L,
  metric = "correlation",
  n_epochs = NULL,
  learning_rate = 1.0,
  min_dist = 0.3,
  spread = 1.0,
  set_op_mix_ratio = 1.0,
  local_connectivity = 1L,
  repulsion_strength = 1,
  negative_sample_rate = 5,
  a = NULL,
  b = NULL,
  seed.use = 42L,
  metric_kwds = NULL,
  angular_rp_forest = FALSE,
  verbose = FALSE){
  if (!reticulate::py_module_available(module = 'umap')) {
    stop("Cannot find UMAP, please install through pip (e.g. pip install umap-learn or reticulate::py_install(packages = 'umap-learn')).")
  }
  set.seed(seed.use)
  reticulate::py_set_seed(seed.use)
  umap_import <- reticulate::import(module = "umap", delay_load = TRUE)
  umap <- umap_import$UMAP(
    n_neighbors = as.integer(n_neighbors),
    n_components = as.integer(n_components),
    metric = metric,
    n_epochs = n_epochs,
    learning_rate = learning_rate,
    min_dist = min_dist,
    spread = spread,
    set_op_mix_ratio = set_op_mix_ratio,
    local_connectivity = local_connectivity,
    repulsion_strength = repulsion_strength,
    negative_sample_rate = negative_sample_rate,
    a = a,
    b = b,
    metric_kwds = metric_kwds,
    angular_rp_forest = angular_rp_forest,
    verbose = verbose
  )
  Rumap <- umap$fit_transform
  umap_output <- Rumap(t(data.use))
  colnames(umap_output) <- paste0('UMAP', 1:ncol(umap_output))
  rownames(umap_output) <- colnames(data.use)
  return(umap_output)
}

.error_if_no_Seurat <- function() {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Seurat installation required for working with Seurat objects")
  }
}


#' Color interpolation
#'
#' This function is modified from https://rdrr.io/cran/circlize/src/R/utils.R
#' Colors are linearly interpolated according to break values and corresponding colors through CIE Lab color space (`colorspace::LAB`) by default.
#' Values exceeding breaks will be assigned with corresponding maximum or minimum colors.
#'
#' @param breaks A vector indicating numeric breaks
#' @param colors A vector of colors which correspond to values in ``breaks``
#' @param transparency A single value in ``[0, 1]``. 0 refers to no transparency and 1 refers to full transparency
#' @param space color space in which colors are interpolated. Value should be one of "RGB", "HSV", "HLS", "LAB", "XYZ", "sRGB", "LUV", see `colorspace::color-class` for detail.
#' @importFrom colorspace coords RGB HSV HLS LAB XYZ sRGB LUV hex
#' @importFrom grDevices col2rgb
#' @return It returns a function which accepts a vector of numeric values and returns interpolated colors.
#' @export
#' @examples
#' \dontrun{
#' col_fun = colorRamp3(c(-1, 0, 1), c("green", "white", "red"))
#' col_fun(c(-2, -1, -0.5, 0, 0.5, 1, 2))
#' }
colorRamp3 = function(breaks, colors, transparency = 0, space = "RGB") {

  if(length(breaks) != length(colors)) {
    stop("Length of `breaks` should be equal to `colors`.\n")
  }

  colors = colors[order(breaks)]
  breaks = sort(breaks)

  l = duplicated(breaks)
  breaks = breaks[!l]
  colors = colors[!l]

  if(length(breaks) == 1) {
    stop("You should have at least two distinct break values.")
  }


  if(! space %in% c("RGB", "HSV", "HLS", "LAB", "XYZ", "sRGB", "LUV")) {
    stop("`space` should be in 'RGB', 'HSV', 'HLS', 'LAB', 'XYZ', 'sRGB', 'LUV'")
  }

  colors = t(grDevices::col2rgb(colors)/255)

  attr = list(breaks = breaks, colors = colors, transparency = transparency, space = space)

  if(space == "LUV") {
    i = which(apply(colors, 1, function(x) all(x == 0)))
    colors[i, ] = 1e-5
  }

  transparency = 1-ifelse(transparency > 1, 1, ifelse(transparency < 0, 0, transparency))[1]
  transparency_str = sprintf("%X", round(transparency*255))
  if(nchar(transparency_str) == 1) transparency_str = paste0("0", transparency_str)

  fun = function(x = NULL, return_rgb = FALSE, max_value = 1) {
    if(is.null(x)) {
      stop("Please specify `x`\n")
    }

    att = attributes(x)
    if(is.data.frame(x)) x = as.matrix(x)

    l_na = is.na(x)
    if(all(l_na)) {
      return(rep(NA, length(l_na)))
    }

    x2 = x[!l_na]

    x2 = ifelse(x2 < breaks[1], breaks[1],
                ifelse(x2 > breaks[length(breaks)], breaks[length(breaks)],
                       x2
                ))
    ibin = .bincode(x2, breaks, right = TRUE, include.lowest = TRUE)
    res_col = character(length(x2))
    for(i in unique(ibin)) {
      l = ibin == i
      res_col[l] = .get_color(x2[l], breaks[i], breaks[i+1], colors[i, ], colors[i+1, ], space = space)
    }
    res_col = paste(res_col, transparency_str[1], sep = "")

    if(return_rgb) {
      res_col = t(grDevices::col2rgb(as.vector(res_col), alpha = TRUE)/255)
      return(res_col)
    } else {
      res_col2 = character(length(x))
      res_col2[l_na] = NA
      res_col2[!l_na] = res_col

      attributes(res_col2) = att
      return(res_col2)
    }
  }

  attributes(fun) = attr
  return(fun)
}

.restrict_in = function(x, lower, upper) {
  x[x > upper] = upper
  x[x < lower] = lower
  x
}

# x: vector
# break1 single value
# break2 single value
# rgb1 vector with 3 elements
# rgb2 vector with 3 elements
.get_color = function(x, break1, break2, col1, col2, space) {

  col1 = colorspace::coords(as(colorspace::sRGB(col1[1], col1[2], col1[3]), space))
  col2 = colorspace::coords(as(colorspace::sRGB(col2[1], col2[2], col2[3]), space))

  res_col = matrix(ncol = 3, nrow = length(x))
  for(j in 1:3) {
    xx = (x - break2)*(col2[j] - col1[j]) / (break2 - break1) + col2[j]
    res_col[, j] = xx
  }

  res_col = get(space)(res_col)
  res_col = colorspace::coords(as(res_col, "sRGB"))
  res_col[, 1] = .restrict_in(res_col[,1], 0, 1)
  res_col[, 2] = .restrict_in(res_col[,2], 0, 1)
  res_col[, 3] = .restrict_in(res_col[,3], 0, 1)
  colorspace::hex(colorspace::sRGB(res_col))
}
