#' computeCommunProb
#'
#' @description
#' Compute the communication probability/strength between any interacting individual cells
#'
#' @param object CellChat object
#' @param LR.use A subset of ligand-receptor interactions used in inferring communication network
#' @param raw.use Whether use the raw data (i.e., `object@data.signaling`) or the projected data (i.e., `object@data.project`).
#' Set raw.use = FALSE to use the projected data when analyzing single-cell data with shallow sequencing depth because the projected data could help to reduce the dropout effects of signaling genes, in particular for possible zero expression of subunits of ligands/receptors.
#' @param Kh Parameter in Hill function
#' @param n Parameter in Hill function
#' @param distance.use Whether to use distance constraints to compute communication probability. Setting `distance.use = TRUE` indicates that the cell-cell communication probability is inversely proportional to the computed distance.
#' @param tol NULL or Numeric. set a distance tolerance when computing cell-cell distances and contact adjacent matrix. By default `tol = NULL` means distance tolerance equals to spot.size/2 or cell.diameter/2.
#' @param interaction.range The maximum interaction/diffusion length of ligands (Unit: microns). This hard threshold is used to filter out the connections between spatially distant individual cells
#' @param scale.distance A scale or normalization factor for the spatial distances when setting `distance.use = TRUE`. For example, scale.distance equals 1, 0.1, 0.01, 0.001, 0.11, or 0.011. We choose this values such that the minimum value of the scaled distances is in [1,2]. This value is not necessary when setting `distance.use = FALSE`.
#' @param use.AGAN Boolean. Whether to take agonist (AG) and antagonist (AN) into consideration when calculating the intercellular communication probability
#' @param contact.range Numeric. The interaction range (Unit: microns) to restrict the contact-dependent signaling.
#' For spatial transcriptomics in a single-cell resolution, `contact.range` is approximately equal to the estimated cell diameter (i.e., the cell center-to-center distance), which means that contact-dependent and juxtacrine signaling can only happens when the two cells are contact to each other.
#' Typically, `contact.range = 10`, which is a typical human cell size. However, for low-resolution spatial data such as 10X visium, it should be the cell center-to-center distance (i.e., `contact.range = 100` for visium data).  The function \link{computeCellDistance} can compute the center-to-center distance.
#' @param contact.dependent Boolean. Whether determining spatially proximal cell groups based on the `contact.range`. By default `contact.dependent = TRUE` when inferring contact-dependent and juxtacrine signaling (including ECM-Receptor and Cell-Cell Contact signaling classified in CellChatDB$interaction$annotation).
#' If only focusing on `Secreted Signaling`, the `contact.dependent` will be automatically set as FALSE except for `contact.dependent.forced = TRUE`.
#' @param contact.dependent.forced Boolean. Whether forcing to determine spatially proximal cell regions based on the `contact.range` for all cells/spots in the CellChat Object.
#' Users can set `contact.dependent.forced = TRUE` to turn `Secreted Signaling` interactions into a contact manner.
#'
#' @return CellChat object
#' @export
#'
#' @examples
computeCommunProb <- function (
    object,
    LR.use = NULL,
    raw.use = TRUE,
    Kh = 0.5,
    n = 1,
    distance.use = TRUE,
    tol = NULL, # will be removed in the future
    interaction.range = 250,
    scale.distance = 0.01,
    use.AGAN = T,
    contact.dependent = TRUE,
    contact.range = 10,
    contact.dependent.forced = FALSE
) {

 if (raw.use) {
    data <- object@data.signaling
    # scale the elements
    data@x <- data@x/max(data@x)
    data.use <- as.matrix(data)

  } else {
    data <- object@data.project
    # scale
    data.use <- data/max(data)
  }

  if (is.null(LR.use)) {
    pairLR.use <- object@LR$LRsig
  } else {
    if (length(unique(LR.use$annotation)) > 1) {
      LR.use$annotation <-
        factor(
          LR.use$annotation,
          levels = c(
            "Secreted Signaling",
            "ECM-Receptor",
            "Non-protein Signaling",
            "Cell-Cell Contact"
          )
        )
      LR.use <- LR.use[order(LR.use$annotation), , drop = FALSE]
      LR.use$annotation <- as.character(LR.use$annotation)
    }
    pairLR.use <- LR.use
  }

  complex_input <- object@DB$complex
  cofactor_input <- object@DB$cofactor

  ptm = Sys.time()

  pairLRsig <- pairLR.use
  group <- object@idents
  geneL <- as.character(pairLRsig$ligand)
  geneR <- as.character(pairLRsig$receptor)
  nLR <- nrow(pairLRsig)
  numCluster <- nlevels(group)
  if (numCluster != length(unique(group))) {
    stop(cli.symbol(2),"Please check `unique(object@idents)` and ensure that the factor levels are correct!\n         You may need to drop unused levels using 'droplevels' function. e.g.,\n         `meta$labels = droplevels(meta$labels, exclude = setdiff(levels(meta$labels),unique(meta$labels)))`")
  }

  nC <- ncol(data.use)

  # working on spatial transcriptomic data and preferring to infer interactions between individual cell
  cat(cli.symbol(),"Analyzing spatial transcriptomic data and preferring to infer interactions between individual cells...\n")


  data.spatial <- BiocGenerics::as.data.frame(object@images$coordinates)

  ### previous ###
  # spot.size.fullres <- object@images$scale.factors$spot
  # spot.size <- object@images$scale.factors$spot.diameter

  ratio <- object@images$spatial.factors[["ratio"]]
  if(is.null(tol)) tol <- object@images$spatial.factors[["tol"]] else NULL

  # compute the cell-to-cell distances
  # res is a list object, containing `d.spatial` matrix and `adj.contact` matrix!
  res <- computeCellDistance(coordinates = data.spatial,
                             ratio = ratio,
                             interaction.range = interaction.range,
                             contact.range = contact.range,
                             tol = tol)
  # long-range distance
  d.spatial <- res$d.spatial
  # short-range distance adjacent matrix for contact-dependent and juxtacrine signaling
  adj.contact <- res$adj.contact
  gc()

  if (distance.use) {
    cat(paste0(cli.symbol(),"Run CellChat on spatial transcriptomic data using distances as constraints <<< [",
               Sys.time(), "]\n"))
    d.spatial@x <- d.spatial@x * scale.distance
    d.min <- min(d.spatial@x, na.rm = F) # d.spatial@x has no NA
    if (d.min < 1) {
      cat(cli.symbol(),"The suggested minimum value of scaled distances is in [1,2], and the calculated value here is ", d.min,"\n")
      stop(cli.symbol(2),"Please increase the value of `scale.distance` and use a value that is slighly smaller than ", format(1/d.min, digits = 2) ,"\n")
    }
    P.spatial <- createPspatialFrom_dspatial(d.spatial,distance.use = T)
    d.spatial@x <- d.spatial@x / scale.distance
  }  else {
    cat(paste0(cli.symbol(),"Run CellChat on transcriptomic imaging data without distances as constraints <<< [",
               Sys.time(), "]\n"))
    P.spatial <- createPspatialFrom_dspatial(d.spatial,distance.use = F)

  }
  rm(d.spatial);gc()

  # set a flag for contact.dependent signaling
  all.contact.dependent <- FALSE
  all.diffusible <- FALSE
  if (contact.dependent.forced == TRUE) {
    # cat(cli.symbol(),"Run with `contact.dependent.forced = T` \n")
    cat(cli.symbol(),"Force to run CellChat in a `contact-dependent` manner for all L-R pairs including secreted signaling.\n")
    P.spatial <- P.spatial * adj.contact
    nLR1 <- nLR
    all.contact.dependent <- TRUE
  } else{ # contact.dependent.forced == F
    # cat(cli.symbol(),"Run with `contact.dependent.forced = F` \n")
    if (contact.dependent == TRUE && length(unique(pairLRsig$annotation))>0 ) {
      if (all(unique(pairLRsig$annotation) == c("Cell-Cell Contact") )) {
        # all interactions in `pairLRsig` are contact-dependent signaling
        cat(cli.symbol(),"All the input L-R pairs are `Cell-Cell Contact` signaling. Run CellChat in a contact-dependent manner. \n")
        P.spatial <- P.spatial * adj.contact
        nLR1 <- nLR
        all.contact.dependent <- TRUE
      } else if (all(unique(pairLRsig$annotation) %in% c("Secreted Signaling", "ECM-Receptor", "Non-protein Signaling"))) {
        # all interactions in `pairLRsig` are not contact-dependent signaling
        cat(cli.symbol(),"Molecules of all the input L-R pairs are diffusible. Run CellChat in a diffusion manner based on the `interaction.range`.\n")
        nLR1 <- nLR
        all.diffusible <- TRUE
      } else {
        # Interactions in `pairLRsig` have both contact-dependent signaling and secreted signaling
        cat(cli.symbol(),"The input L-R pairs have both secreted signaling and contact-dependent signaling. Run CellChat in a contact-dependent manner for `Cell-Cell Contact` signaling, and in a diffusion manner based on the `interaction.range` for other L-R pairs. \n")
        nLR1 <- max(which(pairLRsig$annotation %in% c("Secreted Signaling", "ECM-Receptor", "Non-protein Signaling")))
      }
    } else { # contact.dependent == F or `object@LR$LRsig` does not have `annotation` column, take all interactions as `Secreted Signaling` interactions
      cat(cli.symbol(),"Run CellChat in a diffusion manner based on the `interaction.range` for all L-R pairs. \n")
      cat(cli.symbol(3),"Set`contact.dependent = TRUE` if preferring a contact-dependent manner for `Cell-Cell Contact` signaling. \n")
      nLR1 <- nLR
    }
  }


  # compute the expression of ligand or receptor
  dataLavg <- computeExpr_LR(geneL, data.use, complex_input)
  dataRavg <- computeExpr_LR(geneR, data.use, complex_input)

  # take account into the effect of co-activation and co-inhibition receptors
  dataRavg.co.A.receptor <- computeExpr_coreceptor(cofactor_input,
                                                   data.use, pairLRsig, type = "A")
  dataRavg.co.I.receptor <- computeExpr_coreceptor(cofactor_input,
                                                   data.use, pairLRsig, type = "I")
  dataRavg <- dataRavg * dataRavg.co.A.receptor/dataRavg.co.I.receptor
  rm(dataRavg.co.A.receptor,dataRavg.co.I.receptor);gc();


  # compute the expression of agonist and antagonist
  # index.agonist <- which(!is.na(pairLRsig$agonist) & pairLRsig$agonist != "")
  # index.antagonist <- which(!is.na(pairLRsig$antagonist) & pairLRsig$antagonist != "")

  # Compute the communication probability/strength between any interacting individual cells for each LR pair
  Prob.cell_ <- my_future_lapply(
    X = 1:nLR,
    FUN = function(i) {

      x_ <- as(matrix(dataLavg[i,],nrow = 1),"TsparseMatrix")
      y_ <- as(matrix(dataRavg[i,],nrow = 1),"TsparseMatrix")
      # dataLR <- Matrix::crossprod(matrix(dataLavg[i,],nrow = 1),
      #                             matrix(dataRavg[i,],nrow = 1))
      # dataLR <- as(dataLR, Class = "dgCMatrix");gc();
      dataLR <- Matrix::crossprod(x_,y_)


      # dataLR is a sparse Matrix => P1 will be a sparse Matrix too
      P1 <- HillFunctionFordataLR(dataLR = dataLR, Kh = Kh, n = n)

      # P1_Pspatial is a dgCMatrix
      P1_Pspatial <- P1 * P.spatial;rm(P1);gc()

      if (!use.AGAN) { # use.AGAN = F

        # ####################################
        # Codes below are necessary because we still need to take some contact-dependent signalings
        # into consideration when `use.AGAN=F`
        if (i > nLR1) {
          P1_Pspatial <- P1_Pspatial * adj.contact
        }
        # ####################################

        # Pnull.cell is a sparse matrix
        Pnull.cell = P1_Pspatial
        dimnames(Pnull.cell) <- list(NULL,NULL)
        return(Pnull.cell)

      } else { # use.AGAN = T

        if (i > nLR1) {
          P1_Pspatial <- P1_Pspatial * adj.contact
        }

        # if P1_Pspatial is all-zero matrix, iteration can be terminated in advance
        if (sum(P1_Pspatial) == 0){
          # For test: sum(P1_Pspatial) == 0
          Pnull.cell = P1_Pspatial
          dimnames(Pnull.cell) <- list(NULL,NULL)
          return(Pnull.cell)
        } else {

          # For test: is.element(i, index.agonist)
          # data.agonist => P2
          data.agonist <- computeExpr_agonist(
            data.use = data.use,
            pairLRsig,
            cofactor_input,
            index.agonist = i,
            Kh = Kh,
            n = n
          )

          # P_ = P1*Pspatial*P2
          P_ <- myElementwiseProduct(P1_Pspatial,data.agonist)


          # data.antagonist => P3
          data.antagonist <- computeExpr_antagonist(
            data.use = data.use,
            pairLRsig,
            cofactor_input,
            index.antagonist = i,
            Kh = Kh,
            n = n
          )

          # Note: Pnull.cell = P1*Pspatial*P2*P3
          Pnull.cell <- myElementwiseProduct(P_,data.antagonist)
          rm(P_)
          # Pnull.cell is "dgCMatrix"
          dimnames(Pnull.cell) <- list(NULL,NULL)
          return(Pnull.cell)
        }
      }
    },
    simplify = F # return a list
  )

  # bind the Prob.cell `list` => a `sparse3Darray`
  # then the Prob.cell's shape will be (nC,nC,nLR)
  Prob.cell <- my_as_sparse3Darray(Prob.cell_)
  cat(cli.symbol(),"The number of cells and L-R pairs in Dim(Prob.cell):",dim(Prob.cell),"\n")

  # set `Prob.cell`'s names
  dimnames(Prob.cell) <- list(colnames(data.use), colnames(data.use),
                              rownames(pairLRsig))
  names(Prob.cell_) <- rownames(pairLRsig)

  Tmp <- list(prob.cell = Prob.cell_,Lavg=dataLavg,Ravg=dataRavg) # !important, for parallel iteration
  net <- list(prob.cell = Prob.cell, tmp = Tmp)
  execution.time = Sys.time() - ptm
  object@options$run.time <- as.numeric(execution.time,
                                        units = "secs")
  object@images[["result.computeCellDistance"]] <- res
  object@options$parameter <- list(
    raw.use = raw.use,
    # spot.size = spot.size,
    # spot.size.fullres = spot.size.fullres,
    ratio = ratio,
    tol = tol,
    Kh = Kh,
    n = n,
    nLR = nLR,
    nLR1 = nLR1,
    scale.distance = scale.distance,
    use.AGAN = use.AGAN,
    distance.use = distance.use,
    interaction.range = interaction.range,
    contact.dependent = contact.dependent,
    contact.range = contact.range,
    contact.dependent.forced = contact.dependent.forced,
    all.contact.dependent = all.contact.dependent,
    all.diffusible = all.diffusible
  )

  object@net <- net
  cat(paste0(cli.symbol(symbol = "success")," CellChat inference is done. Parameter values are stored in `object@options$parameter` <<< [",
             Sys.time(), "]", "\n"))

  return(object)
}


#' create P.spatial from d.spatial object
#'
#' @description
#' A function defined for \link{computeCommunProb}. To model ligands' diffusion attenuation, we take the reciprocal of cell-cell distance
#' as a factor (named as P.spatial) determining the cell-cell communication probability. Here we design
#' the function to help take reciprocal of a `d.spatial` matrix.
#'
#' @param d.spatial CsparseMatrix (dsC Matrix) object. A cell-cell distance matrix which has been scaled, see in \link{computeCommunProb}
#' @param distance.use Boolean. Whether to use distance constraints to compute communication probability, passed from \link{computeCommunProb}'s same parameter
#' Setting `distance.use = TRUE` indicates that the cell-cell communication probability is inversely proportional to the computed distance.
#' Otherwise, the cell-cell communication probability share the same computation factor 1 within the maximum interaction/diffusion length of ligands
#' @return CsparseMatrix object in `Matrix` Package
#' @export
#'
#' @examples
createPspatialFrom_dspatial <- function(d.spatial,distance.use=T){
  # d.spatial is a sparse matrix, note if dividing it directly, many 1/0 => inf will arise
  P.spatial <- d.spatial
  P.spatial@x <- 1/P.spatial@x
  # Fill up the diagonal: autocrine has the largest spatial CCC probability
  Matrix::diag(P.spatial) <- max(P.spatial@x)

  # make P.spatial binary
  if (distance.use) {
    return(P.spatial)
  } else {
    P.spatial@x <- rep.int(1,times = length(P.spatial@x))
    return(P.spatial)
  }
}



#' HillFunction for dataLR object
#' @description
#' A function defined for \link{computeCommunProb}. A Hill function was used to model the interactions
#' between L and R with an EC50 parameter Kh whose default value was set to be 0.5 as the input data has a normalized range from 0 to 1.
#' Previous hill function is `dataLR^n/(Kh^n + dataLR^n)`, notice dataLR is a sparse matrix, so we design the function to do element-wise hill function upon each non-zero value in the `dataLR` matrix
#'
#' @param dataLR A dgCMatrix object in `Matrix` Package.
#' @param Kh A hyper parameter in hill function, passed from \link{computeCommunProb}
#' @param n A hyper parameter in hill function, passed from \link{computeCommunProb}
#'
#' @return a normalized dataLR
#' @export
#'
#' @examples
HillFunctionFordataLR <- function(dataLR,Kh,n){
  # dataLR is a dgCMatrix, we map the Hill function on each non-zero in the Matrix
  x_dataLR_ <- purrr::map_dbl(dataLR@x,function(x){x^n/(Kh^n + x^n)})
  # update the values
  dataLR@x <- x_dataLR_
  return(dataLR)
}

#' element-wise product for a sparse matrix and a dense matrix
#'
#' @description
#' A function defined for \link{computeCommunProb}. Define an element-wise product operation for a sparse matrix and a dense matrix,
#' to speed up the product operation
#'
#' @param SparseMat A dgCMatrix object in `Matrix` Package.
#' @param DenseMat A dense matrix shares the same dimensions with `SparseMat`
#'
#' @return element-wise product result. A dgCMatrix object in `Matrix` Package
#' @export
#'
#' @examples
elementwiseProductForSparseMatandDenseMat <- function(SparseMat,DenseMat){
  # SparseMat must be a "dgCMatrix"
  SparseMat <- as(SparseMat,Class = "TsparseMatrix")
  nonZeroNum <- length(SparseMat@i)

  xFromDenseMat <- vector('double',length = nonZeroNum)
  for (k in 1:nonZeroNum) {
    # row indices of non-zero entries in 0-base, so add 1
    # col indices of non-zero entries in 0-base, so add 1
    xFromDenseMat[[k]] <- DenseMat[(SparseMat@i[[k]]+1), (SparseMat@j[[k]]+1)]
  }
  # update the x in SparseMat
  SparseMat@x <- SparseMat@x * xFromDenseMat
  # still return a "dgCMatrix"
  SparseMat <- as(SparseMat,Class = "CsparseMatrix")
  return(SparseMat)
}

#' element-wise product for a sparse matrix and a dense vector's crossprod
#'
#' @description
#' A function defined for \link{computeCommunProb}. Define an element-wise product operation for a sparse matrix and a dense vector's crossprod,
#' to speed up the product operation
#'
#' @param SparseMat A dgCMatrix object in `Matrix` Package.
#' @param DenseVec A numeric vector shares the same dimensions with `SparseMat`
#'
#' @return element-wise product result. A dgCMatrix object in `Matrix` Package
#' @export
myElementwiseProduct <- function(SparseMat,DenseVec){
  # SparseMat must be a "dgCMatrix"
  SparseMat <- as(SparseMat,Class = "TsparseMatrix")
  nonZeroNum <- length(SparseMat@i)

  xFromCrossProduct <- vector('double',length = nonZeroNum)
  for (k in 1:nonZeroNum) {
    # row indices of non-zero entries in 0-base, so add 1
    # col indices of non-zero entries in 0-base, so add 1
    # DenseVec[[(SparseMat@i[[k]]+1)]] * DenseVec[[(SparseMat@j[[k]]+1)]] = crossprod[i,j]
    xFromCrossProduct[[k]] <- DenseVec[[(SparseMat@i[[k]]+1)]] * DenseVec[[(SparseMat@j[[k]]+1)]]
  }
  # update the x in SparseMat
  SparseMat@x <- SparseMat@x * xFromCrossProduct
  # still return a "dgCMatrix"
  SparseMat <- as(SparseMat,Class = "CsparseMatrix")
  return(SparseMat)
}

#' @title relabelSpatialCellChat
#' @description
#' Merge some cell groups into a new cell group, and recalculate cell-cell communication
#' probability/strength at the cell-group level.
#'
#' @param object SpatialCellChat object.
#' @param labelSet Named character vector. Make sure the names are cell group labels in
#' present `object@idents`, and values are new cell group labels respectively.
#' For example,`labelSet = c("Group1"="NewGroup1","Group2"="NewGroup1","Group3"="NewGroup1")`.
#' @param group.by Set a column name in `object@meta` if not use `object@idents` as cell group labels
#' to relabel.
#' @param nboot Integer. See the same parameter in \link{computeAvgCommunProb} or \link{computeAvgCommunProb_Visium}.
#' @param min.cells.sr Integer. See the same parameter in \link{computeAvgCommunProb}.
#' @param min.percent Numeric. See the same parameter in \link{computeAvgCommunProb}.
#' @param do.permutation Boolean. Integer. See the same parameter in \link{computeAvgCommunProb} or \link{computeAvgCommunProb_Visium}.
#' @param min.cells Integer. See the same parameter in \link{filterCommunication}.
#' @param thresh Numeric. See the same parameter in \link{computeCommunProbPathway}.
#'
#' @return SpatialCellChat object.
#' @export
relabelSpatialCellChat <- function(
    object,
    labelSet,
    group.by=NULL,
    nboot = 100,
    min.cells.sr = 10,
    min.percent = 0.1,
    do.permutation = T,
    min.cells = 10,
    thresh = 0.05
){
  # object <- Chat
  #
  # labelSet <- c(
  #   "Basal" = "Keratinocyte",
  #   "Spinous" = "Keratinocyte",
  #   "Supraspinous" = "Keratinocyte",
  #   "TC" = "Immune",
  #   "MYL" = "Immune"
  # )
  if(is.null(group.by)){
    group.lab <- as.character(object@idents)
    group.level <- levels(object@idents)
  } else {
    if (!(group.by %in% colnames(object@meta))) {
      stop("The 'group.by' is not a column name in the `meta` slot, which will be used for cell grouping.")
    } else {
      group.lab <- as.character(object@meta[[group.by]])
      # group.lab <- factor(group.lab,levels = unique(group.lab))
      # object@meta[[group.by]] <- group.lab
      group.level <- unique(group.lab)
    }
  }

  if(is.null(names(labelSet))){
    stop("Please provide the names of cell groups which will be dropped
    by setting the `names` of labelSet.\n")
  }
  group.drop <- names(labelSet)
  if(!all(group.drop%in%group.level)) {
    stop("Please check the `labelSet`, make sure the cell groups
    to be dropped are in the `group.by` column.\n")
  }

  tmp.df <- data.frame(
    pre = group.lab,
    new = labelSet[group.lab]
  )
  tmp.df <- tmp.df %>% mutate(
    new = if_else(is.na(new),pre,new)
  )

  object@meta[["new.ident"]] <- tmp.df[["new"]]
  # spatialDimPlot(object,group.by = "new.ident")
  # spatialDimPlot(object,group.by = "ident")

  if(!is.null(object@net[["tmp"]]$cell.type.decomposition)){
    cell.type.decomposition <- as.matrix(object@net[["tmp"]]$cell.type.decomposition)
    cell.type.decomposition.left <- cell.type.decomposition[,base::setdiff(group.level,group.drop),drop=F]

    cell.type.decomposition.new <- lapply(
      X = unique(labelSet),
      FUN = function(i){
        group.drop.i <- names(labelSet[labelSet==i])
        cell.type.decomposition.i <-
          Matrix::rowSums(cell.type.decomposition[,group.drop.i,drop=F])
        return(cell.type.decomposition.i)
      }
    )
    names(cell.type.decomposition.new) <- unique(labelSet)
    cell.type.decomposition.new <- as.data.frame(cell.type.decomposition.new)
    cell.type.decomposition.new <- cbind(cell.type.decomposition.left,as.matrix(cell.type.decomposition.new))
    object@meta[["new.ident"]] <- factor(object@meta[["new.ident"]],levels = colnames(cell.type.decomposition.new))
    object@idents <- object@meta[["new.ident"]]


    object <- computeAvgCommunProb_Visium(
      object,
      cell.type.decomposition = cell.type.decomposition.new,
      nboot = nboot,
      # min.cells.sr = min.cells.sr,
      # min.percent = min.percent,
      do.permutation = do.permutation
    )
  } else {
    object@meta[["new.ident"]] <- factor(object@meta[["new.ident"]],levels = unique(object@meta[["new.ident"]]))
    object@idents <- object@meta[["new.ident"]]
    object <- computeAvgCommunProb(
      object,
      nboot = nboot,
      min.cells.sr = min.cells.sr,
      min.percent = min.percent,
      do.permutation = do.permutation
    )
  }


  object <- filterCommunication(
    object,
    min.cells = min.cells,
    min.links = NULL,
    min.cells.sr = NULL
  )

  object <- computeCommunProbPathway(
    object,
    thresh = thresh
  )

  object <- aggregateNet(object)

  return(object)
}



#' Compute group-level cell-cell communication for Visium Low Definition(LD)
#'
#' @description
#' For Visium HD ST datasets, please use \link{computeAvgCommunProb} instead.
#'
#' @param object SpatialCellChat object with communication probabilities for pairwise individual cells
#' @param cell.type.decomposition NCells x NTypes Matrix, which comes from any kind of cell type decomposition.
#' @param group.by cell group information used for computing average communication probabilities
# @param min.percent minimum percentage of expressed ligands or receptors per cell group to require for computing the group-level signaling
# @param min.cells.sr minimum number of cells required as senders or receivers per cell group for computing the group-level signaling
#' @param avg.type methods for integrating communication probabilities per cell group
#' @param do.permutation whether performing permutation test
#' @param nboot the number of permutations
#' @param seed.use set a random seed. By default, set the seed to 1.
#' @param colocalization.use whether filtering out spatially distant cell groups based on colocalization analysis between any cell groups
#' @param thresh.colo removal of cell-cell communication with no significant colocalizations (fdr < 0.05)
# #' @inheritParams computeAvgCommunProb_LR_Visium_Avg
# #' @inheritParams computeAvgCommunProb_LR_Visium_Sum
#'
#' @return A CellChat object with updated slot 'net':
#'
#' object@net$prob is the inferred group-level communication probability (strength) array, where the first, second and third dimensions represent a source group, target group and ligand-receptor pair, respectively.
#'
#' USER can access all the inferred cell-cell communications using the function 'subsetCommunication(object)', which returns a data frame.
#'
#' object@net$pval is the corresponding p-values of each interaction
#'
#' @return
#' @export
computeAvgCommunProb_Visium <- function(
    object,
    cell.type.decomposition,
    avg.type = c("avg","sum"),
    group.by = NULL,
    do.permutation = T,
    nboot = 100,
    seed.use = 1L,
    colocalization.use = F,
    thresh.colo = 0.05
){

  if (is.null(group.by)) {
    group <- object@idents
  } else {
    if (!(group.by %in% colnames(object@meta))) {
      stop("The 'group.by' is not a column name in the `object@meta`, which will be used for cell grouping.")
    } else {
      group <- object@meta[[group.by]]
    }
    if (!is.factor(group)) {
      group <- factor(group)
    }
  }
  cat(cli.symbol(),"The cell groups used for averaging cell-cell communication are ", cli::col_red(levels(group)), '\n')
  if (!inherits(x = cell.type.decomposition, what = c("matrix", "Matrix"))) {
    stop(cli.symbol("fail"),"`cell.type.decomposition` must be a matrix!")
  }
  if(!all(colnames(cell.type.decomposition)==levels(group))){
    stop(cli.symbol("fail"),"Please check your `cell.type.decomposition`'s colnames, and make sure they are the same as the cell groups' names in the SpatialCellChat.")
  }
  if(!all(rownames(cell.type.decomposition)==rownames(object@images$coordinates))){
    stop(cli.symbol("fail"),"Please check your `cell.type.decomposition`'s rownames, and make sure they are the same as the cells' names in the SpatialCellChat.")
  }
  cell.type.onehot <- 1*(cell.type.decomposition>0)

  numCluster <- nlevels(group)
  # if (numCluster != length(unique(group))) {
  #   stop("Please check `unique(object@idents)` and ensure that the factor levels are correct!
  #        You may need to drop unused levels using 'droplevels' function. e.g.,
  #        `meta$labels = droplevels(meta$labels, exclude = setdiff(levels(meta$labels),unique(meta$labels)))`")
  # }

  if ( is.null(object@net$prob.cell) ) {
    stop(cli.symbol(2),"Please run `computeCommunProb` to compute the communication probability/strength between any interacting individual cells! ")
  } else {
    prob.cell <- object@net$prob.cell # nC x nC x nLR
    prob.cell_ <- object@net$tmp$prob.cell # a list
    nC <- dim(prob.cell)[[1]]
  }

  LRsig <- dimnames(prob.cell)[[3]]
  nLR <- length(LRsig)
  interaction_input <- object@DB$interaction
  pairLRsig <- interaction_input[LRsig, , drop = FALSE]

  avg.type <- match.arg(avg.type)
  if(avg.type=="avg"){
    computeAvgCommunProb_LR_Visium <- computeAvgCommunProb_LR_Visium_Avg
  } else if (avg.type=="sum"){
    computeAvgCommunProb_LR_Visium <- computeAvgCommunProb_LR_Visium_Sum
  }
  gc()

  if (colocalization.use) {
    data.spatial <- object@images$coordinates
    pval.colo = computeColocalization(coordinates = data.spatial, group = group, nboot = nboot, seed.use = seed.use)
  } else {
    pval.colo <- matrix(0, nrow = numCluster, ncol = numCluster)
  }

  cat(paste0(cli.symbol(),'Compute group-level cell-cell communication... <<< [', Sys.time(),']'),'\n')

  Prob <- array(0, dim = c(numCluster,numCluster,nLR))
  Pval <- array(1, dim = c(numCluster,numCluster,nLR))
  dimnames(Prob) <- list(levels(group), levels(group), rownames(pairLRsig))
  dimnames(Pval) <- dimnames(Prob)

  set.seed(seed.use)

  # retain dim-3, sum up dim-1 && dim-2, `prob.sum` stores each LR's number of cell-level links/interactions
  prob.sum <- purrr::map_dbl(
    .x = prob.cell_,
    .f = function(Mat){
      return(length(Mat@x))
    }
  )
  names(prob.sum) <- LRsig
  object@net$tmp$LRsig.CCC.counts <- prob.sum

  LRsig.use.idx <- which(prob.sum > 0)

  object@net$tmp$LRsig.use.idx <- LRsig.use.idx
  gc()

  if(length(LRsig.use.idx) < 1){
    stop("Each LR pair does not have any spot-level links/interactions.")
  }

  cat(cli.symbol(),"compute the average signaling per cell group...\n")
  Prob.avg_ <- my_future_sapply(
    X = seq_len(length(LRsig.use.idx)),
    FUN = function(x) {
      i <- LRsig.use.idx[[x]] # one LR pair index

      # compute the average signaling per cell group
      prob.cell.i <- prob.cell_[[i]]

      Prob.avg <-
        computeAvgCommunProb_LR_Visium(
          prob.cell.i,
          cell.type.proportion = cell.type.decomposition,
          cell.type.onehot = cell.type.onehot
        )

      Prob.avg[pval.colo > thresh.colo] <- 0
      # Prob: array(0, dim = c(numCluster,numCluster,nLR))
      gc()
      return(Prob.avg)
    },
    simplify = F # return a list
  )

  for (x in seq_len( length(LRsig.use.idx) ) ) {
    i <- LRsig.use.idx[[x]]
    Prob[ , , i] <- Prob.avg_[[x]]
  }

  # update `prob.sum` & `LRsig.use.idx` to do permutation
  # retain dim-3, sum up dim-1 && dim-2, `prob.sum` stores each LR's number of group-level links/interactions
  prob.sum <- apply(Prob > 0, 3, sum) # return a named vector
  # each LR's number of group-level links/interactions
  object@net$tmp$LRsig.GGC.counts <- prob.sum

  LRsig.use.idx <- which(prob.sum > 0)
  if (do.permutation) {
    cat(paste0(cli.symbol(),'Perform permutation test for group-level communication... <<< [', Sys.time(),']'),'\n')
    permutation <- replicate(nboot, sample.int(numCluster, size = numCluster)) # numCluster x nboot
    Pval_ <- my_future_lapply(
      # LRsig.use.idx is a numeric vector
      X = seq_len(length(LRsig.use.idx)),
      FUN = function(x){
        i <- LRsig.use.idx[[x]] # one LR pair index

        # compute the average signaling per cell group after permutation
        prob.cell.i <- prob.cell_[[i]]
        Pnull <- as.vector(Prob[ , , i])

        Pboot <- sapply(
          X = 1:nboot,
          FUN = function(nE) {
            # permutation is nGroups x nboot
            decomposition.boot <- cell.type.decomposition[ ,permutation[, nE]] # group labels order unchanged, change ratio order
            onehot.boot <- 1*(decomposition.boot>0)
            Pboot.avg <- computeAvgCommunProb_LR_Visium(
              prob = prob.cell.i,
              cell.type.proportion = decomposition.boot,
              cell.type.onehot = onehot.boot
            )
            return(as.vector(Pboot.avg))
          }
        )
        gc()
        Pboot <- matrix(unlist(Pboot), nrow=length(Pnull), ncol = nboot, byrow = FALSE)

        nReject <- rowSums(Pboot - Pnull > 0)
        p = nReject/nboot
        Pval.i <- matrix(p, nrow = numCluster, ncol = numCluster, byrow = FALSE)
        return(Pval.i)
      },
      simplify = F, # return a list
      hint.message = "do permutation..."
    )

    for (x in seq_len(length(LRsig.use.idx))) {
      # get correct index
      i <- LRsig.use.idx[[x]]
      # update the values
      Pval[ , , i] <- Pval_[[x]]
    }

    Pval[Prob == 0] <- 1
  } else {
    Pval <- NULL
  }

  object@net$prob <- Prob
  object@net$pval <- Pval
  # object@options$parameter$min.percent <- min.percent;
  # object@options$parameter$min.cells.sr <- min.cells.sr;
  object@options$parameter$do.permutation <- do.permutation;
  object@options$parameter$nboot <- nboot;
  object@options$parameter$avg.type <- avg.type;
  object@options$parameter$seed.use <- seed.use;
  object@options$parameter$colocalization.use <- colocalization.use;
  object@options$parameter$thresh.colo <- thresh.colo;
  # object@options$parameter$do.filter <- do.filter;
  # object@net$tmp$Lavg <- NULL;object@net$tmp$Ravg <- NULL; # clean the cache
  object@net$tmp$cell.type.decomposition <- cell.type.decomposition

  if (colocalization.use) {
    object@images$colocalization <- pval.colo
  }
  cat(paste0(cli.symbol(symbol = "success"),'Inference of group-level cell-cell communication is done. Parameter values are stored in `object@options$parameter` <<< [', Sys.time(),']'))
  return(object)
}


#' @title computeAvgCommunProb_LR_Visium_Sum
#'
#' @description
#' Compute average communication probabilities of pairwise cell groups for one particular ligand-receptor pair/signaling pathway
#' for 10X Visium Low Definition.
#'
#' @param prob a matrix of communication probabilities for pairwise individual cells for one particular ligand-receptor pair/signaling pathway
#' @param cell.type.proportion Ncells x NTypes Matrix, which comes from any kind of cell type decomposition method.
#' @param cell.type.onehot NULL. Placeholder.
#'
#' @return A NTypes x NTypes matrix of communication probabilities for pairwise individual cell types for one particular ligand-receptor pair/signaling pathway
#' @export
computeAvgCommunProb_LR_Visium_Sum  <- function (
    prob,
    cell.type.proportion,
    cell.type.onehot = NULL
){
  # calculate communication strength across cell type pairs
  Prob.avg <- Matrix::crossprod(x = cell.type.proportion,
                                y = prob %*% cell.type.proportion)
  Prob.avg <- as.matrix(Prob.avg)
  return(Prob.avg)
}


#' @title computeAvgCommunProb_LR_Visium_Avg
#'
#' @description
#' Compute average communication probabilities of pairwise cell groups for one particular ligand-receptor pair/signaling pathway
#' for 10X Visium Low Definition.
#'
#' @param prob a matrix of communication probabilities for pairwise individual cells for one particular ligand-receptor pair/signaling pathway
#' @param cell.type.proportion Ncells x NTypes Matrix, which comes from any kind of cell type decomposition method.
#' @param cell.type.onehot Ncells x NTypes Matrix, which differs from \code{cell.type.proportion}, can get through \code{cell.type.onehot <- 1*(cell.type.proportion>0)}
#'
#' @return A NTypes x NTypes matrix of communication probabilities for pairwise individual cell types for one particular ligand-receptor pair/signaling pathway
#' @export
computeAvgCommunProb_LR_Visium_Avg  <- function (
    prob,
    cell.type.proportion,
    cell.type.onehot
){
  # code for test
  # prob <- prob.cell.i
  # group <- group
  # dataLR = dataLR_temp
  # min.percent = 0.1
  # min.cells.sr = 5
  # cell.type.proportion <- cell.type.decomposition

  # calculate communication strength across cell type pairs
  Prob.avg <- Matrix::crossprod(x = cell.type.proportion,
                                y = prob %*% cell.type.proportion)

  # count communication links  across cell type pairs
  prob@x <- rep.int(1, times = length(prob@x))
  Prob.scale.factor <- Matrix::crossprod(x = cell.type.onehot,
                                         y = prob %*% cell.type.onehot)
  Prob.avg <- Prob.avg/Prob.scale.factor
  Prob.avg <- as.matrix(Prob.avg)
  Prob.avg[is.nan(Prob.avg)] <- 0
  return(Prob.avg)
}



# #' @title computeAvgCommunProb_LR_Visium
# #'
# #' @description
# #' Compute average communication probabilities of pairwise cell groups for one particular ligand-receptor pair/signaling pathway
# #' for 10X Visium Low Definition.
# #'
# #' @param prob a matrix of communication probabilities for pairwise individual cells for one particular ligand-receptor pair/signaling pathway
# #' @param group cell group information used for computing averaged communication probabilities
# #' @param cell.type.decomposition Ncells x NTypes Matrix, which comes from any kind of cell type decomposition.
# #' @param dataLR a nCell*2 expression matrix of a given pair of ligand-receptor.
# #' @param do.filter Boolean. Whether to filter out some interactions of a given pair of ligand-receptor
# #' between several cell groups according to the parameters `min.percent` and `min.cells.sr`
# #' when computing the average communication probability of a given pair of ligand-receptor between cell groups.
# #' Default is TRUE.
# #' @param min.percent Numeric from 0 to 1. Minimum percentage of expressed ligands or receptors per cell group to require for computing the group-level signaling
# #' Default is 0.1.
# #' @param min.cells.sr Integer greater than 0. Minimum number of cells required as senders or receivers per cell group for computing the group-level signaling
# #' Default is 5.
# #'
# #' @return
# #' @export
# computeAvgCommunProb_LR_Visium  <- function (
#     prob,
#     group,
#     cell.type.decomposition,
#     dataLR = NULL,
#     do.filter = T,
#     min.percent = 0.1,
#     min.cells.sr = 5
# ){
#   # code for test
#   # prob <- prob.cell.i
#   # group <- group
#   # dataLR = dataLR_temp
#   # min.percent = 0.1
#   # min.cells.sr = 5
#
#
#   # group is vector/factor
#   if (!is.factor(group)) {
#     group <- factor(group)
#   }
#   cell.type.mat <- cell.type.decomposition
#
#   if(do.filter){
#     # min.percent
#     dataLR_temp <- 1 * (dataLR > 0)
#     dataLR_temp <- aggregate(dataLR_temp, list(group), FUN = mean)
#     dataLR_percent <- 1 * (format(dataLR_temp[, -1], digits = 1) >= min.percent)
#     Prob_percent <- Matrix::crossprod(matrix(dataLR_percent[, 1], nrow = 1), matrix(dataLR_percent[, 2], nrow = 1))
#
#     if(sum(Prob_percent)==0){
#       Prob.avg <- Prob_percent # NGroup x NGroup All-zero Matrix
#     } else {
#
#       Prob.avg <- Matrix::crossprod(x = cell.type.mat,
#                                     y = prob %*% cell.type.mat)
#
#       # min.cells.sr
#       prob@x <- rep.int(1, times = length(prob@x))
#       sender.counts <- Matrix::rowSums(prob)
#       receptor.counts <- Matrix::colSums(prob)
#       cells.sr <- cbind(sender.counts, receptor.counts)
#       cells.sr <- aggregate(cells.sr, list(group), FUN = sum)
#       cells.sr <- 1 * (cells.sr[,-1] >= min.cells.sr)
#       cells.sr <- Matrix::crossprod(matrix(cells.sr[, 1], nrow = 1), matrix(cells.sr[, 2], nrow = 1))
#
#       Prob.avg <- Prob.avg * Prob_percent * cells.sr
#     }
#   } else {
#     Prob.avg <- Matrix::crossprod(x = cell.type.mat,
#                                   y = prob %*% cell.type.mat)
#   }
#
#   # return
#   Prob.avg <- as.matrix(Prob.avg)
#   return(Prob.avg)
# }




#' Compute group-level cell-cell communication
#'
#' @param object SpatialCellChat object with communication probabilities for pairwise individual cells
#' @param group.by cell group information used for computing average communication probabilities
#' @param avg.type methods for integrating communication probabilities per cell group
# @param type methods for computing the average gene expression per cell group.
#
# By default = "triMean", defined as a weighted average of the distribution's median and its two quartiles (https://en.wikipedia.org/wiki/Trimean);
#
# When setting `type = "truncatedMean"`, a value should be assigned to 'trim'. See the function `base::mean`.
#
# @param trim the fraction (0 to 0.25) of observations to be trimmed from each end of x before the mean is computed.
#' @param do.permutation whether performing permutation test
#' @param nboot the number of permutations
#' @param seed.use set a random seed. By default, set the seed to 1.
#' @param colocalization.use whether filtering out spatially distant cell groups based on colocalization analysis between any cell groups
#' @param thresh.colo removal of cell-cell communication with no significant colocalizations (fdr < 0.05)
#' @inheritParams computeAvgCommunProb_LR_Avg
#' @inheritParams computeAvgCommunProb_LR_Sum
#'
#' @return A CellChat object with updated slot 'net':
#'
#' object@net$prob is the inferred group-level communication probability (strength) array, where the first, second and third dimensions represent a source group, target group and ligand-receptor pair, respectively.
#'
#' USER can access all the inferred cell-cell communications using the function 'subsetCommunication(object)', which returns a data frame.
#'
#' object@net$pval is the corresponding p-values of each interaction

#' @export
#'
computeAvgCommunProb <- function(
    object,
    group.by = NULL,
    avg.type = c("avg","sum"),
    min.percent = 0.1,
    min.cells.sr = 5,
    do.permutation = T,
    nboot = 100,
    seed.use = 1L,
    colocalization.use = F,
    thresh.colo = 0.05
){
  if (is.null(group.by)) {
    group <- object@idents
  } else {
    if (!(group.by %in% colnames(object@meta))) {
      stop("The 'group.by' is not a column name in the `object@meta`, which will be used for cell grouping.")
    } else {
      group <- object@meta[[group.by]]
    }
    if (!is.factor(group)) {
      group <- factor(group)
    }
  }
  cat(cli.symbol(),"The cell groups used for averaging cell-cell communication are ", cli::col_red(levels(group)), '\n')
  numCluster <- nlevels(group)
  if (numCluster != length(unique(group))) {
    stop("Please check `unique(object@idents)` and ensure that the factor levels are correct!
         You may need to drop unused levels using 'droplevels' function. e.g.,
         `meta$labels = droplevels(meta$labels, exclude = setdiff(levels(meta$labels),unique(meta$labels)))`")
  }

  if (object@options$parameter$raw.use) {
    data <- object@data.signaling
    # scale the elements
    data@x <- data@x/max(data@x)
    data.use <- as.matrix(data)

  } else {
    data <- object@data.project
    # scale
    data.use <- data/max(data)
  }

  nC <- ncol(data.use)

  if ( is.null(object@net$prob.cell) ) {
    stop(cli.symbol(2),"Please run `computeCommunProb` to compute the communication probability/strength between any interacting individual cells! ")
  } else {
    prob.cell <- object@net$prob.cell
    prob.cell_ <- object@net$tmp$prob.cell # a list
  }

  LRsig <- dimnames(prob.cell)[[3]]
  nLR <- length(LRsig)

  interaction_input <- object@DB$interaction
  complex_input <- object@DB$complex
  cofactor_input <- object@DB$cofactor

  pairLRsig <- interaction_input[LRsig, , drop = FALSE]
  # geneL <- as.character(pairLRsig$ligand)
  # geneR <- as.character(pairLRsig$receptor)

  # compute the expression of ligand or receptor
  # cat(cli.symbol(),"compute the expression of ligand or receptor...\n")
  # dataLavg <- computeExpr_LR(geneL, data.use, complex_input)
  # dataRavg <- computeExpr_LR(geneR, data.use, complex_input)
  #
  # cat(cli.symbol(),"take account into the effect of co-activation and co-inhibition receptors...\n")
  # dataRavg.co.A.receptor <- computeExpr_coreceptor(cofactor_input, data.use, pairLRsig, type = "A")
  # dataRavg.co.I.receptor <- computeExpr_coreceptor(cofactor_input, data.use, pairLRsig, type = "I")
  # dataRavg <- dataRavg * dataRavg.co.A.receptor/dataRavg.co.I.receptor
  dataLavg <- object@net$tmp$Lavg
  dataRavg <- object@net$tmp$Ravg

  avg.type <- match.arg(avg.type)
  if(avg.type=="avg"){
    computeAvgCommunProb_LR <- computeAvgCommunProb_LR_Avg
  } else if (avg.type=="sum"){
    computeAvgCommunProb_LR <- computeAvgCommunProb_LR_Sum
  }

  gc()

  if (colocalization.use) {
    data.spatial <- object@images$coordinates
    pval.colo = computeColocalization(coordinates = data.spatial, group = group, nboot = nboot, seed.use = seed.use)
  } else {
    pval.colo <- matrix(0, nrow = numCluster, ncol = numCluster)
  }

  cat(paste0(cli.symbol(),'Compute group-level cell-cell communication... <<< [', Sys.time(),']'),'\n')

  Prob <- array(0, dim = c(numCluster,numCluster,nLR))
  Pval <- array(1, dim = c(numCluster,numCluster,nLR))
  dimnames(Prob) <- list(levels(group), levels(group), rownames(pairLRsig))
  dimnames(Pval) <- dimnames(Prob)

  set.seed(seed.use)

  # retain dim-3, sum up dim-1 && dim-2, `prob.sum` stores each LR's number of cell-level links/interactions
  prob.sum <- purrr::map_dbl(
    .x = prob.cell_,
    .f = function(Mat){
      return(length(Mat@x))
    }
  )
  names(prob.sum) <- LRsig
  object@net$tmp$LRsig.CCC.counts <- prob.sum

  # Previous code:
  # prob.sum <- c()
  # for (i in 1:dim(prob.cell)[[3]]) {
  #   prob.sum[i] <- sum(prob.cell[,,i] > 0)
  # }

  LRsig.use.idx <- which(prob.sum > 0)

  object@net$tmp$LRsig.use.idx <- LRsig.use.idx
  gc()

  if(length(LRsig.use.idx) < 1){
    stop("Each LR pair does not have any cell-level links/interactions.")
  }

  cat(cli.symbol(),"compute the average signaling per cell group...\n")
  Prob.avg_ <- my_future_sapply(
    X = seq_len(length(LRsig.use.idx)),
    FUN = function(x) {
      i <- LRsig.use.idx[[x]] # one LR pair index
      # compute the average signaling per cell group
      prob.cell.i <- prob.cell_[[i]]
      dataLR_temp <- cbind(dataLavg[i, ], dataRavg[i, ])
      Prob.avg <-
        computeAvgCommunProb_LR(
          prob.cell.i,
          group = group,
          dataLR = dataLR_temp,
          min.percent = min.percent,
          min.cells.sr = min.cells.sr
        )
      # Pnull <- as.vector(Prob.avg)
      Prob.avg[pval.colo > thresh.colo] <- 0
      # Prob: array(0, dim = c(numCluster,numCluster,nLR))
      gc()
      return(Prob.avg)
    },
    simplify = F # return a list
  )

  for (x in seq_len( length(LRsig.use.idx) ) ) {
    i <- LRsig.use.idx[[x]]
    Prob[ , , i] <- Prob.avg_[[x]]
  }

  # update `prob.sum` & `LRsig.use.idx` to do permutation
  # retain dim-3, sum up dim-1 && dim-2, `prob.sum` stores each LR's number of group-level links/interactions before permutation
  prob.sum <- apply(Prob > 0, 3, sum) # return a named vector
  # each LR's number of group-level links/interactions
  object@net$tmp$LRsig.GGC.counts <- prob.sum

  LRsig.use.idx <- which(prob.sum > 0)
  if (do.permutation) {
    cat(paste0(cli.symbol(),'Perform permutation test for group-level communication... <<< [', Sys.time(),']'),'\n')
    permutation <- replicate(nboot, sample.int(nC, size = nC))
    Pval_ <- my_future_lapply(
      # LRsig.use.idx is a numeric vector
      X = seq_len(length(LRsig.use.idx)),
      FUN = function(x){
        i <- LRsig.use.idx[[x]]
        # compute the average signaling per cell group after permutation
        prob.cell.i <- prob.cell_[[i]]
        dataLR_temp <- cbind(dataLavg[i,], dataRavg[i,])
        Pnull <- as.vector(Prob[ , , i])

        Pboot <- sapply(
          X = 1:nboot,
          FUN = function(nE) {
            groupboot <- group[permutation[, nE]]
            Pboot.avg <- computeAvgCommunProb_LR(
              prob.cell.i,
              group = groupboot,
              dataLR = dataLR_temp,
              min.percent = min.percent,
              min.cells.sr = min.cells.sr
            )
            return(as.vector(Pboot.avg))
          }
        )
        gc()
        Pboot <- matrix(unlist(Pboot), nrow=length(Pnull), ncol = nboot, byrow = FALSE)
        nReject <- rowSums(Pboot - Pnull > 0)
        p = nReject/nboot
        Pval.i <- matrix(p, nrow = numCluster, ncol = numCluster, byrow = FALSE)
        return(Pval.i)
      },
      simplify = F, # return a list
      hint.message = "do permutation..."
    )

    for (x in seq_len(length(LRsig.use.idx))) {
      # get correct index
      i <- LRsig.use.idx[[x]]
      # update the values
      Pval[ , , i] <- Pval_[[x]]
    }

    Pval[Prob == 0] <- 1

  } else {
    Pval <- NULL
  }

  # Pval[Prob == 0] <- 1
  # dimnames(Prob) <- list(levels(group), levels(group), rownames(pairLRsig))
  # dimnames(Pval) <- dimnames(Prob)
  object@net$prob <- Prob
  object@net$pval <- Pval
  object@options$parameter$min.percent <- min.percent
  object@options$parameter$min.cells.sr <- min.cells.sr
  object@options$parameter$do.permutation <- do.permutation

  object@options$parameter$nboot <- nboot
  object@options$parameter$avg.type <- avg.type
  object@options$parameter$seed.use <- seed.use
  object@options$parameter$colocalization.use <- colocalization.use

  object@options$parameter$thresh.colo <- thresh.colo
  # object@net$tmp$Lavg <- NULL;object@net$tmp$Ravg <- NULL; # clean the cache

  if (colocalization.use) {
    object@images$colocalization <- pval.colo
  }
  cat(paste0(cli.symbol(symbol = "success"),'Inference of group-level cell-cell communication is done. Parameter values are stored in `object@options$parameter` <<< [', Sys.time(),']'))
  return(object)
}


#' Compute average communication probabilities of pairwise cell groups for one particular ligand-receptor pair/signaling pathway
#'
#' @param prob a matrix of communication probabilities for pairwise individual cells for one particular ligand-receptor pair/signaling pathway
#' @param dataLR a nCell*2 data matrix of a given pair of ligand-receptor.
#' @param group Character vector. Cell group information used for computing averaged communication probabilities
# #' @param type methods for computing the average gene expression per cell group.
# #'
# #' By default = "triMean", defined as a weighted average of the distribution's median and its two quartiles (https://en.wikipedia.org/wiki/Trimean);
# #'
# #' When setting `type = "truncatedMean"`, a value should be assigned to 'trim'. See the function `base::mean`.
# #'
# #' @param trim the fraction (0 to 0.25) of observations to be trimmed from each end of x before the mean is computed.
#' @param min.percent Numeric from 0 to 1. Minimum percentage of expressed ligands or receptors per cell group to require for computing the group-level signaling
#' Default is 0.1.
#' @param min.cells.sr Integer greater than 0. Minimum number of cells required as senders or receivers per cell group for computing the group-level signaling
#' Default is 5.
#'
#' @return Returns a matrix containing the interaction weights between any two cell groups.
#' @export
computeAvgCommunProb_LR_Sum <- function (
    prob,
    group,
    dataLR = NULL,
    min.percent = 0.1,
    min.cells.sr = 5
){
  # code for test
  # prob <- prob.cell.i
  # group <- group
  # dataLR = dataLR_temp
  # min.percent = 0.1
  # min.cells.sr = 5

  # group is vector/factor
  if (!is.factor(group)) {
    group <- factor(group)
  }
  cell.type.mat <- model.matrix(~group-1)


  # min.percent
  dataLR_temp <- 1 * (dataLR > 0)
  dataLR_temp <- aggregate(dataLR_temp, list(group), FUN = mean)
  dataLR_percent <- 1 * (format(dataLR_temp[, -1], digits = 1) >= min.percent)
  Prob_percent <- Matrix::crossprod(matrix(dataLR_percent[, 1], nrow = 1), matrix(dataLR_percent[, 2], nrow = 1))

  if(sum(Prob_percent)==0){
    Prob.avg <- Prob_percent # NGroup x NGroup All-zero Matrix
  } else {

    Prob.avg <- Matrix::crossprod(x = cell.type.mat,
                                  y = prob %*% cell.type.mat)

    # min.cells.sr
    prob@x <- rep.int(1, times = length(prob@x))
    sender.counts <- Matrix::rowSums(prob)
    receptor.counts <- Matrix::colSums(prob)
    cells.sr <- cbind(sender.counts, receptor.counts)
    cells.sr <- aggregate(cells.sr, list(group), FUN = sum)
    cells.sr <- 1 * (cells.sr[,-1] >= min.cells.sr)
    cells.sr <- Matrix::crossprod(matrix(cells.sr[, 1], nrow = 1), matrix(cells.sr[, 2], nrow = 1))

    Prob.avg <- Prob.avg * Prob_percent * cells.sr
  }

  # return
  Prob.avg <- as.matrix(Prob.avg)
  return(Prob.avg)
}


#' Compute average communication probabilities of pairwise cell groups for one particular ligand-receptor pair/signaling pathway
#'
#' @param prob a matrix of communication probabilities for pairwise individual cells for one particular ligand-receptor pair/signaling pathway
#' @param dataLR a nCell*2 data matrix of a given pair of ligand-receptor.
#' @param group Character vector. Cell group information used for computing averaged communication probabilities
# #' @param type methods for computing the average gene expression per cell group.
# #'
# #' By default = "triMean", defined as a weighted average of the distribution's median and its two quartiles (https://en.wikipedia.org/wiki/Trimean);
# #'
# #' When setting `type = "truncatedMean"`, a value should be assigned to 'trim'. See the function `base::mean`.
# #'
# #' @param trim the fraction (0 to 0.25) of observations to be trimmed from each end of x before the mean is computed.
#' @param min.percent Numeric from 0 to 1. Minimum percentage of expressed ligands or receptors per cell group to require for computing the group-level signaling
#' Default is 0.1.
#' @param min.cells.sr Integer greater than 0. Minimum number of cells required as senders or receivers per cell group for computing the group-level signaling
#' Default is 5.
#'
#' @return Returns a matrix containing the interaction weights between any two cell groups.
#' @export
computeAvgCommunProb_LR_Avg <- function (
    prob,
    group,
    dataLR = NULL,
    min.percent = 0.1,
    min.cells.sr = 5
){
  # code for test
  # prob <- prob.cell.i
  # group <- group
  # dataLR = dataLR_temp
  # min.percent = 0.1
  # min.cells.sr = 5

  # group is vector/factor
  if (!is.factor(group)) {
    group <- factor(group)
  }
  cell.type.mat <- model.matrix(~group-1)


  # min.percent
  dataLR_temp <- 1 * (dataLR > 0)
  dataLR_temp <- aggregate(dataLR_temp, list(group), FUN = mean)
  dataLR_percent <- 1 * (format(dataLR_temp[, -1], digits = 1) >= min.percent)
  Prob_percent <- Matrix::crossprod(matrix(dataLR_percent[, 1], nrow = 1), matrix(dataLR_percent[, 2], nrow = 1))

  if(sum(Prob_percent)==0){
    Prob.avg <- Prob_percent # NGroup x NGroup All-zero Matrix
  } else {
    Prob.avg <- Matrix::crossprod(x = cell.type.mat,
                                  y = prob %*% cell.type.mat)

    # min.cells.sr
    prob@x <- rep.int(1, times = length(prob@x))
    Prob.scale.factor <- Matrix::crossprod(x = cell.type.mat,
                                           y = prob %*% cell.type.mat)
    Prob.avg <- Prob.avg/Prob.scale.factor
    Prob.avg[is.nan(Prob.avg)] <- 0

    sender.counts <- Matrix::rowSums(prob)
    receptor.counts <- Matrix::colSums(prob)
    cells.sr <- cbind(sender.counts, receptor.counts)
    cells.sr <- aggregate(cells.sr, list(group), FUN = sum)
    cells.sr <- 1 * (cells.sr[,-1] >= min.cells.sr)
    cells.sr <- Matrix::crossprod(matrix(cells.sr[, 1], nrow = 1), matrix(cells.sr[, 2], nrow = 1))

    Prob.avg <- Prob.avg * Prob_percent * cells.sr
  }
  Prob.avg <- as.matrix(Prob.avg)

  return(Prob.avg)
}



# #' Compute average communication probabilities of pairwise cell groups for one particular ligand-receptor pair/signaling pathway
# #'
# #' @param prob a matrix of communication probabilities for pairwise individual cells for one particular ligand-receptor pair/signaling pathway
# #' @param dataLR a nCell*2 data matrix of a given pair of ligand-receptor.
# #' @param group Character vector. Cell group information used for computing averaged communication probabilities
# #' @param type methods for computing the average gene expression per cell group.
# #'
# #' By default = "triMean", defined as a weighted average of the distribution's median and its two quartiles (https://en.wikipedia.org/wiki/Trimean);
# #'
# #' When setting `type = "truncatedMean"`, a value should be assigned to 'trim'. See the function `base::mean`.
# #'
# #' @param trim the fraction (0 to 0.25) of observations to be trimmed from each end of x before the mean is computed.
# #' @param min.percent Numeric from 0 to 1. Minimum percentage of expressed ligands or receptors per cell group to require for computing the group-level signaling
# #' Default is 0.1.
# #' @param min.cells.sr Integer greater than 0. Minimum number of cells required as senders or receivers per cell group for computing the group-level signaling
# #' Default is 5.
# #'
# #' @return Returns a matrix containing the interaction weights between any two cell groups.
# #' @export
# computeAvgCommunProb_LR <- function (
#     prob,
#     group,
#     dataLR = NULL,
#     min.percent = 0.1,
#     min.cells.sr = 5
# ){
#   # code for test
#   # prob <- prob.cell.i
#   # group <- group
#   # dataLR = dataLR_temp
#   # min.percent = 0.1
#   # min.cells.sr = 5
#
#   # group is vector/factor
#   if (!is.factor(group)) {
#     group <- factor(group)
#   }
#   cell.type.mat <- model.matrix(~group-1)
#
#
#   # min.percent
#   dataLR_temp <- 1 * (dataLR > 0)
#   dataLR_temp <- aggregate(dataLR_temp, list(group), FUN = mean)
#   dataLR_percent <- 1 * (format(dataLR_temp[, -1], digits = 1) >= min.percent)
#   Prob_percent <- Matrix::crossprod(matrix(dataLR_percent[, 1], nrow = 1), matrix(dataLR_percent[, 2], nrow = 1))
#
#   if(sum(Prob_percent)==0){
#     Prob.avg <- Prob_percent # NGroup x NGroup All-zero Matrix
#   } else {
#
#     Prob.avg <- Matrix::crossprod(x = cell.type.mat,
#                                   y = prob %*% cell.type.mat)
#
#     # min.cells.sr
#     prob@x <- rep.int(1, times = length(prob@x))
#     sender.counts <- Matrix::rowSums(prob)
#     receptor.counts <- Matrix::colSums(prob)
#     cells.sr <- cbind(sender.counts, receptor.counts)
#     cells.sr <- aggregate(cells.sr, list(group), FUN = sum)
#     cells.sr <- 1 * (cells.sr[,-1] >= min.cells.sr)
#     cells.sr <- Matrix::crossprod(matrix(cells.sr[, 1], nrow = 1), matrix(cells.sr[, 2], nrow = 1))
#
#     Prob.avg <- Prob.avg * Prob_percent * cells.sr
#   }
#
#   # return
#   Prob.avg <- as.matrix(Prob.avg)
#   return(Prob.avg)
# }
#


# computeAvgCommunProb_LR <- function (prob, group, dataLR = NULL, min.percent = 0.1, min.cells.sr = 5)
# {
#   if (!is.factor(group)) {
#     group <- factor(group)
#   }
#   level.use <- levels(group)
#   level.use <- level.use[level.use %in% unique(group)]
#   numCluster <- length(level.use)
#
#   # min.percent
#   dataLR_temp <- 1 * (dataLR > 0)
#
#   dataLR_temp <- aggregate(dataLR_temp, list(group), FUN = mean)
#   dataLR_percent <- 1 * (format(dataLR_temp[, -1], digits = 1) >= min.percent)
#   Prob_percent <- Matrix::crossprod(matrix(dataLR_percent[, 1], nrow = 1), matrix(dataLR_percent[, 2], nrow = 1))
#   melt_Prob_percent <- reshape2::melt(Prob_percent,value.name = "Binary.Prob")
#
#   Prob.avg <- matrix(0, nrow = numCluster, ncol = numCluster)
#   # min.cells.sr
#   if(NROW(melt_Prob_percent[melt_Prob_percent$Binary.Prob==1, ,drop=F]) == 0){
#     return(Prob.avg)
#   } else {
#     df_nonZero_Prob <- melt_Prob_percent[melt_Prob_percent$Binary.Prob==1, ,drop=F]
#     for (x in seq_len(NROW(df_nonZero_Prob))) {
#
#       ii <- df_nonZero_Prob[x,1,drop=T]
#       jj <- df_nonZero_Prob[x,2,drop=T]
#       prob.ij <- prob[group %in% level.use[ii], group %in% level.use[jj], drop = FALSE]
#       if ((sum(Matrix::rowSums(prob.ij > 0) > 0) >= min.cells.sr) &
#           (sum(Matrix::colSums(prob.ij > 0) > 0) >= min.cells.sr)) {
#
#         # Prob.avg[ii, jj] <- sum(prob.ij@x)/(cumprod(prob.ij@Dim)[[2]])
#         Prob.avg[ii, jj] <- sum(prob.ij@x)
#       }
#     } # forloop
#     return(Prob.avg)
#   }
# }


# Prob.avg <- matrix(0, nrow = numCluster, ncol = numCluster)
# for (ii in 1:numCluster) {
#   for (jj in 1:numCluster) {
#     Prob.temp <- as.vector(P.spatial[group %in% level.use[ii], group %in% level.use[jj]]) > 0
#     # Prob.sum <- sum(Prob.temp)
#     # Prob.avg[ii,jj] <- Prob.sum/(sum(group %in% level.use[ii]) * sum(group %in% level.use[jj]))
#     Prob.avg[ii,jj] <- sum(Prob.temp)/(sum(group %in% level.use[ii]) * sum(group %in% level.use[jj]))
#   }
# }


#' Compute the communication probability on signaling pathway level by summarizing all related ligands/receptors
#'
#' @param object CellChat object
#' @param net A list from object@net; If net = NULL, net = object@net
#' @param pairLR.use A dataframe giving the ligand-receptor interactions; If pairLR.use = NULL, pairLR.use = object@LR$LRsig
#' @param thresh threshold of the p-value for determining significant interaction
#' @param do.group whether to compute the group-level signaling based on the cell group information in `object@idents`
#' @param do.cell whether to compute the individual-cell signaling at signaling pathway level. This works when "prob.cell" exists in `object@net`.
#' @return A CellChat object with updated slot 'netP':
#'
#' object@netP$prob is the communication probability array on signaling pathway level; USER can convert this array to a data frame using the function 'reshape2::melt()',
#'
#' e.g., `df.netP <- reshape2::melt(object@netP$prob, value.name = "prob"); colnames(df.netP)[1:3] <- c("source","target","pathway_name")` or access all significant interactions using the function \code{\link{subsetCommunication}}
#'
#' object@netP$pathways list all the signaling pathways with significant communications.
#'
#' From version >= 1.1.0, pathways are ordered based on the total communication probabilities. NB: pathways with small total communication probabilities might be also very important since they might be specifically activated between only few cell types.
#'
#' @export

computeCommunProbPathway <- function(
    object = NULL,
    net = NULL,
    pairLR.use = NULL,
    thresh = 0.05,
    do.group = TRUE,
    do.cell = TRUE
) {
  if (is.null(net)) {
    net <- object@net
  }
  if (is.null(pairLR.use)) {
    pairLR.use <- object@LR$LRsig
  }
  if (do.group) {
    if ( is.null(net$prob) ) {
      stop("Please run `computeAvgCommunProb` to compute the group-level signaling! ")
    }
    cat(cli.symbol(),"Compute the communication probability between cell groups at signaling pathway level by summarizing all related ligands/receptors...\n")
    prob <- net$prob
    prob[net$pval >= thresh] <- 0
    pathways <- unique(pairLR.use$pathway_name)

    group <- factor(pairLR.use$pathway_name, levels = pathways)

    prob.pathways <- aperm(apply(prob, c(1, 2), by, group, sum), c(2, 3, 1))
    # `apply(prob.pathways, 3, sum) != 0` means group-group signaling exists in the specific pathway
    # so `pathways.sig` stores the pathways whose group-group signaling number is no zero
    pathways.sig <- pathways[apply(prob.pathways, 3, sum) != 0]

    # subset the `prob.pathways`
    prob.pathways.sig <- prob.pathways[,,pathways.sig, drop = FALSE]
    # sort `prob.pathways.sig` according to group-group communication probability sum of `prob.pathways.sig`
    idx <- sort(apply(prob.pathways.sig, 3, sum), decreasing=TRUE, index.return = TRUE)$ix
    pathways.sig <- pathways.sig[idx]
    prob.pathways.sig <- prob.pathways.sig[, , idx, drop = FALSE]
  } else {
    pathways.sig <- NULL
    prob.pathways.sig <- NULL
  }
  netP = list(pathways = pathways.sig, prob = prob.pathways.sig)


  if (do.cell) {
    if ("prob.cell" %in% names(net)) {
      prob.cell <- net$prob.cell
      prob.cell_ <- net$tmp$prob.cell # a list

      # unique
      pathways <- unique(pairLR.use$pathway_name)

      nC <- dim(prob.cell)[[1]]
      dns <- dimnames(prob.cell)
      # nrun <- length(pathways)

      cat(cli.symbol(),"Compute the communication probability between individual cells at signaling pathway level by summarizing all related ligands/receptors...\n")
      gc()
      prob.all <- pbapply::pbsapply(
        X = pathways,
        FUN = function(one_pathway) {
          # not unique, so the `prob.cell.i` may be 3 dims
          # Older Code: prob.cell.i <- prob.cell[,,pairLR.use$pathway_name == one_pathway, drop = FALSE]
          prob.cell.i <-
            prob.cell_[pairLR.use$pathway_name == one_pathway] # a named list, can receive logical vector as indexes
          prob.cell.i <-
            my_as_sparse3Darray(prob.cell.i)
          # retain dim1 & dim2, sum up dim3
          res <-
            spatstat.sparse::marginSumsSparse(prob.cell.i, MARGIN = c(1, 2))

          return(res)
        },
        simplify = F # return a list
      )
      names(prob.all) <- pathways

      # so the shape of prob.cell.pathways will be [nC,nC,length(pathways)]
      prob.cell.pathways <- my_as_sparse3Darray(prob.all)
      # cat(cli.symbol(),"Dim(prob.cell.pathways):",dim(prob.cell.pathways),"\n")


      # retain dim-3, sum up dim-1 & dim-2, `prob.sum` stores each pathway's number of cell-level links/interactions
      prob.sum <- as.vector(spatstat.sparse::marginSumsSparse(prob.cell.pathways,MARGIN = c(3)))

      # non zero index
      # After `filterCommunication`, there will be many allZeroMat, prob.sum = 0 is possible
      # `prob.sum != 0` means cell-cell signaling exists in the specific pathway
      PathwaySig.use.idx <- which(prob.sum > 0)


      # sort the `PathwaySig.use.idx` according to its corresponding prob.sum value
      idx <- sort(as.array(prob.sum[PathwaySig.use.idx]),
                  decreasing=TRUE,
                  index.return = TRUE)$ix
      PathwaySig.sort.idx <- PathwaySig.use.idx[idx]

      # Notice: len(pathways) = len(prob.sum) = dim3(prob.all)
      pathways.sig.cell <- pathways[PathwaySig.sort.idx]

      # sort `prob.cell.pathways`
      cat(cli.symbol(),"Subset the pathways with non-zero communication probability and arrange them in a decreasing order based on the total communication probabilities ...\n")
      prob.cell.pathways.sig_ <- prob.all[PathwaySig.sort.idx] # a list
      prob.cell.pathways.sig <- my_as_sparse3Darray(prob.cell.pathways.sig_)
      gc();cat(cli.symbol(),"The number of cells and pathways in Dim(prob.cell.pathways) are :",dim(prob.cell.pathways.sig),"\n")
      # Older Code:
      # pathways.sig.cell <- pathways[idx]
      # prob.cell.pathways.sig <- prob.cell.pathways[, , idx, drop = FALSE]
      # prob.cell.pathways.sig <- spatstat.sparse::as.sparse3Darray(prob.cell.pathways.sig)

      dimnames(prob.cell.pathways.sig) <- list(dns[[1]], dns[[2]], pathways.sig.cell)
      names(prob.cell.pathways.sig_) <- pathways.sig.cell
      Tmp <- list(prob.cell = prob.cell.pathways.sig_) # !important, for parallel iteration
      netP$pathways.cell = pathways.sig.cell; netP$prob.cell = prob.cell.pathways.sig
      netP$tmp <- Tmp

    }
  }
  # group-level: pathways;prob
  # individual cell-level: pathways.cell;prob.cell
  if (is.null(object)) {
    # netP = list(pathways = pathways.sig, prob = prob.pathways.sig, pathways.cell = pathways.sig.cell, prob.cell = prob.cell.pathways.sig)
    cat(cli.symbol(1),"Computing the communication probability on signaling pathway level is done. \n")
    return(netP)
  } else {
    object@netP <- netP
    cat(cli.symbol(1),"Computing the communication probability on signaling pathway level is done. \n")
    return(object)
  }
}

# computeAvgCommunProb_LR <- function (prob, group, dataLR = NULL, min.percent = 0.1, min.cells.sr = 5)
# {
#   if (!is.factor(group)) {
#     group <- factor(group)
#   }
#   level.use <- levels(group)
#   level.use <- level.use[level.use %in% unique(group)]
#   numCluster <- length(level.use)
#
#   # min.percent
#   dataLR_temp <- 1 * (dataLR > 0)
#
#   dataLR_temp <- aggregate(dataLR_temp, list(group), FUN = mean)
#   dataLR_percent <- 1 * (format(dataLR_temp[, -1], digits = 1) >= min.percent)
#   Prob_percent <- Matrix::crossprod(matrix(dataLR_percent[, 1], nrow = 1), matrix(dataLR_percent[, 2], nrow = 1))
#   melt_Prob_percent <- reshape2::melt(Prob_percent,value.name = "Binary.Prob")
#
#   Prob.avg <- matrix(0, nrow = numCluster, ncol = numCluster)
#   # min.cells.sr
#   if(NROW(melt_Prob_percent[melt_Prob_percent$Binary.Prob==1, ,drop=F]) == 0){
#     return(Prob.avg)
#   } else {
#     df_nonZero_Prob <- melt_Prob_percent[melt_Prob_percent$Binary.Prob==1, ,drop=F]
#     for (x in seq_len(NROW(df_nonZero_Prob))) {
#
#       ii <- df_nonZero_Prob[x,1,drop=T]
#       jj <- df_nonZero_Prob[x,2,drop=T]
#       prob.ij <- prob[group %in% level.use[ii], group %in% level.use[jj], drop = FALSE]
#       if ((sum(Matrix::rowSums(prob.ij > 0) > 0) >= min.cells.sr) &
#           (sum(Matrix::colSums(prob.ij > 0) > 0) >= min.cells.sr)) {
#
#         # Prob.avg[ii, jj] <- sum(prob.ij@x)/(cumprod(prob.ij@Dim)[[2]])
#         Prob.avg[ii, jj] <- sum(prob.ij@x)
#       }
#     } # forloop
#     return(Prob.avg)
#   }
# }


# Prob.avg <- matrix(0, nrow = numCluster, ncol = numCluster)
# for (ii in 1:numCluster) {
#   for (jj in 1:numCluster) {
#     Prob.temp <- as.vector(P.spatial[group %in% level.use[ii], group %in% level.use[jj]]) > 0
#     # Prob.sum <- sum(Prob.temp)
#     # Prob.avg[ii,jj] <- Prob.sum/(sum(group %in% level.use[ii]) * sum(group %in% level.use[jj]))
#     Prob.avg[ii,jj] <- sum(Prob.temp)/(sum(group %in% level.use[ii]) * sum(group %in% level.use[jj]))
#   }
# }




#' Filter cell-cell communication if there are only few number of cells in certain cell groups or only few interactions
#'
#' @param object CellChat object
#' @param min.cells the minimum number of cells required in each cell group for filtering cell group-level communication
#' @param min.links the minimum number of links/interactions required in the ligand-receptor pair for filtering individual cell-level communication
#' @param min.cells.sr the minimum number of cells required as senders or receivers for filtering individual cell-level communication
#' @return CellChat object with an updated slot net
#' @export
#'
filterCommunication <- function(object, min.cells = 10, min.links = 5, min.cells.sr = 5) {

  if (!is.null(min.cells)) {
    message("Filter cell-group level communication...",'\n')
    net <- object@net
    cell.excludes <- which(as.numeric(table(object@idents)) < min.cells)
    if (length(cell.excludes) > 0) {
      cat(cli.symbol(),"The cell-cell communication related with the following cell groups are excluded due to the few number of cells: ", levels(object@idents)[cell.excludes],'\n')
      # dim(net$prob) = nCellGroup x nCellGroup x nPairLRsig
      net$prob[cell.excludes,,] <- 0
      net$prob[,cell.excludes,] <- 0
      if (!is.null(net$pval)) {
        net$pval[net$prob == 0] <- 1
      }
      object@net <- net
    }
    rm(net)
    gc()
  }

  if ("prob.cell" %in% names(object@net)) {
    net <- object@net
    prob.cell <- net$prob.cell
    prob.cell_ <- net$tmp$prob.cell # a list

    if (!is.null(min.links) | !is.null(min.cells.sr)) {
      message("Filter individual cell-level communication...",'\n')

      # binary.prob.cell <- modifySparse3Darray(prob.cell,cutoff = 0,remain.cutoff.v = F,do.binary = T)
      # retain dim-3, sum up dim-1 && dim-2, `prob.sum` stores each LR's number of cell-level links/interactions
      # as.vector(spatstat.sparse::marginSumsSparse(binary.prob.cell, MARGIN = 3))

      prob.sum <- purrr::map_dbl(
        .x = prob.cell_,
        .f = function(Mat){
          return(length(Mat@x))
        }
      )

      dimArr <- dim(prob.cell)

      # define a allzero matrix (CsparseMatrix)
      AllzeroMat <- Matrix::sparseMatrix(
        i = integer(0),
        j = integer(0),
        x = numeric(0),
        repr = "C", # default repr in CellChat
        dims = dimArr[c(1, 2)],
        # dimnames = dns[c(1,2)] # too large, not use!
        dimnames = list(NULL,NULL),
        index1 = T # i and j are interpreted as 1-based indices, following the R convention
      )

      # filter communication according to min.links
      gc()
      if (!is.null(min.links)) {
        cat(cli.symbol(),"Filter communication according to min.links...\n")
        idx.signaling.excludes <- which((prob.sum < min.links) & (prob.sum > 0))
        if (length(idx.signaling.excludes) > 0) {
          cat("The cell-cell communication related with #", length(idx.signaling.excludes),'L-R pairs are excluded due to the few number of interactions.','\n')

          # prob.cell[,,idx.signaling.excludes] <- 0
          pb <- utils::txtProgressBar(min = 0, max = length(idx.signaling.excludes), style = 3, file = stderr(), width = 80);i=0
          for (x in idx.signaling.excludes) {
            prob.cell_[[x]] <- AllzeroMat
            utils::setTxtProgressBar(pb = pb, value = (i=i+1))
          } # forloop
          close(con = pb)

        }
      }

      # filter communication according to min.cell.sr
      gc()
      if (!is.null(min.cells.sr)) {
        cat(cli.symbol(),"Filter communication according to min.cell.sr... \n")
        pathways0 <- dimnames(prob.cell)[[3]] # L-R pairs' names
        if (is.null(min.links)) {
          pathways <- pathways0[prob.sum > 0]
        } else {
          pathways <- pathways0[prob.sum >= max(1, min.links)] # Prevent `min.links` from being smaller than 1
        }

        if (length(pathways)<1){
          NULL # not filter
        } else {
          # Check:
          # > SparseLogMat <- matrix(c(T,F,F,F,F,T),2,3)
          # > SparseLogMat
          # [,1]  [,2]  [,3]
          # [1,]  TRUE FALSE FALSE
          # [2,] FALSE FALSE  TRUE
          # > SparseLogMat <- as(SparseLogMat,"CsparseMatrix")
          # > SparseLogMat
          # 2 x 3 sparse Matrix of class "lgCMatrix"
          #
          # [1,] | . .
          # [2,] . . |
          #   > rowSums(SparseLogMat)
          # [1] 1 1
          # > colSums(SparseLogMat)
          # [1] 1 0 1

          pathways.remove <- pbapply::pblapply(
            X = seq_len(length(pathways)),
            FUN = function(x) {
              prob.cell.i <- prob.cell_[[ pathways[[x]] ]] > 0 # `prob.cell.i` is a logical sparse matrix
              if ((sum(Matrix::rowSums(prob.cell.i) > 0) < min.cells.sr) | (sum(Matrix::colSums(prob.cell.i) > 0) < min.cells.sr)) {
                gc()
                return(pathways[[x]])
              }
            }
            # simplify = T,
            # # pathways.remove is a vec
            # hint.message = "Filtering...",
            # future.label = "pathways.remove_my_future_sapply-%d"
          )
          pathways.remove <- unlist(pathways.remove) # vec2vec
          pathways.remove.idx <- which(pathways0 %in% pathways.remove)

          if (length(pathways.remove) > 0) {
            cat(
              "The cell-cell communication related with #",
              length(pathways.remove),
              'L-R pairs are excluded due to the few number of sending/receiving cells.',
              '\n'
            )

            # prob.cell[, , pathways.remove.idx] <- 0
            pb <- utils::txtProgressBar(min = 0, max = length(pathways.remove.idx), style = 3, file = stderr(),width = 80);i=0
            for (x in pathways.remove.idx) {
              prob.cell_[[x]] <- AllzeroMat
              utils::setTxtProgressBar(pb = pb, value = (i=i+1))
            } # forloop
            close(con = pb)

            # pbapply::pblapply(
            #   X = pathways.remove.idx,
            #   FUN = function(x) {
            #     prob.cell_[[x]] <- AllzeroMat
            #   }
            # )
          }
        }
      }

      # update the obj
      net$tmp$prob.cell <- prob.cell_ # a list
      dns <- dimnames(prob.cell) # dimnames
      net$prob.cell <- my_as_sparse3Darray(prob.cell_)
      dimnames(net$prob.cell) <- dns
      object@net <- net
      cat(paste0(cli.symbol(1),'Filtering cell-cell communication is done.<<< [', Sys.time(),']', '\n'))

    } # !is.null(min.links) | !is.null(min.cells.sr)
  } else {
    stop(
      cli.symbol(2),
      "Please run `computeCommunProb` to compute the communication probability/strength between any interacting individual cells! "
    )
  }
  return(object)
}

#' @title modifySparse3Darray
#'
#' @param array3d sparse3Darray. See details in \code{\link[spatstat.sparse]{sparse3Darray}}
#' @param do.binary Boolean. Use 1 to replace the non-zero values in the 3Darray.
#' @param cutoff Numeric. Cut off the sparse3Darray's values to make
#' array' values >= cutoff.
#' @param remain.cutoff.v Boolean. Use `>=` or `>`.
#'
#' @return sparse3Darray.
#' @export
modifySparse3Darray <- function(
    array3d,
    do.binary=T,
    cutoff=NULL,
    remain.cutoff.v=T
){
  if (!inherits(array3d, "sparse3Darray")) {
    stop("Please use sparse3Darray as input!")
  }

  if((!is.null(cutoff)) & is.numeric(cutoff)){
    if(length(array3d$x) > 0){

      if(remain.cutoff.v){
        remain.idx <- which(array3d$x >= cutoff)
      } else {
        remain.idx <- which(array3d$x > cutoff)
      }

      # update
      array3d$i <- array3d$i[remain.idx]
      array3d$j <- array3d$j[remain.idx]
      array3d$k <- array3d$k[remain.idx]
      array3d$x <- array3d$x[remain.idx]
    }
  }

  if(do.binary){
    array3d$x <- rep.int(1,times = length(array3d$x))
  }

  return(array3d)
}


#' @title scMatrixTruncation
#' @description
#' Single cell data is usually saved as "dgCMatrix" or "dgTMatrix" in R,
#' use this function to cut off Expression Matrix's values to make
#' Mat'v >= cutoff.
#'
#' @param Mat dgCMatrix or dgTMatrix
#' @param cutoff Numeric. Cut off Expression Matrix's values to make
#' Mat' values >= cutoff.
#' @param remain.cutoff.v Boolean. Use `>=` or `>`.
#' @param repr Character string. One of "C", "T", specifying the representation of the sparse matrix result.
#'
#' @return dgCMatrix or dgTMatrix
#' @export
scMatrixTruncation <- function(Mat,cutoff=NULL,remain.cutoff.v=T,repr = c("C", "T")){
  if (!inherits(x = Mat, what = c("dgCMatrix","dgTMatrix"))) {
    stop("Please use `dgCMatrix` or `dgTMatrix` as a Matrix input!")
  }

  if((!is.null(cutoff)) & is.numeric(cutoff)){
    if(length(Mat@x) > 0){
      Mat <- as(Mat,Class = "TsparseMatrix")

      if(remain.cutoff.v){
        remain.idx <- which(Mat@x >= cutoff)
      } else{
        remain.idx <- which(Mat@x > cutoff)
      }

      # update
      Mat@i <- Mat@i[remain.idx]
      Mat@j <- Mat@j[remain.idx]
      Mat@x <- Mat@x[remain.idx]
    }
  }

  repr <- match.arg(repr)
  if(repr=="C"){
    Mat <- as(Mat,Class = "CsparseMatrix")
  } else if (repr=="T"){
    Mat <- as(Mat,Class = "TsparseMatrix")
  }
  return(Mat)
}



#' @title filterProbability
#' @description
#' Filter out statistically non-significant communication probability at the level of individual cells after running \link{computeCommunProb}
#'
#' @param object CellChat object
#' @param nboot Numeric. The number of bootstrap samples, 100 by default.
#' @param seed.use Integer. The random seed used when taking a sample from the communication probabilities
#' @param thresh Numeric. The threshold for defining significant individual cell-cell communication at (1-thresh) of a shuffled distribution of each L-R pair.
#'
#' @return CellChat object
#' @export
filterProbability <- function (
    object,
    nboot = 100,
    seed.use = 666L,
    thresh = 0.05
){
  quantile.prob <- 1-thresh
  if(quantile.prob==0){
    cat(cli.symbol(1),"Do not filter any CCC probability!")
    return(object)
  } else { # quantile.prob<1
    if (is.null(object@net$tmp$prob.cell)) {
      stop(
        cli.symbol(2),
        "Please run `computeCommunProb` to compute the communication probability/strength between any interacting individual cells! "
      )
    } else {
      cat(paste0(cli.symbol(), "Filter out non-significant communication with a probability quantile being ",
                 quantile.prob, " for each L-R pair... \n"))
      prob.cell_ <- object@net$tmp$prob.cell
      pair.LR.use <- names(prob.cell_)
      cell.names <- spatstat.sparse::dimnames.sparse3Darray(object@net$prob.cell)[[1]]
    }

    d.spatial <- object@images$result.computeCellDistance$d.spatial
    Matrix::diag(d.spatial) <- 1
    adj.contact <- object@images$result.computeCellDistance$adj.contact

    nLR <- object@options$parameter$nLR
    nLR1 <- object@options$parameter$nLR1
    nC <- NROW(d.spatial)

    set.seed(seed.use)
    if (object@options$parameter[["all.contact.dependent"]] == TRUE) {
      d.spatial <- adj.contact
    }

    # dim(permutation) = nboot x nLR
    permutation <- replicate(nLR, base::sample(x = 1:nC, size = nboot, replace = F))

    prob.cell_ <- my_future_lapply(X = 1:nLR, FUN = function(i) {

      if (i <= nLR1) {
        d_spatial <- d.spatial
      } else {
        d_spatial <- adj.contact
      }
      sample.cells <- permutation[ ,i,drop=T]
      Prob.cell.i <- prob.cell_[[i]]

      sample.prob.cell.i <- purrr::map(.x = sample.cells,
                                       .f = function(cell.index) {
                                         prob.index <- which(d_spatial[cell.index, ,drop = T] > 0)
                                         nboot.prob.cell.i <- Prob.cell.i[cell.index,prob.index, drop = T] # get a dense vec
                                         return(nboot.prob.cell.i)
                                       }) %>% unlist()
      nboot.quantile <- quantile(sample.prob.cell.i, probs = quantile.prob)
      gc()

      if(nboot.quantile == 0){
        # sparse enough, return directly
        return(Prob.cell.i)
      } else if (nboot.quantile > 0) {
        # filter `Prob.cell.i` to make it sparse enough
        Prob.cell.i <- scMatrixTruncation(Prob.cell.i,cutoff = nboot.quantile,remain.cutoff.v = T,repr = "C")
      }

      return(Prob.cell.i)
    }, simplify = F, hint.message = "filtering...")
    prob.cell <- my_as_sparse3Darray(prob.cell_)
    dimnames(prob.cell) <- list(cell.names, cell.names, pair.LR.use)
    names(prob.cell_) <- pair.LR.use
    object@net$prob.cell <- prob.cell
    object@net$tmp$prob.cell <- prob.cell_

    cat(cli.symbol(1), "Filtering is done.\n")
    return(object)
  } # whether to filter out
}


#' Calculate the aggregated network by counting the number of links or summarizing the communication probability
#'
#' @param object CellChat object
#' @param sources.use,targets.use,signaling,pairLR.use Please check the description in function \code{\link{subsetCommunication}}
#' @param remove.isolate whether removing the isolate cell groups without any interactions when applying \code{\link{subsetCommunication}}
#' @param thresh threshold of the p-value for determining significant interaction
#' @param return.object whether return an updated CellChat object
#' @importFrom  dplyr group_by summarize groups
#' @importFrom stringr str_split
#'
#' @return Return an updated CellChat object:
#'
#' `object@net$count` is a matrix: rows and columns are sources and targets respectively, and elements are the number of interactions between any two cell groups. USER can convert a matrix to a data frame using the function `reshape2::melt()`
#'
#' `object@net$weight` is also a matrix containing the interaction weights between any two cell groups
#'
#' `object@net$sum` is deprecated. Use `object@net$weight`
#'
#' @export
#'
aggregateNet <- function(object, sources.use = NULL, targets.use = NULL, signaling = NULL, pairLR.use = NULL, remove.isolate = TRUE, thresh = 0.05, return.object = TRUE) {
  net <- object@net
  if (is.null(sources.use) & is.null(targets.use) & is.null(signaling) & is.null(pairLR.use)) {
    prob <- net$prob
    pval <- net$pval
    pval[prob == 0] <- 1
    prob[pval >= thresh] <- 0
    net$count <- apply(prob > 0, c(1,2), sum)
    net$weight <- apply(prob, c(1,2), sum)
    net$weight[is.na(net$weight)] <- 0
    net$count[is.na(net$count)] <- 0
    net$LR.sig <- dimnames(prob)[[3]][apply(prob, 3, sum) > 0]
  } else {
    df.net <- subsetCommunication(object, slot.name = "net",
                                  sources.use = sources.use, targets.use = targets.use,
                                  signaling = signaling,
                                  pairLR.use = pairLR.use,
                                  thresh = thresh)
    df.net$source_target <- paste(df.net$source, df.net$target, sep = "_")
    df.net2 <- df.net %>% group_by(source_target) %>% summarize(count = n(), .groups = 'drop')
    df.net3 <- df.net %>% group_by(source_target) %>% summarize(prob = sum(prob), .groups = 'drop')
    df.net2$prob <- df.net3$prob
    a <- stringr::str_split(df.net2$source_target, "_", simplify = T)
    df.net2$source <- as.character(a[, 1])
    df.net2$target <- as.character(a[, 2])
    cells.level <- levels(object@idents)
    if (remove.isolate) {
      message("Isolate cell groups without any interactions are removed. To block it, set `remove.isolate = FALSE`")
      df.net2$source <- factor(df.net2$source, levels = cells.level[cells.level %in% unique(df.net2$source)])
      df.net2$target <- factor(df.net2$target, levels = cells.level[cells.level %in% unique(df.net2$target)])
    } else {
      df.net2$source <- factor(df.net2$source, levels = cells.level)
      df.net2$target <- factor(df.net2$target, levels = cells.level)
    }

    count <- tapply(df.net2[["count"]], list(df.net2[["source"]], df.net2[["target"]]), sum)
    prob <- tapply(df.net2[["prob"]], list(df.net2[["source"]], df.net2[["target"]]), sum)
    net$count <- count
    net$weight <- prob
    net$weight[is.na(net$weight)] <- 0
    net$count[is.na(net$count)] <- 0
  }

  if ("prob.cell" %in% names(object@net)) {
    prob.cell <- net$prob.cell

    # net$count.cell <-  as(apply(prob.cell > 0, c(1,2), sum), "dgCMatrix")
    # prob.cell > 0 generate a sparse logical array, but cannot sum directly
    prob.cell.positive <- prob.cell > 0
    prob.cell.positive$x <- rep.int(1,length(prob.cell.positive$x))
    net$count.cell <-  spatstat.sparse::marginSumsSparse(prob.cell.positive,MARGIN = c(1,2))

    # net$weight.cell <- as(apply(prob.cell, c(1,2), sum), "dgCMatrix")
    net$weight.cell <- spatstat.sparse::marginSumsSparse(prob.cell,MARGIN = c(1,2))
    prob.cell.sum <- spatstat.sparse::marginSumsSparse(prob.cell,MARGIN = c(3))
    net$LR.sig.cell <- dimnames(prob.cell)[[3]][prob.cell.sum@i]
    if (!is.null(sources.use) | !is.null(targets.use) | !is.null(signaling) | !is.null(pairLR.use)) {
      message("Subsetting cells or signaling is not applicable to individual cell-based `prob.cell`!", '\n')
    }
  }


  if (return.object) {
    object@net <- net
    # object@net$tmp <- NULL
    return(object)
  } else {
    return(net)
  }
}




#' Compute averaged expression values for each cell group
#'
#' @param object CellChat object
#' @param features a char vector giving the used features. default use all features
#' @param group.by cell group information; default is `object@idents` when input is a single object and `object@idents$joint` when input is a merged object; otherwise it should be one of the column names of the meta slot
#' @param type methods for computing the average gene expression per cell group.
#'
#' By default = "triMean", defined as a weighted average of the distribution's median and its two quartiles (https://en.wikipedia.org/wiki/Trimean);
#'
#' When setting `type = "truncatedMean"`, a value should be assigned to 'trim'. See the function \code{\link[base]{mean}}.
#'
#' @param trim the fraction (0 to 0.25) of observations to be trimmed from each end of x before the mean is computed.
#' @param slot.name the data in the slot.name to use
#' @param data.use a customed data matrix. Default: data.use = NULL and the expression matrix in the 'slot.name' is used
#'
#' @return Returns a matrix with genes as rows, cell groups as columns.

#' @export
#'
computeAveExpr <- function(object, features = NULL, group.by = NULL, type = c("triMean", "truncatedMean", "median"), trim = NULL,
                           slot.name = c("data.signaling", "data"), data.use = NULL) {
  type <- match.arg(type)
  slot.name <- match.arg(slot.name)
  FunMean <- switch(type,
                    triMean = triMean,
                    truncatedMean = function(x) mean(x, trim = trim, na.rm = TRUE),
                    median = function(x) median(x, na.rm = TRUE))

  if (is.null(data.use)) {
    data.use <- slot(object, slot.name)
  }
  if (is.null(features)) {
    features.use <- row.names(data.use)
  } else {
    features.use <- intersect(features, row.names(data.use))
  }
  data.use <- data.use[features.use, , drop = FALSE]
  data.use <- as.matrix(data.use)

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
  # compute the average expression per group
  data.use.avg <- aggregate(t(data.use), list(labels), FUN = FunMean)
  data.use.avg <- t(data.use.avg[,-1])
  rownames(data.use.avg) <- features.use
  colnames(data.use.avg) <- levels(labels)
  return(data.use.avg)
}



#' Compute the expression of complex in individual cells using geometric mean
#' @param complex_input the complex_input from CellChatDB
#' @param data.use data matrix (row are genes and columns are cells or cell groups)
#' @param complex the names of complex
#' @return
#' @importFrom dplyr select starts_with
#' @importFrom future nbrOfWorkers
#' @importFrom future.apply future_sapply
#' @importFrom pbapply pbsapply
#' @export
computeExpr_complex <- function(complex_input, data.use, complex) {
  Rsubunits <- complex_input[complex,] %>% dplyr::select(starts_with("subunit"))

  data.complex = my_future_sapply(
    X = 1:nrow(Rsubunits),
    FUN = function(x) {
      RsubunitsV <- unlist(Rsubunits[x,], use.names = F)
      RsubunitsV <- RsubunitsV[RsubunitsV != ""]
      return(geometricMean(data.use[RsubunitsV,,drop=F]))
    },
    hint.message = "Computing complex...",
    strategy.message = FALSE
  )
  data.complex <- t(data.complex)
  return(data.complex)
}

# Compute the average expression of complex per cell group using geometric mean
# @param complex_input the complex_input from CellChatDB
# @param data.use data matrix (rows are genes and columns are cells)
# @param complex the names of complex
# @param group a factor defining the cell groups
# @param FunMean the function for computing mean expression per group
# @return
# @importFrom dplyr select starts_with
# @importFrom future nbrOfWorkers
# @importFrom future.apply future_sapply
# @importFrom pbapply pbsapply
# #' @export
.computeExprGroup_complex <- function(complex_input, data.use, complex, group, FunMean) {
  Rsubunits <- complex_input[complex,] %>% dplyr::select(starts_with("subunit"))
  my.sapply <- ifelse(
    test = future::nbrOfWorkers() == 1,
    yes = sapply,
    no = future.apply::future_sapply
  )
  data.complex = my.sapply(
    X = 1:nrow(Rsubunits),
    FUN = function(x) {
      RsubunitsV <- unlist(Rsubunits[x,], use.names = F)
      RsubunitsV <- RsubunitsV[RsubunitsV != ""]
      RsubunitsV <- intersect(RsubunitsV, rownames(data.use))
      if (length(RsubunitsV) > 1) {
        data.avg <- aggregate(t(data.use[RsubunitsV,]), list(group), FUN = FunMean)
        data.avg <- t(data.avg[,-1])
      } else if (length(RsubunitsV) == 1) {
        data.avg <- aggregate(matrix(data.use[RsubunitsV,], ncol = 1), list(group), FUN = FunMean)
        data.avg <- t(data.avg[,-1])
      } else {
        data.avg = matrix(0, nrow = 1, ncol = length(unique(group)))
      }
      return(geometricMean(data.avg))
    }
  )
  data.complex <- t(data.complex)
  return(data.complex)
}

#' Compute the expression of ligands or receptors using geometric mean
#' @param geneLR a char vector giving a set of ligands or receptors
#' @param data.use data matrix (row are genes and columns are cells or cell groups)
#' @param complex_input the complex_input from CellChatDB
# #' @param group a factor defining the cell groups; If NULL, compute the expression of ligands or receptors in individual cells; otherwise, compute the average expression of ligands or receptors per cell group
# #' @param FunMean the function for computing average expression per cell group
#' @return
#' @export
computeExpr_LR <- function(geneLR, data.use, complex_input){
  nLR <- length(geneLR)
  numCluster <- ncol(data.use)
  index.singleL <- which(geneLR %in% rownames(data.use))
  dataL1avg <- data.use[geneLR[index.singleL],]
  dataLavg <- matrix(nrow = nLR, ncol = numCluster)
  dataLavg[index.singleL,] <- dataL1avg
  index.complexL <- setdiff(1:nLR, index.singleL)
  if (length(index.complexL) > 0) {
    complex <- geneLR[index.complexL]
    data.complex <- computeExpr_complex(complex_input, data.use, complex)
    dataLavg[index.complexL,] <- data.complex
  }
  return(dataLavg)
}


#' Modeling the effect of coreceptor on the ligand-receptor interaction
#'
#' @param data.use data matrix
#' @param cofactor_input the cofactor_input from CellChatDB
#' @param pairLRsig a data frame giving ligand-receptor interactions
#' @param type when type == "A", computing expression of co-activation receptor; when type == "I", computing expression of co-inhibition receptor.
#' @return
#' @importFrom future nbrOfWorkers
#' @importFrom future.apply future_sapply
#' @importFrom pbapply pbsapply
#' @export
computeExpr_coreceptor <- function(cofactor_input, data.use, pairLRsig, type = c("A", "I")) {
  type <- match.arg(type)
  if (type == "A") {
    coreceptor.all = pairLRsig$co_A_receptor
  } else if (type == "I"){
    coreceptor.all = pairLRsig$co_I_receptor
  }
  index.coreceptor <- which(!is.na(coreceptor.all) & coreceptor.all != "")
  if (length(index.coreceptor) > 0) {

    coreceptor <- coreceptor.all[index.coreceptor]
    coreceptor.ind <- cofactor_input[coreceptor, grepl("cofactor" , colnames(cofactor_input) )]
    data.coreceptor.ind = my_future_sapply(
      X = 1:nrow(coreceptor.ind),
      FUN = function(x) {
        coreceptor.indV <- unlist(coreceptor.ind[x,], use.names = F)
        coreceptor.indV <- coreceptor.indV[coreceptor.indV != ""]
        coreceptor.indV <- intersect(coreceptor.indV, rownames(data.use))
        if (length(coreceptor.indV) == 1) {
          return(1 + data.use[coreceptor.indV, ])
        } else if (length(coreceptor.indV) > 1) {
          return(apply(1 + data.use[coreceptor.indV, ], 2, prod))
        } else {
          return(matrix(1, nrow = 1, ncol = ncol(data.use)))
        }
      },
      hint.message = "Computing coreceptor...",
      strategy.message = FALSE
    )
    data.coreceptor.ind <- t(data.coreceptor.ind)
    data.coreceptor <- matrix(1, nrow = length(coreceptor.all), ncol = ncol(data.use))
    data.coreceptor[index.coreceptor,] <- data.coreceptor.ind
  } else {
    data.coreceptor <- matrix(1, nrow = length(coreceptor.all), ncol = ncol(data.use))
  }
  return(data.coreceptor)
}

# Modeling the effect of coreceptor on the ligand-receptor interaction
#
# @param data.use data matrix
# @param cofactor_input the cofactor_input from CellChatDB
# @param pairLRsig a data frame giving ligand-receptor interactions
# @param type when type == "A", computing expression of co-activation receptor; when type == "I", computing expression of co-inhibition receptor.
# @param group a factor defining the cell groups
# @param FunMean the function for computing mean expression per group
# @return
# @importFrom future nbrOfWorkers
# @importFrom future.apply future_sapply
# @importFrom pbapply pbsapply
# #' @export
.computeExprGroup_coreceptor <- function(cofactor_input, data.use, pairLRsig, type = c("A", "I"), group, FunMean) {
  type <- match.arg(type)
  if (type == "A") {
    coreceptor.all = pairLRsig$co_A_receptor
  } else if (type == "I"){
    coreceptor.all = pairLRsig$co_I_receptor
  }
  index.coreceptor <- which(!is.na(coreceptor.all) & coreceptor.all != "")
  if (length(index.coreceptor) > 0) {
    my.sapply <- ifelse(
      test = future::nbrOfWorkers() == 1,
      yes = pbapply::pbsapply,
      no = future.apply::future_sapply
    )
    coreceptor <- coreceptor.all[index.coreceptor]
    coreceptor.ind <- cofactor_input[coreceptor, grepl("cofactor" , colnames(cofactor_input) )]
    data.coreceptor.ind = my.sapply(
      X = 1:nrow(coreceptor.ind),
      FUN = function(x) {
        coreceptor.indV <- unlist(coreceptor.ind[x,], use.names = F)
        coreceptor.indV <- coreceptor.indV[coreceptor.indV != ""]
        coreceptor.indV <- intersect(coreceptor.indV, rownames(data.use))
        if (length(coreceptor.indV) > 1) {
          data.avg <- aggregate(t(data.use[coreceptor.indV,]), list(group), FUN = FunMean)
          data.avg <- t(data.avg[,-1])
          return(apply(1 + data.avg, 2, prod))
          # return(1 + apply(data.avg, 2, mean))
        } else if (length(coreceptor.indV) == 1) {
          data.avg <- aggregate(matrix(data.use[coreceptor.indV,], ncol = 1), list(group), FUN = FunMean)
          data.avg <- t(data.avg[,-1])
          return(1 + data.avg)
        } else {
          return(matrix(1, nrow = 1, ncol = length(unique(group))))
        }
      }
    )
    data.coreceptor.ind <- t(data.coreceptor.ind)
    data.coreceptor <- matrix(1, nrow = length(coreceptor.all), ncol = length(unique(group)))
    data.coreceptor[index.coreceptor,] <- data.coreceptor.ind
  } else {
    data.coreceptor <- matrix(1, nrow = length(coreceptor.all), ncol = length(unique(group)))
  }

  return(data.coreceptor)
}

#' Modeling the effect of agonist on the ligand-receptor interaction
#' @param data.use data matrix
#' @param cofactor_input the cofactor_input from CellChatDB
#' @param pairLRsig the L-R interactions
#' @param group a factor defining the cell groups
#' @param index.agonist the index of agonist in the database
#' @param Kh a parameter in Hill function
#' @param FunMean the function for computing mean expression per group
#' @param n Hill coefficient
#' @return
#' @export
#' @importFrom stats aggregate
computeExprGroup_agonist <- function(data.use, pairLRsig, cofactor_input, group, index.agonist, Kh, FunMean, n) {
  agonist <- pairLRsig$agonist[index.agonist]
  agonist.ind <- cofactor_input[agonist, grepl("cofactor" , colnames(cofactor_input))]
  agonist.indV <- unlist(agonist.ind, use.names = F)
  agonist.indV <- agonist.indV[agonist.indV != ""]
  agonist.indV <- intersect(agonist.indV, rownames(data.use))
  if (length(agonist.indV) == 1) {
    data.avg <- aggregate(matrix(data.use[agonist.indV,], ncol = 1), list(group), FUN = FunMean)
    data.avg <- t(data.avg[,-1])
    data.agonist <- 1 + data.avg^n/(Kh^n + data.avg^n)
  } else if (length(agonist.indV) > 1) {
    data.avg <- aggregate(t(data.use[agonist.indV,]), list(group), FUN = FunMean)
    data.avg <- t(data.avg[,-1])
    data.agonist <- apply(1 + data.avg^n/(Kh^n + data.avg^n), 2, prod)
  } else {
    data.agonist = matrix(1, nrow = 1, ncol = length(unique(group)))
  }
  return(data.agonist)
}

#' Modeling the effect of antagonist on the ligand-receptor interaction
#'
#' @param data.use data matrix
#' @param cofactor_input the cofactor_input from CellChatDB
#' @param pairLRsig the L-R interactions
#' @param group a factor defining the cell groups
#' @param index.antagonist the index of antagonist in the database
#' @param Kh a parameter in Hill function
#' @param n Hill coefficient
#' @param FunMean the function for computing mean expression per group
#' @return
#' @export
#' @importFrom stats aggregate
computeExprGroup_antagonist <- function(data.use, pairLRsig, cofactor_input, group, index.antagonist, Kh, FunMean, n) {
  antagonist <- pairLRsig$antagonist[index.antagonist]
  antagonist.ind <- cofactor_input[antagonist, grepl( "cofactor" , colnames(cofactor_input) )]
  antagonist.indV <- unlist(antagonist.ind, use.names = F)
  antagonist.indV <- antagonist.indV[antagonist.indV != ""]
  antagonist.indV <- intersect(antagonist.indV, rownames(data.use))
  if (length(antagonist.indV) == 1) {
    data.avg <- aggregate(matrix(data.use[antagonist.indV,], ncol = 1), list(group), FUN = FunMean)
    data.avg <- t(data.avg[,-1])
    data.antagonist <- Kh^n/(Kh^n + data.avg^n)
  } else if (length(antagonist.indV) > 1) {
    data.avg <- aggregate(t(data.use[antagonist.indV,]), list(group), FUN = FunMean)
    data.avg <- t(data.avg[,-1])
    data.antagonist <- apply(Kh^n/(Kh^n + data.avg^n), 2, prod)
  } else {
    data.antagonist = matrix(1, nrow = 1, ncol = length(unique(group)))
  }
  return(data.antagonist)
}


#' Modeling the effect of agonist on the ligand-receptor interaction
#' @param data.use data matrix
#' @param cofactor_input the cofactor_input from CellChatDB
#' @param pairLRsig the L-R interactions
# #' @param group a factor defining the cell groups
#' @param index.agonist the index of agonist in the database
#' @param Kh a parameter in Hill function
# #' @param FunMean the function for computing mean expression per group
#' @param n Hill coefficient
#' @return
#' @export
#' @importFrom stats aggregate
computeExpr_agonist <- function(data.use, pairLRsig, cofactor_input, index.agonist, Kh,  n) {
  agonist <- pairLRsig$agonist[index.agonist]
  agonist.ind <- cofactor_input[agonist, grepl("cofactor" , colnames(cofactor_input))]
  agonist.indV <- unlist(agonist.ind, use.names = F)
  agonist.indV <- agonist.indV[agonist.indV != ""]
  agonist.indV <- intersect(agonist.indV, rownames(data.use))
  if (length(agonist.indV) == 1) {
    # data.avg <- aggregate(matrix(data.use[agonist.indV,], ncol = 1), list(group), FUN = FunMean)
    # data.avg <- t(data.avg[,-1])
    data.avg <- data.use[agonist.indV,, drop = FALSE]
    data.agonist <- 1 + data.avg^n/(Kh^n + data.avg^n)
  } else if (length(agonist.indV) > 1) {
    # data.avg <- aggregate(t(data.use[agonist.indV,]), list(group), FUN = FunMean)
    # data.avg <- t(data.avg[,-1])
    data.avg <- data.use[agonist.indV,, drop = FALSE]
    data.agonist <- apply(1 + data.avg^n/(Kh^n + data.avg^n), 2, prod)
  } else {
    # data.agonist = matrix(1, nrow = 1, ncol = length(unique(group)))
    data.agonist = matrix(1, nrow = 1, ncol = ncol(data.use))
  }
  return(data.agonist)
}

#' Modeling the effect of antagonist on the ligand-receptor interaction
#'
#' @param data.use data matrix
#' @param cofactor_input the cofactor_input from CellChatDB
#' @param pairLRsig the L-R interactions
# #' @param group a factor defining the cell groups
#' @param index.antagonist the index of antagonist in the database
#' @param Kh a parameter in Hill function
#' @param n Hill coefficient
# #' @param FunMean the function for computing mean expression per group
#' @return
#' @export
#' @importFrom stats aggregate
computeExpr_antagonist <- function(data.use, pairLRsig, cofactor_input, index.antagonist, Kh, n) {
  antagonist <- pairLRsig$antagonist[index.antagonist]
  antagonist.ind <- cofactor_input[antagonist, grepl( "cofactor" , colnames(cofactor_input) )]
  antagonist.indV <- unlist(antagonist.ind, use.names = F)
  antagonist.indV <- antagonist.indV[antagonist.indV != ""]
  antagonist.indV <- intersect(antagonist.indV, rownames(data.use))
  if (length(antagonist.indV) == 1) {
    # data.avg <- aggregate(matrix(data.use[antagonist.indV,], ncol = 1), list(group), FUN = FunMean)
    # data.avg <- t(data.avg[,-1])
    data.avg <- data.use[antagonist.indV,, drop = FALSE]
    data.antagonist <- Kh^n/(Kh^n + data.avg^n)
  } else if (length(antagonist.indV) > 1) {
    # data.avg <- aggregate(t(data.use[antagonist.indV,]), list(group), FUN = FunMean)
    # data.avg <- t(data.avg[,-1])
    data.avg <- data.use[antagonist.indV,, drop = FALSE]
    data.antagonist <- apply(Kh^n/(Kh^n + data.avg^n), 2, prod)
  } else {
    # data.antagonist = matrix(1, nrow = 1, ncol = length(unique(group)))
    data.antagonist = matrix(1, nrow = 1, ncol = ncol(data.use))
  }
  return(data.antagonist)
}


#' Compute the geometric mean
#' @param x a numeric vector
#' @param na.rm whether remove na
#' @return Numeric.
#' @export
geometricMean <- function(x,na.rm=TRUE){
  if (is.null(nrow(x))) {
    exp(mean(log(x),na.rm=na.rm))
  } else {
    exp(apply(log(x),2L,mean,na.rm=na.rm))
  }
}


# arithmeticMean <- function(x,na.rm=TRUE){
#   if (is.null(nrow(x))) {
#     BiocGenerics::mean(x,na.rm = na.rm)
#   } else {
#     apply(x,2, BiocGenerics::mean ,na.rm=na.rm)
#   }
# }


#' Compute the Tukey's trimean
#' @param x a numeric vector
#' @param na.rm whether remove na
#' @return
#' @importFrom stats quantile
#' @export
triMean <- function(x, na.rm = TRUE) {
  mean(stats::quantile(x, probs = c(0.25, 0.50, 0.50, 0.75), na.rm = na.rm))
}

#' Compute the average expression per cell group when the percent of expressing cells per cell group larger than a threshold
#' @param x a numeric vector
#' @param trim the percent of expressing cells per cell group to be considered as zero
#' @param na.rm whether remove na
#' @return
#' @importFrom Matrix nnzero
#' @export
thresholdedMean <- function(x, trim = 0.1, na.rm = TRUE) {
  percent <- Matrix::nnzero(x)/length(x)
  if (percent < trim) {
    return(0)
  } else {
    return(mean(x, na.rm = na.rm))
  }
}

#' Identify all the significant interactions (L-R pairs) from some cell groups to other cell groups
#'
#' @param object CellChat object
#' @param from a vector giving the index or the name of source cell groups
#' @param to a corresponding vector giving the index or the name of target cell groups. Note: The length of 'from' and 'to' must be the same, giving the corresponding pair of cell groups for communication.
#' @param bidirection whether show the bidirectional communication, i.e., both 'from'->'to' and 'to'->'from'.
#' @param pair.only whether only return ligand-receptor pairs without pathway names and communication strength
#' @param pairLR.use0 ligand-receptor pairs to use; default is all the significant interactions
#' @param thresh threshold of the p-value for determining significant interaction
#'
#' @return
#' @export
#'
identifyEnrichedInteractions <- function(object, from, to, bidirection = FALSE, pair.only = TRUE, pairLR.use0 = NULL, thresh = 0.05){
  pairwiseLR <- object@net$pairwiseRank
  if (is.null(pairwiseLR)) {
    stop("The interactions between pairwise cell groups have not been extracted!
         Please first run `object <- rankNetPairwise(object)`")
  }
  group.names.all <- names(pairwiseLR)
  if (!is.numeric(from)) {
    from <- match(from, group.names.all)
    if (sum(is.na(from)) > 0) {
      message("Some input cell group names in 'from' do not exist!")
      from <- from[!is.na(from)]
    }
  }
  if (!is.numeric(to)) {
    to <- match(to, group.names.all)
    if (sum(is.na(to)) > 0) {
      message("Some input cell group names in 'to' do not exist!")
      to <- to[!is.na(to)]
    }
  }
  if (length(from) != length(to)) {
    stop("The length of 'from' and 'to' must be the same!")
  }
  if (bidirection) {
    from2 <- c(from, to)
    to <- c(to, from)
    from <- from2
  }
  if (is.null(pairLR.use0)) {
    k <- 0
    pairLR.use0 <- list()
    for (i in 1:length(from)){
      pairwiseLR_ij <- pairwiseLR[[from[i]]][[to[i]]]
      idx <- pairwiseLR_ij$pval < thresh
      if (length(idx) > 0) {
        k <- k +1
        pairLR.use0[[k]] <- pairwiseLR_ij[idx,]
      }
    }
    pairLR.use0 <- do.call(rbind, pairLR.use0)
  }

  k <- 0
  pval <- matrix(nrow = length(rownames(pairLR.use0)), ncol = length(from))
  prob <- pval
  group.names <- c()
  for (i in 1:length(from)) {
    k <- k+1
    pairwiseLR_ij <- pairwiseLR[[from[i]]][[to[i]]]
    pairwiseLR_ij <- pairwiseLR_ij[rownames(pairLR.use0),]
    pval_ij <- pairwiseLR_ij$pval
    prob_ij <- pairwiseLR_ij$prob
    pval_ij[pval_ij > 0.05] = 1
    pval_ij[pval_ij > 0.01 & pval_ij <= 0.05] = 2
    pval_ij[pval_ij <= 0.01] = 3
    prob_ij[pval_ij ==1] <- 0
    pval[,k] <- pval_ij
    prob[,k] <- prob_ij
    group.names <- c(group.names, paste(group.names.all[from[i]], group.names.all[to[i]], sep = " - "))
  }
  prob[which(prob == 0)] <- NA
  # remove rows that are entirely NA
  pval <- pval[rowSums(is.na(prob)) != ncol(prob), ,drop = FALSE]
  pairLR.use0 <- pairLR.use0[rowSums(is.na(prob)) != ncol(prob), ,drop = FALSE]
  prob <- prob[rowSums(is.na(prob)) != ncol(prob), ,drop = FALSE]
  if (pair.only) {
    pairLR.use0 <- dplyr::select(pairLR.use0, ligand, receptor)
  }
  return(pairLR.use0)
}


#' Compute the region distance based on the spatial locations of each splot/cell of the spatial transcriptomics
#'
#' @param coordinates a data matrix in which each row gives the spatial locations/coordinates of each cell/spot
#' @param group a factor vector defining the regions/labels of each cell/spot
#' @param trim the fraction (0 to 0.25) of observations to be trimmed from each end of x before computing the average distance per cell group.
#' @param interaction.length The maximum interaction/diffusion length of ligands. This hard threshold is used to filter out the connections between spatially distant regions
#' @param spot.size theoretical spot size; e.g., 10x Visium (spot.size = 65 microns)
#' @param spot.size.fullres The number of pixels that span the diameter of a theoretical spot size in the original,full-resolution image.
#' @param k.min the minimum number of interacting cell pairs required for defining adjacent cell groups
# #' @param k.spatial Number of neighbors in a knn graph, which is used to filter out the connections between spatially distant regions that do not share many neighbor spots/cells
#' @importFrom BiocNeighbors queryKNN KmknnParam
#' @return A square matrix giving the pairwise region distance
#'
#' @export
computeRegionDistance <- function(coordinates, group, trim = 0.1,
                                  interaction.length = NULL, spot.size = NULL, spot.size.fullres = NULL, k.min = 10
) {
  if (ncol(coordinates) != 2) {
    stop("Please check the input 'coordinates' and make sure it is a two column matrix.")
  }
  if (!is.factor(group)) {
    stop("Please input the `group` as a factor!")
  }
  # type <- match.arg(type)
  type <- "truncatedMean"
  FunMean <- switch(type,
                    triMean = triMean,
                    truncatedMean = function(x) mean(x, trim = trim, na.rm = TRUE),
                    thresholdedMean = function(x) thresholdedMean(x, trim = trim, na.rm = TRUE),
                    median = function(x) median(x, na.rm = TRUE))

  numCluster <- nlevels(group)
  level.use <- levels(group)
  level.use <- level.use[level.use %in% unique(group)]
  d.spatial <- matrix(NaN, nrow = numCluster, ncol = numCluster)
  adj.spatial <- matrix(0, nrow = numCluster, ncol = numCluster)
  for (i in 1:numCluster) {
    for (j in 1:numCluster) {
      data.spatial.i <- coordinates[group %in% level.use[i], , drop = FALSE]
      data.spatial.j <- coordinates[group %in% level.use[j], , drop = FALSE]
      qout <- suppressWarnings(BiocNeighbors::queryKNN(data.spatial.j, data.spatial.i, k = 1, BNPARAM = BiocNeighbors::KmknnParam(), get.index = TRUE))
      if (!is.null(spot.size) & !is.null(spot.size.fullres)) {
        qout$distance <- qout$distance*spot.size/spot.size.fullres
        idx <- qout$distance - interaction.length < spot.size/2
        adj.spatial[i,j] <- (length(unique(qout$index[idx])) >= k.min) * 1
      }
      d.spatial[i,j] <- FunMean(qout$distance) # since distances are positive values, different ways for computing the mean have little effects.

    }
  }
  d.spatial <- (d.spatial + t(d.spatial))/2
  if (!is.null(spot.size) & !is.null(spot.size.fullres)) {
    adj.spatial <- adj.spatial * t(adj.spatial) # if one is zero, then both are zeros.
    adj.spatial[adj.spatial == 0] <- NaN
    d.spatial <- d.spatial * adj.spatial
  }

  rownames(d.spatial) <- levels(group); colnames(d.spatial) <- levels(group)
  return(d.spatial)

}

