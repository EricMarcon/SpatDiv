#' SpatDiv
#'
#' Spatially Explicit Measures of Diversity
#' 
#' This package extends the \emph{entropart} package \insertCite{Marcon2014c}{SpatDiv}.
#' It provides spatially explicit measures of diversity such as the mixing index.
#'
#' @name SpatDiv
#' @docType package
#' @import Rcpp
#' @importFrom RcppParallel RcppParallelLibs
#' @importFrom Rdpack reprompt
#' @useDynLib SpatDiv, .registration = TRUE
#' @references
#' \insertAllCited{}
NULL

#' Paracou plot 6
#'
#' This dataset is from Paracou field station, French Guiana, managed by \href{http://www.cirad.fr}{Cirad}.
#' It contains the position, species and basal area of all trees above 10 cm diameter at breast height (DBH) in a 6.25ha plot.
#'
#' @format An object of class \code{\link{wmppp}}.
#' @source Permanent data census of Paracou: \url{http://paracou.cirad.fr/}
"Paracou6"


#' CheckSpatDivArguments
#'
#' Checks the arguments of a function of the package SpatDiv
#'
#' The function compares the arguments passed to its parent function to the type they should be and performs some extra tests (\emph{e.g.} probabilities must be positive and sum to 1). It stops if an argument is not correct.
#'
#' @return Returns \code{TRUE} or stops if a problem is detected.
#' 
#' @export
#'
#' @author Eric Marcon <Eric.Marcon@ecofog.gf>
#'
#' @keywords internal
CheckSpatDivArguments <-
function() {

  # Verify that the package is attached
  if (! "SpatDiv" %in% .packages()) {
    warning("Function arguments cannot be checked because the SpatDiv package is not attached. Add CheckArguments=FALSE to suppress this warning or run library('SpatDiv')")
    return (TRUE)
  }
  # Get the list of arguments of the parent function
  ParentFunction <- sys.call(-1)[[1]]
  # If apply() or similar was used, the function name is not in ParentFunction: sys.call(-1)[[1]] returns "FUN"
  if (ParentFunction == "FUN") {
    warning("Function arguments cannot be checked, probably because you used apply(). Add CheckArguments=FALSE to suppress this warning.")
    return (TRUE)
  }

  ErrorFunction <- paste("Error in ", ParentFunction, ":")

  # Find the arguments. match.fun does not work with SpatDiv::function
  ParentFunctionNoNS <- as.name(gsub("SpatDiv::", "", as.character(ParentFunction)))
  Args <- formals(match.fun(ParentFunctionNoNS))

  ErrorMessage <- function(Message, Argument) {
    cat(deparse(substitute(Argument)), "cannot be:\n")
    print(utils::head(Argument))
    cat(paste(ErrorFunction, Message))
    stop("Check the function arguments.", call. = FALSE)
  }

  # MC
  if (!is.na(names(Args["MC"]))) {
    MC <- eval(expression(MC), parent.frame())
    if (!entropart::is.MetaCommunity(MC))
      ErrorMessage("MC must be a MetaCommunity object.", MC)
  }

  # q
  if (!is.na(names(Args["q"]))) {
    q <- eval(expression(q), parent.frame())
    if (!is.numeric(q) | length(q)!=1)
      ErrorMessage("q must be a number.", q)
  }
  # q.seq
  if (!is.na(names(Args["q.seq"]))) {
    q.seq <- eval(expression(q.seq), parent.frame())
    if (!is.vector(q.seq))
      ErrorMessage("q.seq must be a numeric vector.", q.seq)
  }

  # BootstrapMethod
  if (!is.na(names(Args["BootstrapMethod"]))) {
    BootstrapMethod <- eval(expression(BootstrapMethod), parent.frame())
    if (!is.character(BootstrapMethod))
      ErrorMessage("BootstrapMethod must be a string.", BootstrapMethod)
  }

  # mean
  if (!is.na(names(Args["mean"]))) {
    mean <- eval(expression(mean), parent.frame())
    if (!is.numeric(mean) | length(mean)!=1)
      ErrorMessage("mean must be a number.", mean)
  }

  # n
  if (!is.na(names(Args["n"]))) {
    n <- eval(expression(n), parent.frame())
    if (!is.numeric(n) | length(n)!=1)
      ErrorMessage("n must be a number.", n)
    if (any(n < 1))
      ErrorMessage("n must be at least 1.", n)
    if (as.integer(n) != n)
      ErrorMessage("n must be an integer.", n)
  }

  # NumberOfSimulations
  if (!is.na(names(Args["NumberOfSimulations"]))) {
    NumberOfSimulations <- eval(expression(NumberOfSimulations), parent.frame())
    if (!is.numeric(NumberOfSimulations))
      ErrorMessage("NumberOfSimulations must be a number.", NumberOfSimulations)
    if (NumberOfSimulations < 0)
      ErrorMessage("NumberOfSimulations must be positive.", NumberOfSimulations)
  }

  # NorP
  if (!is.na(names(Args["NorP"]))) {
    NorP <- eval(expression(NorP), parent.frame())
    if (!is.numeric(NorP))
      ErrorMessage("NorP must be numeric.", NorP)
    if (any(NorP < 0))
      ErrorMessage("All NorP values must be positive.", NorP)
    if (!is.vector(NorP) & !entropart::is.SpeciesDistribution(NorP)) {
      # NorP may be a true vector or a SpeciesDistribution. Then dim(NorP) is NULL, and nothing more has to be checked
      # or a "named vector" whose attributes are not "names". Then dim() returns the vector's length.
      if (length(dim(NorP)) != 1) {
        # or a 2D numeric object
        if (length(dim(NorP)) == 2) {
          if (dim(NorP)[2] > 2) {
            # then it must have 1 or 2 columns
            ErrorMessage("NorP may be a vector or a two-column matrix.", NorP)
          }
        } else {
          ErrorMessage("NorP may be a vector or a two-column matrix.", NorP)
        }
      }
    }
  }

  # sd
  if (!is.na(names(Args["sd"]))) {
    sd <- eval(expression(sd), parent.frame())
    if (!is.numeric(sd) | length(sd)!=1)
      ErrorMessage("sd must be a number.", sd)
    if (any(sd < 0))
      ErrorMessage("sd must be positive.", sd)
  }

  # Tree
  if (!is.na(names(Args["Tree"]))) {
    Tree <- eval(expression(Tree), parent.frame())
    if (!is.null(Tree)) {
      if (!inherits(Tree, "phylo") & !inherits(Tree, "phylog") & !inherits(Tree, "hclust") & !inherits(Tree, "PPtree"))
        ErrorMessage("Tree may be NULL or an object of class hclust or phylo or phylog or PPtree.", Tree)
      if (inherits(Tree, "phylog")) {
        if (is.null(Tree$Wdist))
          ErrorMessage("phylog Tree must contain a distance matrix (use add.tools=TRUE when creating it).", Tree)
      }
    }
  }
  # PhyloTree
  if (!is.na(names(Args["PhyloTree"]))) {
    PhyloTree <- eval(expression(PhyloTree), parent.frame())
    if (!is.null(PhyloTree)) {
      if (!inherits(PhyloTree, "phylo") & !inherits(PhyloTree, "phylog") & !inherits(PhyloTree, "hclust") & !inherits(PhyloTree, "PPtree"))
        ErrorMessage("PhyloTree may be NULL or an object of class hclust or phylo or phylog or PPtree", PhyloTree)
      if (inherits(PhyloTree, "phylog")) {
        if (is.null(PhyloTree$Wdist))
          ErrorMessage("phylog PhyloTree must contain a distance matrix (use add.tools=TRUE when creating it).", PhyloTree)
      }
    }
  }

  # prob
  if (!is.na(names(Args["prob"]))) {
    prob <- eval(expression(prob), parent.frame())
    if (!is.numeric(prob) | length(prob)!=1)
      ErrorMessage("prob must be a number.", prob)
    if (any(prob < 0) | any(prob > 1))
      ErrorMessage("prob must be between 0 and 1.", prob)
  }

  # r.seq
  if (!is.na(names(Args["r.seq"]))) {
    r.seq <- eval(expression(r.seq), parent.frame())
    if (!is.null(r.seq)) {
      if (!is.numeric(r.seq) && !is.vector(r.seq))
        stop(paste(ErrorFunction, "r.seq must be a numeric vector"))
      if (length(r.seq) < 2)
        stop(paste(ErrorFunction, "r.seq has length", length(r.seq), "- must be at least 2"))
      if (r.seq[1] != 0)
        stop(paste(ErrorFunction, "First r.seq value must be 0"))
      if (any(diff(r.seq) <= 0))
        stop(paste(ErrorFunction, "successive values of r.seq must be increasing"))
    }
  }

    # S
  if (!is.na(names(Args["S"]))) {
    S <- eval(expression(S), parent.frame())
    if (!is.numeric(S) | length(S)!=1)
      ErrorMessage("S must be a number.", S)
    if (any(S < 1))
      ErrorMessage("S must be at least 1.", S)
    if (as.integer(S) != S)
      ErrorMessage("S must be an integer.", S)
  }

  # size
  if (!is.na(names(Args["size"]))) {
    size <- eval(expression(size), parent.frame())
    if (!is.numeric(size) | length(size)!=1)
      ErrorMessage("size must be a number.", size)
    if (any(size < 1))
      ErrorMessage("size must be at least 1.", size)
    if (as.integer(size) != size)
      ErrorMessage("size must be an integer.", size)
  }

  # Z
  if (!is.na(names(Args["Z"]))) {
    Z <- eval(expression(Z), parent.frame())
    if (!is.null(Z)) {
      if (!is.matrix(Z)) {
        ErrorMessage("Z must be a square matrix.", Z)
      } else {
        if (dim(Z)[1] != dim(Z)[2])
          ErrorMessage("Z must be a square matrix.", Z)
        if (!is.null(colnames(Z)) | !is.null(rownames(Z))) {
          # If the matrix is named, rows and columns must have the same names
          if (!identical(colnames(Z), rownames(Z)))
            ErrorMessage("Z row and column names must be identical.", Z)
        }
        # Must be a relatedness matrix
        if (any(Z<0))
          ErrorMessage("All terms of the relatedness matrix Z must be positive.", Z)
        if (any(diag(Z)<0))
          ErrorMessage("All terms of the relatedness matrix Z diagonal must be strictly positive.", Z)
      }
    }
  }

  return (TRUE)
}
