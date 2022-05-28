#' Common constraints
#'
#' Build a constraint matrix from common simple constraints.
#' 
#' @param p The number of variables.
#' @param first Indicates sign constraint for first coefficient.
#' Recommended for identifiability purposes but overriden by any 
#' other constraint. 
#' @param sign Constrained sign for each coefficient. \code{0}: no constraint, 
#' @param monotone Monotonicity constraint. \code{0}: no constraint, \code{-1}:
#'  decreasing coefficients and \code{1}: increasing coefficients.
#' @param convex Convexity constraint. \code{0}: no constraint, \code{-1}:
#'  convex coefficients and \code{1}: concave coefficients.
#'
#' @details 
#' For monotonicity and convexity / concavity, the function assumes the 
#' coefficients are ordered. For instance, for increasing monotone coefficients,
#' the the first one will be lower than the second, which be lower than the
#' third and so on.
#' 
#' If a different value than -1, 0 or 1 is passed to an argument, 
#' the function uses its sign.
#' 
#' @export
build_constraints <- function(p, first = 1, sign = 0, monotone = 0, convex = 0)
{
  # Initialize list of constraint matrix and constraints
  cmatlist <- vector("list", 4)
  allcons <- sign(c(first, sign, monotone, convex))
  
  # Identifiability constraint
  # Only if there is no overall sign constraint
  if (allcons[1] & all(allcons[-1] == 0)){
    cmatlist[[1]] <- c(1, rep_len(0, p - 1))
  }
  
  # Sign constraint
  if (allcons[2]){
    cmatlist[[2]] <- diag(p)
  }
  
  # Monotone and convexity constraints
  for (i in 3:4){
    if (allcons[i]){
      # Matrix
      cmatlist[[i]] <- allcons[i] * diff(diag(p), diff = i - 2)
      
      # Remove redundant sign constraints
      if (allcons[i - 1]){
        ind <- if (allcons[i - 1] == allcons[i]) 1 else p
        cmatlist[[i - 1]] <- cmatlist[[i - 1]][ind,]
      }
    }
  }
  
  # Put everything together and return
  do.call(rbind, cmatlist)
}
