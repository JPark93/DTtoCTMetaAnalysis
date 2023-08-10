#' Convert DT VAR(1) matrices to an approximate CT drift matrix
#'
#' @author Jonathan Park
#'
#' This function takes in a series of `lists`; one for the discrete-time VAR(1) transition matrices and a separate list of the corresponding delta-t values for each VAR(1) model. The values are then transformed and a CT drift coefficients matrix is constructed via approximation with the `nloptr` package.
#'
#' @return Returns a pxp drift matrix
#'
#' @references Add references to [dynr.var()] here.
#'
#' @param matrices a list of matrices that are of the same dimension, pxp.
#' @param deltas a vector or list of values that describe the delta-t by which the corresponding discrete-time matrix was gathered (i.e., delta-t = 1.00 is 1 day then 2-measures/day = delta-t = 0.50).
#' @param tol A numeric value that describes the degree of numeric precision desired during optimization. Set to `1e-6` by default.
#' @param eval A numeric value which designates the number of iterations to perform during optimization. By default, this value is set to 2000.
#' @param algorithm The optimization algorithm. Drawn from `nloptr`. By default, `NLOPT_LN_NELDERMEAD` is used.
#' @param lb The lower-bound estimate allowed during estimation. If null, the value of -2 will be used for all variables.
#' @param ub The upper-bound estimate allowed during estimation. If null, the value of 2 will be used for all variables.
#' @import dynr
#' @export
dt2ct = function(matrices = NULL, deltas = NULL, 
                 tol = 1e-6, eval = 2000, algorithm = 'NLOPT_LN_NELDERMEAD',
                 lb = NULL, ub = NULL){
  nv = ncol(matrices[[1]])
  if(is.null(lb)){
    lb = rep(-2, nv)
  }
  if(is.null(ub)){
    ub = rep(2, nv)
  }
# Initial values for optimization
  params = rnorm(nv^2, 0, 0.1)
# Objective Function to be Minimized
  eval_f0 <- function(params, matrices, deltas) {
    # Reshape params into a matrix
    #  target_matrix = drift matrix 
      target_matrix <- matrix(params, nrow = length(params)/2, ncol = length(params)/2, byrow = TRUE)
    # Compute the sum of absolute differences between matrices
      distances = matrix(0, length(matrices), 1)
      for(i in 1:length(matrices)){
        distances[i,] = sum(abs(expm(target_matrix*deltas[[i]]) - matrices[[i]]))
      }
      distance = sum(distances)
      return(distance)
    }
# Options for nloptr
  opts = list("algorithm" = algorithm,
              "xtol_rel" = tol,
              "maxeval" = eval)
# Running the Optimization
  output = nloptr(x0 = params,  #starting values
                  eval_f = eval_f0, 
                  matrices = matrices, 
                  deltas = deltas, 
                  opts = opts,
                  lb = rep(-2, nv^2), # the range can affect the quality
                  ub = rep(2, nv^2))

# Comparing estimated to true matrix, A
  est_A = matrix(output$solution, nv, nv, byrow = TRUE)
  return(est_A)
}
