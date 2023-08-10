# DT to CT Meta Analysis
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
  params = rnorm(nv, 0, 0.1)
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
  opts = list("algorithm" = "algorithm",
              "xtol_rel" = tol,
              "maxeval" = eval)
# Running the Optimization
  output = nloptr(x0 = params,  #starting values
                  eval_f = eval_f0, 
                  matrices = matrices, 
                  deltas = deltas, 
                  opts = opts,
                  lb = rep(-2, nv), # the range can affect the quality
                  ub = rep(2, nv))

# Comparing estimated to true matrix, A
  est_A = matrix(output$solution, nv, nv, byrow = TRUE)
  return(est_A)
}
