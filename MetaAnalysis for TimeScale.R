library(nloptr)
library(Matrix)
set.seed(03021993)
# Number of variables
  ne = 2
# The true, data-generating matrix
  themat = -Matrix(c(0.20,-0.15, -0.23, 0.25), ne, ne, byrow = TRUE)
  # Zita's OU model:
    # themat = -Matrix(c(0.02,-0.05, -0.05, 0.02), ne, ne, byrow = TRUE)
# Plotting OU by Delta-T
  delts = seq(0.0, 20, 0.5)
  transforms = matrix(NA, length(delts), ne^2)
  index = 1
  for(i in delts){
    transforms[index,] = as.numeric(expm(themat*i))
    index = index + 1
  }
  plot(delts, transforms[,1], type = 'l', ylim = c(0,1.1), lwd = 1.5, lty = 2)
    lines(delts, transforms[,2], type = 'l', col = 'red', lwd = 1.5)
    lines(delts, transforms[,3], type = 'l', col = 'blue', lwd = 1.5)
    lines(delts, transforms[,4], type = 'l', col = 'gray', lwd = 1.5, lty = 2)

# "Sampling" 3 studies
  (mat1 = expm(themat*0.5))
  (mat2 = expm(themat*1))
  (mat3 = expm(themat*10))
# Making a list and adding some small noise
  matrices = list(mat1 = mat1+rnorm(4,0,0.05),
                  mat2 = mat2+rnorm(4,0,0.05),
                  mat3 = mat3+rnorm(4,0,0.05))
# Matching Matrices to Corresponding Delta T's
  deltas = list(del1 = 0.5,
                del2 = 1.0,
                del3 = 10)
# Initial values for optimiation
  params = rnorm(4, 0, 0.1)
# Objective Function to be Minimized
  eval_f0 <- function(params, matrices, deltas) {
    # Reshape params into a matrix
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
  opts = list("algorithm" = "NLOPT_LN_NELDERMEAD",
              "xtol_rel" = 1e-6,
              "maxeval" = 2000)
# Running the Optimization
  output = nloptr(x0 = rnorm(4, 0, 0.1), eval_f = eval_f0, 
                  matrices = matrices, deltas = deltas, 
                  opts = opts,
                  lb = rep(-0.9, 4),
                  ub = rep(.9, 4))
# Comparing estimated to true matrix
  matrix(output$solution, 2, 2, byrow = TRUE)
  themat

  