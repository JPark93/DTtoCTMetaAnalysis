
#Steps
# 1. Set the true OU drift matrix (A) >> simulate a data set.
# 2. Sample data at different time scales >> fit VAR(1) models to get estimated AR/CR transition matrices.
# 3. Use the estimated transition matrices as input for the optimization function to estimate the drift matrix >> compare it with the true drift matrix 
# (4. Get the optimal time scale based on the estimated OU process)

library(OpenMx)
library(Matrix)
source('./StateSpaceFunctions.R')
source('./JPmx.R')

library(dynr)
source('./dynrVAR.R') # Auto VAR

library(nloptr)


# 1. Simulate a bivariate OU process ----


##Data Generation ----
set.seed(30214117)
# Large EF condition
y = NULL
maxT <- 10000 # Total number of time points
Deltat = .1 #Resolution of Data Generation
ne = 2 #Number of latent variables
ny = 2 #Number of observed variables
TimeSeq  <- seq(0, (maxT-1)*Deltat,by=Deltat) #Time indices, with, Delta t = 0.1
mu = rep(0, ne) #Home base values
# a0 = c(2,-1, .5,-.5) #c(0,0,0,0)
a0 = rep(0, ne)
P0 = diag(1,ne) #Initial condition co-variance matrix
# Drift
A = -matrix(c(0.20,-0.15, 
              -0.23, 0.25), ne, ne, byrow = TRUE)
# A = -matrix(c(0.02,-0.05, 
#               -0.05, 0.02), ne, ne, byrow = TRUE)

# Mu
b = -A %*% matrix(mu,ncol=1)
# G - Matrix square-root of process noise covariance matrix
G = diag(1,ne)
# Process noise covariance matrix
Q = G%*%t(G) 
# Components related to the discrete-time solution of the OU model
A_del = expm(A*Deltat)
b_del= solve(A)%*%(expm(A*Deltat)-diag(1,ne))%*%b
# Residual covariance matrix, Psi
A_hashtag = kronecker(A,diag(rep(1,ne))) + kronecker(diag(rep(1,ne)),A)
Qtorow = matrix(Q,ncol=1)
Psi = matrix(solve(A_hashtag)%*%(expm(A_hashtag*Deltat)-diag(rep(1,ne^2)))%*%Qtorow,ncol=ne)
# Measurement-related parameters
Lambda = diag(1,ne) #Factor loading matrix
R = diag(1e-5,ne) #Measurement error covariance matrix
thedat = simStateSpaceData(a0 = a0, P0=P0, Q=Psi, R=R, Phi=A_del, Lambda = Lambda,
                           alpha=b_del, tau = rep(0,ny),nt = maxT,np = 1,ne = ne,
                           ny = ny,nx=0,npad=50)
y = data.frame(ID = rep(1,each=maxT), Time = TimeSeq,
               y1 = thedat$obsData[,1], y2 = thedat$obsData[,2])

plot(1:maxT, y$y1, type = 'l')
lines(1:maxT, y$y2, type = 'l', col = 'red')



# 2. Sample data at different time scales & fit VAR(1) models to get estimated transition matrices ----

## Sampling ----
intervals <- c(1, 10, 15)
data <- list()
for (i in intervals) {
  data[[as.character(i)]] <- y[y$Time %% i == 0, ]
}
names(data) <- paste0("dt", intervals) # dt5, dt10, dt100

par(mfrow=c(3,1))
plot(data$dt1$y1, type='l');lines(data$dt1$y2, col='red')
plot(data$dt10$y1, type='l');lines(data$dt10$y2, col='red')
plot(data$dt15$y1, type='l');lines(data$dt15$y2, col='red')

## Fit VAR(1) model to data ----

# Auto VAR (dynr)
DTmatrices = list()
for (i in 1:length(intervals)){
  fullresult = dynr.VAR(data = data[[i]], varnames = paste0('y', 1:2), 
                        id = 'ID', time = 'Time', MLVAR = FALSE)
  DTmatrices[[i]] <- matrix(fullresult[[1]][[1]]$Res$fitted.parameters[1:4], ncol=ne)
  #DTmatrices[[i]] <- fullresult[[1]][[1]]$Betas
}



# Plot OU
par(mfrow=c(1,1))
delts = seq(0.0, 20, 0.5)
transforms = matrix(NA, length(delts), ne^2)
index = 1
for(i in delts){
  transforms[index,] = as.numeric(expm(A*i)) # dynamic factors in DT-framework
  index = index + 1
}
plot(delts, transforms[,1], type = 'l', ylim = c(0,1.1), lwd = 1.5, lty = 2)
lines(delts, transforms[,2], type = 'l', col = 'red', lwd = 1.5)
lines(delts, transforms[,3], type = 'l', col = 'blue', lwd = 1.5)
lines(delts, transforms[,4], type = 'l', col = 'darkgray', lwd = 1.5, lty = 2)

# Add samples (3 studies; AR & CR parameters)
for (i in 1:length(intervals)){
  points(y=DTmatrices[[i]][1], x=intervals[[i]])
  points(y=DTmatrices[[i]][2], x=intervals[[i]], col = 'red')
  points(y=DTmatrices[[i]][3], x=intervals[[i]], col = 'blue')
  points(y=DTmatrices[[i]][4], x=intervals[[i]], col = 'darkgray')
}



# 3. Use the DTmatrices as input for the optimization function to estimate the drift matrix ----

# Initial values for optimization
params = rnorm(4, 0, 0.1)
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
opts = list("algorithm" = "NLOPT_LN_NELDERMEAD",
            "xtol_rel" = 1e-6,
            "maxeval" = 2000)
# Running the Optimization
output = nloptr(x0 = rnorm(4, 0, 0.1),  #starting values
                eval_f = eval_f0, 
                matrices = DTmatrices, 
                deltas = intervals, 
                opts = opts,
                lb = rep(-2, 4), # the range can affect the quality
                ub = rep(2, 4))

# Comparing estimated to true matrix, A
(est_A = matrix(output$solution, 2, 2, byrow = TRUE))
A






# 4. So... what is the optimal time scale? ----

# delta_t where the CR is max
(index_maxCR1 <- which.max(transforms[,2]))
(index_maxCR2 <- which.max(transforms[,3]))
(opt_dt <- delts[(index_maxCR1+index_maxCR2)/2])


# Plot OU by the estimated matrix
delts = seq(0.0, 20, 0.1)
transforms = matrix(NA, length(delts), ne^2)
index = 1
for(i in delts){
  transforms[index,] = as.numeric(expm(est_A*i)) # dynamic factors in DT-framework
  index = index + 1
}
plot(delts, transforms[,1], type = 'l', ylim = c(0,1.1), lwd = 1.5, lty = 2)
lines(delts, transforms[,2], type = 'l', col = 'red', lwd = 1.5)
lines(delts, transforms[,3], type = 'l', col = 'blue', lwd = 1.5)
lines(delts, transforms[,4], type = 'l', col = 'darkgray', lwd = 1.5, lty = 2)
abline(v=opt_dt, lty=2)




