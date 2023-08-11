rm(list=ls())
args=(commandArgs(TRUE))
if(length(args)==0){
  print("No arguments supplied.")
  ##supply default values
}else{
  for(i1 in 1:length(args)){
    eval(parse(text=args[[i1]]))
  }
}

library(OpenMx)
library(Matrix)
library(dynr)
library(nloptr)
source('~/work/StateSpaceFunctions.R')
source('~/work/dynrVAR.R') # Auto VAR
# source(paste0(getwd(),"/../dynrautoVAR/dynrautoVAR/dynrVAR.R"))

thedir = '~/scratch/'
# Just making relevant directories for everything to save in one place
dir.create(file.path(paste0(thedir,'/SC/', sep = '')), showWarnings = F)
dir.create(file.path(paste0(thedir,'/SC/Rep/', sep = '')), showWarnings = F)


# 1. Simulate a bivariate OU process ----
bias = matrix(NA, nrep, 4)
biasB = matrix(NA, nrep, 4)
delb = matrix(NA, nrep, 1)
delB = matrix(NA, nrep, 1)
##Data Generation ----
  set.seed(03021993+x)
  y = NULL
  maxT <- 10000 # Total number of time points
  Deltat = .1 #Resolution of Data Generation
  ne = 2 #Number of latent variables
  ny = 2 #Number of observed variables
  TimeSeq  <- seq(0, (maxT-1)*Deltat,by=Deltat) #Time indices, with, Delta t = 0.1
  mu = rep(0, ne) #Home base values
  a0 = rep(0, ne)
  P0 = diag(1,ne) #Initial condition co-variance matrix
  # Drift
  A = -matrix(c(0.20,-0.15, 
                -0.23, 0.25), ne, ne, byrow = TRUE)
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
  # R = diag(.3,ne) #Measurement error covariance matrix
  thedat = simStateSpaceData(a0 = a0, P0=P0, Q=Psi, R=R, Phi=A_del, Lambda = Lambda,
                             alpha=b_del, tau = rep(0,ny),nt = maxT,np = 1,ne = ne,
                             ny = ny,nx=0,npad=50)
  y = data.frame(ID = rep(1,each=maxT), Time = TimeSeq,
                 y1 = thedat$obsData[,1], y2 = thedat$obsData[,2])
# Auto VAR (dynr)
  DTmatrices = list()
  DTmatricesB = list()
  for (i in 1:length(intervals)){
    fullresult = dynr.VAR(data = data[[i]], varnames = paste0('y', 1:2), 
                          id = 'ID', time = 'Time', MLVAR = FALSE)
    DTmatrices[[i]] <- matrix(fullresult[[1]][[1]]$Res$fitted.parameters[1:4], ncol=ne)
    DTmatricesB[[i]] <- fullresult[[1]][[1]]$Betas
  }
delts = seq(0.0, 20, 0.5)
transforms = matrix(NA, length(delts), ne^2)
index = 1
for(i in delts){
  transforms[index,] = as.numeric(expm(A*i)) # dynamic factors in DT-framework
  index = index + 1
}
opt_dtTRUE = delts[which.max(rowMeans(transforms[,2:3]))]

est_A = dt2ct(matrices = DTmatrices, deltas = intervals, 
              tol = 1e-6, eval = 2000, algorithm = 'NLOPT_LN_NELDERMEAD',
              lb = rep(-1, 4), ub = rep(1, 4))

est_A1 = dt2ct(matrices = DTmatricesB, deltas = intervals, 
              tol = 1e-6, eval = 2000, algorithm = 'NLOPT_LN_NELDERMEAD',
              lb = rep(-1, 4), ub = rep(1, 4))

# Comparing estimated to true matrix, A
delts = seq(0.0, 20, 0.1)
transforms = matrix(NA, length(delts), ne^2)
index = 1
for(i in delts){
  transforms[index,] = as.numeric(expm(est_A*i)) # dynamic factors in DT-framework
  index = index + 1
}
 opt_dt = delts[which.max(rowMeans(transforms[,2:3]))]

transforms = matrix(NA, length(delts), ne^2)
index = 1
for(i in delts){
  transforms[index,] = as.numeric(expm(est_A*i)) # dynamic factors in DT-framework
  index = index + 1
}
 opt_dt1 = delts[which.max(rowMeans(transforms[,2:3]))]
 
 bias[x,] = est_A - A
 biasB[x,] = est_A1 - A
 delb[x,] = opt_dt - opt_dtTRUE
 delbB[x,] = opt_dt1 - opt_dtTRUE

 write.table(bias,file=paste0(thedir,"/SC/Rep/",x,"Bias.txt"),sep=",")
 write.table(biasB,file=paste0(thedir,"/SC/Rep/",x,"BiasB.txt"),sep=",")
 write.table(delb,file=paste0(thedir,"/SC/Rep/",x,"delb.txt"),sep=",")
 write.table(delbB,file=paste0(thedir,"/SC/Rep/",x,"delbB.txt"),sep=",")

q('no')