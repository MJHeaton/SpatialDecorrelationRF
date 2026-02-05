##
## Functions to estimate nugget & range
##

# fxDir <- getSrcDirectory(function(x) {x})
# source(paste0(fxDir,"/mkNNIndx.R"))
# source(paste0(fxDir,"/fitMaternGP.R"))
source("./mkNNIndx.R")
source("./fitMaternGP.R")

##
## Function to transform to independence
##

transform_to_ind <- function(formula,
                           trainData,
                           trainLocs,
                           testData, #Don't include response
                           testLocs,
                           MaternParams=NULL, #Either null or a 2 vector of (rng, nug)
                           smoothness=1/2,
                           M = 30, #num neighbors
                           ncores=detectCores()-1){
  
  ######################################
  ## Figure out the nearest neighbors ##
  ######################################
  nnList <- mkNNIndx(trainLocs, m=M)
  
  #################################################
  ## Estimate a range and nugget if not provided ##
  #################################################
  if(is.null(MaternParams)){
    mFit <- fit.NN.Matern(formula, data=trainData, locs=trainLocs, nu=smoothness,
                          NearNeighs=nnList, num.cores=ncores)
    range <- 1/mFit$decay
    nugget <- mFit$nugget
  } else {
    range <- MaternParams[1]
    nugget <- MaternParams[2]
  }
  
  
  ################################
  ## Transform the Training Set ##
  ################################
  
  ## Define X and y matrices
  Xtrain <- model.matrix(formula, data=trainData)
  ytrain <- as.matrix(trainData[,all.vars(formula)[1]], ncol=1)
  
  ## Apply decorrelating transform by location
  indData <- mclapply(1:nrow(Xtrain), FUN=function(idx){
    if(idx==1){
      y <- ytrain[idx]
      w <- 1
      X <- Xtrain[idx,] / sqrt(w)
    } else if(idx==2){
      D <- rdist(trainLocs[1:idx,])
      R <- (1-nugget)*Matern(D, nu=1/2, range=range) + 
        nugget*diag(nrow(D))
      w <- as.numeric(1-R[1,-1]%*%solve(R[-1,-1])%*%R[-1,1])
      X <- (t(Xtrain[idx,]) - (R[1,-1]%*%solve(R[-1,-1])%*%Xtrain[nnList[[idx]],])) / sqrt(w)
      y <- (ytrain[idx]-R[1,-1]%*%solve(R[-1,-1])%*%(ytrain[nnList[[idx]]]))/sqrt(w)
    } else {
      D <- rdist(trainLocs[c(idx,nnList[[idx]]),])
      R <- (1-nugget)*Matern(D, nu=1/2, range=range) + 
        nugget*diag(nrow(D))
      w <- as.numeric(1-R[1,-1]%*%solve(R[-1,-1])%*%R[-1,1])
      X <- (t(Xtrain[idx,]) - (R[1,-1]%*%solve(R[-1,-1])%*%Xtrain[nnList[[idx]],])) / sqrt(w)
      y <- (ytrain[idx]-R[1,-1]%*%solve(R[-1,-1])%*%(ytrain[nnList[[idx]]]))/sqrt(w)
    }
    
    return(list(y=y, X=X, w=w))
  }, mc.cores=ncores) # End mclapply()
  
  ## Apply decorrelating transform to test data
  Xtest <- model.matrix(formula[-2], data=testData)
  indTestData <- mclapply(1:nrow(Xtest), FUN=function(idx){
    D <- rdist(matrix(testLocs[idx,], nrow=1), trainLocs)
    theNeighbors <- order(D)[1:M]
    R <- rdist(rbind(testLocs[idx,],trainLocs[theNeighbors,]))
    R <- nugget*diag(M+1)+(1-nugget)*Matern(R, range=range, smoothness=smoothness)
    R12 <- R[1,-1]%*%chol2inv(chol(R[-1,-1]))
    w <- as.numeric(1-R12%*%R[-1,1])
    X <- (t(Xtest[idx,])-R12%*%Xtrain[theNeighbors,])/sqrt(w)
    return(list(backTrans=R12%*%matrix(ytrain[theNeighbors,], ncol=1), X=X, 
                w=w))
  }, mc.cores=ncores)
  
  ## Return transformed data
  outList <- list(trainData=data.frame(y=do.call(rbind, lapply(indData, function(x){x$y})),
                                       do.call(rbind, lapply(indData, function(x){x$X}))),
                  testData=data.frame(do.call(rbind, lapply(indTestData, function(x){x$X}))),
                  range=range,
                  nugget=nugget,
                  M=M,
                  formula=formula,
                  backTransformInfo=lapply(indTestData,function(x){x$X<-NULL
                  return(x)}))
  return(outList)
  
  
} # End spatial_to_ind function

back_transform_to_spatial <- function(preds, transformObj){
  
  spatialPreds <- preds*sapply(transformObj$backTransformInfo, function(x){x$w})+
    sapply(transformObj$backTransformInfo, function(x){x$backTrans})
  return(spatialPreds)
  
}

# load("../Linear Simulated Data/LinSimDataSet17.RData")
# transformTest <- transform_to_ind(formula=y~.,
#                                   trainData=trainData,
#                                   trainLocs=trainLocs,
#                                   testData=testData[,-1],
#                                   testLocs=testLocs)
# 
# back_transform_to_spatial(rnorm(nrow(testData)), transformTest)






