#'QSPRpredictor class
#' @description QSPR model construction using the Bayesian linear regression model. 
#' 
#' @field ymin vector representing minimal value in each target property
#' @field ymax vector representing maximum value in each target property
#' @field descriptor function transforming SMILES strings into numeric vector (matrix) 
#' @field fpnames name of fingerprint used in a predictor
#' 
#' @examples data(qspr.data)
#' smis <- paste(qspr.data[,1])
#' ty <- qspr.data[,c(2,5)]
#' trainidx <- sample(1:nrow(qspr.data), 5000)
#' testidx <- (1:nrow(qspr.data))[-trainidx][1:100]
#' 
#' data(qsprpred_EG_5k) 
#' ## same as run => qsprpred_EG_5k <- 
#  ##   QSPRpred$new(smis=smis[trainidx], y=as.matrix(ty[trainidx, ]), v_fpnames="graph")
#' 
#' #-----arguments
#' #smis: SMILES string set (character vector) for training
#' #y: property sets (matrix) for training  
#' #v_ymin: minimum value of target properties
#' #v_ymax: maximum value of target properties
#' #v_descriptor: function transforming SMILES strings to feature matrix 
#' #v_fpnames: character vectors indicating finger print used in the descriptor
#' #w0: matrix representing prior mean of coeffients in linear regression model
#' #V0_inv: matrix representing the prior variance of coeffients in linear regression model
#' #a0: numeric value representing the location parameter in gamma prior 
#' #b0: numeric value representing the shape parameter in gamma prior 
#' #------
#' 
#' predictions <- qsprpred_EG_5k$qspr_predx(smis[testidx])  
#' par(mfrow=c(1,2))
#' plot(predictions[[1]][1,], ty[testidx,1])
#' plot(predictions[[1]][2,], ty[testidx,2])
#' 
#' #computing the probability which the properties of test structures is in the target range  
#' # set the minimal values in 2-d target properties
#' qsprpred_EG_5k$ymin <- c(100, 4)
#' # set the maximum values in 2-d target properties
#' qsprpred_EG_5k$ymin <- c(200, 5.5) 
#' # method inverse_predx returns the probability that input SMILES has property in target range
#' qsprpred_EG_5k$inverse_predx("c1ccccc1O")
#' @import methods
#' @importFrom methods new
#' @import rcdk 
#'
#' @export QSPRpred
#' @exportClass QSPRpred

QSPRpred <- setRefClass(
  Class = "QSPRpred",
  
  fields = list(
    ymin = "numeric",
    ymax = "numeric",
    descriptor = "function",
    fpnames = "character",
    XX = "matrix",
    TVn = "matrix",
    Twn = "matrix",
    TaN = "numeric",
    TbN = "numeric"
  ),
  
  methods = list(
    initialize = function(smis=NULL, y=NULL, v_ymin=NULL, v_ymax=NULL, v_descriptor=NULL, v_fpnames=NULL, w0=NULL, V0_inv=NULL, a0=NULL, b0=NULL){
      "initialize the QSPR predictor"
      if(is.null(smis) | is.null(y)){
        cat("Need SMILES strings and properties\n")
        return(NULL)
      }
      if(is.null(v_fpnames)){
        fpnames <<- "graph"
      }else{
        fpnames <<- v_fpnames
      }
      
      if(is.null(v_descriptor)){
        descriptor <<- function(smis, fpnames){
          fp <- NULL
          fp_temp <- NULL       
          for(fpname in fpnames){
            if(fpname==fpnames[length(fpnames)]){
              cat("\rget fingerprint", fpname ,"from SMILES strings          \n")  
            }else{
              cat("\rget fingerprint", fpname ,"from SMILES strings          ")   
            }
            mols <- parse.smiles(smis, kekulise=F) ## temporaly FALSE for kekulise option
            fp.obj <- lapply(mols, get.fingerprint, type=fpname)
            nullidx <- which(sapply(fp.obj, is.null)==T)
            if(length(nullidx)>0){
              idx <- (1:length(fp.obj))[-nullidx]          
              fp_temp <- matrix(0, length(fp.obj), length(fp.to.matrix(fp.obj[idx[1]])))
              fp_temp[idx,] <- fp.to.matrix(fp.obj[idx])
            }else{
              fp_temp <- fp.to.matrix(fp.obj)
            }      
            fp <- cbind(fp, fp_temp)   
          }  
          return(fp)
        }
      }else{
        descriptor <<- v_descriptor
      }
      
      if(is.null(v_ymin)){
        ymin <<- colMeans(y)
      }
      
      if(is.null(v_ymax)){
        ymax <<- colMeans(y)
      }
      
      if(is.null(a0)){
        a0 <- 0
      }
      if(is.null(b0)){
        b0 <- numeric(ncol(y))
      }
      
      X <- descriptor(smis, fpnames)
      
      if(is.null(V0_inv)){
        V0_inv <- diag(rep(1^2, ncol(X)))  ## almost uniform prior
      }
      if(is.null(w0)){
        w0 <- matrix(0, ncol(X), ncol(y))
      }
      cat("getting posterior distribution of parameters")
      XX <<- t(X) %*% X
      TVn <<- solve(V0_inv+XX)
      Twn <<- TVn %*% (V0_inv %*% w0 + t(X) %*% y)
      TaN <<- a0 + nrow(X)/2
      TbN <<- b0 + diag(0.5*(diag(t(y) %*% y) - t(Twn) %*% (V0_inv + XX) %*% Twn )) 
      cat("\rgetting posterior distribution of parameters ==> done \n\n")
    },
    
    set_targety = function(v_ymin, v_ymax){
      "set the target propety range"
      ymin <<- v_ymin
      ymax <<- v_ymax
    },

    qspr_predx = function(smis=NULL){
      "QSPR prediction of input SMILES strings smis"
      newx <- descriptor(smis, fpnames)
      predy <- t(Twn) %*% t(newx)
      temp1 <- TbN/TaN
      temp2 <- diag((1 + newx %*% TVn %*% t(newx)))
      predvar <- temp1 %*% t(temp2)
      
      return(c(list(predy), list(predvar)))
    },
    
    inverse_predx = function(smis=NULL,  temp=1){   
      "Inverse QSPR prediction of input SMILES strings smis for target property range"
      if(is.null(smis) | is.null(ymin) | is.null(ymax)){
        cat("smiles strings and target range of y is necessary\n\n")
        return(NULL)
      }else{
        newx <- descriptor(smis, fpnames)
        predy <- t(Twn) %*% t(newx)
        
        temp1 <- TbN/TaN
        temp2 <- diag((1 + newx %*% TVn %*% t(newx)))
        predvar <- temp1 %*% t(temp2) * temp^2
        
        res <- rep(1, length(smis))
        for(i in 1:length(ymin)){
          lx <- (ymin[i]-predy[i,])/sqrt(predvar[i,])
          ux <- (ymax[i]-predy[i,])/sqrt(predvar[i,])
          res <- res * (sapply(ux, function(x) pt(x, TaN))
                        - sapply(lx, function(x) pt(x, TaN)))  
        }
        return(res) 
      }
    }
  )
)