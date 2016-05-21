#' iqspr
#' 
#' Generate chemical structures from Inverse-QSPR model
#' 
#' @name iqspr
#' @docType package
#' @import rcdk
#' @examples
#'#sample data
#'data(qspr.data)
#'idx <- sample(nrow(qspr.data), 5000)
#'smis <- paste(qspr.data[idx,1])
#'y <- qspr.data[idx,c(2,5)]
#'
#'#learning a pattern of chemical strings
#'data(trainedSMI)
#'data(engram_5k)  #same as run => engram <- ENgram$new(trainedSMI, order=10)
#'
#'#learning QSPR model
#'data(qsprpred_EG_5k)
#'#same as run => qsprpred <- QSPRpred$new(smis=smis, y=as.matrix(y), v_fpnames="graph")
#'
#'#set target range 
#'qsprpred_EG_5k$ymin <- c(200, 1.5)
#'qsprpred_EG_5k$ymax <- c(350, 2.5)
#'
#'#getting chemical strings from the Inverse-QSPR model
#'smchem <- SmcChem$new(smis = rep("c1ccccc1O", 25), v_qsprpred=qsprpred_EG_5k,
#'                      v_engram=engram_5k,temp=3, decay=0.95)
#'smchem$smcexec(niter=5, preorder=0, nview=4)
#'#if OpenBabel (>= 2.3.1) is installed, you can use reordering for better mixing as 
#'#smchem$smcexec(niter=100, preorder=0.2, nview=4)
#'#see http://openbabel.org
#'
#'#check
#'smiles <- get_smiles(smchem)
#'predict(qsprpred_EG_5k, smiles[1:5]) 
NULL 