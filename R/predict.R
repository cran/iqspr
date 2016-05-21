#'predict properties with forward QSPR model
#' @description predict properties with forward QSPR model
#' @param qsprpred QSPRpred class object
#' @param smis SMILES strings
#' @examples
#' data(qsprpred_EG_5k)
#' predict(qsprpred_EG_5k, c("c1ccccc1O", "c1ccccc1F"))
#'
#' @export  predict
#' 
predict <- function(qsprpred, smis){
  res <- t(qsprpred$qspr_predx(smis)[[1]])
  rownames(res) <- smis
  res
}