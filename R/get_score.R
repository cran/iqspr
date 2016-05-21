#'get QSPR scores for structures in the SmcChem object
#' @description get QSPR scores for structures in the SmcChem object
#' @param smchem SmcChem object
#'
#' @examples 
#' data(engram_5k) 
#' data(qsprpred_EG_5k)
#' smchem <- SmcChem$new(smis = rep("c1ccccc1O", 25), v_qsprpred=qsprpred_EG_5k,
#'                      v_engram=engram_5k,temp=3)
#' scores <- matrix(0, 25, 5)
#' for(i in 1:5){
#' smchem$smcexec(1)
#' scores[,i] <- get_score(smchem)
#' boxplot(scores)
#' }
#' 
#' @export get_score
get_score <- function(smchem){
  csmi <- sapply(smchem$esmi, function(x) x$get_validsmi())
  idx <- which(smchem$isvalid==T)
  
  score <- numeric(length(csmi))
  score[idx] <- smchem$qsprpred$inverse_predx(csmi[idx], smchem$temp)
  score[-idx] <- 0
  score
}