#'get chemical structures with high QSPR score from SmcChem object
#' @description get chemical structures with high QSPR scores from the SmcChem object after excluding similar structures 
#' @param smchem SmcChem class object
#' @param nsmi maximum number of SMILES strings obtained
#' @param exsim excluding structures that are similar with already chosen structures in that the Tanimoto coefficient >= exsim
#' @examples
#' data(engram_5k) 
#' data(qsprpred_EG_5k)
#' smchem <- SmcChem$new(smis = rep("c1ccccc1O", 25), v_qsprpred=qsprpred_EG_5k,
#'                      v_engram=engram_5k,temp=3)
#' smcexec(smchem,10)
#' res <- get_hiscores(smchem, exsim=0.8)
#' viewstr(res[1:4, 1])
#' 
#' @export  get_hiscores
get_hiscores <- function(smchem, nsmi=50, exsim=0.8){
  smchem$get_hiscores(nsmi, exsim)
}