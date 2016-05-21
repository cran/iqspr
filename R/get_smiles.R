#'get SMILES strings from the SmcChem object
#' @description get SMILES strings from the SmcChem object
#' @param smchem SmcChem class object
#' @examples 
#'data(engram_5k) 
#'data(qsprpred_EG_5k)
#'smchem <- SmcChem$new(smis = rep("c1ccccc1O", 25), v_qsprpred=qsprpred_EG_5k,
#'                      v_engram=engram_5k,temp=3)
#'smcexec(smchem, niter=5, preorder=0, nview=4)
#'get_smiles(smchem)
#'
#' @export  get_smiles
get_smiles <- function(smchem){
  smchem$get_smiles()
}