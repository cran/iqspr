#'modify chemical structures with SMC
#' @description modify chemical structures with SMC
#' @param smchem SmcChem object
#' @param niter number of iterations
#' @param nsteps total number of letters extending and contracting in a iteration
#' @param preorder probability for reordering SMILES
#' @param nview number of structures shown after each iteration (if zero, structures are not shown)
#' @examples 
#'data(engram_5k) 
#'data(qsprpred_EG_5k)
#'smchem <- SmcChem$new(smis = rep("c1ccccc1O", 25), v_qsprpred=qsprpred_EG_5k,
#'                      v_engram=engram_5k,temp=3)
#'smcexec(smchem, niter=5, preorder=0, nview=4)
#'#if OpenBabel (>= 2.3.1) is installed, you can use reordering for better mixing as 
#'#smcexec(smchem, niter=100, preorder=0.2, nview=4)
#'#see http://openbabel.org
#' @export  smcexec
smcexec <- function(smchem, niter, nsteps=5, preorder=0, nview=0){
  smchem$smcexec(niter, nsteps, preorder, nview)
}