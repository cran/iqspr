#'view 2D structures from SMILES string vector
#' @description view 2D structures from SMILES string vector
#' @param smis SMILES string vector
#'
#' @examples viewstr(c("c1ccc2ccc3c(NCCN(C)C)cc(nc3c2c1)", "c1ccc2ccc3c(NCCN(CC)CCCl)cc(nc3c2c1)",
#'  "c1ccc2ccc3c(NC(CC)CC)cc(nc3c2c1)", "c1ccc2ccc3c(c2c1)ncc(c3NCCNCC=CCCCC)"))
#' 
#' @export viewstr
viewstr <- function(smis){
  mols <- parse.smiles(smis, kekulise=F)
  img <- lapply(mols, function(x) view.image.2d(x, width = 250, height = 250))
  n <- length(img)
  nn <- trunc(sqrt(n-0.1^5),0)+1
  plot(1:5, xlim=c(0,nn), ylim=c(0,nn), type="n", axes = F, ann=F)
  for(i in 0:nn){
    abline(h=i, lty=3)
    abline(v=i, lty=3)
  }
  for(i in 1:length(img)){
    d1 <- (i-1) %% nn
    d2 <- nn- trunc((i-1)/nn, 0) - 1
    rasterImage(img[[i]], d1+0.05 ,d2+0.05 ,d1+0.95,d2+0.95)
  }
}