#'generate SMILES strings from extended N-gram model
#' @description generate SMILES strings from extended N-gram model
#' @param nsmis number of generating SMILES strings
#' @param engram ENgram object
#' @param order n in ENgram model
#' @param gentype Back-off procedure with "ML" option and Neaser-Nay smoothing with "KN" option
#' @param crange range in the length of output SMILES strings
#' @examples data(engram_5k)
#' smiles <- genENgram(4, engram_5k, 10)
#' viewstr(smiles)
#'
#' @export  genENgram 
genENgram <- function(nsmis, engram, order, gentype="ML", crange=c(10, 100)){
  res <- character(nsmis)
  for(i in 1:nsmis){
    cat("\r", i, "th molecules generated")
    flag <- F
    while(!flag){
      tsmi <- Esmi$new("C", m=order, engram, type=gentype)
      while(tsmi$prev!="#term#"){
        tsmi$chem_local(engram, 0, 1)
      }
      if(class(try(parse.smiles(tsmi$get_validsmi(), kekulise=T), silent=T))=="list"){
        smilength <- nchar(tsmi$get_validsmi())
        if( (smilength >= crange[1]) & (smilength <= crange[2]) ){
          res[i] <- tsmi$get_validsmi()
          flag <- T
        }
      }
    }
  }
  cat("\n")
  return(res)
}