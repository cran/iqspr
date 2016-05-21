#'Extended N-gram model for leraning SMILES strings
#' @description Extended N-gram model for leraning SMILES strings
#'
#' @examples data(trainedSMI)
#' data(engram_5k)  #same as run => engram <- ENgram$new(trainedSMI, order=10)
#' #-----arguments
#' #smis: character vector, SMILES string set for training an extended N-gram model
#' #order: numeric, value representing the maximum order in modified N-gram model
#' #-----arguments
#' 
#' @import methods
#' @importFrom methods new
#' @export ENgram
#' @exportClass ENgram

ENgram <- setRefClass(
  Class = "ENgram",
  
  fields = list(
    mat = "list",
    N1p = "list",
    denom = "list",
    mat_D = "matrix",
    order = "numeric"
  ),
  
  methods = list(
    initialize = function(smis=NULL, order=NULL){
      "Initialize the extend N-gram model with SMILES strings smis and numeric value order"
      if(is.null(order)){
        order <<- 10
      }else{
        order <<- order
      }
      if(is.null(smis)){
        mat <<- list()
      }else{
        update_mat(smis)
      }
      
    },
    
    update_mat = function(data){
      "update the model with additional SMILES strings data"
      ###function
      enlarge_row <- function(tmat, letter){
        tmat <- rbind(tmat, numeric(dim(tmat)[2]))
        rownames(tmat)[dim(tmat)[1]] <- letter
        return(tmat)
      }
      
      enlarge_col <- function(tmat, letter){
        tmat <- cbind(tmat, numeric(dim(tmat)[1]))
        colnames(tmat)[dim(tmat)[2]] <- letter
        return(tmat)
      }
      
      parsenum <- function(x){
        if(length(grep("^[0-9]$", x))==1){
          return(as.numeric(x))
        }else{
          return(x)
        }
      }
      
      comb_les <- function(vec){
        if(length(vec)==1){return(vec)}
        res <- vec[1]
        for(i in 2:length(vec)){
          res <- paste(res, vec[i], sep="")
        }
        return(res)
      }
      
      #######main proc
      #cand <- which(sapply(as.character(data), function(x) length(grep("\\[[A-Z]", x)))==0)
      cand <- 1:length(data)
      cand <- intersect(cand, which(sapply(as.character(data), function(x) length(grep("\\[s", x)))==0))
      cand <- intersect(cand, which(sapply(as.character(data), function(x) length(grep("\\[c", x)))==0))
      #cand <- intersect(cand, which(sapply(as.character(data), function(x) length(grep("\\+", x)))==0))
      #cand <- intersect(cand, which(sapply(as.character(data), function(x) length(grep("\\-", x)))==0))
      #cand <- intersect(cand, which(sapply(as.character(data), function(x) length(grep("\\.", x)))==0))
      #cand <- intersect(cand, which(sapply(as.character(data), function(x) length(grep("#", x)))==0))
      cand <- intersect(cand, which(sapply(as.character(data), function(x) length(grep("%", x)))==0))
      #cand <- intersect(cand, which(sapply(as.character(data), function(x) length(grep("\\/", x)))==0))
      #cand <- intersect(cand, which(sapply(as.character(data), function(x) length(grep("\\\\", x)))==0))
      cand <- intersect(cand, which(sapply(as.character(data), function(x) nchar(x)>=(order+4))))
      
      res <- c()
      if(!((is.list(mat)==T) & (length(mat)==order))){
        mat <<- list()
        for(m in 1:order){
          p_inmat <- c()
          for(k in 1:20){
            tmp <- matrix(numeric(1), 1, 1)
            rownames(tmp) <- comb_les(rep("C", m))
            colnames(tmp) <- rep("C")
            p_inmat <- c(p_inmat, list(tmp))
          }
          mat <<- c(mat, list(p_inmat))
        }
      }
      for(m in 1:order){
        tempmat <- mat[[m]]
        cat("\n")
        for(i in cand){
          cat("\rdim=", m, " no of molecule=", i)
          skip_flag <- F
          str <- as.character(data[[i]])
          numchecker <- numeric(100)
          numstk <- c()
          start <- 2
          numn <- 1
          numbr <- 0
          vstr <- sapply(1:nchar(str), function(x) substr(str, x, x))
          vstr <- c(vstr, "\t")
          while(1){
            p1 <- which(vstr=="[")
            p2 <- which(vstr=="]")
            if(length(p1)>0){
              vstr[p1[1]] <- paste(vstr[p1[1]:p2[1]], collapse="")
              vstr <- vstr[-((p1[1]+1):p2[1])]
            }else{
              break
            }
          }
          p <- which(vstr=="=")
          if(length(p)>0){
            vstr[p] <- paste(vstr[p], vstr[p+1], sep="")
            vstr <- vstr[-(p+1)]
          }
          
          ini_flag <- F; tm <- m; startbr <- c(); next_flag <- F; brstart <- F; brsteps <- 0
          brpos <- c()
          
          brflag <- F; short_flag <- F
          b <- vstr[1]
          n <- 2
          while(1){
            if(skip_flag==F){
              a <- b
              if(length(b) >= m){
                b <- c(b[-1], vstr[n])
              }else{
                b <- c(b, vstr[n])
              }
              if(length(a)<m){
                short_flag <- T
              }else{
                short_flag <- F
              }
              
              lastc <- parsenum(b[length(b)])
              ### double bond is conbinded into a next atom.
              if(lastc=="("){
                if((brstart==T) & (brsteps==2)){
                  if(m==1){
                    tempa <- "#brterm#"
                  }else{
                    tempa <- c(a[2:length(a)], "#brterm#")
                    #cat(tempa, "\n")
                  }
                  startbr <- c(startbr, list(tempa))
                }
                brstart <- T
                brsteps <- 0
              }
              #cat(lastc, " ", brstart, " ", brsteps, " ",  numn, " ", b, "\n")
              if(brstart==T){
                if(brsteps==2){
                  if(m==1){
                    tempa <- "#brterm#"
                  }else{
                    tempa <- c(a[2:length(a)], "#brterm#")
                  }
                  startbr <- c(startbr, list(tempa))
                  brstart <- F
                  brsteps <- 0
                }else{
                  brsteps <- brsteps + 1 
                }
              }
              if(is.numeric(lastc)){
                if(sum(lastc==numstk)>0){        
                  numdiff <- length(numstk)-which(lastc==numstk)
                  numchecker[lastc] <- 0
                  numstk <- numstk[-which(lastc==numstk)]
                  lastc <- paste("#numt_", numdiff, sep="")
                }else{
                  numchecker[lastc] <- 1 
                  numstk <- c(numstk, lastc)
                  lastc <- "#num#"
                }
              }
              if(lastc==")"){
                lastc <- "#brterm#"
                brflag <- T
              }
              if(lastc=="\t"){
                matidx <- 10*(numbr>0)+numn
                lastc <- "#term#"
                b[length(b)] <- lastc
                if(short_flag == F){
                  if(sum(rownames(tempmat[[matidx]])==comb_les(a))==0){
                    tempmat[[matidx]] <- enlarge_row(tempmat[[matidx]], comb_les(a))
                  }
                  if(sum(colnames(tempmat[[matidx]])==lastc)==0){
                    tempmat[[matidx]] <- enlarge_col(tempmat[[matidx]], lastc)
                  }
                  tempmat[[matidx]][comb_les(a), lastc] <- tempmat[[matidx]][comb_les(a), lastc]+1
                }
                break
              }     
              b[length(b)] <- lastc
              
              matidx <- 10*(numbr>0)+numn
              if(short_flag == F){   
                if(sum(rownames(tempmat[[matidx]])==comb_les(a))==0){
                  tempmat[[matidx]] <- enlarge_row(tempmat[[matidx]], comb_les(a))
                }
                if(sum(colnames(tempmat[[matidx]])==lastc)==0){
                  tempmat[[matidx]] <- enlarge_col(tempmat[[matidx]], lastc)
                }
                tempmat[[matidx]][comb_les(a), lastc] <- tempmat[[matidx]][comb_les(a), lastc]+1
              }
              if(brflag==T){
                matidx <- 10*(numbr>0)+numn
                b <- startbr[[length(startbr)]]
                startbr <- startbr[-length(startbr)]
                
                if(short_flag==F){
                  if(sum(rownames(tempmat[[matidx]])==comb_les(b))==0){
                    tempmat[[matidx]] <- enlarge_row(tempmat[[matidx]], comb_les(b))
                  }
                }
                brflag <- F
              }
            }else{
              skip_flag <- F
            }
            
            if(lastc=="("){
              numbr <- numbr + 1
            }else if(lastc=="#brterm#"){
              numbr <- numbr - 1
            }else if(lastc=="#num#"){
              numn <- numn + 1
            }else if(length(grep("numt_", lastc))==1){
              numn <- numn - 1
            }
      
            n <- n+1
          }        
        }
        for(k in 1:length(tempmat)){
          tempmat[[k]] <- tempmat[[k]][order(rownames(tempmat[[k]])),] 
        }
        res <- c(res, list(tempmat))
      }
      mat <<- res
      
      mat_D <<- matrix(0, order, 20)
      for(i in 1:nrow(mat_D)){
        for(j in 1:ncol(mat_D)){
          n1 <- sum(mat[[i]][[j]]==1)
          n2 <- sum(mat[[i]][[j]]==2)
          mat_D[i,j] <<- n1/(n1+2*n2)
        }
      }
      
      denom <<- list()
      N1p <<- list()
      for(i in 1:length(mat)){
        p_denom <- list()
        p_N1p <- list()
        for(j in 1:length(mat[[i]])){
          if(!is.na(mat_D[i,j])){
            p_denom <- c(p_denom, list(rowSums(mat[[i]][[j]])))
            p_N1p <- c(p_N1p, list(apply(mat[[i]][[j]], 1, function(x) sum(x==1))))
          }else{
            p_denom <- c(p_denom, list(NA))
            p_N1p <- c(p_N1p, list(NA))
          }
        }
        denom <<- c(denom, list(p_denom))
        N1p <<- c(N1p, list(p_N1p))
      }
    }
  )
)
