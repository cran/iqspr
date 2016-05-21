#'SMC chemical generator class
#' @description SMILES generator with Sequence Monte Calro sampler
#'
#' @field qsprpred QSPRpred object 
#' @field engram ENgram object
#' @field m numeric value representing the order of extended N-gram model
#' @field v_ESSth numeric value representing the threshold which resample is done when ESS is less than 100 * v_ESSth
#' @field v_decay numeric value decaying temparature for the target distribution in SMC sampler
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
#'                      v_engram=engram_5k,temp=3)
#' #-----arguments
#' #smis: initial SMILES strings in SMC sampler
#' #v_qsprpred: QSPRpred object 
#' #v_engram: ENgram object
#' #v_m: numeric value representing the order of modified N-gram model
#' #ESSth: numeric resampling threshold at 100 * v_ESSth
#' #temp: numeric, annealing parameter in SMC sampler
#' #lambda: numeric 
#' #decay: numeric, decaying rate for temparature 
#' #-----arguments
#'                 
#'smchem$smcexec(niter=5, preorder=0, nview=4)
#'#if OpenBabel (>= 2.3.1) is installed, you can use reordering for better mixing as 
#'#smchem$smcexec(niter=100, preorder=0.2, nview=4)
#'#see http://openbabel.org
#'
#'#check
#'smiles <- get_smiles(smchem)
#'predict(qsprpred_EG_5k, smiles[1:5])
#'
#' @import methods
#' @importFrom graphics abline plot rasterImage
#' @importFrom methods new
#' @export SmcChem
#' @exportClass SmcChem

SmcChem <- setRefClass(
  Class = "SmcChem",
  
  fields = list(
    esmi = "list",
    oldesmi = "list",
    m = "numeric",
    engram = "ENgram",
    qsprpred = "QSPRpred",
    weights = "numeric",
    oldscore = "numeric",
    fkernel = "numeric",
    bkernel = "numeric",
    isvalid = "logical",
    temp = "numeric",
    v_ESSth = "numeric",
    v_decay = "numeric",
    v_maxstock = "numeric",
    smistock = "character",
    scorestock = "numeric"
  ),
  
  methods = list(
    initialize = function(smis=NULL, v_engram=NULL, v_qsprpred=NULL, v_m=NULL, temp=1, ESSth=0.5, lambda=0, decay=0.95, gentype="ML", maxstock=2000){
      "Initialize the SMC chemical generator with initial SMILES strings smis, ENgram class object v_engram and QSPRpred class object v_qsprpred"
      if(is.null(v_engram)){
        cat("Need ENgram object\n")
        return(NULL)
      }else{
        engram <<- v_engram
      }
      
      if(is.null(v_qsprpred)){
        cat("Need QSPRpred object\n")
        return(NULL)
      }else{
        qsprpred <<- v_qsprpred
      }
      
      if(is.null(v_m)){
        m <<- length(engram$mat)
      }else{
        m <<- v_m
      }
      
      v_ESSth <<- ESSth
      
      if(is.null(smis)){
        esmi <<- lapply(rep("c1ccccc1O", 100), function(x) Esmi$new(x, m=m, engram, type=gentype))
      }else{
        esmi <<- lapply(smis, function(x) Esmi$new(x, m=m, engram, type=gentype))
      } 
      
      oldesmi <<- sapply(esmi, function(x) x$copy())
      
      csmi <- sapply(esmi, function(x) x$get_validsmi())
      oldscore <<- qsprpred$inverse_predx(csmi, temp)
      weights <<- rep(1, length(esmi))
      fkernel <<- rep(1, length(esmi))
      bkernel <<- rep(1, length(esmi))
      isvalid <<- rep(TRUE, length(esmi))
      v_decay <<- decay
      v_maxstock <<- maxstock
      temp <<- temp
      smistock <<- character(0)
      scorestock <<- numeric(0)
    },
    
    smcexec = function(niter, nsteps=5, preorder=0, nview=0){
      "modify chemical structures with niter SMC updates"
      for(i in 1:niter){
        localmove(nsteps)
        update_particles()
        reordering(preorder)
        if(nview>0){
          idx <- which(isvalid==T)
          viewstr(idx[1:nview])
        }
      }
    },
    
    localmove = function(nstep){
      oldesmi <<- sapply(esmi, function(x) x$copy())
      for(i in 1:length(esmi)){
        if(i<length(esmi)){
          cat("\rlocal move for ", i, "th molecules")
        }else{
          cat("\rlocal move for ", i, "th molecules\n")
        }
        u <- rbinom(1, nstep, 0.5)
        steps <- c(rep(1, min(u,length(esmi[[i]]$vstr)-2)), rep(0, nstep-u))
        suppressWarnings(tryres <- try(tempk <- esmi[[i]]$chem_local(engram, steps, 1), silent=T))
        if(class(tryres)=="try-error"){
          fkernel[i] <<- 1
          bkernel[i] <<- 1
          isvalid[i] <<- F
        }else{
          fkernel[i] <<- 1
          bkernel[i] <<- 1
          isvalid[i] <<- TRUE
        }

      }
    },
    
    reordering = function(prob){
      idxs <- which(sapply(esmi, function(x) nchar(x$get_validsmi()))>1)
      idxs <- intersect(idxs, which(sapply(esmi, function(x) x$numn+x$bcount)==0))
      idxs <- intersect(idxs, which(isvalid==T))
      idxs <- sample(idxs, round(length(idxs)*prob, 0))
      for(idx in idxs){
        natoms <- esmi[[idx]]$get_natoms()
        prevsmi <- esmi[[idx]]$get_validsmi()
        if(.Platform$OS=="windows"){
          suppressWarnings(newsmi <- system(paste('obabel -:"', prevsmi ,'" -osmi -xf ', sample(natoms, 1) , sep=""), intern=T, show.output.on.console = F, ignore.stderr=T))
        }else{
          suppressWarnings(newsmi <- system(paste('obabel -:"', prevsmi ,'" -osmi -xf ', sample(natoms, 1) , sep=""), intern=T, ignore.stderr=T))
        }
        newsmi <- substr(newsmi[1], 1, nchar(newsmi)-1)
        if(idx!=idxs[length(idxs)]){
          cat("\rreordering prev:", prevsmi, "=> new:", newsmi, paste(rep(" ", 40), collapse=" "))
          flush.console()
        }else{
          cat("\rreordering prev:", prevsmi, "=> new:", newsmi, paste(rep(" ", 40), collapse=" "), "\n\n")
          flush.console()
        }
        suppressWarnings(reschk <- try(tempesmi<- Esmi$new(smi=newsmi, m=m, engram), silent=T))
        if(class(reschk)!="try-error"){ # exception for reordering        
          temp_vstr <- tempesmi$vstr
          idxc <- grep("=\\[S", temp_vstr)
          if(length(idxc)>0) temp_vstr[idxc] <- sapply(idxc, function(x) paste(substr(temp_vstr[x], 1, 1), substr(temp_vstr[x], 3, 3), sep=""))
          idxc <- grep("=\\[C", temp_vstr)
          if(length(idxc)>0) temp_vstr[idxc] <- sapply(idxc, function(x) paste(substr(temp_vstr[x], 1, 1), substr(temp_vstr[x], 3, 3), sep=""))
          idxc <- grep("=\\[N", temp_vstr)
          if(length(idxc)>0) temp_vstr[idxc] <- sapply(idxc, function(x) paste(substr(temp_vstr[x], 1, 1), substr(temp_vstr[x], 3, 3), sep=""))
          idxc <- grep("=\\[O", temp_vstr)
          if(length(idxc)>0) temp_vstr[idxc] <- sapply(idxc, function(x) paste(substr(temp_vstr[x], 1, 1), substr(temp_vstr[x], 3, 3), sep=""))
          idxc <- grep("\\[S", temp_vstr)
          if(length(idxc)>0) temp_vstr[idxc] <- sapply(idxc, function(x) substr(temp_vstr[x], 2, 2))
          idxc <- grep("\\[C", temp_vstr)
          if(length(idxc)>0) temp_vstr[idxc] <- sapply(idxc, function(x) substr(temp_vstr[x], 2, 2))
          idxc <- grep("\\[N", temp_vstr)
          if(length(idxc)>0) temp_vstr[idxc] <- sapply(idxc, function(x) substr(temp_vstr[x], 2, 2))
          idxc <- grep("\\[O", temp_vstr)
          if(length(idxc)>0) temp_vstr[idxc] <- sapply(idxc, function(x) substr(temp_vstr[x], 2, 2))
          tempesmi$vstr <- temp_vstr
          tempesmi$vstr2smi()
          
          
          suppressWarnings(reschk <- try(tempesmi <- Esmi$new(smi=tempesmi$get_validsmi(), m=m, engram), silent=T))
          if(class(reschk)!="try-error"){ # exception for reordering     
            ar <- 1
            if(runif(1, 0, 1) < ar){
              esmi[[idx]] <<- tempesmi 
              weights[idx] <<- weights[idx]*ar 
            }
          }
        }
      }
      weights <<- weights/sum(weights)
    },
    
    update_particles = function(){
      csmi <- sapply(esmi, function(x) x$get_validsmi())
      isvalid[which(sapply(csmi, function(x) class(try(parse.smiles(x, kekulise=T), silent=T))!="list"))] <<- F
      idx <- which(isvalid==T)
      mols <- parse.smiles(csmi[idx], kekulise=T)
      csmi[idx] <-  sapply(mols, function(x) get.smiles(x, aromatic=T))
      
      newscore <- numeric(length(csmi))
      newscore[idx] <- qsprpred$inverse_predx(csmi[idx], temp)
      newscore[-idx] <- 0

      smistock <<- c(smistock, csmi[idx])
      scorestock <<- c(scorestock, newscore[idx]) 
      or <- order(scorestock, decreasing=T)
      smistock <<- smistock[or[1:min(length(or),v_maxstock)]]
      scorestock <<- scorestock[or[1:min(length(or),v_maxstock)]]
      
      weights <<- weights*(newscore/oldscore)
      weights[which(oldscore==0)] <<- 0
      weights <<- weights/sum(weights)
      
      ESS <- 1/sum(weights^2)
      if(ESS < v_ESSth*length(esmi)){
        cat("weight update done, ESS=", ESS, "\n")
        idx <- sample(1:length(esmi), prob=weights, replace=T)
        esmi <<- lapply(esmi[idx], function(x) x$copy())
        weights <<- rep(1/length(esmi), length(esmi))
        oldscore <<- newscore[idx]
        cat("resampling done\n\n")
      }else{
        cat("weight update done, ESS=", ESS, "\n\n")
        oldscore <<- newscore
      }
      temp <<- temp^v_decay
    },
    
    get_smiles = function(){
      "get SMILES strings from the SmcChem object (same as get_smiles function) "
      csmi <- sapply(esmi, function(x) x$get_validsmi())
      idx <- which(isvalid==T)
      mols <- parse.smiles(csmi[idx], kekulise=T)
      csmi[idx] <-  sapply(mols, function(x) get.smiles(x, aromatic=T))
      csmi
    },
    
    get_hiscores = function(nsmi=100, exsim=0.8){
      "get chemical structures with high QSPR score from SmcChem object (same as get_hiscores function) "
      SMILES <- character(0)
      QSPRScore <- numeric(0)
      tsmi <- smistock
      tscore <- scorestock
      fp <- qsprpred$descriptor(paste(tsmi), fpnames=qsprpred$fpnames)
      cat("\n")
      j <- 1
      while((j<=nsmi) & (length(tsmi)>2)){
        SMILES <- c(SMILES, tsmi[1])
        QSPRScore <- c(QSPRScore, tscore[1])
        sim <- apply(fp[-1,], 1, function(x) sum(x*fp[1,])/sum((x+fp[1,])>=1))
        fp <- fp[-1,]
        tsmi <- tsmi[-1]
        if(length(tsmi)>1){
          if(sum(sim>=exsim)>0){
            fp <- fp[-which(sim>=exsim),]
            tsmi <- tsmi[-which(sim>=exsim)]
            tscore <- tscore[-which(sim>=exsim)]
          }
        }
        cat("\r", j, "th molecules is chosen")
        j <- j + 1
      }
      cat("\n")
      return(cbind(SMILES, QSPRScore))
    },
    
    viewstr = function(idx){
      "view 2D structures from SMILES string vector with index idx (same as viewstr function) "
      tsmi <- sapply(smchem$esmi, function(x) x$get_validsmi())[idx]
      mols <- parse.smiles(tsmi, kekulise=F)
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
  )
)