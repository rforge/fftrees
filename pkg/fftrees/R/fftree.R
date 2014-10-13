######################################
# Class for Fftree, builds upon fftcue 
# 28.04.2014; Marc Giesmann
######################################

# A fftcue includes
# -> cue with criterion
# -> list of cues

#' Class fftree
#'
#' Class \code{fftree} is a tree, consisting of \code{cues} and a \code{criterion}.
#'
#' @name fftree-class
#' @rdname fftree-class
#' @exportClass fftree
setClass(Class = "fftree", 
         slots = c( 
                    name            = "character", # name
                    call            = "ANY",
                    call.formula    = "ANY",
                    criterion       = "cue",       # criterion
                    fftcues         = "list",      # list of fftcues
                    result          = "logical",   # predictions
                    
                    resHi           = "numeric",   # vector: hits at cue index
                    resFa           = "numeric",   # vector: false alarms at cue index
                    resMi           = "numeric",   # vector: misses at cue index
                    resCr           = "numeric",   # vector: correctRej at cue index
                    
                    endHi           = "numeric",   # integer: Final tree hits
                    endFa           = "numeric",   # integer: Final tree false alarms
                    endMi           = "numeric",   # integer: Final tree misses
                    endCr           = "numeric",   # integer: Final tree correct rejections
                    
                    restHi          = "numeric",   # integer: Hits of "REST"
                    restFa          = "numeric",   # integer: false alarms of "REST"
                    restMi          = "numeric",   # integer: misses of "REST"
                    restCr          = "numeric",   # integer: correct rejections of "REST"
                    
                    classcounter    = "numeric",   # vector: how many have been calculated by cue?
                    rest            = "numeric",   # vector: how many have still to be predicted after
                    potential       = "numeric",   # vector: how many have still to be predicted before
                    restCumulated   = "numeric",
                    predictedBy     = "numeric",
                    
                    crossvalidationdata = "ANY"
         ),contains = c("logical")
)

#' Creates a Fast- and Frugal tree
#'
#' A \code{fftree} is an S4 object, that represents a Fast- and Frugal tree.
#' Fast- and Frugal trees use \code{cues} to predict data. It's success is measured
#' at a \code{criterion}, which is also defined as a cue.
#'
#' @name fftree
#' @rdname fftree-class
#'
#' @param criterion cue
#' @param fftcues list with cues. The order of the list represents the hierarchy of the \code{cues}
#' @param treename charactervector, name of the tree (just for description, without function).
#' 
#' @return fftree S4 object
#' @seealso \code{\link{cue}}
Fftree <- function(criterion, fftcues, treename = ""){
  logi <- vector(mode="logical",length=length(criterion@result))
  tree <- new("fftree",criterion = criterion, fftcues = fftcues, name = treename,result = logi)
  return ( updateFftree(tree) )
}

#' Checks if a specific \code{cue} in \code{fftree} is trivial
#'
#' The last cue is defined as trivial, if the last cue could replace the
#' second to last, without changing tree- outcomes.
#'
#' @param tree \code{\link{fftree}}
#' @param cl Optional list with cues. Default is the cue list of the given tree
#' @param n Optional, which cue of tree should be analysed? Default is last cue
#'  
#' @return logical; TRUE if last cue is trivial, else FALSE
#' @seealso \code{\link{isCueLazy}}, \code{\link{isCuePrunable}}
isCueTrivial <- function(tree, cl=tree@fftcues, n = length(cl)){  
  b <- n - 1
    
  if(n <= 1){
    return(FALSE)
  }
  
  
  if(cl[[n]]@pred != cl[[b]]@pred){
    return(FALSE)
  }
  
  bin <- (tree@predictedBy >= b & tree@predictedBy <= n)
  npred <- cl[[n]]@predicts
  
  return(all(npred[bin]))
  
}

#' Checks if a specific \code{cue} in \code{fftree} is lazy
#'
#' A cue is defined as \code{lazy}, if it's \code{cPot} is zero. It can be \code{pruned}
#' without changing tree- outcomes.
#'
#' @param tree \code{\link{fftree}}
#' @param cl Optional list with cues. Default is the cue list of the given tree
#' @param n Optional, which cue of tree should be analysed? Default is last cue
#'  
#' @return logical; TRUE if last cue is lazy, else FALSE
#' @seealso \code{\link{getCueEfficiency}}, \code{\link{isCueTrivial}}, \code{\link{isCuePrunable}}
isCueLazy <- function(tree, cl=tree@fftcues, n = length(cl)){
  return(tree@classcounter[n] == 0)
}

#' Checks if a specific \code{cue} in \code{fftree} is prunable
#'
#' A cue is defined as \code{prunable}, if it's \code{lazy} or it's follower is \code{trivial} 
#'
#' @param tree \code{\link{fftree}}
#' @param asMatrix optional logical. If set to TRUE a data.frame is returned, with all information about possible prunings
#' @param n vector or single value
#'  
#' @return logical vector/data.frame
#' @seealso \code{\link{getCueEfficiency}}, \code{\link{isCueTrivial}}, \code{\link{isCueLazy}}
isCuePrunable <- function(tree, asMatrix = FALSE, n=1:length(tree@fftcues)){
  #n <- length(tree@fftcues)
  trivial <- sapply(n, isCueTrivial, tree=tree, cl=tree@fftcues)
  lazy    <- sapply(n, isCueLazy   , tree=tree, cl=tree@fftcues)
  
  #Create logical array, with FALSE.
  prunable <- replicate(length(n), FALSE)
  
  #If a cue is trivial, the cue before the trivial cue is prunable.
  #To accomplish that, get the indexes of TRUE and set them one index lower. Because the first
  #index cannot be trivial, a shift from 1 to 0 isn't possible
  prunable[ which(trivial) - 1 ] <- TRUE

  #Now check for cPotential == 0.
  prunable <- (prunable | lazy)
  
  if(asMatrix){
    prunable <- data.frame(triv. = trivial, lazy = lazy, prunable = prunable, check.rows = F, check.names = F,stringsAsFactors = F)
  }
  
  return(prunable) 
}

######################
#  Extended generic methods from sdt, needs sdt.R
######################

#' Converts \code{fftree} to a \code{SDT} vector
#'
#' The Sdt- function needs \code{hits}, \code{false alarms}, \code{correct rejections} and \code{misses}.
#' By comparing the \code{prediciton} aka the \code{result} of the tree with the criterion,
#' these values can calculated.
#' @name Sdt.fftree

#' @param tree fftree
#' 
#' @return sdt vector
#' @seealso \code{\link{Sdt}}
Sdt.fftree <- function(tree){
  return(Sdt.default(tree@endHi, tree@endFa, tree@endMi, tree@endCr))
}

##### Internal convenience methods

#' Converts \code{fftree} to a \code{data.frame}
#'
#' @name getTreeDevelopment
#' @param tree fftree
#' 
#' @return data.frame
#' @seealso \code{\link{Sdt}}, \code{\link{fftree}}
getTreeDevelopment <- function(tree){
  cuelist <- tree@fftcues
  names   <- vector(mode = "numeric", length(tree@fftcues))
  
  sdt.matrix <- matrix(nrow = length(tree@fftcues), 
                       ncol = length((Sdt(1,0,0,0))),
                       data = NA)
  
  colnames(sdt.matrix) <- attr(Sdt(1,2,3,4) , "names")
  
  #Show development, by showing SDT-Data for all
  #cue- steps
  for(i in 1:length(tree@fftcues)){
    tr <-Fftree(tree@criterion,cuelist[1:i])
    
    sdt.matrix[i,] <- (Sdt(tr))
    names[i]       <- toString( cuelist[[i]] )
  }
  
  sdt.matrix <- cbind(names = names, as.data.frame(sdt.matrix))
  return(sdt.matrix)
}

#' Returns a vector with the cue efficiency values of a specific cue
#'
#' @name getSingleCueEfficiency
#' @param tree fftree
#' @param id integer. Specifies which \code{cue} within \code{fftree} is meant. Default: Last cue.
#' @param includeRest logical: Should "rest" be added to last cue (because last cue automatically predicts rest)? Default: TRUE
#' @return Vector with "cEff", "triv.", "lazy", "prunable"
#' @details All values are cue-related. 
#' @seealso \code{\link{Sdt},\link{fftree}, \link{getCueEfficiency}}
getSingleCueEfficiency <- function(tree, id = length(tree@fftcues), includeRest = TRUE){
  
  if(includeRest && id == length(tree@fftcues)){
    cHi <- tree@restHi + tree@resHi[id]
    cCr <- tree@restCr + tree@resCr[id]
    
    counter <- tree@restCumulated + tree@classcounter[id]
  }else{
    cHi <- tree@resHi[id]
    cCr <- tree@resCr[id]
    
    counter <- tree@classcounter[id]
  }
  
  
  cEff <-  (cHi + cCr)/counter
  lazy <- isCueLazy(tree,n = id)
  triv <- isCueTrivial(tree,n = id)
  prun <- (lazy || triv)
  
  
  return(c(cEff = cEff, triv. = triv, lazy = lazy, prunable = prun))
}


#' Converts \code{fftree} to a \code{data.frame} with it's efficiency values
#'
#' @name getCueEfficiency
#' @param tree fftree
#' @param hideRest logical: Should "rest" be hidden? Default FALSE.
#' @param includeRest logical: Should "rest" be added to last cue (because last cue automatically predicts rest)
#' @return data.frame
#' @details All values are cue-related. \code{fPot} is the number of unpredicted criterions, \code{cPot} is the quantity of unpredicted criterions WITHIN the spectrum of the cue, \code{boldness} represents the cPot/fPot ratio, \code{cEff} the correct predictions/cPot ratio. All columns can be used by \code{\link{chase}}.
#' @seealso \code{\link{Sdt},\link{fftree}}
getCueEfficiency <- function(tree, hideRest = FALSE, includeRest = FALSE){
  sdt.matrix <- data.frame(stringsAsFactors=FALSE,check.names = F,check.rows = F,
    names   = sapply(tree@fftcues, toString),
    cHi     = tree@resHi,                               #column "potential"
    cFa     = tree@resFa,                               #column "potential"
    cMi     = tree@resMi,                               #column "potential"
    cCr     = tree@resCr,                               #column "potential"
    
    fPot     = tree@potential,                          #column "potential"
    cPot     = tree@classcounter,                       #column "predicted"
    boldness = tree@classcounter/tree@potential,        #column "boldness"
      
    cEff = (tree@resHi + tree@resCr)/tree@classcounter  #column "efficiency"
  )
  
  #Add "prunable" columns
  sdt.matrix <- cbind(sdt.matrix,  isCuePrunable(tree, TRUE))
  
  if(includeRest){
    nr <- nrow(sdt.matrix)

    sdt.matrix[nr, 'cHi']  <- sdt.matrix[nr, 'cHi'] + tree@restHi
    sdt.matrix[nr, 'cFa']  <- sdt.matrix[nr, 'cFa'] + tree@restFa
    sdt.matrix[nr, 'cMi']  <- sdt.matrix[nr, 'cMi'] + tree@restMi
    sdt.matrix[nr, 'cCr']  <- sdt.matrix[nr, 'cCr'] + tree@restCr
    sdt.matrix[nr, 'boldness'] <- 1
    sdt.matrix[nr, 'cEff'] <- (sdt.matrix[nr, 'cHi'] + sdt.matrix[nr, 'cCr']) / sdt.matrix[nr, 'fPot']
  }
  
  if(!hideRest){
    #Now some dirty stuff. Last row should be "rest", so the last cues efficiency won't "suffer" from the rest-predictions.
    sdt.matrix[nrow(sdt.matrix)+1,2:(ncol(sdt.matrix)-3)] <- c(tree@restHi,tree@restFa,tree@restMi,tree@restCr,tree@restCumulated,tree@restCumulated,1,(tree@restHi+tree@restCr)/tree@restCumulated)
    sdt.matrix[nrow(sdt.matrix),'names'] <- "rest"
  }

  return(sdt.matrix)
}

prune <- function(x,...) UseMethod("prune")

#' Prunes a \code{cue} from given \code{fftree}
#'
#' @name prune
#' @param tree fftree.
#' @param n integer/integer vector. Index of \code{cue}, which should be pruned, If set to "NULL", all prunable cues will be pruned automaticly. See \code{\link{isCuePrunable}} for details.
#' @return fftree.
#' @seealso \code{\link{cue}},\code{\link{fftree}}, \code{\link{isCuePrunable}}, \code{\link{isCueLazy}}, \code{\link{isCueTrivial}}
prune.fftree <- function(tree, n = NULL){
  if(is.null(n)){
    n <- which(isCuePrunable(tree))
  }
  
  if(length(n) == 0){
    return(tree)
  }
  
  tree@fftcues <- tree@fftcues[-n]
  
  return(update(tree))
}

#' Crossvalidates \code{fftree} with another data- sample
#'
#' @name cv
#' @param tree fftree.
#' @param data dataframe or list with dataframes  similar to data the tree was fitted on
#' @return list with 2 items. "MAE.leafs": Errors of leafs (exit structures of the fftree) crossvalidated predictions and "mean absolute error". "MAE.overall": Errors of fftree crossvalidated predictions and "mean absolute error".
#' @seealso \code{\link{cue}},\code{\link{fftree}}, \code{\link{isCuePrunable}}, \code{\link{isCueLazy}}, \code{\link{isCueTrivial}}
crossvalidate <- function(tree, data=tree@crossvalidationdata){
  
  if(!is.list(data)){
    data <- list(data)
  }
  
  n.cues <- length(tree@fftcues)
  n.data <- length(data)
  
  #get cue names
  cue.names <- getCueEfficiency(tree)[,'names']
  
  #temporary list for cues
  clist <- vector(mode="list", length=n.cues)
  
  #temporary list for leaf- errors
  leaf.error    <- vector(mode="list", length=n.data)
  overall.error <- vector(mode="list", length=n.data)
  
  for(i in 1:n.data){
    current.data <- data[[i]]
    
    #build temporary tree
    temptree <- update(tree, current.data)
    
    leaf.error[[i]]    <- 1 - as.vector(getCueEfficiency(temptree,hideRest=F,includeRest=F)[,'cEff'])
    overall.error[[i]] <- 1 - as.vector(Sdt(temptree)['percCorr'])
  }
  
  leaf.error.table    <- do.call(cbind, leaf.error)
  leaf.mae.table      <- apply(leaf.error.table,MARGIN=1, sum, na.rm = TRUE) / n.data
  leaf.error.table    <- data.frame(names=cue.names, leaf.error.table, MAE=leaf.mae.table)
  colnames(leaf.error.table) <- c("names", paste0("err.",1:n.data, sep=""), "MAE")
  
  overall.error.table <- do.call(cbind, overall.error)
  overall.error.table <- cbind(overall.error.table, apply(overall.error.table,MARGIN=1, sum, na.rm = TRUE) / n.data) 
  colnames(overall.error.table) <- c(paste0("err.",1:n.data, sep=""), "MAE")
  
  ret <- list(MAE.leafs=leaf.error.table, MAE.overall=overall.error.table)
  
  return( ret )
  
}


#' Returns a SDT table for all splitvalues of a specific independentend variable
#'
#' @name splitvalanalysis
#' @param criterion \code{\link{cue}}, which indicates the criterion to be used
#' @param data dataframe or vector. If data is a data.frame, you have to specify the \code{name} of the column provided in the dataframe of \code{data}
#' @param name of independent variable in column
#' @return dataframe with SDT values for all possible splitvalues
#' @seealso \code{\link{cue}}, code{\link{Sdt}}
splitvalanalysis <- function(criterion, data, name = NULL){
  
  if(!is.vector(data)){
    data <- unlist(rawdata[,name])
  }
  
  un   <- sort(unique(data))
  n    <- length(un)
  
  sdtlist <- vector(mode="list",length=n)
  for(i in 1:n){
    currentcue   <- Cue(inputvector=data, name=name, split=un[i])
    sdtlist[[i]] <- c(toString(currentcue), round(Sdt.fftree( Fftree(criterion=criterion, list(  currentcue  ))), digits=5))
  }
  
  sdtlist <- data.frame(do.call(rbind, sdtlist))
  colnames(sdtlist) <- c("cue", names(Sdt(0,0,0,0)))
  
  return(sdtlist)
}

# stabilitytable <- function(tree,steps=1){
#   
#   n <- length(tree@criterion@data)
#   ncriterions <- seq(from=n,by=-steps)
#   sdt.vec <- Sdt(0,0,0,0)
#   sdt.names <- names(sdt.vec)
#   
#   res <- matrix(nrow = length(ncriterions), ncol= (1 + length(sdt.vec)))
#   colnames(res) <- c("criterions",sdt.names)
#     
#   i <- 1
#   while(n >= 1){
#     tree@criterion <- update(tree@criterion, tree@criterion@data[-steps])
#     tree <- updateFftree2(tree, lapply(tree@fftcues, function(x) update(x, x@data[-steps])))
#     
#     res[i,1] <- n
#     res[i,2:ncol(res)] <- Sdt(tree)
#     
#     n <- length(tree@criterion@data)
#     i <- i + 1
#   }
#   
#   return(res)
# }

######################
#  Extended generic methods from base
######################

#' Converts \code{fftree} to a \code{data.frame}
#'
#' @aliases as.data.frame
#' @param x Fast and frugal tree object
#' @param row.names not used
#' @param optional string Can either be "treedevelopment" or "cueefficiency". Default is "treedevelopment"
#' @param ... not used
#' 
#' @return data.frame
#' @seealso \code{\link{Sdt},\link{fftree}}
setMethod("as.data.frame", "fftree", function(x, row.names = NULL, optional = NULL, ...) {
  
  if(is.null(optional)){
    return(getTreeDevelopment(x))
  }
  
  if(optional=="treedevelopment"){
    return(getTreeDevelopment(x))
  }
  
  if(optional=="cueefficiency"){
    return(getCueEfficiency(x))
  }
  
  stop("Unknown dataframe format '",optional,"'.")
})


#' Shows \code{fftree}
#'
#' @name show
#' 
#' @seealso \code{\link{Sdt},\link{fftree}, \link{getCueEfficiency}, \link{getTreeDevelopment}}
setMethod("show", "fftree", function(object) {
  cat("Cue efficiency: \n\n")
  show(getCueEfficiency(object))
  
  if(!is.null(object@crossvalidationdata)){
    cat("\n\nMean average error (crossvalidation): \n\n")
    show(crossvalidate(object))
    
  }
  
  cat("\n\nDevelopment (last entry ~ overall tree performance) \n\n")
  show(getTreeDevelopment(object))
})

#' Predicts data with \code{fftree}
#'
#' @param object \code{\link{fftree}}
#' @param newdata OPTIONAL data.frame. Columns of \code{newdata} have to match the \code{names} of the \code{\link{cue}}s in the \code{\link{fftree}}. If \code{newdata} is \code{NULL}, the fitting values of the \code{fftree} are shown. 
#' @param type OPTIONAL string. Has to be "response" or "class". Default is "response".
#' 
#' @seealso \code{\link{fftree}}
 setMethod("predict", "fftree", function(object, newdata=NULL,  type="response", ...) {
   if(!is.null(newdata)){
     newdata <- data.frame(newdata)
    
     cnames  <- sapply(object@fftcues, function(x) x@name)
     ucnames <- unique(cnames)
     columns <- ucnames %in% colnames(newdata)
     
     if(ncol(newdata) < length(ucnames)){
       stop("Insufficient parameters! For this tree you need ", length(ucnames), " parameters:",paste0(ucnames, sep=', '))
     }
     
     if(!all(columns)){
       stop("All cues have to have a name that matches one of the columns of newdata! Not found: ",paste0(ucnames[!columns], sep=', '))
     }
     
     result    <- vector(length=nrow(newdata))
     predicted <- replicate(nrow(newdata),FALSE)
     effitbl   <- getCueEfficiency(object,hideRest=T,includeRest=T)
     effi      <- matrix(ncol=2,nrow=nrow(newdata))
     colnames(effi) <- c(0, 1)
     
     for(i in 1:nrow(newdata)){
       j <- 0
       for(jcue in object@fftcues){
         j <- j + 1
         p <- predict(jcue, unlist(newdata[i, jcue@name]))
         
         if(p == jcue@pred){
           result[i]    <- p
           predicted[i] <- T
           
           if(jcue@pred){
             effi[i, 2] <- effitbl[j, 'cEff']
             effi[i, 1] <- 1 - effi[i, 2]
           }else{
             effi[i, 1] <- effitbl[j, 'cEff']
             effi[i, 2] <- 1 - effi[i, 1]
           }
           
           break
         }
       }
     }
     
     #set all not yet predicted results to the opposite of last cues prediction
     n <- length(object@fftcues)
     p <- !object@fftcues[[ n ]]@pred
     e <- effitbl[n, 'cEff']
     
     result[!predicted] <- p
     effi[!predicted, ] <- ifelse(p, c(1-e, e), c(e, 1-e))
     
   }
   if(type=="response"){
    return(effi)
   }
   
   if(type=="class"){
     return(result)
   }
   
   stop("Unknown type '", type, "'. Use 'class' or 'response'.")
   
 })


######################
#  Extended generic methods from stats
######################

#' Updates/Refreshes \code{\link{fftree}} 
#' 
#' Updates/Refreshes \code{\link{fftree}} internal caches and results, like predictions, hit-/false alarm- counters, rest-data etc.
#'
#' @name update
#' @param object fftree
#' @param data optional dataframe to update the model on
#' @seealso \code{\link{fftree}}
setMethod("update", "fftree", function(object,data = NULL) {
  
  if(!is.null(data)){
    n.cues <- length(object@fftcues)
    
    object@criterion <- update(object=object@criterion, data=unlist(data[, object@criterion@name]))
    
    for(i in 1:n.cues){
      cue  <- object@fftcues[[i]]
      name <- cue@name
      
      object@fftcues[[i]] <- update(object=cue, data=unlist(data[, name])) 
    }
  
  }
  return(updateFftree(object))
})


#Plots the tree
#' Plots \code{fftree}
#'
#' @name plot.fftree
#' @aliases plot
#' @param object fftree
#' @param Further parameter.
#' 
#' @seealso \code{\link{fftree}}
plot.fftree <- function(object,...) {   
  
  fftree <- object
  
  names               <- sapply(fftree@fftcues, function(x) toString(x,suppressNameBrackets = T,suppressThen = T, suppressName = F, suppressTest = F))
  names.maxlen.index  <- which.max(nchar(names))
  
  binarys        <- sapply(fftree@fftcues, function(x) x@pred)
  binarys.string <- sapply(binarys, toString)
  n              <- length(fftree@fftcues)
  
  #------------------------------------------
  # config of stuff
  #------------------------------------------+
  offset.x <- 1
  offset.y <- 1
  
  min.y <- 0
  max.y <- (-1-n) * offset.y
  
  min.x     <- 0
  max.x     <- +offset.x
  
  vec.binarys.x <- NULL
  vec.binarys.y <- NULL
  vec.cues.x    <- NULL
  vec.cues.y    <- NULL
  
  current.x <- min.x
  current.y <- min.y
  
  for(b in binarys){
    vec.cues.x <- c(vec.cues.x, current.x)
    vec.cues.y <- c(vec.cues.y, current.y)
    
    #The next binary will be one step deeper than the current cue
    vec.binarys.y <- c(vec.binarys.y, current.y - offset.y)
    
    if(b){
      #Set the binary to the opposite: in this case to the left
      vec.binarys.x <- c(vec.binarys.x, current.x - offset.x)
      
      #And wander to the right
      current.x <- current.x + offset.x      
    }else{
      #Set the binary to the opposite, in this case to the right
      vec.binarys.x <- c(vec.binarys.x, current.x + offset.x)
      
      #And wander to the left
      current.x <- current.x - offset.x
    }
    
    #After every step, we need to decrease the current.y value
    current.y <- current.y - offset.y
    
  }
  
  #The last cue needs to have 2 exits.
  if(b){
    vec.binarys.x <- c(vec.binarys.x, current.x)
    binarys.string <- c(binarys.string, "FALSE")
    binarys        <- c(binarys,FALSE)
  }else{
    vec.binarys.x <- c(vec.binarys.x, current.x)
    binarys.string <- c(binarys.string, "TRUE")
    binarys        <- c(binarys,TRUE)
  }
  
  #And the last binary
  vec.binarys.y <- c(vec.binarys.y, current.y)
  
  
  #Calculate Boxes height and length:
  cue.max.width   <- strwidth (  names[names.maxlen.index]  ,units = "inches") - .5
  cue.max.height  <- strheight(  names[names.maxlen.index]  ,units = "inches") + .35
  bin.max.width   <- strwidth (  "FALSE"  ,units = "inches")
  bin.max.height  <- strheight(  "FALSE"  ,units = "inches") + .2
  
  plo <- qplot(geom="blank", x = vec.cues.x, y = vec.cues.y)
  plo <- plo + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
                     axis.line=element_blank(),
                     axis.text.x=element_blank(),
                     axis.text.y=element_blank(),
                     axis.ticks=element_blank(),
                     axis.title.x=element_blank(),
                     axis.title.y=element_blank())
  
  for(i in 1:length(vec.cues.x)){
    from.x <- vec.cues.x[i]
    from.y <- vec.cues.y[i]
    
    to.x <- vec.binarys.x[i]
    to.y <- vec.binarys.y[i]
    
    #Connect cue with binary
    plo <- plo + annotate("line", x=c(from.x,to.x), y=c(from.y,to.y), lty=2, col="black") 
    
    if(i != length(vec.cues.x)){
      plo <- plo + annotate("line", x=c(from.x,vec.cues.x[i+1]), y=c(from.y,vec.cues.y[i+1]), lty=1, col="black") 
    }else{
      plo <- plo + annotate("line", x=c(from.x,vec.binarys.x[i+1]), y=c(from.y,vec.binarys.y[i+1]), lty=1, col="black")     
    }
  }
  plo <- plo +   geom_rect(data = data.frame(names), aes(xmin = vec.cues.x - cue.max.width/2 , xmax = vec.cues.x + cue.max.width/2, 
                                            ymin = vec.cues.y - cue.max.height/2, ymax = vec.cues.y + cue.max.height/2), fill = "grey80") + 
    geom_text(data = data.frame(names), aes(x = vec.cues.x, y = vec.cues.y, label = names), size = 4) +  
    
    
    
    annotate("rect", xmin = vec.binarys.x - bin.max.width/2 , xmax = vec.binarys.x + bin.max.width/2, 
                     ymin = vec.binarys.y - bin.max.height/2, ymax = vec.binarys.y + bin.max.height/2, fill = "snow") + 
    
    annotate("text",x = vec.binarys.x, y = vec.binarys.y, label = binarys.string, size = 4)
    
  return(plo)   
}


