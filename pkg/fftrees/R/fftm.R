######################################

# Fast- and frugal- tree model object
#.Data <- ""

# if(!require(Combinations)){
#   install.packages("Combinations", repos = "http://www.omegahat.org/R")
#   require(Combinations)
# }

# source("./permutation.R",chdir = TRUE)
# source("./sdt.R"        ,chdir = TRUE)
# source("./cue.R"        ,chdir = TRUE)
# source("./fftree.R"     ,chdir = TRUE)


#communicate new class with r
# setClass(Class = "fftm", 
#          slots = c( form       = "formula",
#                     criterion  = "cue",
#                     tree       = "fftree",
#                     results    = "list")
#)



######################
#  Generic methods
######################

#' Wrapper function fftm
#'
#' @name Fftm
#' @param formula an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted. The details of model specification are given under ‘Details’.
#' @seealso \code{\link{Fftm.default}}
Fftm <- function(formula,...) UseMethod("Fftm")

#' Fitting  Fast- and Frugal tree models
#'
#' \code{Fftm} is used to fit Fast- and Frugal tree models, specified by giving a 
#' symbolic description of the predictors and the criterion.
#'
#' @name Fftm
#'
#' @param method stringvector which specifies the method. \link{bruteforcemaximize}, \link{maximize}, \link{montecarlo} or \link{chase}
#' @param data Dataframe containing the variables in the model.
#' @param criterion \code{\link{cue}}; OPTIONAL: If criterion in \code{data} isn't binary, use a predefined criterion cue.
#' @param autoprune logical: If set to true \code{\link{prune}} is called before returning the tree or all trees of the list. Default: TRUE
#' @param xvalidate.parts integer, how many party of \code{data} should be generated (odd-even method) to crossvalidate the model
#' @param sdtfun optional function with 4 parameters (\code{hits}, \code{falsealarms}, \code{misses}, \code{correctrejections}). Has to return a named vector. The returned value can be maximized.
#' @param cuelist optional list with cues. The \code{name}s of the cues have to match the columns provided in \code{data}.
#' @param ... Additional parameter which should be passed to called \code{method}.
#' @return ANY; Depending on used method, return value can either be a \code{\link{list}}, or a \code{\link{fftree}}.
#' @seealso \code{\link{bruteforcemaximize}}, \code{\link{maximize}}, \code{\link{chase}}, \code{\link{montecarlo}}
Fftm.default <-  function(formula, method, data, criterion=NULL, autoprune=TRUE, xvalidate.parts = 2,sdtfun=NULL, cuelist = NULL,...){
  
  call.full    <- match.call()
  call.formula <- toString(formula)
  
  
  
  #If custom function for Sdt is provided, override Sdt.default and Sdt.fftree manually.
  if(!is.null(sdtfun)){
    Sdt.old <- Sdt.default 
    Sdt.default <- function(hi,fa,mi,cr){
      return(c(Sdt.old(hi,fa,mi,cr), sdtfun(hi,fa,mi,cr)))
    }
    
    Sdt.fftree <- function(x) Sdt.default(x@endHi,x@endFa,x@endMi,x@endCr)
  }
  
  ##### INPUT
  
  #Extract everything
  allv  <- rownames(attr(terms(formula),"factors"))
  terms <- allv[2:length(allv)]
  crit  <- allv[1]
  
  
  ##### INPUT VALIDATION
  #Look, if data hast colnames
  if(length(which(allv %in% colnames(data))) != length(allv)){
    
    faultyTerms <- terms[!(allv %in% colnames(data))]
    err <- ""

    for(t in faultyTerms){
      t <- paste("'", t, "'", sep = "")
      
      if(err != ""){
        sep <- " or "
      }else{
        sep <- ""
      }
        err <- paste(err,t, sep = sep)
      
    }
    
    stop("No column(s) with name ", err  , " in data found!")
  }
  
  ####### CREATE CROSSVALIDATION DATAPARTS
  data.original <- data
  if(xvalidate.parts > 1){
    cv.list <- vector(mode="list", length=xvalidate.parts)
    for(i in 1:xvalidate.parts){
      data.rows    <- seq(i, nrow(data), xvalidate.parts)
      cv.list[[i]] <- data.frame(data[data.rows, ])
    }
    data <- cv.list[[1]]
  }else{
   cv.list <- NULL
  }
  
  
  ####### CREATE CRITERION
  if(is.null(criterion)){
    cdat <- as.vector(unlist(data[,crit]))
    splval <- max(cdat)
    
    criterion <- Cue(inputvector=cdat,test="==",split=splval,pred=TRUE,name=crit)
    
    if(length(unique(cdat)) != 2) {
      warning("More than 2 unique values of '",crit, "' found! Check, if your generated criterion is valid!",immediate.=TRUE)
    }
    
    message("Automatic generated criterion: ", toString(criterion))
    
  }else if(class(criterion) != "cue"){
    stop("Criterion has to be of class 'cue'! ")
  }

  
  ####### CREATE CUELIST
  if(is.null(cuelist)){
    cuelist   <- vector(mode="list",length=0)
    
    for(t in terms){
      #Append all variants of each cue to cuelist
      cuelist <- c(cuelist, cue.getAllVariants(inputvector = as.vector(unlist(data[,t])), name = t) )
    }
  }
  
  ####### CREATE TREE
  tree <- Fftree(criterion=criterion,fftcues=cuelist,treename=toString(formula))
  tree@crossvalidationdata <- cv.list
  
  res <- NULL
  
  #### START methods
  if(method == "bruteforcemaximize"){
    environment(fftm.bruteforce.maximize) <- environment()
    res <- fftm.bruteforce.maximize(tree, ...)
  }
  
  if(method == "maximize"){
    environment(fftm.maximize) <- environment()
    res <- fftm.maximize(tree, ...)
  }
  
  if(method == "chase"){
    environment(fftm.chase) <- environment()
    res <- fftm.chase(tree, ...)
  }
  
  if(method == "montecarlo"){
    environment(fftm.montecarlo) <- environment()
    res <- fftm.montecarlo(tree, ...)
  }
  
  #Method check
  if(is.null(res)){
    knownmethods <- c("bruteforcemaximize","maximize","chase","montecarlo")
    stop("Don't know method '", method  , "'!' Possible methods: ", toString(knownmethods))
  }
  
  #Autoprune if requested
  if(autoprune){
    if(is.list(res)){
      res <- lapply(res, function(x){
          tree <- prune(x)
          tree@call         <- call.full
          tree@call.formula <- call.formula
          return(tree)
      })
    }else{
      res <- prune(res)
      res@call          <- call.full
      res@call.formula  <- call.formula
    }
  }
  
  #Change back data to original, unsplitted data
  if(is.list(res)){
    res <- lapply(res, FUN = update, data=data.original)
  }else{
    res <- update(object=res, data=data.original)
  }
  
  return(res)
  
}

#' Maximizes a \code{\link{Sdt}} value by bruteforcing the tree permutation
#'
#' Takes given \code{\link{fftree}}, and permutates given \code{\link{cue}s} in all possible variants
#' to find a permutation of cues, which has the highest value of \code{whatToMaximize} in \code{\link{Sdt}}. Careful: This may take a LONG while!
#'
#' @name bruteforcemaximize
#'
#' @param tree \code{\link{fftree}}
#' @param whatToMaximize string; Which parameter of the tree should be maximized? See \link{Sdt} for possible values.
#' @return list with best cues
#' @details This is the most inefficient method, but the most secure. For a more performant approach with similar results, use \link{maximize} instead. Everything from \link{Sdt} can be maximized.
#' @seealso \code{\link{Fftm}}, \code{\link{maximize}}, \code{\link{chase}}, , \code{\link{montecarlo}}
fftm.bruteforce.maximize <- function(tree, whatToMaximize){
  cnames <- names(Sdt(0,0,0,0))
  if(!(whatToMaximize %in% cnames)){
    stop("Can't maximize ",whatToMaximize, ". Use values provided by function Sdt().")
  }
  
  #put all cues in a list to have fastest access to them
  cuelist <- tree@fftcues

  #build workertree with length(cue.list) fftcues
  tree.worker         <- tree
  
  ###CONFIG
  permutation          <- length(cuelist)
  permutation.length   <- gamma(permutation + 1)
  permutation.iterator <- 0
  
  chunk.length   <- 100000
  chunk.iterator <- 0
  
  #if permutation is smaller as requested chunk
  if(permutation.length < chunk.length ){
    chunk.length <- permutation.length
  }
  
  chunk.buffer <- matrix(nrow = chunk.length, ncol = permutation, data=NA)
  
  #Now prepare result-storage
  storage.length <- 10 # ATTENTION: HAS TO BE EVEN!
  storage.matrix <- matrix(ncol=permutation, nrow=storage.length*2)
  storage.vector <- vector(mode="numeric", length = storage.length*2)
  
  #Signal Detection Object
  rating.sdt <- Sdt(0,0,0,0)
  
  #Variables for statistics to show while waiting
  stat.starttime         <- Sys.time()
  stat.time              <- format(Sys.time(), "%X")
  stat.timeelapsed       <- 0
  stat.averagecalcpersec <- 0
  stat.timeneeded        <- 0
  stat.perc              <- 0
  
  #call self-written, inperformant permutation algorithm 
  permu.old(1:permutation, function(current){
    #Iterate
    chunk.iterator       <<- chunk.iterator + 1
    permutation.iterator <<- permutation.iterator + 1
    
    #Fill permutationbuffer
    chunk.buffer[chunk.iterator, ] <<- current
    if(chunk.iterator == chunk.length || permutation.iterator == permutation.length){
      # DO CALCULATIONS
      all <-        apply(X      = chunk.buffer, 
                          MARGIN = 1, 
                          FUN    =  function(p) {
                                      #Set tree with current permutation of cues
                                      #slot(tree.##worker,"fftcues",check=F) <-  cuelist[ p ]
                                      tree.worker <- updateFftree2(tree.worker,cuelist[ p ])
                                                                            
                                      #rate, extract and return requested value (e.g. dPr) from sdt-rating
                                      return(Sdt(tree.worker)[whatToMaximize])
                                    }
                          )
      #Append storage
      vector.from <- storage.length + 1
      vector.to   <- storage.length * 2
      
      ##Get the indexes of the top 50 (storage.length)
      best <- order(all, na.last =TRUE, decreasing=TRUE) [1:storage.length]
      
      #Append together.
      storage.vector[vector.from : vector.to]   <<- all[best]
      storage.matrix[vector.from : vector.to, ] <<- chunk.buffer[best,]
      
      #Sort again and cleanup second buffer part
      best <- order(storage.vector, na.last =TRUE, decreasing=TRUE) [1:storage.length]
      
      storage.vector[1:storage.length]         <<- storage.vector[best]
      storage.vector[vector.from : vector.to]  <<- NA
      
      storage.matrix[1:storage.length, ]       <<- storage.matrix[best, ]
      storage.matrix[vector.from : vector.to,] <<- NA
      
      #when no end (break didn't happen), display statistics about calculation
      stat.time              <<- format(Sys.time(), "%X")
      stat.timeelapsed       <<- as.integer(Sys.time()-stat.starttime)
      stat.averagecalcpersec <<- round((chunk.length)       /stat.timeelapsed, 1)
      stat.timeneeded        <<- round((permutation.length - permutation.iterator  )/stat.averagecalcpersec/60, 1)
      stat.perc              <<- round((permutation.iterator / permutation.length)                  *100,1)
      stat.starttime         <<- Sys.time()
      .logger(stat.perc,"% (",permutation.iterator," of ",permutation.length,") Estimated time left: ",stat.timeneeded," min @",stat.averagecalcpersec," calcs/sec)")
      # END CALCULATIONS
      
      #If next chunk would be bigger than needed (no more permutations left),
      #reconfigure chunk.length
      if((permutation.iterator + chunk.length) > permutation.length){
        chunk.length <<- permutation.length - permutation.iterator
      }
      
      #Reallocate chunk.buffer, reset iterator
      chunk.buffer   <<- matrix(nrow = chunk.length, ncol = permutation, data=NA)
      chunk.iterator <<- 0
    }
    
    
  })
  
  #Build trees
  res <- list()
  i <- 0
  for(i in 1:storage.length){
    permutation <- storage.matrix[i,]
    
    #Filter NA- permutations. This happens, if
    #permutation length was smaller than storage
    if(!length(which(is.na(permutation)))){
      tree.worker@fftcues <- cuelist[ permutation ]
      tree.worker <- updateFftree(tree.worker)
      
      res[[i]] <- tree.worker
    }
  }
  
  return(res)
}


#' Maximizes a \code{\link{Sdt}} value
#'
#' Takes given \code{\link{fftree}}, and permutates given \code{\link{cue}s} in all possible variants
#' to find a permutation of cues, which has the highest value of \code{whatToMaximize}. 
#' Maximize is similar to bruteforce, but much smarter. Still: Be careful: This may take a LONG while!
#'
#' @name maximize
#'
#' @param tree \code{\link{fftree}}
#' @param whatToMaximize string; Which parameter of the tree should be maximized? See \link{Sdt} for possible values.
#' @param savebest numeric; How big should the returning list be? If set to '100', the 'Top 100' of maximized cues will be returned.
#' @param cluster If you have a HPC cluster, specifiy configuration here. See \code{\link{makeCluster}} from \code{snow} package.
#' @param cores If multithreading is wanted (recommended for bigger computations!), you can specify your cores here.
#' @return List with best cues
#' @details This is a bruteforce- variant, to maximize one of \link{Sdt}s values. Normally one should use this method instead of \link{maximize}.
#' @seealso \code{\link{Fftm}}, \code{\link{bruteforcemaximize}}, \code{\link{chase}}, \code{\link{montecarlo}}
#' @examples 
#' rawdata <- fftm.Titanic.data()
#' tree <- Fftm(Survived ~ Age + Sex + Class, "maximize",rawdata, whatToMaximize = "percCorr", cores=4)
fftm.maximize <- function(tree, whatToMaximize, savebest = 1, cluster = NULL, cores = 1){
  
  #put all cues in a list to have fastest access to them;
  #here they are global for our callback functions to see
  cuelist     <- tree@fftcues
  tree.worker <- Fftree(tree@criterion,list())
  #Signal Detection Object
  rating.sdt <- Sdt(0,0,0,0)
  
  if(!(whatToMaximize %in% names(rating.sdt))){
    stop("Can't maximize ",whatToMaximize, ". Use values provided by function Sdt().")
  }
  
  #create callbacks for permutation function
  perm.callback <- function(x,resultobject){
    #Set cues in current permutation, update tree
    #tree.worker <- setCuelist(tree.worker,cuelist[ x ])
    #slot(tree.worker,"fftcues",check=F) <- cuelist[ x ]
    
    #Update tree
    #tree.worker <- updateFftree(tree.worker)
    
    #Version 2
    clist <- cuelist[ x ]
    tree.worker <- updateFftree2(tree.worker, clist)
    
    #Most important: if this cue didn't predict (Cue-Boldbess = 0) anything, then just skip, and don't save.
    #Nothing was done by this one
    if(!tree.worker@classcounter[length(x)]){
      resultobject[CONST_SKIP] <- TRUE
      resultobject[CONST_SAVE] <- FALSE
      return(resultobject)
    }
    
    if(isCueTrivial(tree.worker, clist)){
       resultobject[CONST_SKIP] <- TRUE
       resultobject[CONST_SAVE] <- FALSE
       return(resultobject)
    }
    
    
    
    #If we got this far, this means, this cue did something! Judge it!
    #Get SDT rating
    rating.sdt <- Sdt(tree.worker)
    
    #Get requested Value
    requestedValue <- rating.sdt[whatToMaximize]
    
    #Save, if requested Value is bigger than calculated values before
    if(requestedValue >= resultobject[CONST_VAL]){
      resultobject[CONST_VAL]  <- requestedValue 
      resultobject[CONST_SAVE] <- TRUE
    }
    
    #Skip, if there is nothing left to predict (FFT-Potential = 0)
    resultobject[CONST_SKIP] <- (tree.worker@restCumulated == 0)
    
    return(resultobject)
  }
  
  ##Callback for saving
  perm.callback.save <- function(resultlist){
    res <- resultlist[['resultobject']]
    lis <- resultlist[['results']]
    
    #GET VALUES
    allvals <- sapply(lis, function(x) x[['values']] ) # get all values
    lis     <- lis    [  order( allvals,decreasing=T )  ]      # order list by values
    allvals <- allvals[  order( allvals,decreasing=T )  ]      # order values
    
    if(length(lis) > savebest){
      lis     <- lis[1:savebest]                                 # only save Top 100 (savebest)
      allvals <- allvals[1:savebest]
    }
    
    #The next best "save-worthy" value has to be better than the smallest of the top 100
    res[CONST_VAL] <- min(allvals)
    
    .logger("Step complete. Current best ", whatToMaximize,": ",max(allvals)," Minimum: ", res[CONST_VAL])
    
    return(list(results=lis, resultobject= res))
  }
  
  result <- permu.new(perm    = 1:length(cuelist), 
                      fun     = perm.callback,
                      funsave = perm.callback.save,
                      values  = 1, 
                      cluster = cluster,
                      cores   = cores,
                      savemax = 1000000)
  
  #NOW CREATE TREES
  p        <- lapply(result, function(x) x[['perm']] ) #extract permutations
  
  treelist <- lapply(p, function(x) {
                        tree@fftcues <- cuelist[x]
                        tree <- update(tree)
                        return(tree)
                        })
  
  if(savebest==1){
    treelist <- treelist[[1]]
  }
  
  return(treelist)
}


#' Chases a \code{\link{Sdt}} value
#'
#' Takes all \code{\link{cue}s} in \code{\link{fftree}} and builds a fftree from scratch,
#' by appending the 'next best' or prepending the 'previous best' cue according to given parameters. You can chase multiple parameters by
#' specifying an array of \code{whatToMaximize} and \code{chaseHow} like this: \code{whatToMaximize = c("dPr","percCorrect","crit")}, \code{chaseHow=c("max","max","min"))} and \code{forward=c(TRUE,TRUE,FALSE))}.
#' Use \code{\link{Fftm}} instead of calling this method directly.
#'
#' @name chase
#'
#' @param tree \code{\link{fftree}}
#' @param whatToMaximize stringvector; Which parameters of the tree should be chased? The index of the vector represents the sequence.  See \link{Sdt} for possible values.
#' @param chaseHow stringvector; Every item has to be either \code{max} or \code{min}.
#' @param forward OPTIONAL logicalvector: Should the next test be appended or prepended to the current tree? Default is \code{forward = TRUE}
#' @return \code{\link{Fftree}}
#' @details You can chase all values of \code{\link{Sdt}} and \code{\link{getCueEfficiency}}
#' @seealso \code{\link{Fftm}}, \code{\link{bruteforcemaximize}}, \code{\link{chase}}, \code{\link{montecarlo}}, \code{\link{getCueEfficiency}}, \code{\link{Sdt}}
#' @examples 
#' rawdata <- fftm.Titanic.data()
#' tree <- Fftm(Survived ~ Age, "chase",rawdata, whatToMaximize = "percCorr", chaseHow = "max", forward=TRUE)
#' tree2 <- Fftm(Survived ~ Age, "chase",rawdata, whatToMaximize = c("percCorr","dPr"), chaseHow = c("max","min"), forward=c(TRUE,FALSE))
fftm.chase <- function(tree, whatToMaximize, chaseHow = "max", forward = TRUE){

  #Prepare cluster
  #multithreading wanted
#   if(!is.null(cluster) && cores > 1){
#     stop("You can't define parameter 'cluster' AND 'cores'. If you want local multithreading use 'cores' (> 0), else specify cluster (see 'makeCluster()').")
#   }
#   
#   if(is.null(cluster) && !is.null(cores) ){
#     cluster <- makeCluster(cores, type = "SOCK")
#   }
  
  
  if(length(whatToMaximize) != length(chaseHow) || length(whatToMaximize) != length(forward)){
    stop("Vector chase, chaseHow and forward have to be the same length!")
  }
  
  cnames <- c(names(Sdt(0,0,0,0)), names(getCueEfficiency(tree)))
  if(!(whatToMaximize %in% cnames)){
    stop("Can't chase ",whatToMaximize, ". Use values provided by function Sdt() or getCueEfficiency().")
  }
  
  
  if(whatToMaximize %in% names(Sdt(0,0,0,0))){
      rating <- function(tree, forward){ c(Sdt.fftree(tree),lazy=isCueLazy(tree))  }  
      
    }else{
      rating <- function(tree, forward){ 
        #ceff <- getCueEfficiency(tree,hideRest = T, includeRest = T)
        if(forward){
          
          return(getSingleCueEfficiency(tree))
          
          
          #return(unlist(ceff[nrow(ceff),-1]))
        }else{
          #return(unlist(ceff[1,-1]))
          return(getSingleCueEfficiency(tree,id = 1))
        }
      }
      
    }
  
  #Reset bestlist
  #example.rating    <- rating(tree,F)
  #treesdtlist.cols  <- length(example.rating)
  #treesdtlist.names <- names(example.rating)
  
  cue.list     <- tree@fftcues
  tree@fftcues <- list()
  fftcuelist   <- list()
  
  chaseItem  <- 0 #index of chase and chaseHow
  chaseIndex <- 0
  
  #repeat while list isn't empty and there is still tree potential
  while(length(cue.list > 0) && (tree@restCumulated > 0 || chaseIndex == 0) ){
    chaseItem  <- chaseItem  + 1
    chaseIndex <- chaseIndex + 1
    
    if(chaseItem > length(whatToMaximize)){
      chaseItem <- 1
    }
    
    #Wrap max or min function in "internaltest"
    if(chaseHow[chaseItem] == "min"){
      internaltest <- which.min
    }else if(chaseHow[chaseItem] == "max"){
      internaltest <- which.max
    }else{
      stop(chaseHow[chaseItem], " can't be chased.")
    }
    
    forw <- forward[chaseItem]
    if(forw){
      insert.index             <- chaseIndex
      fftcuelist               <- tree@fftcues
      fftcuelist               <- append(fftcuelist, NA)
    }else{
      insert.index             <- 1
      fftcuelist               <- tree@fftcues
      fftcuelist               <- append(fftcuelist, NA, 0)
    }
    
    ratinglist <- lapply(cue.list, function(current.cue){
                                      fftcuelist[[insert.index]]   <- current.cue
                                      slot(tree,"fftcues",check=F) <- fftcuelist
                                      
                                      return( rating(updateFftree(tree), forw) )
                                    })
    
    ratinglist <- do.call(rbind, ratinglist)
    
    id.best  <- internaltest(ratinglist[,whatToMaximize[chaseItem]])
    #.logger("ID Best:",id.best,"\n")
    ids.lazy <- which(ratinglist[,"lazy"] == TRUE)
    
    #Create new "best" tree
    fftcuelist[insert.index] <- cue.list[id.best]
    slot(tree,"fftcues",check=F) <- fftcuelist
    tree <- updateFftree(tree)
    
    #Remove lazy cues and used cue
    #.logger("Old number of Cues: ", length(cue.list),".\n")
    cue.list <- cue.list[  -c(id.best,ids.lazy)  ]
    
    #.logger("Remaining: ", length(cue.list), ". Current:",Sdt(tree)['percCorr'],"\n")
  }
   
  #stopCluster(cluster)
  
  return(tree)
    
#     res<-foreach(i=1:length(cue.list), .combine=c, .inorder=F,.packages="FFTEST") %dopar% {
#       fftcuelist[insert.index]     <- cue.list[[i]]
#       slot(tree,"fftcues",check=F) <- fftcuelist
#       
#       current.tree <- updateFftree(tree)
#     }
    
#     treesdtlist <- matrix(nrow=length(cue.list), ncol = treesdtlist.cols )
#     colnames(treesdtlist) <- treesdtlist.names
#     
#     fftcuelist <- tree@fftcues
#     fftcuelist[length(fftcuelist) + 1] <- NA
#     
#     #build every possible tree
#     for(i in 1:length(cue.list)){
#       
#       if((i %% 100)==0){
#         cat(Sys.time(), i, "of ", length(cue.list), "\n")
#       }
#       
#       if(forward[chaseItem]){
#         #forward: appending
#         fftcuelist[chaseIndex] <- cue.list[i]
#         slot(tree,"fftcues",check=F) <- fftcuelist
#         #tree@fftcues[chaseIndex] <- cue.list[i]
#       }else{
#         #backward: prepending
#         
#         #If it's the first iteration, prepend item. Else overwrite first item
#         if(i == 1){
#           #tree@fftcues <- as.list(c(cue.list[i], tree@fftcues))
#           slot(tree,"fftcues",check=F) <-  as.list(c(cue.list[i], tree@fftcues))
#         }else{
#           slot(tree,"fftcues",check=F) <-  cue.list[i]          
#         }
#       }
#       tree <- updateFftree(tree)
#       sdtRating <- Sdt(tree)
#       effRating <- getCueEfficiency(tree,hideRest=TRUE)[,-1]
#       
#       if(forward[chaseItem]){
#         effRating <- effRating[nrow(effRating),] # forward, we need the last row, of the appended
#       }else{
#         effRating <- effRating[1,] #prepending, we need the first row
#       }
#       
#       #Save current values of tree in new row
#       #treesdtlist <- rbind(treesdtlist, c(sdtRating, effRating))
#       treesdtlist[i,] <- unlist(  c(sdtRating, unlist(effRating) ))
#     }
#     
#     #find out which tree was "best or"max" or "min"
#     id   <- internaltest(treesdtlist[,whatToMaximize[chaseItem]])
#     lazy <- which(treesdtlist[,"lazy"] == TRUE)   
# 
#     #Append/Prepend best cue
#     if(forward[chaseItem]){
#       #forward: appending
#       tree@fftcues[chaseIndex] <- cue.list[id]
#     }else{
#       #backward: prepending to item "1" because the previous for-loop already created the
#       #first item in the fírst iteration.
#       tree@fftcues[1] <- cue.list[id]
#     }
#     tree <- updateFftree(tree) 
#     
#     #Delete used cue from cuelist
#     cue.list <- cue.list[-c(id,lazy)]
#     
#     cat("Remaining: ", length(cue.list), ". Current:",Sdt(tree)['percCorr'],"\n")
#  }
  
  return(tree)
}

#' Maximizes a \code{\link{Sdt}} value by \code{monte carlo} method
#'
#' Takes all \code{\link{cue}s} in \code{\link{fftree}} and builds \code{blocksize} random generated fftrees.
#' If no better fftree is found, or the the value didn't improve greater than \code{threshold}, the algorithm finishes.
#' Use \code{\link{Fftm}} instead of calling this method directly.
#'
#' @name montecarlo
#'
#' @param tree \code{\link{fftree}}
#' @param whatToMaximize stringvector; Which parameters of the tree should be chased? The index of the vector represents the sequence.  See \link{Sdt} for possible values.
#' @param blocksize integer; How many random trees should be generated, before \code{threshold} is checked?
#' @param threshold double; How much of an improvement in \code{whatToMaximize} should be found to generate another block?
#' @param cluster If you have a HPC cluster, specifiy configuration here. See \code{\link{makeCluster}} from \code{snow} package.
#' @param cores Integer; If multithreading is wanted (recommended for bigger computations!), you can specify your cores here.
#' @return \code{\link{Fftree}}
#' @details You can chase all values provided by \code{\link{Sdt}}
#' @seealso \code{\link{Fftm}}, \code{\link{bruteforcemaximize}}, \code{\link{chase}}, \code{\link{maximize}} 
#' @examples 
#' rawdata <- fftm.Titanic.data()
#' tree <- Fftm(Survived ~ Age + Sex + Class, "montecarlo",rawdata, whatToMaximize = "percCorr", blocksize = 1000, threshold=.001, cores=4)
fftm.montecarlo <- function(tree, whatToMaximize, blocksize, threshold, cluster = NULL, cores = 1){
  
  cnames <- names(Sdt(0,0,0,0))
  if(!(whatToMaximize %in% cnames)){
    stop("Can't maximize ",whatToMaximize, ". Use values provided by function Sdt().")
  }
  
  mc.worker <- function(tree, clist, whatToMaximize, blocksize){
    n <- length(clist)
    randomvector <- 1:n
    res <- matrix(nrow=blocksize, ncol =n+1 )
    
    i<- NULL
    for(i in 1:blocksize){
      vec <- sample(randomvector,n,replace=F)
      ftr <- updateFftree2(tree, clist[ vec ])
      res[i, ] <- c(vec, Sdt.fftree(ftr)[whatToMaximize])
    }          
    
    id <- which.max(  res[, ncol(res)]  )
    return(as.vector(res[id, ]))
  }
  
  #Prepare cluster
  #multithreading wanted
  if(!is.null(cluster) && cores > 1){
    stop("You can't define parameter 'cluster' AND 'cores'. If you want local multithreading use 'cores' (> 0), else specify cluster (see 'makeCluster()').")
  }
  
  if(is.null(cluster) && !is.null(cores) ){
    cluster <- makeCluster(cores, type = "SOCK")
  }
  
  clist <- tree@fftcues
  .logger("Possible variants: ", gamma(1+length(clist))," (",length(clist),"!)")
  diff <- Inf
  best.tree <- NULL
  best.val  <- 0
  breaker <- FALSE
  
  while(breaker == FALSE){
    time <- proc.time()
    res <- clusterCall(cluster, mc.worker, tree=tree, clist=tree@fftcues, whatToMaximize=whatToMaximize, blocksize=blocksize)
    res <- do.call(rbind, res)
    
    best.intermediate.id  <- which.max( res[, ncol(res)] )
    best.intermediate.vec <- res[best.intermediate.id, 1:ncol(res)-1]
    best.intermediate.val <- res[best.intermediate.id, ncol(res)]
    
    time  <- as.integer((proc.time() - time)['elapsed'])
    speed <- (cores*blocksize)/time
    .logger("Current: ", best.val,"; Found: ", best.intermediate.val,"; Difference: ", (best.intermediate.val - best.val)," @ ",speed," calcs/s")
    
    breaker <- (best.intermediate.val - best.val) <= threshold
    
   if(best.intermediate.val > best.val){
     best.val  <- best.intermediate.val
     best.tree <- tree
     best.tree@fftcues <- clist[best.intermediate.vec]
     best.tree <- updateFftree(best.tree)
   }
  }
  
  stopCluster(cluster)
  
  return(best.tree)
}