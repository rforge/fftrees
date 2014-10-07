#' CONST_SKIP
#' @details internal constant for \code{\link{permu.new}}. Represents the vector-index of the SKIP flag.
CONST_SKIP <- 1

#' CONST_SAVE
#' @details internal constant for \code{\link{permu.new}} Represents the vector-index of the SAVE flag.
CONST_SAVE <- 2

#' CONST_VAL
#' @details internal constant for \code{\link{permu.new}} Represents the vector-index of the user defined VAL-parameter.
CONST_VAL  <- 3

#' permu.old
#'
#' Generates permutations of \code{perm} and calls a callback \code{fun} for every permutation, without preallocating
#' memory. This helps, when permutating very large vectors.
#'
#' @param perm numeric vector; numbers which should be permuted
#' @param fun function; callback function with 1 parameter
#' @param steps logical; Generate only full permutations, or generate steps between them? Default: \code{FALSE}
#' @param current numeric vector; Prepends given vector to permutation before calling \code{fun}
#' 
#' @return nothing. If you want to save calculations, you need to write them into global workspace via \code{fun}
permu.old <- function(perm,fun,steps=FALSE,current=NULL){
  for(i in 1: length(perm)){
    
    fix  <- c(current,perm[i])   # calculated elements; fix at this point
    rest <- perm[-i]  # elements yet to permutate
    
    #Call callback.
    if(steps || !length(rest)){
      result <- fun(fix)
    }
    
    if(length(rest)){
      result <- permu.old(rest, fun, steps, fix)
    }
  }#end for
}

#' permu.new
#'
#' INTERNAL. Generates permutations of \code{perm} and calls a callback \code{fun} for every permutation, without preallocating
#' memory. It uses SNOW-package for parallel computing. This helps, when permutating very large vectors.
#'
#' @param perm numeric vector; numbers which should be permuted
#' @param fun function; callback function with 1 parameter. See examples for more information.
#' @param funsave function; callback function with 2 parameters: \code{x} for current permutations and \code{resultobject} for saving purposes. See examples for more information.
#' @param values numeric; Cache for \code{resultobject}
#' @param savemax numeric; Cache of results calculated by \code{fun}
#' @param cores numeric; Number of preocessors should be used
#' @param cluster cluster-variable created by \code{\link{makeCluster}} if more computers should be used.
#' 
#' @return list of all return values of \code{funsave}.
#' @examples
#'#callback "fun"
#'perm.callback.FUN.STUMP <- function(x,resultobject){
#'  resultobject[CONST_SAVE] <-TRUE
#'  resultobject[CONST_SKIP] <-FALSE
#'  
#'  return(resultobject)
#'}
#'
#'#callback "funsave"
#'perm.callback.FUNSAVE.STUMP <- function(resultlist){
#'  res <- resultlist[['resultobject']]
#'  lis <- resultlist[['results']]
#'  
#'  return(list(results=lis, resultobject= res))
#'}
#'
#'# Example:
#'#----------------------------------------------------------------
#'# EXAMPLE CALLBACK
#'#----------------------------------------------------------------
#'# Prunes, if 3 is last number in permutation
#'# Saves only, if sum() of permutation is the highest found yet.
#'# IMPORTANT: return has to be a "resultobject", which is provided
#'# through the parameters. 
#'# Use 
#'# resultobject[CONST_SKIP] <- TRUE/FALSE (prune after this permutation T/F)
#'# resultobject[CONST_SAVE] <- TRUE/FALSE (return this permutation, save it T/F)
#'# resultobject[CONST_VAL]  <- NUMERIC (use this to save something for the process)
#'#-----------------------------------------------------------------
#'perm.callback <- function(x,resultobject){
#'  
#'  #CALCULATE STUFF HERE;
#'  
#'  #SKIP EXAMPLE
#'  #Skip this one? skip next permutations if the last number is 3
#'  resultobject[CONST_SKIP] <- (x[length(x)] == 3)
#'   
#'  #SAVE EXAMPLE
#'  #Should we save this permutation?
#'  #Save only, if sum of permutation is bigger than the ones we already saved. 
#'  s <- sum(x)
#'  if(s > resultobject[CONST_VAL]){
#'    resultobject[CONST_VAL]  <- s
#'    resultobject[CONST_SAVE] <-TRUE
#'  }else{
#'    resultobject[CONST_SAVE] <-FALSE
#'  }
#'  
#'  return(resultobject)
#'}
#'
#'#----------------------------------------------------------------
#'# EXAMPLE CALLBACK FOR SAVING
#'#-----------------------------------------------------------------
#'# Orders resultlist, saves only "TOP 50"
#'#
#'# INPUT: List with 2 items: 'resultobject' and 'results'
#'# res <- resultlist[['resultobject']]
#'# lis <- resultlist[['results']]
#'#
#'# OUTPUT: (list with item 'results'(resultlist)) and resultobject,
#'#         to pass on to further calculations
#'# return(list(results=lis, resultobject= res))
#'#-----------------------------------------------------------------
#'perm.callback.save <- function(resultlist){
#'  res <- resultlist[['resultobject']]
#'  lis <- resultlist[['results']]
#'  
#'  # ORDER RESULTLIST
#'  if(length(lis) > 0){
#'    
#'    #For Example "TOP 50": Only save the top 50!
#'    allvals <- sapply(lis, function(x) x[['values']], simplify = TRUE ) # get all values
#'    
#'    lis <- lis[  order( allvals )  ]                   # order list by values
#'   
#'    if(length(lis) > 50)
#'      lis <- lis[1:50]                                 # only save Top 50
#'    
#'    res[CONST_VAL] <- max(allvals)
#'  }
#'  
#'  return(list(results=lis, resultobject= res))
#'}
#'
#'#Execution
#'result <- permu.new(perm=1:10, fun=perm.callback,funsave=perm.callback.save,values=1, cores = 4)
permu.new <- function(perm,fun, funsave, values = 0, savemax = 100000, cores = NULL, cluster = NULL){
  
  #Prepare cluster
  #multithreading wanted
  if(!is.null(cluster) && cores > 1){
    stop("You can't define parameter 'cluster' AND 'cores'. If you want local multithreading use 'cores' (> 0), else specify cluster (see 'makeCluster()').")
  }
  
  if(is.null(cluster) && !is.null(cores) ){
    cluster <- makeCluster(cores, type = "SOCK")
    
    #Export Everythin
    #clusterExport(cluster,ls(envir = .GlobalEnv),.GlobalEnv)
    #clusterExport(cluster,ls(), environment())
    
    registerDoSNOW(cluster)
    
  }
  
  if(is.null(cores)){
    cores <- 1
  }
  
  
  #DEFINE INTERNAL FUNCTIONS

  
  #CREATES RESULTOBJECT
  robj <- function(vals){
    ob<- vector(mode="numeric",length=2+vals)
    ob[CONST_VAL] <- -Inf
    return(ob)
  }
  

  #DEFINE INTERNAL END
  
  
  #BEGIN FUNCTION----------!
  resultobject <- robj(values) #resultobject for
  endresult <- list()
  
  all.calculations            <- sum(gamma(1 + length(perm)) /  gamma(length(perm) - perm + 1) )
  all.calculations.percore    <- all.calculations / cores
  
  calc.bruteforce <- gamma(1 + length(perm))
  
  timer <- proc.time() #set timer
  
  for(from in seq(1,length(perm),cores )){
    to <- from + cores- 1
    if(to > length(perm)){
      to <- length(perm)
    }
    
    calculations.per.step <- all.calculations/(cores)
    calculated <- calculations.per.step * (from-1)
    t         <- as.integer((proc.time() - timer)['elapsed'])
    speed     <- calculated/t
    ETA.MIN   <- round( ((all.calculations-calculated) / speed) /60 ,1)
    ETA.H     <- round(ETA.MIN /60, 1)
    if(from > 1){
      .logger(calculated, " / ",all.calculations," (",calc.bruteforce,")", " calculated. @",speed,"/s. Estimated:",ETA.MIN,"min (~",ETA.H,"h )", digits=4)
    }

    .logger("Calculating step ", from, " to ", to, " of ", length(perm)," (",calculations.per.step, " permutations)", digits=4)
    
    res<-foreach(i=from:to, .combine=c, .inorder=F,.packages="fftrees") %dopar% {      
        #WORKERBEE. Does the funpart of recursion and calling the callbacks
        permu.worker <- function(perm, current, resultobject, fun){
          resultobject[1:2] <- 0 #reset internal values.
          
          for(i in 1: length(perm)){
            
            fix  <- c(current,perm[i])   # calculated elements; fix at this point
            rest <- perm[-i]  # elements yet to permutate
            
            #Call callback.
            resultobject <- fun(x=fix, resultobject = resultobject)
            
            #Save permutation?
            if(resultobject[CONST_SAVE]){
              permu.worker.save(fix, resultobject[CONST_VAL])
            }
            
            #if this is the call with the last
            #value (the deepest,recursive call) or object wanted
            #to skip next iterations stop recursion
            if(length(rest) && !resultobject[CONST_SKIP]){
              resultobject <- permu.worker(rest, fix, resultobject, fun)
            }
          }#end for
          
          return(resultobject)
        }
      
        permu.worker.save.max        <- savemax
        permu.worker.save.count      <- 1
        permu.worker.global.savelist <- vector(mode="list",length = permu.worker.save.max)
        
        #Saves permutation. If there are more to save than in savemax defined,
        #it primitlively appends a entry to the list
        permu.worker.save <- function(permutation, values){
          if(permu.worker.save.count > permu.worker.save.max){
            .logger("Worker says: MEMORY LIMIT EXCEEDED! Trying to allocate one more. \n")
            permu.worker.global.savelist[[length(permu.worker.global.savelist)+1]] <<- list(perm=permutation,values=values)
          }else{
            permu.worker.global.savelist[[permu.worker.save.count]] <<- list(perm=permutation,values=values)
          }
          permu.worker.save.count <<- permu.worker.save.count + 1 
        }
        
        #calculate the first permutation manually
        resultobject <- permu.worker(perm[i], NULL, resultobject, fun)
      
        if(!resultobject[CONST_SKIP]){
          #now do the funny, recursive stuff
          resultobject <- permu.worker(perm[-i], perm[i], resultobject, fun)
        }
  
        # Now we're ready for the next permutation.
        # Save all the things we need
        return(permu.worker.global.savelist[1:permu.worker.save.count-1])
        
    }#end foreach
    
    #Only save stuff, if found stuff
    if(length(res)){
      save.function.result <- funsave(list(results = c(endresult,res), resultobject = resultobject))
      resultobject <- save.function.result[['resultobject']]
      endresult    <- save.function.result[['results']]
    }
    #CLEANUP RES
    
  }#end for
    
  stopCluster(cluster)
  
  
return(endresult) 
}


#----------------------------------------------------------------
#EXAMPLE CALLBACK
# Prunes, if 3 is last number in permutation
# Saves only, if sum() of permutation is the highes found yet.
# IMPORTANT: return has to be a "resultobject", which is provided
# through the parameters. 
# Use 
# resultobject[CONST_SKIP] <- TRUE/FALSE (prune after this permutation T/F)
# resultobject[CONST_SAVE] <- TRUE/FALSE (return this permutation, save it T/F)
# resultobject[CONST_VAL]  <- NUMERIC (use this to save something for the process)
#-----------------------------------------------------------------
# perm.callback <- function(x,resultobject){
#   
#   #CALCULATE STUFF HERE;
#   #Example a global counter;(works only singlethreaded)
#   #counter <<- counter + 1
#   
#   #SKIP EXAMPLE
#   #Skip this one? skip next permutations if the last number is 3
#   resultobject[CONST_SKIP] <- (x[length(x)] == 3)
#   
#   #if(resultobject[CONST_SKIP]){
#     #another global counter (works only singlethreaded)
#     #skipped <<- skipped + 1 
#   #}
#   
#   #SAVE EXAMPLE
#   #Should we save this permutation?
#   #Save only, if sum of permutation is bigger than own value 
#   s <- sum(x)
#   if(s > resultobject[CONST_VAL]){
#     resultobject[CONST_VAL]  <- s
#     resultobject[CONST_SAVE] <-TRUE
#     
#     #yet another example-counter. (works only singlethreaded)
#     #saved <<- saved + 1 
#   }else{
#     resultobject[CONST_SAVE] <-FALSE
#   }
#   
#   return(resultobject)
# }

#----------------------------------------------------------------
#EXAMPLE CALLBACK FOR SAVING
# Orders resultlist, saves only "TOP 50"
#
#INPUT: List with 2 items: 'resultobject' and 'results'
#res <- resultlist[['resultobject']]
#lis <- resultlist[['results']]
#
#OUTPUT: (list with item 'results'(resultlist)) and resultobject,
#         to pass on to further calculations
#return(list(results=lis, resultobject= res))
#-----------------------------------------------------------------
# perm.callback.save <- function(resultlist){
#   res <- resultlist[['resultobject']]
#   lis <- resultlist[['results']]
#   
#   # ORDER RESULTLIST
#   if(length(lis) > 0){
#     
#     #For Example "TOP 50": Only save the top 50!
#     allvals <- sapply(lis, function(x) x[['values']], simplify = TRUE ) # get all values
#     
#     lis <- lis[  order( allvals )  ]                   # order list by values
#     
#     if(length(lis) > 50)
#       lis <- lis[1:50]                                 # only save Top 50
#     
#     res[CONST_VAL] <- max(allvals)
#   }
#     
#   return(list(results=lis, resultobject= res))
# }


# Callback STUMP
# perm.callback.STUMP <- function(x,resultobject){
#   resultobject[CONST_SAVE] <-TRUE
#   resultobject[CONST_SKIP] <-FALSE
#     
#   return(resultobject)
# }
# 
# perm.callback.save.STUMP <- function(resultlist){
#   res <- resultlist[['resultobject']]
#   lis <- resultlist[['results']]
# 
#   return(list(results=lis, resultobject= res))
# }


#---------- MAINEXAMPLECODE
# #counter/skipped/saved are working in singlethreading mode,
# #See usage in perm.callback().
# #
#Variables show, how many...
#counter <- 0 # ...permutations have been calculated 
#skipped <- 0 # ... have been skipped (last digit was 3)
#saved   <- 0 # ... were saved and returned

#registerDoMC(4) #uncomment for multithreading
# Rprof(filename = "Rprof.out")
# stime <- system.time(gcFirst = TRUE, expr ={
# result <- permu.new(perm=1:10, fun=perm.callback,funsave=perm.callback.save,values=1, cores = 4)
# })
# Rprof(NULL)
# summaryRprof(filename = "Rprof.out")
# cat(as.double(stime[3]), "seconds; ~", (2959228 / as.double(stime[3])), " calculations/second")
