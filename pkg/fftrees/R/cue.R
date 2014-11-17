######################################
# Base- Class for cuetest.
#
# Contains cue-variables and cached
# results
# 27.04.2014; Marc Giesmann
######################################

.fftest.global.env        <- new.env()
.fftest.global.env$rawcue <- NULL

#' Class cue
#'
#' Class \code{cue} is a node in a tree
#' with parameters to binarize the data
#'
#' @name cue-class
#' @rdname cue-class
#' @exportClass cue
setClass(Class = "cue", 
         slots = c( data           = "numeric",   #inputvector
                    test           = "character", #test-case  parameter
                    split          = "numeric",   #splitvalue parameter
                    pred           = "logical",   #predict- parameter
                    result         = "logical",   #internal result vector
                    predicts       = "logical",   #vector, where cue predicted stuff
                    name           = "character", #optional: name
                    uniquepredicts = "logical"    #named vector with uniques of given data and prediction
                    #checksum = "integer"    #crc32 as integer checksum.
                  ),
         contains = c("logical")
         )

######################
#  Extended generic methods from base
######################

######################

#' Converts \code{cue} to a \code{data.frame}
#'
#' @aliases as.data.frame
#' @param x cue
#' @param row.names not used.
#' @param optional not used
#' @param ... not used
#' 
#' @return data.frame
#' @seealso \code{\link{cue}}
setMethod("as.data.frame", "cue", function(x, row.names = NULL, optional = NULL, ...) {
  return(cbind(data.frame(Data = x@data), Prediction = x@result))
})

# #' Converts \code{cue} in named vector.
# #'
# #' @param object cue
# #' 
# #' @return numeric vector
# #' @seealso \code{\link{Cue}} \code{\link{cue-class}}
# #' ...
# getDataPart.cue <- function(object) {
#   return(
#     c(
#       test  = object@test,
#       split = object@split,
#       pred  = object@pred
#     )
#   )
# }

#' Converts \code{cue} in readable string format.
#'
#' @param x cue
#' @param suppressNameBrackets logical
#' @param suppressThen logical
#' @param suppressName logical
#' @param suppressTest logical
#' @return string
#' @seealso \code{\link{Cue}} \code{\link{cue-class}}
setMethod("toString", "cue", function(x,suppressNameBrackets= FALSE, suppressThen = FALSE, suppressName = FALSE, suppressTest = FALSE) {
  s <- ""
  
  if(!suppressName){
    s <- x@name
    if(!suppressNameBrackets){
      s <- paste("[",s, "]",sep="")
    }
  }
  
  if(!suppressTest){
    s <- paste(s, x@test, x@split)
  }
  
  if(!suppressThen){
    s <- paste(s, "THEN",x@pred)
  }
  return( s )
})

#' Shows \code{cue}
#'
#' @name show
#' @aliases toString
#' 
#' @seealso \code{\link{cue}}
#' ...
setMethod("show", "cue", function(object) {
  cat(toString(object),"\n")
  #show(
  #  as.data.frame(object)
  #  )
})

######################
#  Extended generic methods from stats
######################
#' Predicts value with parameter of cue
#' 
#' If parameters of cue have been changed manually, this has to be called. Internal.
#'
#' @name predict
#' @param object cue
#' @param newdata numerical vector with data which should be predictes
#' 
#' @return logical. TRUE/FALSE
#' @seealso \code{\link{cue}} \code{\link{cue-class}}
setMethod("predict", "cue", function(object, newdata = NULL) {
  if(is.null(newdata)){
   newdata <- cue@data 
  }
  
  if (object@test == "<") {
    return(object@pred == (newdata < object@split))
    
  }else if (object@test == ">"){
    (object@pred == (newdata > object@split))
    
  }else if (object@test == "=="){
    return(object@pred == (newdata == object@split))
    
  }else if (object@test == ">="){
    return(object@pred == (newdata >= object@split))
    
  }else if (object@test == "<="){
    return(object@pred == (newdata <= object@split))
    
  }else {
    #Error
    stop("Predict: Can't translate operator ",object@test,"!")
  }
  
  
})


#' Refreshes/updated cue- prediction- cache.
#' 
#' If parameters of cue have been changed manually, this has to be called. Internal.
#'
#' @param object cue
#' @param data numerical vector with new prediciton data
#' 
#' @return list with S4 cue objects
#' @seealso \code{\link{Cue}} \code{\link{cue-class}}
setMethod("update", "cue", function(object, data) {
  if (object@test == "<") {
    slot(object,"result",check=F)  <- (object@pred == (data < object@split))
    
  }else if (object@test == ">"){
    slot(object,"result",check=F)  <- (object@pred == (data > object@split))
    
  }else if (object@test == "=="){
    slot(object,"result",check=F)  <- (object@pred == (data == object@split))
    
  }else if (object@test == ">="){
    slot(object,"result",check=F)  <- (object@pred == (data >= object@split))
    
  }else if (object@test == "<="){
    slot(object,"result",check=F)  <- (object@pred == (data <= object@split))
    
  }else {
    #Error
    stop("Predict: Can't translate operator ",object@test,"!")
  }
  
  #save data
  slot(object,"data",check=F) <- data
  
  #cue predicted there effectivly, where pred == result
  slot(object,"predicts",check=F) <- (object@pred == object@result)
  
  pred.unique                           <- (object@pred == predict(object, unique(data)))
  names(pred.unique)                    <- unique(data)
  slot(object,"uniquepredicts",check=F) <- pred.unique
  
  return(object)
})

################# INTERNAL METHODS

#' Creates list of all possible variants of cue- configurations.
#' 
#' By defining the optional parameters, these variants will be fix.
#' Parameters with NULL- value will be varied
#'
#' @param inputvector Data, which should be predicted. 
#' @param allSplitValues Numeric value; at which the prediction switches
#' @param allTests character-/ stringvector, which describes the binarization- operator. ('==', '<=', '<', etc.). \code{NULL} means "Autodetect". Autodetect will use c("<=",">") for linear data (more than 2 \code{unique} values in \code{inputvector}) and c("==") for dichotome data (2  \code{unique} values in \code{inputvector})
#' @param allPredictions Logical value; what the cue will predict, if the \code{test} at the particular \code{split} succeeds
#' @param name Charactervector; how the cues should be named
#' 
#' @return list with S4 cue objects
#' @seealso \code{\link{Cue}} \code{\link{Fftree}}
cue.getAllVariants <- function(inputvector, allSplitValues = NULL, allTests = NULL, allPredictions = NULL, name = "Cue"){
    
    if(is.null(allSplitValues)){
      #Autodetect Splitvalues. Uniques from input.
      allSplitValues <- sort(unique(inputvector) , method = c("quick"))
    }
    
    # Autodetect?
    if(is.null(allTests)){ 
      
      #if allSplitvalues are dichotom/binary, force test to be "=="
      if(length(allSplitValues) == 2){
        allTests <- "=="
      }else{
        allTests <- c("<=",">") #all possible tests
      }
    }
    
    #Autodetect?
    if(is.null(allPredictions)){
      allPredictions <- c(TRUE,FALSE) #"Autodetect"
    }
    
    #how many variants will there be?
    variants <- length(allPredictions) * length(allTests) * length(allSplitValues)
    
    #create a list, that can contain all possible cues
    cues <- vector("list", variants)
    
    #Now cache all possible variants: build cue and append it to list
    i = 1
    
    for(p in allPredictions){
      for(t in allTests){
        for(s in allSplitValues){
          cues[[i]] <- Cue(inputvector,t,s,p,name) # call constructor of cue
          
          #TODO: Entferne unnütze Datensätze.
          #le <- length(which( cues[[i]]@result))
          #  
          #  if(le == length(cues[[i]]@result) || le == 0){
          #    j = j +1
          #  }
          
          i = i+1
        }
      }
    }
  
  
  return(cues)
}

.cuelist.normalize <- function(clist){
  if(!is.list(clist)){
    stop("Parameter clist has to be of type 'list'")
  }
  
  
  clist <- lapply(clist, function(cue){
    if(cue@test == "<"){
      if(cue@split == min(cue@data)){
        return(NULL)
      }else{
        uniquedata <- sort(unique(cue@data), decreasing = F)
        cue@split  <- uniquedata[which(uniquedata == cue@split) - 1]
        cue@test   <- "<="
        return(cue)
      }
    }
    
    
    if(cue@test == ">"){
      if(cue@split == max(cue@data)){
        return(NULL)
      }else{
        uniquedata <- sort(unique(cue@data), decreasing = F)
        cue@split  <- uniquedata[which(uniquedata == cue@split) + 1]
        cue@test   <- ">="
        return(cue)
      }
    }
    
    
    return(cue)
  })
  
  return(clist[! sapply(clist,is.null)])
}

#' Creates a cue.
#'
#' A \code{cue} is an S4 object, that represents a node in a tree.
#' Fast- and Frugal trees use binary predicitons, but data
#' normally is numeric. By defining the cue- parameters
#' the cue automatically binarizes given data by a set of parameters.
#'
#' @name cue
#' @rdname cue-class
#'
#' @param inputvector Data, which should be predicted. If you just want to build a tree without measuring it's performance, set to c(0)
#' @param test Charactervector, which describes the binarization- operator. Possible values \code{'==', '<=', '<', '>', '>='}
#' @param split Numeric value; at which the prediction switches
#' @param pred Logical value; what the cue will predict, if the \code{test} at the particular \code{split} succeeds
#' @param name Charactervector; how the cue is named (just for description, without function)
#' 
#' @return cue S4 object
#' @seealso \code{\link{Fftree}}, \code{\link{cue.getAllVariants}}
Cue<-function(inputvector, test = "<=", split = 0, pred = TRUE, name = ""){
  
  # return instance of cue
  if(is.null(.fftest.global.env$rawcue)){
    .fftest.global.env$rawcue <- new("cue",test = "<=", split = 0, pred = TRUE, name = "RAWCUE")
  }
  
  C <- .fftest.global.env$rawcue
  slot(C,"test",check=F)  <- test
  slot(C,"split",check=F) <- split
  slot(C,"pred",check=F)  <- pred
  slot(C,"name",check=F)  <- name
   
  #C <- new("cue",test = test, split = split, pred = pred, name = name)
  
  #update (calculate prediciton data) and return
  return ( update(C, inputvector)  )
}