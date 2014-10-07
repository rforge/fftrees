######################################

# Signal-Detection Vector
######################
#  Generic methods
######################

# Generic Constructor for sdt
#' Wrapper function sdt
#'
#' @name Sdt
#' @param ... further parameter
#' @export
Sdt <- function(hi,...) UseMethod("Sdt")

#' Creates a 'Signal Detection Theory' vector
#'
#' @name Sdt
#'
#' @param hi numeric; hits / true positives
#' @param fa numeric; false alarms / false positives
#' @param mi numeric; misses / false negatives
#' @param cr numeric; correct rejection / true negatives
#' 
#' @return numeric vector with signal-detection values
#' @seealso \code{\link{Sdt.fftree}}
#' @references \url{http://kangleelab.com/signal detection theory.html} \url{http://en.wikipedia.org/wiki/Matthews_correlation_coefficient}
#' @details This function returns: hitrate (sensitivity/TPR), specifity (true negative rate/SPC), false alarm rate (fall-out/FPR), false discovery rate (FDR), 
#' an estimated d' (qnorm(hitrate)-qnorm(false alarm rate)) and the MCC, the
#' "Matthews correlation efficient".
#'  
#' Some results are adjusted, to make them calculatable. If one of the contingency-values \code{hi}, \code{fa}, \code{mi} or \code{cr} equals zero, 
#' all of them will gain .25: \code{Sdt(1, 0, 2, 4)} equals \code{Sdt(1.25, .25, 2.25, 4.25)}. 
#' The denominator of the Matthews correlation coefficient is adjusted to 1 if \code{(hi + fa) == 0}, \code{(hi + mi) == 0}, \code{(fa + cr) == 0} or \code{(cr + mi) == 0}.
Sdt.default <-  function(hi,fa,mi,cr){ 

  #Basics
  sum_all <- hi+fa+mi+cr
  per_correct <- (hi+cr)/sum_all
  
  #Corrections
  if(hi == 0 || fa == 0 || mi == 0 || cr == 0){
    hi.corrected <- hi + .25
    fa.corrected <- fa + .25
    mi.corrected <- mi + .25
    cr.corrected <- cr + .25
  }else{
    hi.corrected <- hi
    fa.corrected <- fa
    mi.corrected <- mi
    cr.corrected <- cr
  }

  
  #standard calculation
  sum_fa_cr <- fa.corrected+cr.corrected
  sum_hi_mi <- hi.corrected+mi.corrected
  
  hitrate            <- hi.corrected/sum_hi_mi
  falsealarmrate     <- fa.corrected/sum_fa_cr 
  falsediscoveryrate <- fa.corrected/(fa.corrected+hi.corrected)
  specifity          <- 1-falsealarmrate
  
  qh <- suppressWarnings(qnorm(hitrate))
  qf <- suppressWarnings(qnorm(falsealarmrate))
  
  #d' aka dPrime
  dPrime <- qh-qf
  
  
  #Matthews correlation coefficient
  
  #"In this equation, TP is the number of true positives, 
  #TN the number of true negatives, FP the number of false positives and 
  #FN the number of false negatives. 
  #If any of the four sums in the denominator is zero, 
  #the denominator can be arbitrarily set to one; 
  #this results in a Matthews correlation coefficient of zero, 
  #which can be shown to be the correct limiting value."
  # Source: http://en.wikipedia.org/wiki/Matthews_correlation_coefficient 
  
  hifa <- hi + fa
  himi <- hi + mi
  facr <- fa + cr
  crmi <- cr + mi
  
  if(hifa == 0 || himi==0 || facr == 0 || crmi==0){
    mcc.denominator <- 1
  }else{
    mcc.denominator <- sqrt( prod(hifa,himi,facr,crmi))
  }
  
  MCC <- ((hi * cr) - (fa * mi)) /  mcc.denominator  
  
  #BETA / natural log
  #beta_natLog   <-  -dPrime*0.5*(qh + qf)
  
  #BETA / ratio
  #beta_ratio <- exp(beta_natLog)
  
  #criterion c
  #C <- -.5*(qh + qf)
  
  #normalized c (c')
  #C_norm <- C/dPrime
  
  
  return(  c(hi = hi, fa = fa, mi = mi, cr = cr, 
                    hiRate = hitrate, spec = specifity, faRate = falsealarmrate, fdRate = falsediscoveryrate,
                    dPr = dPrime, MCC = MCC,
                    #betaNl = beta_natLog, betaRa = beta_ratio,
                    #crit = C, critNrm = C_norm ,
                    percCorr = per_correct)) 
  
}

#' Creates a 'Signal Detection Theory' vector
#'
#' @name Sdt
#'
#' @param criterion logical vector
#' @param prediction logical vector
#' 
#' @return numeric vector with signal-detection values
#' @seealso \code{\link{Sdt.fftree}}
Sdt.logical <- function(criterion, prediction){
  hi <- length(which(criterion[ prediction ] == TRUE ))
  fa <- length(which(criterion[ prediction ] == FALSE))
  
  mi <- length(which(prediction[ criterion ]  == FALSE ))
  cr <- length(which(prediction[ !criterion ] == FALSE  ))
  
  
  #call default constructor
  return(Sdt(hi,fa,mi,cr))
}

######################
#  Extended generic methods from base
######################

#GET- method
# PRINT- Method, std-instance
# setMethod("as.vector", "sdt", function(x) {
#   object <- x  
#   return( c(
#       hits  = (object@hi),
#       fa    = (object@fa),
#       mi    = (object@mi),
#       cr    = (object@cr),
#       
#       #now all variants,
#       hitRate  = object@hiRate,  #hitrate
#       faRate   = object@faRate, #false alarm rate
#       dPr      = object@dPr,     #dPrime
#       betaNl   = object@betaNl,  #beta nat log
#       betaRa   = object@betaRa,  #beta ratio
#       crit     = object@crit,    #crit
#       critNrm  = object@critNrm, #crit Norm
#       percCorr = object@percCorr #percent correct
#       )
#     
#   )
#   
# } )

# # PRINT- Method, std-instance
# setMethod("show", "sdt", function(object) {
#   
#   print(as.vector(object), quote = FALSE)
#   
# } )