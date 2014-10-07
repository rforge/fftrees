#include <Rcpp.h>
using namespace Rcpp;


//// [[Rcpp::export]]
//S4 setCuelist(S4 x, Rcpp::List cl){
//     //S4 fftree(x) ;
//     
//     x.slot("fftcues") = cl;
//     return x;
//}

//' @name update
//' @param fftree fftree
//' @return fftree
//' @details This a somewhat internal function. Use it with caution. Normally one should use the wrapped version \code{update} instead of the direct call \code{updateFftree}.
//' @seealso \code{\link{update},\link{fftree}}
//' ...
// [[Rcpp::export]]
S4 updateFftree(S4 fftree){
     //S4 fftree(x) ;
     S4 crit(fftree.slot( "criterion" ));      
     
     //Cuelist size
     Rcpp::List cuelist(fftree.slot( "fftcues" )); 
     int cLsize = cuelist.size();

     //Criterion size
     Rcpp::LogicalVector criterion(crit.slot( "result"));
     int cRsize = criterion.size();
     
     //Tree result
     Rcpp::LogicalVector treePrediction(cRsize);

     //SDT Counter, for every possible cue
     Rcpp::NumericVector hi(cLsize);
     Rcpp::NumericVector fa(cLsize);
     Rcpp::NumericVector mi(cLsize);
     Rcpp::NumericVector cr(cLsize);
     
     //Logical: Has criterion been predicted
     Rcpp::NumericVector predicted(cRsize);
     
     //cuePotential, fftreePotential,rest, restCumulated
     Rcpp::NumericVector cPot(cLsize);
     Rcpp::NumericVector fPot(cLsize);
     Rcpp::NumericVector rest(cLsize);
     int restCumulated = cRsize;
     
     //Cue prediction
     int pred = 0;
     
     //TREE-stats
     int endHi=0;
     int endFa=0;
     int endMi=0;
     int endCr=0;
     
     int restHi = 0;
     int restFa = 0;
     int restMi = 0;
     int restCr = 0;
     
     //iterate through all cues
     for(int i = 0; i < cLsize; i++){
       //Read from list
        S4 currentCue(cuelist[i]);
        pred = currentCue.slot("pred");
        Rcpp::IntegerVector cuePredicts(currentCue.slot("predicts"));
        
        //iterate through all criterion-booleans
        for(int j = 0; j < cRsize; j++){
          
          if(predicted[j] == 0 && cuePredicts[j]){ //If criterion j hasn't been predicted yet
            
            //Mark as predicted
            predicted[j] = i+1;
            
            //save Treeprediction
            treePrediction[j] = pred;
            
            //Count
            if(criterion[j] == 1 && pred == 1)  {hi[i]++;} 
            if(criterion[j] == 0 && pred == 1)  {fa[i]++;}
            if(criterion[j] == 1 && pred == 0)  {mi[i]++;}
            if(criterion[j] == 0 && pred == 0)  {cr[i]++;}
          }
          
        }//end FOR cRsize
        
        //Update cue-stats
        cPot[i] = hi[i] + fa[i] + mi[i] + cr[i];
        fPot[i] = restCumulated;
        restCumulated-=cPot[i];
        
        //Update tree stats
        endHi+=hi[i];
        endFa+=fa[i];
        endMi+=mi[i];
        endCr+=cr[i];
    }//end FOR cLsize
 
    //Now another for loop for all unpredicted criterions. The last cues prediction is still stored in "pred". 
    //iterate through all criterion-booleans
    pred = !pred;
    for(int j = 0; j < cRsize; j++){
      
      if(predicted[j] == 0){ //If criterion j hasn't been predicted yet
          //Count
          if(criterion[j] == 1 && pred == 1)  {restHi++;} 
          if(criterion[j] == 0 && pred == 1)  {restFa++;}
          if(criterion[j] == 1 && pred == 0)  {restMi++;}
          if(criterion[j] == 0 && pred == 0)  {restCr++;}
      }
    }
 
    //Update tree stats
    endHi+=restHi;
    endFa+=restFa;
    endMi+=restMi;
    endCr+=restCr;
 
    
    fftree.slot("result")        = treePrediction;
    fftree.slot("classcounter")  = cPot;
    fftree.slot("rest")          = fPot;
    fftree.slot("potential")     = fPot;
    fftree.slot("restCumulated") = restCumulated;
    fftree.slot("predictedBy")    = predicted;
    
    fftree.slot("resHi")    = hi;
    fftree.slot("resFa")    = fa;
    fftree.slot("resMi")    = mi;
    fftree.slot("resCr")    = cr;
    
    fftree.slot("endHi")    = endHi;
    fftree.slot("endFa")    = endFa;
    fftree.slot("endMi")    = endMi;
    fftree.slot("endCr")    = endCr;
    
    fftree.slot("restHi")    = restHi;
    fftree.slot("restFa")    = restFa;
    fftree.slot("restMi")    = restMi;
    fftree.slot("restCr")    = restCr;
    return fftree;
 
}

//' Updates tree cache
//'
//' Updates trees internal caches and results without writing new cues into trees internal storage.
//'
//' @name updateFftree2
//' @param fftree fftree
//' @param cuelist list. List with cues.
//' @return fftree
//' @details Use it with caution. Normally one should use the wrapped version \code{\link{update}}. This specific variant avoids writing the used cues back into the trees internal storage to avoid the list-bottlenecks.
//' @seealso \code{\link{update},\link{fftree}}
//' ...
// [[Rcpp::export]]
S4 updateFftree2(S4 fftree, Rcpp::List cuelist){
     //S4 fftree(x) ;
     S4 crit(fftree.slot( "criterion" ));      
     
     //Cuelist size
     //Rcpp::List cuelist(fftree.slot( "fftcues" )); //Rcpp::List cuelist(fftree.slot( "fftcues" )); 
     int cLsize = cuelist.size();

     //Criterion size
     Rcpp::LogicalVector criterion(crit.slot( "result"));
     int cRsize = criterion.size();
     
     //Tree result
     Rcpp::LogicalVector treePrediction(cRsize);

     //SDT Counter, for every possible cue
     Rcpp::NumericVector hi(cLsize);
     Rcpp::NumericVector fa(cLsize);
     Rcpp::NumericVector mi(cLsize);
     Rcpp::NumericVector cr(cLsize);
     
     //Logical: Has criterion been predicted
     Rcpp::NumericVector predicted(cRsize);
     
     //cuePotential, fftreePotential,rest, restCumulated
     Rcpp::NumericVector cPot(cLsize);
     Rcpp::NumericVector fPot(cLsize);
     Rcpp::NumericVector rest(cLsize);
     int restCumulated = cRsize;
     
     //Cue prediction
     int pred = 0;
     
     //TREE-stats
     int endHi=0;
     int endFa=0;
     int endMi=0;
     int endCr=0;
     
     int restHi = 0;
     int restFa = 0;
     int restMi = 0;
     int restCr = 0;
     
     //iterate through all cues
     for(int i = 0; i < cLsize; i++){
       //Read from list
        S4 currentCue(cuelist[i]);
        pred = currentCue.slot("pred");
        Rcpp::IntegerVector cuePredicts(currentCue.slot("predicts"));
        
        //iterate through all criterion-booleans
        for(int j = 0; j < cRsize; j++){
          
          if(predicted[j] == 0 && cuePredicts[j]){ //If criterion j hasn't been predicted yet
            
            //Mark as predicted
            predicted[j] = i+1;
            
            //save Treeprediction
            treePrediction[j] = pred;
            
            //Count
            if(criterion[j] == 1 && pred == 1)  {hi[i]++;} 
            if(criterion[j] == 0 && pred == 1)  {fa[i]++;}
            if(criterion[j] == 1 && pred == 0)  {mi[i]++;}
            if(criterion[j] == 0 && pred == 0)  {cr[i]++;}
          }
          
        }//end FOR cRsize
        
        //Update cue-stats
        cPot[i] = hi[i] + fa[i] + mi[i] + cr[i];
        fPot[i] = restCumulated;
        restCumulated-=cPot[i];
        
        //Update tree stats
        endHi+=hi[i];
        endFa+=fa[i];
        endMi+=mi[i];
        endCr+=cr[i];
    }//end FOR cLsize
 
    //Now another for loop for all unpredicted criterions. The last cues prediction is still stored in "pred". 
    //iterate through all criterion-booleans
    pred = !pred;
    for(int j = 0; j < cRsize; j++){
      
      if(predicted[j] == 0){ //If criterion j hasn't been predicted yet
          //Count
          if(criterion[j] == 1 && pred == 1)  {restHi++;} 
          if(criterion[j] == 0 && pred == 1)  {restFa++;}
          if(criterion[j] == 1 && pred == 0)  {restMi++;}
          if(criterion[j] == 0 && pred == 0)  {restCr++;}
      }
    }
 
    //Update tree stats
    endHi+=restHi;
    endFa+=restFa;
    endMi+=restMi;
    endCr+=restCr;
 
    //fftree.slot("fftcues")       = cuelist;
    fftree.slot("result")        = treePrediction;
    fftree.slot("classcounter")  = cPot;
    fftree.slot("rest")          = fPot;
    fftree.slot("potential")     = fPot;
    fftree.slot("restCumulated") = restCumulated;
    fftree.slot("predictedBy")    = predicted;
    
    fftree.slot("resHi")    = hi;
    fftree.slot("resFa")    = fa;
    fftree.slot("resMi")    = mi;
    fftree.slot("resCr")    = cr;
    
    fftree.slot("endHi")    = endHi;
    fftree.slot("endFa")    = endFa;
    fftree.slot("endMi")    = endMi;
    fftree.slot("endCr")    = endCr;
    
    fftree.slot("restHi")    = restHi;
    fftree.slot("restFa")    = restFa;
    fftree.slot("restMi")    = restMi;
    fftree.slot("restCr")    = restCr;
    return fftree;
 
}
