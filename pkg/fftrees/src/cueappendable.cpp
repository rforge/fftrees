#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export(".isCueAppendable")]]
  bool isCueAppendable(S4 fftree, Rcpp::List cuelist){     
     //Cuelist size
     int cLsize = cuelist.size();
     
     if(cLsize <=1){
       //Rcpp::Rcout << "(Obviously) nice! " << std::endl;
       return TRUE;
     }
     
     S4 cue(cuelist[cLsize-1]);
     
     //Data of new cue
     std::string   newCueName       = cue.slot("name");
     //LogicalVector newCueUniquePred = cue.slot("uniquepredicts");
     float         newCueSplit      = cue.slot("split");
     std::string   newCueTest       = cue.slot("test");
     bool          newCuePred       = cue.slot("pred");
     
     //iterators
     int i = 0;
     
     if(newCueName != ""){
       //Check for overlapping
       //iterate through all cues
       
       for(i = 0; i < cLsize-1; i++){
         ///Read from list
         S4 currentCue(cuelist[i]);
         std::string   currentName        = currentCue.slot("name");
         //LogicalVector currentUniquePred  = currentCue.slot("uniquepredicts");
         float   currentSplit             = currentCue.slot("split");
         std::string   currentTest        = currentCue.slot("test");
         
         //If not same nameskip this iteration
         if(currentName != newCueName){ continue; }
         
         //If both predict one of the same unique values in data, then their spectrum overlaps
         if((currentTest == "<=") && (newCueSplit <= currentSplit)){
             //Rcpp::Rcout << "Trivial: " << currentTest << currentSplit << " + " << newCueTest <<  newCueSplit    <<std::endl;
             return FALSE;
         }
         
         if((currentTest == ">=") && (newCueSplit >= currentSplit)){
             //Rcpp::Rcout << "Trivial: " << currentTest << currentSplit << " + " << newCueTest <<  newCueSplit    <<std::endl;
             return FALSE;
         }
         
         if((currentTest == "==") && (newCueSplit == currentSplit)){
             //Rcpp::Rcout << "Trivial: " << currentTest << currentSplit << " + " << newCueTest <<  newCueSplit    <<std::endl;
             return FALSE;
         }
         
         
       }//for
     }//if
     
     //Check for triviality
     S4 lastCue(cuelist[cLsize - 2]);
     std::string   lastCueName       = lastCue.slot("name");
     float         lastCueSplit      = lastCue.slot("split");
     std::string   lastCueTest       = lastCue.slot("test");
     bool          lastCuePred       = lastCue.slot("pred");
     
     if(newCueName == lastCueName && newCueTest == lastCueTest && newCuePred == lastCuePred && newCueTest != "=="){
       //Rcpp::Rcout << "Trivial: " << lastCueTest << lastCueSplit << " + " << newCueTest <<  newCueSplit    <<std::endl;
       return FALSE;
     }
     
//     //Check for lazyness
//     LogicalVector newCuePredicts    = cue.slot("predicts");
//     NumericVector fftreePredictedBy = fftree.slot( "predictedBy" );
//     for(i=0; i<fftreePredictedBy.size(); i++){
//       if((fftreePredictedBy[i] == 0) && newCuePredicts[i]){
//         return TRUE;
//       } 
//     }
     
     //Check if it is trivial or lazy
     //Rcpp::Rcout << "Lazy! " << std::endl;
     //Rcpp::Rcout << "Nice! " << std::endl;
     return TRUE;
}
