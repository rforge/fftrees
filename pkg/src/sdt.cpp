#include <Rcpp.h>
using namespace Rcpp;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar

//// [[Rcpp::export]]
//NumericVector internalSdtlogical(LogicalVector x, LogicalVector y) {
//      Rcpp::NumericVector out(4);
//        
//      Rcpp::LogicalVector criterion(x);
//      Rcpp::LogicalVector prediction(y);
//
//      int hi = 0;
//      int fa = 0;
//      int mi = 0;
//      int cr = 0;
//      int n = criterion.size();
//          
//      for(int i=0; i<n; i++){
//        if(criterion[i] == 1 && prediction[i] == 1)  { hi++; }
//        if(criterion[i] == 0 && prediction[i] == 1)  { fa++; }
//        if(criterion[i] == 1 && prediction[i] == 0)  { mi++; }
//        if(criterion[i] == 0 && prediction[i] == 0)  { cr++; }
//      }
//
//      out[0] = hi;
//      out[1] = fa;
//      out[2] = mi;
//      out[3] = cr;
//
//       return out;
//}
