// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// updateFftree
S4 updateFftree(S4 fftree);
RcppExport SEXP fftrees_updateFftree(SEXP fftreeSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< S4 >::type fftree(fftreeSEXP );
        S4 __result = updateFftree(fftree);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// updateFftree2
S4 updateFftree2(S4 fftree, Rcpp::List cuelist);
RcppExport SEXP fftrees_updateFftree2(SEXP fftreeSEXP, SEXP cuelistSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< S4 >::type fftree(fftreeSEXP );
        Rcpp::traits::input_parameter< Rcpp::List >::type cuelist(cuelistSEXP );
        S4 __result = updateFftree2(fftree, cuelist);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}