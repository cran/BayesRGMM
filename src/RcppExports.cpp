// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// ProbitMCMCHSD
RcppExport SEXP ProbitMCMCHSD(SEXP i_Num_of_iterations, SEXP list_Data, SEXP logic_Robust, SEXP list_InitialValues, SEXP list_HyperPara, SEXP list_UpdatePara, SEXP list_TuningPara, SEXP logic_Interactive);
RcppExport SEXP _BayesRGMM_ProbitMCMCHSD(SEXP i_Num_of_iterationsSEXP, SEXP list_DataSEXP, SEXP logic_RobustSEXP, SEXP list_InitialValuesSEXP, SEXP list_HyperParaSEXP, SEXP list_UpdateParaSEXP, SEXP list_TuningParaSEXP, SEXP logic_InteractiveSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type i_Num_of_iterations(i_Num_of_iterationsSEXP);
    Rcpp::traits::input_parameter< SEXP >::type list_Data(list_DataSEXP);
    Rcpp::traits::input_parameter< SEXP >::type logic_Robust(logic_RobustSEXP);
    Rcpp::traits::input_parameter< SEXP >::type list_InitialValues(list_InitialValuesSEXP);
    Rcpp::traits::input_parameter< SEXP >::type list_HyperPara(list_HyperParaSEXP);
    Rcpp::traits::input_parameter< SEXP >::type list_UpdatePara(list_UpdateParaSEXP);
    Rcpp::traits::input_parameter< SEXP >::type list_TuningPara(list_TuningParaSEXP);
    Rcpp::traits::input_parameter< SEXP >::type logic_Interactive(logic_InteractiveSEXP);
    rcpp_result_gen = Rcpp::wrap(ProbitMCMCHSD(i_Num_of_iterations, list_Data, logic_Robust, list_InitialValues, list_HyperPara, list_UpdatePara, list_TuningPara, logic_Interactive));
    return rcpp_result_gen;
END_RCPP
}
// ProbitMCMCARMAKB
RcppExport SEXP ProbitMCMCARMAKB(SEXP i_Num_of_iterations, SEXP list_Data, SEXP logic_Robust, SEXP list_InitialValues, SEXP list_HyperPara, SEXP list_UpdatePara, SEXP list_TuningPara, SEXP ARMA_Order, SEXP logic_Interactive);
RcppExport SEXP _BayesRGMM_ProbitMCMCARMAKB(SEXP i_Num_of_iterationsSEXP, SEXP list_DataSEXP, SEXP logic_RobustSEXP, SEXP list_InitialValuesSEXP, SEXP list_HyperParaSEXP, SEXP list_UpdateParaSEXP, SEXP list_TuningParaSEXP, SEXP ARMA_OrderSEXP, SEXP logic_InteractiveSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type i_Num_of_iterations(i_Num_of_iterationsSEXP);
    Rcpp::traits::input_parameter< SEXP >::type list_Data(list_DataSEXP);
    Rcpp::traits::input_parameter< SEXP >::type logic_Robust(logic_RobustSEXP);
    Rcpp::traits::input_parameter< SEXP >::type list_InitialValues(list_InitialValuesSEXP);
    Rcpp::traits::input_parameter< SEXP >::type list_HyperPara(list_HyperParaSEXP);
    Rcpp::traits::input_parameter< SEXP >::type list_UpdatePara(list_UpdateParaSEXP);
    Rcpp::traits::input_parameter< SEXP >::type list_TuningPara(list_TuningParaSEXP);
    Rcpp::traits::input_parameter< SEXP >::type ARMA_Order(ARMA_OrderSEXP);
    Rcpp::traits::input_parameter< SEXP >::type logic_Interactive(logic_InteractiveSEXP);
    rcpp_result_gen = Rcpp::wrap(ProbitMCMCARMAKB(i_Num_of_iterations, list_Data, logic_Robust, list_InitialValues, list_HyperPara, list_UpdatePara, list_TuningPara, ARMA_Order, logic_Interactive));
    return rcpp_result_gen;
END_RCPP
}
// CumulativeProbitMCMC
RcppExport SEXP CumulativeProbitMCMC(SEXP i_Num_of_iterations, SEXP list_Data, SEXP logic_Robust, SEXP list_InitialValues, SEXP list_HyperPara, SEXP list_UpdatePara, SEXP list_TuningPara, SEXP logic_Interactive);
RcppExport SEXP _BayesRGMM_CumulativeProbitMCMC(SEXP i_Num_of_iterationsSEXP, SEXP list_DataSEXP, SEXP logic_RobustSEXP, SEXP list_InitialValuesSEXP, SEXP list_HyperParaSEXP, SEXP list_UpdateParaSEXP, SEXP list_TuningParaSEXP, SEXP logic_InteractiveSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type i_Num_of_iterations(i_Num_of_iterationsSEXP);
    Rcpp::traits::input_parameter< SEXP >::type list_Data(list_DataSEXP);
    Rcpp::traits::input_parameter< SEXP >::type logic_Robust(logic_RobustSEXP);
    Rcpp::traits::input_parameter< SEXP >::type list_InitialValues(list_InitialValuesSEXP);
    Rcpp::traits::input_parameter< SEXP >::type list_HyperPara(list_HyperParaSEXP);
    Rcpp::traits::input_parameter< SEXP >::type list_UpdatePara(list_UpdateParaSEXP);
    Rcpp::traits::input_parameter< SEXP >::type list_TuningPara(list_TuningParaSEXP);
    Rcpp::traits::input_parameter< SEXP >::type logic_Interactive(logic_InteractiveSEXP);
    rcpp_result_gen = Rcpp::wrap(CumulativeProbitMCMC(i_Num_of_iterations, list_Data, logic_Robust, list_InitialValues, list_HyperPara, list_UpdatePara, list_TuningPara, logic_Interactive));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_BayesRGMM_ProbitMCMCHSD", (DL_FUNC) &_BayesRGMM_ProbitMCMCHSD, 8},
    {"_BayesRGMM_ProbitMCMCARMAKB", (DL_FUNC) &_BayesRGMM_ProbitMCMCARMAKB, 9},
    {"_BayesRGMM_CumulativeProbitMCMC", (DL_FUNC) &_BayesRGMM_CumulativeProbitMCMC, 8},
    {NULL, NULL, 0}
};

RcppExport void R_init_BayesRGMM(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
