// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "ProbitML.h"
#include "ProbitML_ARMA_KB.h"
//RNGScope scope;



// [[Rcpp::export]]
RcppExport SEXP ProbitMCMCHSD(SEXP i_Num_of_iterations, SEXP list_Data, SEXP list_InitialValues, SEXP list_HyperPara, SEXP list_UpdatePara, SEXP list_TuningPara)
{
    List lData(list_Data);
    List lInitialValues(list_InitialValues);
    List lHyperPara(list_HyperPara);
    List lUpdatePara(list_UpdatePara);
    List lTuningPara(list_TuningPara);
    
    List PosteriorSamples;
    
    int iNum_of_iterations = Rcpp::as<int> (i_Num_of_iterations);
    
    ProbitMLModelSelection DoMLModelSelectionMCMC(iNum_of_iterations, lData, lInitialValues, lHyperPara, lUpdatePara, lTuningPara);
    
    PosteriorSamples = DoMLModelSelectionMCMC.MCMC_Procedure();
    
    //Rcout << "Check1" << endl;
    return PosteriorSamples;
    
}


// [[Rcpp::export]]
RcppExport SEXP ProbitMCMCARMAKB(SEXP i_Num_of_iterations, SEXP list_Data, SEXP list_InitialValues, SEXP list_HyperPara, SEXP list_UpdatePara, SEXP list_TuningPara, SEXP ARMA_Order)
{
    List lData(list_Data);
    List lInitialValues(list_InitialValues);
    List lHyperPara(list_HyperPara);
    List lUpdatePara(list_UpdatePara);
    List lUpdateTuningPara(list_TuningPara);
    
    List PosteriorSamples;
    
    int iNum_of_iterations = Rcpp::as<int> (i_Num_of_iterations);
    
    vec vARMA_Order = as<vec>(ARMA_Order);
    
    ProbitMLModelSelectionARMAKB DoMLModelSelectionMCMC(iNum_of_iterations, lData, lInitialValues, lHyperPara, lUpdatePara, lUpdateTuningPara, vARMA_Order);
    
    PosteriorSamples = DoMLModelSelectionMCMC.MCMC_Procedure();
    
    //Rcout << "Check1" << endl;
    return PosteriorSamples;
    
}
