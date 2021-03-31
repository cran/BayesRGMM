//
//  ProbitML.cpp
//  
//
//  Created by kuojung on 2020/2/26.
//

// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "ProbitML_ARMA.h"
//#include "tmvrnormGibbs.h"
#include "tmvrnormGibbs_KJLEE.h"
//RNGScope scope;
//#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]


ProbitMLModelSelectionARMA::ProbitMLModelSelectionARMA(int iNum_of_iterations, List list_Data, List list_InitialValues, List list_HyperPara, List list_UpdatePara, List list_TuningPara, vec vARMA_Order)
{
    Num_of_iterations = iNum_of_iterations;
    Data = list_Data;
    InitialValues = list_InitialValues;
    HyperPara = list_HyperPara;
    UpdatePara = list_UpdatePara;
    TuningPara = list_TuningPara;
    
    phi_tune = as<double>(TuningPara["TuningPhi"]);
    psi_tune = as<double>(TuningPara["TuningPsi"]);
    
    ARMA_Order = vARMA_Order;
    
    //Rcout << "ARMA order = " << ARMA_Order << endl;
    
    //Rcout<< "Read Data" << endl;
    Y = as<mat>(Data["Y"]);
    X = as<cube>(Data["X"]);
    Z = as<cube>(Data["Z"]);

    TimePointsAvailable = as<vec>(Data["TimePointsAvailable"]);
    
    
    Num_of_obs = Y.n_cols;
    Num_of_Timepoints = Y.n_rows;
    Num_of_RanEffs = Z.n_cols;
    Num_of_covariates = X.n_cols;
    
    updateystar = as<bool>(UpdatePara["UpdateYstar"]);
    updateb = as<bool>(UpdatePara["UpdateRandomEffect"]);
    updatenu = as<bool>(UpdatePara["UpdateNu"]);
    updatebeta = as<bool>(UpdatePara["UpdateBeta"]);
    updateSigma = as<bool>(UpdatePara["UpdateSigma"]);
    updatephi = as<bool>(UpdatePara["UpdatePhi"]);
    updatepsi = as<bool>(UpdatePara["UpdatePsi"]);
    
    SinglePhiPsi = as<bool>(UpdatePara["SinglePhiPsi"]);
    //Unconstraint = as<bool>(UpdatePara["Unconstraint"]);
    //updateomega = as<bool>(UpdatePara["UpdateOmega"]);
    //SIMPLE = as<bool>(UpdatePara["SIMPLE"]);

    
    b_samples.set_size(Num_of_RanEffs, Num_of_obs, Num_of_iterations);
    nu_samples.set_size(Num_of_obs, Num_of_iterations);
    beta_samples.set_size(Num_of_covariates,Num_of_iterations);
    Sigma_samples.set_size(Num_of_RanEffs, Num_of_RanEffs, Num_of_iterations);
    
    phi_samples.set_size(ARMA_Order(0), Num_of_obs, Num_of_iterations);
    psi_samples.set_size(ARMA_Order(1), Num_of_obs, Num_of_iterations);
    
    //phi_samples.zeros();
    //psi_samples.zeros();
    
    omega_samples.set_size(Num_of_Timepoints, Num_of_Timepoints, Num_of_iterations);
    b_mean.set_size(Num_of_RanEffs, Num_of_obs);
    nu_mean.set_size(Num_of_obs);
    beta_mean.set_size(Num_of_covariates);
    Sigma_mean.set_size(Num_of_RanEffs, Num_of_RanEffs);
    
    phi_mean.set_size(ARMA_Order(0), Num_of_obs);
    psi_mean.set_size(ARMA_Order(1), Num_of_obs);

    //omega_samples.set_size(Num_of_RanEffs, Num_of_RanEffs, Num_of_iterations);
    
    //Rcout<< "Initial Values Y" << endl;
    //Rcout << "size(Y)" << size(Y) << endl;
    //Y_star_sample.set_size(size(Y));
    //Y_star_sample.zeros();
    //Rcout << "Y_star = " << Y_star_sample.n_rows << "\t" << Y_star_sample.n_cols << endl;
    //Rcout<< "Initial Values y.star" << endl;
    Y_star_sample = as<mat>(InitialValues["y.star"]);
    
    //Rcout<< Y_star_sample.submat(0, 0, 4, 4) << endl;
    //Rcout<< "Initial Values b" << endl;
    b_samples.slice(0) = as<mat>(InitialValues["b"]);
    //Rcout<< "Initial Values nu" << endl;
    nu_samples.col(0) = as<vec>(InitialValues["nu"]);
    //Rcout<< "Initial Values beta" << endl;
    beta_samples.col(0) = as<vec>(InitialValues["beta"]);
    //Rcout<< "Initial Values Sigma" << endl;
    Sigma_samples.slice(0) = as<mat>(InitialValues["Sigma"]);
    
    //Rcout<< "Initial Values phi" << endl;
    phi_samples.slice(0) = as<mat>(InitialValues["phi"]);
    //Rcout<< "Initial Values psi" << endl;
    psi_samples.slice(0) = as<mat>(InitialValues["psi"]);

    

    //Rcout<< "Read Hyperparameters." << endl;
    // Hyperparameters
    v_gamma = as<double>(HyperPara["v.gamma"]);
    sigma2_beta = as<double>(HyperPara["sigma2.beta"]);

    Vb = as<double>(HyperPara["InvWishart.df"]);
    Lambda = as<mat>(HyperPara["InvWishart.Lambda"]);
    
    Ib_diag.eye(Num_of_RanEffs, Num_of_RanEffs);

    
    //vec v = { -1, -2, 1, 3}, del = {-0.5, -0.3};
    //vec v_lower = {datum::inf, datum::inf, 0, 0}, v_upper= {0, 0, -datum::inf, -datum::inf};
    //Rcout << rtmvnorm_gibbs_KJLEE(1, v, Ri(1, 4, del), v_lower, v_upper, 0, zeros<vec>(4), 1) << endl;

    acc_phi_rate = 0.;
    acc_psi_rate = 0.;
}


/*
mat ProbitMLModelSelectionARMA::CovARMA(int i, int tp, double phi, double psi)
{
    //Rcout << "ARMA" << endl;
    mat Phi = eye(tp, tp);
    mat Psi = eye(tp, tp);
    mat CovARMA_tmp;
    for(int i=0; i< (tp-1); i++){
        Phi(i+1, i) = -phi;
        Psi(i+1, i) = psi;
    }

        //Rcout << "Phi = " << endl << Phi << endl;
        //Rcout << "Psi = " << endl << Psi << endl;

    
    CovARMA_tmp = Phi.i()* Psi* Psi.t()*(Phi.t().i());
    return CovARMA_tmp;
}
*/

mat ProbitMLModelSelectionARMA::CovARMA(int tp, vec phi, vec psi)
{
    
    mat Phi = eye(tp, tp);
    mat Psi = eye(tp, tp);
    mat CovARMA_tmp;
    for(int t=1; t<tp; t++){
        for(int j=(t-1); (t-j)<=ARMA_Order(0) && j>=0; j--)
            Phi(t, j) = -phi(t-j-1);
        for(int j=(t-1); (t-j)<=ARMA_Order(1) && j>=0; j--)
            Psi(t, j) = psi(t-j-1);
    }
    
    CovARMA_tmp = Phi.i()* Psi* Psi.t()*(Phi.t().i());
    return CovARMA_tmp;
}


void ProbitMLModelSelectionARMA::Update_ystar_b_nu_beta_Sigma(int iter)
{
    //if(iter % 100 == 0)
    //    Rcout << "Update ystar, b, nu, beta, Sigmab simultaneously" << endl;
    vec mu_tmp_b, res_b, mu_tmp_ystar;

    vec res_beta, mu_tmp_beta = zeros<vec>(Num_of_covariates);

    mat Omegai_tmp, Omegai_inv;

    mat Sigma_tmp_beta = zeros<mat>(Num_of_covariates, Num_of_covariates);
    mat Sigma_tmp_b = zeros<mat>(Num_of_RanEffs, Num_of_RanEffs);
    mat Sigma_tmp = zeros<mat>(Num_of_RanEffs, Num_of_RanEffs);

    //vec v = { -1, -2, 1, 3};
    //vec v_lower = {datum::inf, datum::inf, 0, 0}, v_upper= {0, 0, -datum::inf, -datum::inf};
    //rtmvnorm_gibbs_KJLEE(1, v, Ri(1, 4, delta_samples.col(iter)), v_lower, v_upper, 100, zeros<vec>(4), 5);

    int tp;
    //vec y_star_tmp;
    vec lower, upper;
    mat X_tmp, Z_tmp;
    //vec  init_state;
    vec b_vec;
    double beta_tmp;
    //double sigma_tmp, mu_tmp_1;
    double alpha_tmp = 0.5*(v_gamma + Num_of_RanEffs);
    
    //Rcout << "iter = " << iter << endl;
    for(int i=0; i<Num_of_obs; i++){
        //Rcout << "iter = " << iter << "\t i = " << i << endl;
        tp = TimePointsAvailable(i);

        X_tmp = X(span(0, tp-1), span(0, Num_of_covariates-1), span(i));
        Z_tmp = (Z.slice(i).rows(0, tp-1));
        mu_tmp_ystar = X_tmp*beta_samples.col(iter) + Z_tmp*b_samples.slice(iter).col(i);
        if(SinglePhiPsi)
            Omegai_tmp = CovARMA(tp, phi_samples.slice(iter).col(0), psi_samples.slice(iter).col(0));
        else
            Omegai_tmp = CovARMA(tp, phi_samples.slice(iter).col(i), psi_samples.slice(iter).col(i));
        
        //Rcout << "Omega = " << Omegai_tmp << endl;
        Omegai_inv = inv_sympd(Omegai_tmp);
        //Rcout << "Omega_inv " << Omegai_inv << endl;
        
        
        lower.set_size(tp);
        upper.set_size(tp);
        lower.elem(find(Y(span(0, tp-1), i)>0)).zeros();
        lower.elem(find(Y(span(0, tp-1), i)==0)).ones();
        lower.elem(find(Y(span(0, tp-1), i)==0)) *= -datum::inf;
        
        upper.elem(find(Y(span(0, tp-1), i)==0)).zeros();
        upper.elem(find(Y(span(0, tp-1), i)>0)).ones();
        upper.elem(find(Y(span(0, tp-1), i)>0)) *= datum::inf;
        //if(i==368)
            //Rcout << "iter = " << iter << "\t Y = " << Y.col(i).t() << "\t lower = " << lower.t() << "\t upper " << upper.t() << endl;
        
        mu_tmp_ystar.elem( find( ( (Y(span(0, tp-1), i)-1) % mu_tmp_ystar) < 0 )).zeros();
        mu_tmp_ystar.elem( find( ( (Y(span(0, tp-1), i) ) % mu_tmp_ystar) < 0 )).zeros();

        if(tp == 1)
            Y_star_sample(0, i) = rtruncnorm(1, as_scalar(mu_tmp_ystar), as_scalar(Omegai_tmp), as_scalar(lower),  as_scalar(upper))(0);
        else
            Y_star_sample(span(0, tp-1), i) = rtmvnorm_gibbs_KJLEE(1, mu_tmp_ystar, Omegai_tmp, lower, upper, 100, zeros<vec>(tp), 5).t();

        if(Y_star_sample.col(i).has_nan())
            Y_star_sample.col(i).zeros();
        Y_star_sample.col(i) = clamp(Y_star_sample.col(i), -5, 5);
        //Y_star_sample.col(i).replace(find(Y_star_sample.col(i))>5, 5);
        //Y_star_sample.col(i).replace(find(Y_star_sample.col(i))<-5, -5);
        
        if(Y_star_sample(span(0, tp-1), i).has_nan()){
            Rcout << "iter = " << iter << "\t i = " << i << endl;
            Rcout << "mu_tmp_ystar = " << mu_tmp_ystar << endl;
            Rcout << "Omegai_tmp = \n" << Omegai_tmp << endl;
            Rcout << "Y(span(0, tp-1), i) = " << Y(span(0, tp-1), i) << endl;
            Rcout << "Y_star_sample(span(0, tp-1), i)=" << Y_star_sample(span(0, tp-1), i) << endl;
        }
            
        
        //Y_star_sample(span(0, tp-1), i) = rtmvnorm_gibbs_KJLEE(1, mu_tmp_ystar, Omegai_tmp, lower, upper, 100, zeros<vec>(tp), 5).t();
        
        //if(iter==39 || iter==40)
        //Rcout << "Y_star_sample = " << Y_star_sample(span(0, tp-1), i) << endl;
        
        //if(i==368){
        //Rcout << "mu_tmp_ystar = " << mu_tmp_ystar.t() << endl;
        //Rcout << "Y.star = " << Y_star_sample(span(0, tp-1), i).t() << endl;
        //}
        
        //Rcout << "Sigma_tmp_b" <<endl;
        
        Sigma_tmp_b = (Z_tmp.t()*Omegai_inv*Z_tmp+nu_samples(i, iter)*Sigma_samples.slice(iter).i());
        if(!Sigma_tmp_b.is_sympd())
            Sigma_tmp_b.eye();
        else
            Sigma_tmp_b = inv_sympd(Sigma_tmp_b);
  
        res_b = Y_star_sample(span(0, tp-1), i) - X_tmp*beta_samples.col(iter);
        mu_tmp_b = Sigma_tmp_b*Z_tmp.t()*Omegai_inv*res_b;
        
        
        if(Sigma_tmp_b.has_nan()){
            Sigma_tmp_b.eye();
            mu_tmp_b.zeros();
        }

        b_samples.slice(iter+1).col(i) = mvnrnd(mu_tmp_b, Sigma_tmp_b);
        b_vec = b_samples.slice(iter+1).col(i);
        beta_tmp = 0.5*(as_scalar(b_vec.t()*Sigma_samples.slice(iter).i()*b_vec) + v_gamma);
        
        
        nu_samples(i, iter+1) = randg( 1, distr_param(alpha_tmp, 1./beta_tmp))(0);
        //nu_samples(i, iter+1) = randg( 1, distr_param(1., 1.))(0);
        
        
        

        Sigma_tmp_beta += (X_tmp.t()*Omegai_inv*X_tmp);
        
        //Rcout << "Sigma_tmp =" << Sigma_tmp << endl;
        //yi_star =  y_star_samples.slice(iter).col(i);
        //res = Y_star_sample((span(0, tp-1), i))- (Z.slice(i))*b_samples.slice(iter).col(i);
        
        res_beta = Y_star_sample(span(0, tp-1), i)- Z_tmp*b_vec;
        
        //Rcout << "i = " << i << "\t res_beta = " << res_beta << endl;
        mu_tmp_beta += X_tmp.t()*Omegai_inv*res_beta;
        
        //b_vec = b_samples.slice(iter+1).col(i);
        //if(i==368)
            //Rcout << i << "\t" << "b_vec = "<< b_vec << endl;
        Sigma_tmp += nu_samples(i, iter+1)*(b_vec*b_vec.t());
        
      }
    
    //Rcout << "Sigma_tmp_beta = " << sigma2_beta << endl;
    //Rcout << "mu_tmp_beta = " << mu_tmp_beta << endl;
    
    Sigma_tmp_beta.diag() += 1./sigma2_beta;
    Sigma_tmp_beta = inv_sympd(Sigma_tmp_beta);
    mu_tmp_beta = Sigma_tmp_beta * mu_tmp_beta;
    
    beta_samples.col(iter+1) = mvnrnd(mu_tmp_beta, Sigma_tmp_beta);
    
    
    //Rcout << "Sigma_tmp = " << Sigma_tmp << endl;
    
    Sigma_tmp = (Sigma_tmp + Lambda);

    //Rcout << "Sigma_tmp_inv = " << Sigma_tmp << "\t Num_of_obs + Vb =" << Num_of_obs + Vb << endl;
    Sigma_samples.slice(iter+1) =iwishrnd( Sigma_tmp, (Num_of_obs + Vb)); //1./randg( 1, distr_param(gamma_shape, gamma_scale) );// //riwish(df,
}



void ProbitMLModelSelectionARMA::Update_phi(int iter)
{
    if(iter % 100 == 0)
    Rcout << "iter = " << iter << " Update phi" << endl;
    double phi_den = 0., phi_num = 0.;
    vec phi_cand;
    
    vec res;
    mat X_tmp, Z_tmp;
    int tp;
    mat Cov_tmp_inv;
    for(int i=0; i<Num_of_obs; i++){
        //phi_cand = mvnrnd(phi_samples.slice(iter).col(i), 0.1*eye(ARMA_Order(0),ARMA_Order(0)));
        
        if(ARMA_Order(0)==1){
            do{
                phi_cand = mvnrnd(phi_samples.slice(iter).col(i), phi_tune*eye(ARMA_Order(0),ARMA_Order(0)));
                //phi_cand =  phi_samples(i, iter)+0.01*(2*randu()-1);
            }
            while(abs(as_scalar(phi_cand))>1);
        }
        else{
            do{
                //phi_cand =  phi_samples(i, iter)+0.01*(2*randu()-1);
                phi_cand = mvnrnd(phi_samples.slice(iter).col(i), phi_tune*eye(ARMA_Order(0),ARMA_Order(0)));
                //bool_phi = (sum(phi_cand)<1) && (as_scalar(diff(phi_cand)<1)) && (abs(phi_cand(1))<1);
            }
            while((sum(phi_cand)>1) || (as_scalar(diff(phi_cand)>1)) || (abs(phi_cand(1))>1));
        }
        
        tp = TimePointsAvailable(i);
        X_tmp = X.slice(i).rows(0, tp-1);  //X(span(0, tp-1), span(0, Num_of_covariates-1), span(i));
        Z_tmp = Z.slice(i).rows(0, tp-1);
        
        res = Y_star_sample(span(0, tp-1), i) - X_tmp*beta_samples.col(iter+1)-Z_tmp*b_samples.slice(iter+1).col(i);
        
        Cov_tmp_inv = inv_sympd( CovARMA(tp, phi_samples.slice(iter).col(i), psi_samples.slice(iter).col(i)));
        phi_den = 0.5*log(det(Cov_tmp_inv)) - 0.5*as_scalar(res.t()* Cov_tmp_inv*res);
        
        Cov_tmp_inv = inv_sympd( CovARMA(tp, phi_cand, psi_samples.slice(iter).col(i)));
        phi_num = 0.5*log(det(Cov_tmp_inv)) - 0.5*as_scalar(res.t()* Cov_tmp_inv*res);
        
        //Rcout << "phi_num = " << phi_num << "\t" << "phi_den = " << phi_den << "\tphi_num - phi_den = " << phi_num - phi_den << endl;
        
        if(log(randu()) < phi_num - phi_den ){
            phi_samples.slice(iter+1).col(i) = phi_cand;
        }
        else
            phi_samples.slice(iter+1).col(i) = phi_samples.slice(iter).col(i);
        
    }
    
    
}

void ProbitMLModelSelectionARMA::Update_single_phi(int iter)
{
    //if(iter % 100 == 0)
    //    Rcout << "iter = " << iter << " Update single phi" << endl;
    double phi_den = 0., phi_num = 0.;
    vec phi_cand;
    
    vec res;
    mat X_tmp, Z_tmp;
    int tp;
    mat Cov_tmp_inv;
    
    if(ARMA_Order(0)==1){
        do{
            phi_cand = mvnrnd(phi_samples.slice(iter).col(0), phi_tune*eye(ARMA_Order(0),ARMA_Order(0)));
            //phi_cand =  phi_samples(i, iter)+0.01*(2*randu()-1);
        }
        while(abs(as_scalar(phi_cand))>1);
    }
    else{
        do{
            //phi_cand =  phi_samples(i, iter)+0.01*(2*randu()-1);
            phi_cand = mvnrnd(phi_samples.slice(iter).col(0), phi_tune*eye(ARMA_Order(0),ARMA_Order(0)));
            //bool_phi = (sum(phi_cand)<1) && (as_scalar(diff(phi_cand)<1)) && (abs(phi_cand(1))<1);
        }
        while((sum(phi_cand)>1) || (as_scalar(diff(phi_cand)>1)) || (abs(phi_cand(1))>1));
    }

    for(int i=0; i<Num_of_obs; i++){
        //phi_cand = mvnrnd(phi_samples.slice(iter).col(i), 0.1*eye(ARMA_Order(0),ARMA_Order(0)));
        
        
        tp = TimePointsAvailable(i);
        X_tmp = X.slice(i).rows(0, tp-1);  //X(span(0, tp-1), span(0, Num_of_covariates-1), span(i));
        Z_tmp = Z.slice(i).rows(0, tp-1);
        
        res = Y_star_sample(span(0, tp-1), i) - X_tmp*beta_samples.col(iter+1)-Z_tmp*b_samples.slice(iter+1).col(i);
        
        Cov_tmp_inv = inv_sympd( CovARMA(tp, phi_samples.slice(iter).col(0), psi_samples.slice(iter).col(0)));
        phi_den += - 0.5*as_scalar(res.t()* Cov_tmp_inv*res);
        
        Cov_tmp_inv = inv_sympd( CovARMA(tp, phi_cand, psi_samples.slice(iter).col(0)));
        phi_num += - 0.5*as_scalar(res.t()* Cov_tmp_inv*res);
        
        //Rcout << "phi_num = " << phi_num << "\t" << "phi_den = " << phi_den << "\tphi_num - phi_den = " << phi_num - phi_den << endl;
    }
    if(log(randu()) < phi_num - phi_den ){
        phi_samples.slice(iter+1).col(0) = phi_cand;
        acc_phi_rate++;
    }
    else
        phi_samples.slice(iter+1).col(0) = phi_samples.slice(iter).col(0);
    
    if((iter+1)%500 == 0){
        //Rcout << "tuning_delta = " << tuning_delta << endl;
        //Rcout << "acc_rate_delta/iter = " << acc_rate_delta/iter << endl;
        if( acc_phi_rate/iter<0.25 )
            phi_tune = phi_tune/2.;
        if( (1.*acc_phi_rate)/iter>0.50 )
            phi_tune = 2*phi_tune;
        
        //Rcout << "tuning_delta = " << tuning_delta << endl;
        
    }
    //Rcout <<"Done for phi" << endl;
}

void ProbitMLModelSelectionARMA::Update_psi(int iter)
{
    if(iter % 100 == 0)
        Rcout << "iter = " << iter << " Update psi" << endl;
    double psi_den = 0., psi_num = 0.;
    vec psi_cand;
    
    vec res;
    mat X_tmp, Z_tmp;
    int tp;
    mat Cov_tmp_inv;
    
    for(int i=0; i<Num_of_obs; i++){
        
        if(ARMA_Order(1)==1){
            do{
                psi_cand = mvnrnd(psi_samples.slice(iter).col(i), psi_tune*eye(ARMA_Order(1),ARMA_Order(1)));
                //phi_cand =  phi_samples(i, iter)+0.01*(2*randu()-1);
            }
            while(abs(as_scalar(psi_cand))>1);
        }
        else{
            do{
                //phi_cand =  phi_samples(i, iter)+0.01*(2*randu()-1);
                psi_cand = mvnrnd(psi_samples.slice(iter).col(i), psi_tune*eye(ARMA_Order(1),ARMA_Order(1)));
                //bool_phi = (sum(phi_cand)<1) && (as_scalar(diff(phi_cand)<1)) && (abs(phi_cand(1))<1);
            }
            while((sum(psi_cand)>1) || (as_scalar(diff(psi_cand)>1)) || (abs(psi_cand(1))>1));
        }

        
        tp = TimePointsAvailable(i);
        X_tmp = X.slice(i).rows(0, tp-1);  //X(span(0, tp-1), span(0, Num_of_covariates-1), span(i));
        Z_tmp = Z.slice(i).rows(0, tp-1);
        
        res = Y_star_sample(span(0, tp-1), i) - X_tmp*beta_samples.col(iter+1)-Z_tmp*b_samples.slice(iter+1).col(i);
        
        Cov_tmp_inv = inv_sympd( CovARMA(tp, phi_samples.slice(iter+1).col(i), psi_samples.slice(iter).col(i)));
        psi_den = 0.5*log(det(Cov_tmp_inv)) - 0.5*as_scalar(res.t()* Cov_tmp_inv*res);
        
        Cov_tmp_inv = inv_sympd( CovARMA(tp, phi_samples.slice(iter+1).col(i), psi_cand));
        psi_num = 0.5*log(det(Cov_tmp_inv)) - 0.5*as_scalar(res.t()* Cov_tmp_inv*res);
        
        //Rcout << "phi_num = " << phi_num << "\t" << "phi_den = " << phi_den << "\tphi_num - phi_den = " << phi_num - phi_den << endl;
        
        if(log(randu()) < psi_num - psi_den ){
            psi_samples.slice(iter+1).col(i) = psi_cand;
        }
        else
            psi_samples.slice(iter+1).col(i) = psi_samples.slice(iter).col(i);
        
    }
}

void ProbitMLModelSelectionARMA::Update_single_psi(int iter)
{
    //if(iter % 100 == 0)
    //    Rcout << "iter = " << iter << " Update single psi" << endl;
    double psi_den = 0., psi_num = 0.;
    vec psi_cand;
    
    vec res;
    mat X_tmp, Z_tmp;
    int tp;
    mat Cov_tmp_inv;
    
    
    if(ARMA_Order(1)==1){
        do{
            psi_cand = mvnrnd(psi_samples.slice(iter).col(0), psi_tune*eye(ARMA_Order(1),ARMA_Order(1)));
            //phi_cand =  phi_samples(i, iter)+0.01*(2*randu()-1);
        }
        while(abs(as_scalar(psi_cand))>1);
    }
    else{
        do{
            //phi_cand =  phi_samples(i, iter)+0.01*(2*randu()-1);
            psi_cand = mvnrnd(psi_samples.slice(iter).col(0), psi_tune*eye(ARMA_Order(1),ARMA_Order(1)));
            //bool_phi = (sum(phi_cand)<1) && (as_scalar(diff(phi_cand)<1)) && (abs(phi_cand(1))<1);
        }
        while((sum(psi_cand)>1) || (as_scalar(diff(psi_cand)>1)) || (abs(psi_cand(1))>1));
    }
    
    for(int i=0; i<Num_of_obs; i++){
        
        tp = TimePointsAvailable(i);
        X_tmp = X.slice(i).rows(0, tp-1);  //X(span(0, tp-1), span(0, Num_of_covariates-1), span(i));
        Z_tmp = Z.slice(i).rows(0, tp-1);
        
        res = Y_star_sample(span(0, tp-1), i) - X_tmp*beta_samples.col(iter+1)-Z_tmp*b_samples.slice(iter+1).col(0);
        
        Cov_tmp_inv = inv_sympd( CovARMA(tp, phi_samples.slice(iter+1).col(0), psi_samples.slice(iter).col(0)));
        psi_den += - 0.5*as_scalar(res.t()* Cov_tmp_inv*res);
        
        Cov_tmp_inv = inv_sympd( CovARMA(tp, phi_samples.slice(iter+1).col(0), psi_cand));
        psi_num += - 0.5*as_scalar(res.t()* Cov_tmp_inv*res);
   
        //Rcout << "phi_num = " << phi_num << "\t" << "phi_den = " << phi_den << "\tphi_num - phi_den = " << phi_num - phi_den << endl;
    }
    if(log(randu()) < psi_num - psi_den ){
        psi_samples.slice(iter+1).col(0) = psi_cand;
        acc_psi_rate++;
    }
    else
        psi_samples.slice(iter+1).col(0) = psi_samples.slice(iter).col(0);
    
    //Rcout << "Done for psi" << endl;
    if( acc_psi_rate/iter<0.25 )
        psi_tune = psi_tune/2.;
    if( (1.*acc_psi_rate)/iter>0.50 )
        psi_tune = 2*psi_tune;
}


void ProbitMLModelSelectionARMA::ParameterEstimation()
{

    b_mean = mean(b_samples, 2);
    Sigma_mean = mean(Sigma_samples, 2);
    beta_mean = mean(beta_samples, 1);
    nu_mean = mean(nu_samples, 1);
    phi_mean = mean(phi_samples, 2);
    psi_mean = mean(psi_samples, 2);
    
    //Rcout << size(b_mean) << endl;
    //rowvec X_tmp, Z_tmp;
    
    rowvec X_tmp, Z_tmp;
    //vec mu_tmp;
    double pit, CPO_tmp, ESS=0, GP=0, ESS_GP_tmp, RJ1, RJ2;
    logL = 0.;
    
    mat Djt(Num_of_Timepoints, Num_of_covariates, fill::zeros);
    mat Omega_I(Num_of_covariates, Num_of_covariates, fill::zeros), M_LZ(Num_of_covariates, Num_of_covariates, fill::zeros);
    vec mu_it(Num_of_Timepoints), p_it(Num_of_Timepoints);
    mat A_sqrt, Cov_Y, V, V_inv, Omega_I_inv, V_LZ, Gamma_RJ, Djt_sub;

    int tp;
    
    CIC = 0.;
    RJ_R = 0.;
    ACC = 0.;
    
    mat CPO = zeros<mat>(Num_of_obs, TimePointsAvailable.max());


    for(int i=0; i<Num_of_obs; i++){
        //Rcout << "i = " << i << endl;
        tp = TimePointsAvailable(i);
        for(int t=0; t<tp; t++){
            
            X_tmp = X.slice(i).row(t);
            Z_tmp = Z.slice(i).row(t);
            //Rcout << "X_tmp = " << X_tmp << endl;
            //Rcout << "beta_mean = " << beta_mean << endl;
            mu_it(t) = as_scalar( X_tmp*beta_mean + Z_tmp*b_mean.col(i));
            p_it(t) = normcdf(mu_it(t));
            
            ACC += 1.*(1.*(mu_it(t)>0) == Y(t, i));
            
            for(int j=0; j<Num_of_covariates; j++)
                Djt(t, j) = X(t,j,i)/normpdf(Rf_pnorm5(mu_it(t), 0., 1., 1, 0));
        }
        //Rcout << "============== 1 ============"<<endl;
        Djt_sub = Djt.head_rows(tp);
        //Rcout << "============== 2 ============"<<endl;
        Cov_Y = (Y(span(0, tp-1), i)-p_it.head(tp))*(Y(span(0, tp-1), i)-p_it.head(tp)).t();
        //Rcout << "============== 3 ============" << endl;
        A_sqrt = diagmat(sqrt(p_it.head(tp)));
        V = A_sqrt*A_sqrt;
        //Rcout << "============== 4 ============"<<endl;
        V_inv = V.i();
        Omega_I += Djt_sub.t()*V_inv*Djt_sub;
        //Rcout << "============== 5 ============"<<endl;
        M_LZ += Djt_sub.t()*V_inv*Cov_Y*V_inv*Djt_sub;
        //Rcout << "============== 6 ============"<<endl;
        ESS_GP_tmp = as_scalar((Y(span(0, tp-1), i)-mu_it.head(tp)).t()*V_inv*(Y(span(0, tp-1), i)-mu_it.head(tp)));
        //Rcout << "============== 7 ============"<<endl;
        ESS += ESS_GP_tmp;
        GP += -0.5*ESS_GP_tmp + log(det(V));
    }
    Omega_I_inv = Omega_I.i();
    Gamma_RJ = Omega_I_inv*M_LZ;
    V_LZ = Gamma_RJ*Omega_I_inv;
    
    RJ1 = trace(Gamma_RJ)/Num_of_covariates;
    RJ2 = accu((diagvec(Gamma_RJ*Gamma_RJ)))/Num_of_covariates;

    RJ_R = sqrt((1-RJ1)*(1-RJ1)+(1-RJ2)*(1-RJ2));
    //SC = ESS/(N-P-a
    //GP = -0.5*GP

    CIC = trace(Omega_I*V_LZ);
    //Rcout << "RJ_R = " << RJ_R << "\tCIC = " << CIC << endl;
    //cat("RJ.R = ", RJ.R, "\t", "SC = ", SC, "\n")

    
    for(int i=0; i<Num_of_obs; i++){
        //Rcout << "i=" << i << endl;
        //tp = TimePointsAvailable(i);
        //Rcout << "tp = " << tp << endl;
        //tp = TimePointsAvailable(i);
        //X_tmp = X.slice(i).rows(0, tp-1);  //X(span(0, tp-1), span(0, Num_of_covariates-1), span(i));
        //Z_tmp = Z.slice(i).rows(0, tp-1);
        //mu_tmp = X_tmp*beta_mean+Z_tmp*b_mean.col(i);
        //Ri_tmp = Ri(i, tp, delta_mean);

        for(int t=0; t<TimePointsAvailable(i); t++){
            X_tmp = X.slice(i).row(t); //(span(t), span(0, Num_of_covariates-1)), span(i));
            Z_tmp = Z.slice(i).row(t);
            //mu_tmp = X_tmp*beta_mean+Z_tmp*b_mean.col(i);
            //Rcout << X_tmp << endl;
            //Rcout << Z_tmp << endl;
            //Rcout << "beta = " << endl << beta_mean << endl;
            //Rcout << "b = " << endl << b_mean.col(i) << endl;
            //cond_mu = mu_tmp(t) + negSubRow(Ri_tmp.row(t)) * inv_sympd( sub1(Ri_tmp, t) ) * ()
            
            for(int iter = Num_of_iterations/2; iter<Num_of_iterations; iter++){
                pit = normcdf( as_scalar( X_tmp*beta_samples.col(iter) + Z_tmp*b_samples.slice(iter).col(i) ), 0., 1.);
                if(pit == 1 && Y(t, i) == 1){
                    DIC += 0.;
                    CPO_tmp = 0.;
                }
                else if(pit == 0 && Y(t, i) == 0){
                    DIC += 0.;
                //else if(pit == 0 && Y(t, i) == 1)
                //    Likelihood += 0.;
                //else if(pit == 1 && Y(t, i) == 0)
                //    Likelihood += 0.;
                    CPO_tmp = 0.;
                }
                else{
                    //Likelihood *= pow(pit, Y(t, i))*pow( (1-pit), (1-Y(t, i)) );
                    CPO_tmp = Y(t, i)*log(pit) + (1-Y(t, i))*log(1-pit);
                    DIC += CPO_tmp;
                    
                }
                
                CPO(i, t) += exp(-CPO_tmp);
            }
    
            
            
            pit = normcdf(as_scalar( X_tmp*beta_mean + Z_tmp*b_mean.col(i) ), 0., 1.);
            //Rcout << "i=" << i << ", t=" << t << "\tpit=" << pit << endl;
            
            if(pit == 1 && Y(t, i) == 1)
                logL += 0.;
            else if(pit == 0 && Y(t, i) == 0)
                logL += 0.;
            //else if(pit == 0 && Y(t, i) == 1)
            //    Likelihood += 0.;
            //else if(pit == 1 && Y(t, i) == 0)
            //    Likelihood += 0.;
            else
                //Likelihood *= pow(pit, Y(t, i))*pow( (1-pit), (1-Y(t, i)) );
                logL += Y(t, i)*log(pit) + (1-Y(t, i))*log(1-pit);
        }
    }


    CPO = 1./CPO;

    //Rcout << "CPO = " << endl << CPO.submat(0, 0, 9, 3) << endl;

    CPO.elem( find_nonfinite(CPO) ).zeros();

    MPL = accu(CPO);
    DIC = -4*DIC/(Num_of_iterations/2) + 2*logL;
    AIC = -2*logL + 2 * (Num_of_covariates+Num_of_obs*Num_of_RanEffs + accu(ARMA_Order) );
    BIC = -2*logL + log(Num_of_obs) * (Num_of_covariates+Num_of_obs*Num_of_RanEffs+ accu(ARMA_Order) );

}


SEXP ProbitMLModelSelectionARMA::MCMC_Procedure()
{
    Rcout << "============= FMR: MCMC Starts=============="<< endl;
    
    //X1 <- rtmvnorm(n=20000, mean = c(-1, -2, 1, 3), sigma=C, lower=c(-Inf, -Inf, 0, 0), upper=c(0, 0, Inf, Inf), algorithm="gibbs", burn.in.samples=100, thinning=5)
    //cor(X1)
    
    //Y_star_sample(span(0, tp-1), i) = rtmvnorm_gibbs(1, mu_tmp, Ri_tmp, lower, upper, Y_star_sample.col(i) ).t();

    //Rcout << cor(rtmvnorm_gibbs(2000, vec , Ri_tmp, lower, upper, Y_star_sample.col(i) ))
    
    
    List PosteriorSamples;
    List PosteriorEstimates;
    List Posterior;
    
    time_t start = time(NULL);
    
    int iter = 0, percent=0;
    
    while(iter < Num_of_iterations-1){
        //Rcout << "iter = " << iter << endl;
        Update_ystar_b_nu_beta_Sigma(iter);

        //if(updatenu)
        //    Update_nu(iter);
        //else
        //    nu_samples.col(iter+1) = nu_samples.col(iter);

        //if(updatebeta)
          //  Update_beta(iter);
        //else
            //beta_samples.col(iter+1) = beta_samples.col(iter);

//        if(SIMPLE){
//            if(updateomega)
//                Update_omega(iter);
//            else{
//                omega_samples.slice(iter+1) = omega_samples.slice(iter);
//                Sigma_samples.slice(iter+1) = Sigma_samples.slice(iter);
//            }
//        }
//        else{
        
//        }
        //if(updatedelta)
        //    Update_delta(iter);
        //else
        //    delta_samples.col(iter+1) = delta_samples.col(iter);
        
        if(SinglePhiPsi){
            if(updatephi)
                Update_single_phi(iter);
            else
                phi_samples.slice(iter+1) = phi_samples.slice(iter);
            
            if(updatepsi)
                Update_single_psi(iter);
            else
                psi_samples.slice(iter+1) = psi_samples.slice(iter);
        }
        else{
            if(updatephi)
                Update_phi(iter);
            else
                phi_samples.slice(iter+1) = phi_samples.slice(iter);
            
            if(updatepsi)
                Update_psi(iter);
            else
                psi_samples.slice(iter+1) = psi_samples.slice(iter);
        }
        
        percent = (100 *iter) / (Num_of_iterations-2) ;
        iter++;
        
        if(percent%5==0){
            Rcout << "\r" << "[" << std::string(percent / 5, (char)61) << std::string(100 / 5 - percent / 5, ' ') << "]" << "\t" << percent << "%";
            //Rcout << percent << "%" << " [Iteration " << iter + 1 << " of " << Num_of_iterations << "]";
            Rcout.flush();
        }

        //if((iter+1)%100==0)
        //    Rcout << iter+1 << endl;
        //Rcout<<"\t"<< round((iter+1.)/Num_of_iterations*100)<<" %"<<'\r';
    }
    
    ParameterEstimation();
    
    Rcout << endl<< "============= FMR: MCMC: Done =============="<< endl;
    time_t end = time(NULL);
    Rcout<<"Execution Time: "<< (double)(end-start)<<" Seconds"<<std::endl;
    
    PosteriorSamples["ystar.samples"] = Y_star_sample;
    PosteriorSamples["b.samples"] = b_samples;
    PosteriorSamples["nu.samples"] = nu_samples;
    PosteriorSamples["beta.samples"] = beta_samples;
    PosteriorSamples["Sigma.samples"] = Sigma_samples;
    //PosteriorSamples["delta.samples"] = delta_samples;
    PosteriorSamples["phi.samples"] = phi_samples;
    PosteriorSamples["psi.samples"] = psi_samples;
 
    PosteriorEstimates["beta.mean"] = beta_mean;
    PosteriorEstimates["nu.mean"] = nu_mean;
    //PosteriorEstimates["delta.mean"] = delta_mean;
    PosteriorEstimates["b.mean"] = b_mean;
    PosteriorEstimates["Sigma.mean"] = Sigma_mean;
    PosteriorEstimates["phi.mean"] = phi_mean;
    PosteriorEstimates["psi.mean"] = psi_mean;


    PosteriorEstimates["AIC"] = AIC;
    PosteriorEstimates["BIC"] = BIC;
    PosteriorEstimates["CIC"] = CIC;
    PosteriorEstimates["logL"] = logL;

    PosteriorEstimates["DIC"] = DIC;
    PosteriorEstimates["RJR"] = RJ_R;

    PosteriorEstimates["MPL"] = MPL;
    PosteriorEstimates["ACC"] = ACC/accu(TimePointsAvailable);

 
    Posterior["PosteriorEstimates"] = PosteriorEstimates;
    Posterior["PosteriorSamples"] = PosteriorSamples;
 
    Rcout << endl << "=======================================" << endl;
    Rcout << "acceptance rate for phi= " << acc_phi_rate/Num_of_iterations << endl;
    Rcout << "acceptance rate for psi= " << acc_psi_rate/Num_of_iterations << endl;
    Rcout << "=======================================" << endl;
    return (Posterior);
}

