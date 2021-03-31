//
//  ProbitML.cpp
//  
//
//  Created by kuojung on 2020/2/26.
//

// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "ProbitML.h"
//#include "tmvrnormGibbs.h"
#include "tmvrnormGibbs_KJLEE.h"
//RNGScope scope;
//#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]


ProbitMLModelSelection::ProbitMLModelSelection(int iNum_of_iterations, List list_Data, List list_InitialValues, List list_HyperPara, List list_UpdatePara, List list_TuningPara)
{
    //Rcout<< "Read Data 0" << endl;
    Num_of_iterations = iNum_of_iterations;
    Data = list_Data;
    InitialValues = list_InitialValues;
    HyperPara = list_HyperPara;
    UpdatePara = list_UpdatePara;
    TuningPara = list_TuningPara;
    
    
    updateystar = as<bool>(UpdatePara["UpdateYstar"]);
    updateb = as<bool>(UpdatePara["UpdateRandomEffect"]);
    updatenu = as<bool>(UpdatePara["UpdateNu"]);
    updatebeta = as<bool>(UpdatePara["UpdateBeta"]);
    updateSigma = as<bool>(UpdatePara["UpdateSigma"]);
    updatedelta = as<bool>(UpdatePara["UpdateDelta"]);
    
    //Rcout << "updatedelta = " << updatedelta << endl;
    //Unconstraint = as<bool>(UpdatePara["Unconstraint"]);
    //updateomega = as<bool>(UpdatePara["UpdateOmega"]);
    //SIMPLE = as<bool>(UpdatePara["SIMPLE"]);

    
    //Rcout<< "Read Data" << endl;
    Y = as<mat>(Data["Y"]);
    X = as<cube>(Data["X"]);
    Z = as<cube>(Data["Z"]);
    
    Num_of_obs = Y.n_cols;
    Num_of_Timepoints = Y.n_rows;
    Num_of_RanEffs = Z.n_cols;
    Num_of_covariates = X.n_cols;
    
    //Rcout << "Num_of_obs = " << Num_of_obs << endl;
    //Rcout << "Num_of_Timepoints = " << Num_of_Timepoints << endl;
    //Rcout << "Num_of_RanEffs = " << Num_of_RanEffs << endl;
    //Rcout << "Num_of_covariates = " << Num_of_covariates << endl;

    if(updatedelta){
        //Rcout << "updatedelta" << endl;
        U = as<cube>(Data["U"]);
        Num_of_deltas = U.n_slices/Num_of_obs;
        UU.set_size(Num_of_deltas);
        
        for(int delta_index = 0; delta_index<Num_of_deltas; delta_index++){
            //Rcout << "delta_index = " << delta_index << endl;
            UU(delta_index) = U.slices( (Num_of_obs*delta_index), (Num_of_obs*(delta_index+1)-1));
        }
        //Rcout << "Num_of_deltas = " << Num_of_deltas << endl;
        //Rcout << "updatedelta 2" << endl;
        delta_samples.set_size(Num_of_deltas, Num_of_iterations);
        //Rcout << "updatedelta 3" << endl;
        delta_mean.set_size(Num_of_deltas);
        //Rcout << "updatedelta 4" << endl;
        delta_samples.col(0) = as<vec>(InitialValues["delta"]);
        //Rcout << "updatedelta 5" << endl;
        sigma2_delta = as<double>(HyperPara["sigma2.delta"]);
        
        //Rcout<< "Read Tuning parameters." << endl;

        tuning_delta = as<double>(TuningPara["TuningDelta"]);
        Idelta_diag.eye(Num_of_deltas, Num_of_deltas);
        acc_rate_delta = 0;
        
        //Rcout << "Num_of_deltas = " << Num_of_deltas << endl;

    }
    
    //if(!Rf_isNull(Data["HSD.cov"]))
    //    HSD_Cov =as<mat>(Data["HSD.cov"]);
    //if(!Rf_isNull(Data["U"]))
     //   U = as<cube>(Data["U"]);
    
    //UU = as<field<cube>>(Data["U"]);
    //Rcout << "U = " << size(U) << endl;
    //Rcout << "HSD_Cov" << HSD_Cov << endl;
    //if(U.is_empty())
    //    U = cube();
   
    //if(HP.is_empty())
    //    HP_Cov = mat();

    //W = as<cube>(Data["W"]);
    //HP_Model = as<int>(Data["HP.Model"]);
    TimePointsAvailable = as<vec>(Data["TimePointsAvailable"]);
    
    b_samples.set_size(Num_of_RanEffs, Num_of_obs, Num_of_iterations);
    nu_samples.set_size(Num_of_obs, Num_of_iterations);
    beta_samples.set_size(Num_of_covariates,Num_of_iterations);
    Sigma_samples.set_size(Num_of_RanEffs, Num_of_RanEffs, Num_of_iterations);
    
    
    phi_samples.set_size(Num_of_obs, Num_of_iterations);
    psi_samples.set_size(Num_of_obs, Num_of_iterations);
    
    omega_samples.set_size(Num_of_Timepoints, Num_of_Timepoints, Num_of_iterations);

    
    b_mean.set_size(Num_of_RanEffs, Num_of_obs);
    nu_mean.set_size(Num_of_obs);
    beta_mean.set_size(Num_of_covariates);
    Sigma_mean.set_size(Num_of_RanEffs, Num_of_RanEffs);
    
    

    //omega_samples.set_size(Num_of_RanEffs, Num_of_RanEffs, Num_of_iterations);
    
    //Rcout<< "Initial Values" << endl;
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
    //Rcout<< "Initial Values delta" << endl;


    //Rcout<< "Read Hyperparameters." << endl;
    // Hyperparameters
    v_gamma = as<double>(HyperPara["v.gamma"]);
    sigma2_beta = as<double>(HyperPara["sigma2.beta"]);

    Vb = as<double>(HyperPara["InvWishart.df"]);
    Lambda = as<mat>(HyperPara["InvWishart.Lambda"]);

    //Rcout<< "Read Tuning parameters." << endl;
    
    
    Ib_diag.eye(Num_of_RanEffs, Num_of_RanEffs);
    
    //vec v = { -1, -2, 1, 3}, del = {-0.5, -0.3};
    //vec v_lower = {datum::inf, datum::inf, 0, 0}, v_upper= {0, 0, -datum::inf, -datum::inf};
    //Rcout << rtmvnorm_gibbs_KJLEE(1, v, Ri(1, 4, del), v_lower, v_upper, 0, zeros<vec>(4), 1) << endl;

}

mat ProbitMLModelSelection::Ri_Version2(int i, int tp, vec delta)
{
    //Rcout << "delta = " << delta << endl;
    mat F_tmp(tp, tp), F(tp, tp);
    F.zeros();
    F_tmp.zeros();
    //Rcout << "delta = " << delta << endl;
    //Rcout << "Ri" << endl;
    //Rcout << "W = " << endl << W << endl;
    
    
    /*
    for(int l=0; l<tp; l++)
        for(int m=0; m<tp; m++)
            for(int delta_ind_U = 0; delta_ind_U<Num_of_deltas; delta_ind_U++)
                F_tmp(l, m) = delta(delta_ind_U)*UU(delta_ind_U)(l, m, i);
    */
    
    
    for(int delta_ind_U = 0; delta_ind_U<Num_of_deltas; delta_ind_U++){
        if(0){
            Rcout << "UU(delta_ind_U).slice(i)(0, 0, size(tp, tp) )" << endl;
            Rcout << "delta i = " << delta_ind_U << "\t" << delta(delta_ind_U) << endl;
            Rcout << UU(delta_ind_U).slice(i)(0, 0, size(tp, tp) ) << endl;
        }
        F_tmp += delta(delta_ind_U)*UU(delta_ind_U).slice(i)(0, 0, size(tp, tp) );
    }
    
    //Rcout << "F_tmp = " << endl << F_tmp << endl;
    //F_tmp = F_tmp + 10.;
    //for(int delta_ind_HSD_Cov = 0; delta_ind_HSD_Cov<HSD_Cov.n_cols; delta_ind_HSD_Cov++)
        //F_tmp += delta(delta_ind_HSD_Cov+U.n_slices)*U.slice(delta_ind_HSD_Cov)(0, 0, size(tp, tp) );
    
    //if(i == 0)
        //Rcout << "F_tmp = " << F_tmp << endl;
        
    F_tmp = datum::pi*exp(F_tmp)/(1.+exp(F_tmp));
    
    //Rcout << "F_tmp = " << endl << F_tmp << endl;
    F(0, 0) = 1;
    
    //Rcout << "Dit 2" << endl;
    for(int t=1; t<tp; t++)
        F(t, 0) = cos(F_tmp(t, 0));
    for(int j = 1; j<tp-1; j++)
        for(int t = j+1; t<tp; t++)
            F(t, j) = cos(F_tmp(t, j))*prod(sin(F_tmp(t, span(0, j-1) )));
    //Rcout << "Dit 3" << endl;
    for(int t=1; t<tp; t++)
        F(t, t) = prod(sin(F_tmp(t, span(0, t-1) )));
    //Rcout << "Dit 4" << endl;
    mat Ri = F * F.t();
    //Rcout << "Dit 5" << endl;
    return (Ri);
}


mat ProbitMLModelSelection::Ri(int i, int tp, vec delta)
{
    mat F_tmp(tp, tp), F(tp, tp);
    F.zeros();
    
    //Rcout << "delta = " << delta << endl;
    //Rcout << "Ri" << endl;
    //Rcout << "W = " << endl << W << endl;
    for(int t=0; t<tp; t++)
        for(int j=0; j<tp; j++)
            if(HP_Model==1)
                F_tmp(t, j) = (abs(t-j)==1)*delta(0); // omega*delta
            else if(HP_Model==2)
                F_tmp(t, j) = (abs(t-j)==1)*delta(0)+X(1, 1, i)*(abs(t-j)==1)*delta(1); // I|t-j|, Sex
            else if(HP_Model==3)
                F_tmp(t, j) = (abs(t-j)==1)*delta(0)+X(t, 3, i)*(abs(t-j)==1)*delta(1)+X(t, 4, i)*(abs(t-j)==1)*delta(2); // I|t-j|, Drink past, present
            else if(HP_Model==4)
                F_tmp(t, j) = (abs(t-j)==1)*delta(0)+X(t, 5, i)*(abs(t-j)==1)*delta(1)+X(t, 6, i)*(abs(t-j)==1)*delta(2); // I|t-j|, Smoke past, present
            else if(HP_Model==5)
                F_tmp(t, j) = (abs(t-j)==1)*delta(0)+(abs(t-j)==2)*delta(1); // omega*delta
            else if(HP_Model==6)
                F_tmp(t, j) = delta(0)+abs(t-j)*delta(1); // omega*delta
            else if(HP_Model==7)
                F_tmp(t, j) = delta(0)+abs(t-j)*delta(1)+ pow((t-j), 2.)*delta(2); // omega*delta
            else if(HP_Model==8)
                F_tmp(t, j) = (abs(t-j)==1)*delta(0)+(abs(t-j)==2)*delta(1)+(abs(t-j)==3)*delta(2); // omega*delta
            else F_tmp(t, j) = 0; // independence
        
    F_tmp = datum::pi*exp(F_tmp)/(1.+exp(F_tmp));
    
    //Rcout << "F_tmp = " << endl << F_tmp << endl;
    F.diag() += 1.;
    
    //Rcout << "Dit 2" << endl;
    for(int t=1; t<tp; t++)
        F(t, 0) = cos(F_tmp(t, 0));
    for(int j = 1; j<tp-1; j++)
        for(int t = j+1; t<tp; t++)
            F(t, j) = cos(F_tmp(t, j))*prod(sin(F_tmp(t, span(0, j-1) )));
    //Rcout << "Dit 3" << endl;
    for(int t=1; t<tp; t++)
        F(t, t) = prod(sin(F_tmp(t, span(0, t-1) )));
    //Rcout << "Dit 4" << endl;
    mat Ri = F * F.t();
    //Rcout << "Dit 5" << endl;
    return (Ri);
}

mat ProbitMLModelSelection::Ri_Unconstraint(int i, int tp, mat omega_mat)
{
    mat F(tp, tp);
    F.zeros();
    
    //Rcout << "F_tmp = " << endl << F_tmp << endl;
    F.diag() += 1.;
    
    //Rcout << "Dit 2" << endl;
    for(int t=1; t<tp; t++)
        F(t, 0) = cos(omega_mat(t, 0));
    for(int j = 1; j<tp-1; j++)
        for(int t = j+1; t<tp; t++)
            F(t, j) = cos(omega_mat(t, j))*prod(sin(omega_mat(t, span(0, j-1) )));
    //Rcout << "Dit 3" << endl;
    for(int t=1; t<tp; t++)
        F(t, t) = prod(sin(omega_mat(t, span(0, t-1) )));
    //Rcout << "Dit 4" << endl;
    mat Ri = F * F.t();
    //Rcout << "Dit 5" << endl;
    return (Ri);
}


mat ProbitMLModelSelection::Sigma_b(int iter, mat omega)
{
    mat F(Num_of_RanEffs, Num_of_RanEffs), F_tmp(Num_of_RanEffs, Num_of_RanEffs);
    mat Ri_tmp(Num_of_RanEffs, Num_of_RanEffs);
    if(Num_of_RanEffs == 1)
        Ri_tmp(1, 1) = 1.;
    else{
        F.zeros();
        F.diag() += 1.;
        for(int j = 1; j<Num_of_RanEffs-1; j++)
            for(int t = j+1; t<Num_of_RanEffs; t++)
                F_tmp(t, j) = Rf_runif(0., datum::pi);
        
        for(int t=1; t<Num_of_RanEffs; t++)
            F(t, 0) = cos(F_tmp(t, 0));
        for(int j = 1; j<Num_of_RanEffs-1; j++)
            for(int t = j+1; t<Num_of_RanEffs; t++)
                F(t, j) = cos(F_tmp(t, j))*prod(sin(F_tmp(t, span(0, j-1) )));
        //Rcout << "Dit 3" << endl;
        for(int t=1; t<Num_of_RanEffs; t++)
            F(t, t) = prod(sin(F_tmp(t, span(0, t-1) )));
        Ri_tmp = F * F.t();
    }
    
    //Rcout << "Dit 5" << endl;
    return (Ri_tmp);
}


void ProbitMLModelSelection::Update_ystar(int iter)
{
    if(iter % 100 == 0)
        Rcout << "Update ystar" << endl;
    vec mu_tmp, res, res_tmp;
    //uvec q1;
    mat Ri_tmp;
    int tp;
    //vec y_star_tmp;
    vec lower, upper;
    mat X_tmp, Z_tmp;
    //vec  init_state;
    //double sigma_tmp, mu_tmp_1;
    for(int i=0; i<Num_of_obs; i++){
        //Rcout << "i = " << i << endl;
        
        tp = TimePointsAvailable(i);
        lower.zeros(tp);
        upper.zeros(tp);
        //init_state.zeros(tp);
        X_tmp = X(span(0, tp-1), span(0, Num_of_covariates-1), span(i));
        Z_tmp = (Z.slice(i).rows(0, tp-1));
        mu_tmp = X_tmp*beta_samples.col(iter) + Z_tmp*b_samples.slice(iter).col(i);
        Rcout << "mu_tmp =" << mu_tmp << endl;
        //if(Unconstraint)
            //Ri_tmp = Ri(i, tp, omega_samples(iter));
        //else
        if(updatedelta)
            Ri_tmp = Ri_Version2(i, tp, delta_samples.col(iter));
        else
            Ri_tmp = eye(tp, tp);
        //if(i<2)
            //Rcout << "Ri_tmp = " << endl << Ri_tmp << endl;
        
        //Y_star_sample.col(i) = rmvnorm(1, mu_tmp, Ri_tmp).t();
        //Y_star_sample.col(i).elem(q1) = 0.;
        //y_star_tmp = Y_star_sample.col(i);
        
        //q1 = find(y_star_tmp%Y.col(i)<0);
        //y_star_tmp.elem(q1).zeros();
        //X( row_number, span(first_col, last_col) )
        lower.elem(find(Y(span(0, tp-1), i)>0)).zeros();
        lower.elem(find(Y(span(0, tp-1), i)==0)).ones();
        lower.elem(find(Y(span(0, tp-1), i)==0)) *= -datum::inf;
        
        upper.elem(find(Y(span(0, tp-1), i)==0)).zeros();
        upper.elem(find(Y(span(0, tp-1), i)>0)).ones();
        upper.elem(find(Y(span(0, tp-1), i)>0)) *= datum::inf;
        
        
        //Rcout << "Test 1 = " << r_truncnorm(1, 1, -datum::inf, 0) << endl;
        //Rcout << "Test 2 = " << r_truncnorm(1, 1, 0, datum::inf) << endl;
        
        //for(int t=0; t<tp; t++){
        //    Y_star_sample(t, i) = r_truncnorm(mu_tmp(t), 1, lower(t), upper(t));
        //}
        
        //Rcout << "Test 1 = " << endl;
        // double r_truncnorm(const double mu, const double sigma, const double a, const double b)
        
        //q1 = find(y_star_tmp%(Y.col(i)-1)<0);
        //y_star_tmp.elem(q1).zeros();
        
        //Y_star_sample.col(i) = y_star_tmp;
        //Y_star_sample.col(i) = ( find(Y_star_sample.col(i)%Y.col(i))<0 ).zeros();
        //Rcout << join_rows(join_rows(Y.col(i), join_rows(lower, upper)), mu_tmp) << endl;
        //Rcout << mu_tmp << endl;
        //Rcout << Ri_tmp << endl;
        //Rcout << rtmvnorm_gibbs(1, mu_tmp, Ri_tmp, lower, upper, init_state) << endl;
        //rtmvnorm_gibbs(int n, arma::vec mu, arma::mat sigma, arma::vec lower, arma::vec upper, arma::vec init_state)
        //Y_star_sample(span(0, tp-1), i) = rtmvnorm_gibbs(1, mu_tmp, Ri_tmp, lower, upper, Y_star_sample.col(i) ).t();
        
        Y_star_sample(span(0, tp-1), i) = rtmvnorm_gibbs_KJLEE(1, mu_tmp, Ri_tmp, lower, upper, 100, zeros<vec>(tp), 5).t();
        
        //vec v = { -1, -2, 1, 3};
        //vec v_lower = {-datum::inf, -datum::inf, 0, 0}, v_upper= {0, 0, datum::inf, datum::inf};
        //Rcout << "cor = " << endl << cor(rtmvnorm_gibbs(20000, v, Ri_tmp, lower, upper, zeros<vec>(4))) <<endl;
        
        //Rcout << "mean = " << endl << mean(rtmvnorm_gibbs(20000, v, Ri_tmp, lower, upper, zeros<vec>(4)), 0) <<endl;
        
        //mat rtmvnorm_gibbs_KJLEE(int n, vec mean, mat Sigma, vec lower, vec upper, int burn_in, vec start_value, int thinning);

        //Rcout << "cor = " << endl << cor(rtmvnorm_gibbs_KJLEE(20000, v, Ri_tmp, v_lower, v_upper, 100, zeros<vec>(4), 5)) <<endl;
        
        //Rcout << "mean = " << endl << mean(rtmvnorm_gibbs_KJLEE(20000, v, Ri_tmp, v_lower, v_upper, 100, zeros<vec>(4), 5)) <<endl;

        //Rcout << "norm = " << endl << rtmvnorm_gibbs(5, mu_tmp, Ri_tmp, lower, upper, mu_tmp) << endl;
        //Y_star_sample.col(i) = rmvnorm(1, mu_tmp, Ri_tmp).t();
        
        //rmvnorm(const arma::uword n, const arma::vec& mu, const arma::mat& S)
        //Rcout << join_rows(join_rows(join_rows(Y.col(i), join_rows(lower, upper)), mu_tmp), Y_star_sample.col(i)) << endl;
    }
    
    //Rcout << "Cor = " << endl << cor(Y_star_sample.t()) << endl;
        //Rcout << "res = " << res << endl;
        /*
        for(int t=0; t<tp; t++){
            //res = Y_star_sample.col(i) - mu_tmp;
            res_tmp = Y_star_sample.col(i) - mu_tmp;
            //Rcout << "t = " << t << endl;
            res_tmp.shed_row(t);
            Sigma22 = Ri_tmp;
            Sigma12 = Ri_tmp.row(t);
            Sigma12.shed_col(t);
            Sigma22.shed_row(t);
            Sigma22.shed_col(t);
            Sigma22_inv = Sigma22.i();
            
            mu_tmp_1 = mu_tmp(t) + as_scalar(Sigma12*Sigma22_inv*res_tmp);
            sigma_tmp = Ri_tmp(t, t) - as_scalar(Sigma12*Sigma22_inv*Sigma12.t());
            
            //Rcout <<"mu_tmp_1=" << mu_tmp_1 << "\t sigma_tmp=" << sigma_tmp << endl;
            //Rcout << "check 1" << endl;
            if(Y(t, i)>0){
                if(mu_tmp_1<0)
                    mu_tmp_1 = 0.;
                //Rcout << "mu_tmp_1 = " << mu_tmp_1 << endl;
                //Rcout << "Test = " << r_truncnorm(0, 1, 0, datum::inf) <<endl;
                Y_star_sample(t, i) = r_truncnorm(mu_tmp_1, sqrt(sigma_tmp), 0, datum::inf);
                //Rcout << "Y_star_sample(t, i) = " << Y_star_sample(t, i) << endl;
            }
            else{
                if(mu_tmp_1>0)
                    mu_tmp_1 = 0.;
                Y_star_sample(t, i) = r_truncnorm(mu_tmp_1, sqrt(sigma_tmp), -datum::inf, 0);
            }
        }
        
    }
         */
}

void ProbitMLModelSelection::Update_b(int iter)
{
    if(iter % 100 == 0)
        Rcout << "Update b" << endl;
    int tp;
    vec mu_tmp;
    mat Sigma_tmp(Num_of_RanEffs, Num_of_RanEffs), Ri_inv;
    Sigma_tmp.zeros();
    vec yi_star, res;
    mat X_tmp, Z_tmp;
    
    
    
    for(int i=0; i<Num_of_obs; i++){
        //Rcout << "i=" << i << endl;
        tp = TimePointsAvailable(i);
        Ri_inv = inv_sympd( Ri_Version2(i, tp, delta_samples.col(iter)) );
        
        //Rcout << "Ri_inv = " << endl << Ri_inv << endl;
        X_tmp =  X(span(0, tp-1), span(0, Num_of_covariates-1), span(i));
        Z_tmp = (Z.slice(i).rows(0, tp-1));
        Sigma_tmp = inv_sympd(Z_tmp.t()*Ri_inv*Z_tmp+nu_samples(i, iter)*Sigma_samples.slice(iter).i());
        //Rcout << "Sigma_tmp = " << endl << Sigma_tmp << endl;
        //yi_star =  y_star_samples.slice(iter).col(i);
        //Rcout << "X.slice(i)*beta_samples.col(iter) = " << endl << X.slice(i)*beta_samples.col(iter) << endl;
        //Rcout << "Y_star_sample.col(i) =" << endl << Y_star_sample.col(i) << endl;
        res = Y_star_sample(span(0, tp-1), i) - X_tmp*beta_samples.col(iter);
        
        //Rcout << "res = " << endl << res << endl;
        mu_tmp = Sigma_tmp*Z_tmp.t()*Ri_inv*res;
        
        //Rcout << "mu_tmp = " << endl << mu_tmp << endl;
        //Rcout << rmvnorm(1, mu_tmp, Sigma_tmp) << endl;
        b_samples.slice(iter+1).col(i) = mvnrnd(mu_tmp, Sigma_tmp);
    }
    
}




void ProbitMLModelSelection::Update_ystar_b_beta_Sigma(int iter)
{
    //if(iter == 1)
        //Rcout << "Update ystar, b, beta, Sigmab simultaneously" << endl;
    vec mu_tmp_b, res_b, mu_tmp_ystar;
    vec res_beta, mu_tmp_beta = zeros<vec>(Num_of_covariates);

    mat Ri_tmp, Ri_inv;

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
    //double sigma_tmp, mu_tmp_1;
    
    //do {
        for(int i=0; i<Num_of_obs; i++){
            tp = TimePointsAvailable(i);

            X_tmp = X(span(0, tp-1), span(0, Num_of_covariates-1), span(i));
            Z_tmp = (Z.slice(i).rows(0, tp-1));
            mu_tmp_ystar = X_tmp*beta_samples.col(iter) + Z_tmp*b_samples.slice(iter).col(i);
            
            if(updatedelta)
                Ri_tmp = Ri_Version2(i, tp, delta_samples.col(iter));
            else
                Ri_tmp = eye(tp, tp);

            if(!Ri_tmp.is_sympd())
                Ri_tmp.eye();
            Ri_inv = inv_sympd(Ri_tmp);
            
            lower.set_size(tp);
            upper.set_size(tp);
            lower.elem(find(Y(span(0, tp-1), i)>0)).zeros();
            lower.elem(find(Y(span(0, tp-1), i)==0)).ones();
            lower.elem(find(Y(span(0, tp-1), i)==0)) *= -datum::inf;
            
            upper.elem(find(Y(span(0, tp-1), i)==0)).zeros();
            upper.elem(find(Y(span(0, tp-1), i)>0)).ones();
            upper.elem(find(Y(span(0, tp-1), i)>0)) *= datum::inf;


            //mu_tmp_ystar.elem( find( ( (Y(span(0, tp-1), i)-1) % mu_tmp_ystar) < 0 )).zeros();
            //mu_tmp_ystar.elem( find( ( (Y(span(0, tp-1), i) ) % mu_tmp_ystar) < 0 )).zeros();
            //if(mu_tmp_ystar.has_nan()){
            //    Rcout << "mu_tmp_ystar iter = " << iter << "\ti= " << i << endl;
            //}
            /*
            if(iter==1 && (i==1)){
                Rcout << "============== " << i << "============= " << endl;
                Rcout << "X_tmp*beta_samples.col(iter) = " << X_tmp*beta_samples.col(iter) << endl;
                Rcout << "Z_tmp*b_samples.slice(iter).col(i) = " << Z_tmp*b_samples.slice(iter).col(i) << endl;
                Rcout << "lower = " << lower << endl;
                Rcout << "upper = " << upper << endl;
                Rcout << "mu_tmp_ystar = " << mu_tmp_ystar << endl;
                Rcout << "Ri_tmp = " << Ri_tmp << endl;
            }
            */
            
            if(updateystar){
            
                if(tp == 1)
                    Y_star_sample(0, i) = rtruncnorm(1, as_scalar(mu_tmp_ystar), as_scalar(Ri_tmp), as_scalar(lower),  as_scalar(upper))(0);
                else
                    Y_star_sample(span(0, tp-1), i) = rtmvnorm_gibbs_KJLEE(1, mu_tmp_ystar, Ri_tmp, lower, upper, 100, zeros<vec>(tp), 5).t();
                
                if(Y_star_sample.col(i).has_nan())
                    Y_star_sample.col(i).zeros();
                Y_star_sample.col(i) = clamp(Y_star_sample.col(i), -10, 10);
            
            

                if(Y_star_sample(span(0, tp-1), i).has_nan()){
                    Rcout << endl << "iter = " << iter << "\ti= " << i << "\tY_star_sample = " << Y_star_sample(span(0, tp-1), i) <<endl;
                    Rcout << "mu_tmp_ystar = " << mu_tmp_ystar << "\tlower = " << lower << "\t upper = " << upper << endl;
                    Rcout << "Y(span(0, tp-1), i) = " << Y(span(0, tp-1), i) <<  "\tmu_tmp_ystar" <<  mu_tmp_ystar << "\t both = " <<Y(span(0, tp-1), i) % mu_tmp_ystar << endl;
                    //std::exit(1);
                    Rcpp::stop("NaN!\n"); 
                    Rcout << "============== " << i << "============= " << endl;
                    Rcout << "X_tmp*beta_samples.col(iter) = " << X_tmp*beta_samples.col(iter) << endl;
                    Rcout << "Z_tmp*b_samples.slice(iter).col(i) = " << Z_tmp*b_samples.slice(iter).col(i) << endl;
                    Rcout << "lower = " << lower << endl;
                    Rcout << "upper = " << upper << endl;
                    Rcout << "mu_tmp_ystar = " << mu_tmp_ystar << endl;
                    Rcout << "Ri_tmp = " << Ri_tmp << endl;
                    Rcout << "Y_star_sample = " << Y_star_sample.col(i) << endl;
                }

            }
        
            /*
            if(iter==121 && (i==54|| i ==55)){
                Rcout << "============== Before " << i << "============= " << endl;
                Rcout << "Sigma_tmp_b = " << Sigma_tmp_b << endl;
                Rcout << "res_b = " << res_b << endl;
                Rcout << "mu_tmp_b = " << mu_tmp_b << endl;
                Rcout << "b_samples.slice(iter+1).col(i) = " << b_samples.slice(iter+1).col(i) << endl;
                
                if(i==54)
                    b_samples.slice(iter+1).col(i) = 0;
                Rcout << "b_samples.slice(iter+1).col(i) = " << b_samples.slice(iter+1).col(i) << endl;
            }
            */
            
            
            
            Sigma_tmp_b = Z_tmp.t()*Ri_inv*Z_tmp+nu_samples(i, iter)*Sigma_samples.slice(iter).i();
            if(!Sigma_tmp_b.is_sympd())
                Sigma_tmp_b.eye();
            else
                Sigma_tmp_b = inv_sympd(Sigma_tmp_b);
            res_b = Y_star_sample(span(0, tp-1), i) - X_tmp*beta_samples.col(iter);
            mu_tmp_b = Sigma_tmp_b*Z_tmp.t()*Ri_inv*res_b;
            
            if(Sigma_tmp_b.has_nan()){
                Sigma_tmp_b.eye();
                mu_tmp_b.zeros();
            }
            
            if(updateb)
                b_samples.slice(iter+1).col(i) = mvnrnd(mu_tmp_b, Sigma_tmp_b);
            else
                b_samples.slice(iter+1).col(i) = b_samples.slice(iter).col(i);

            /*
            if(iter==121 && (i==54|| i ==55)){
                Rcout << "============== " << i << "============= " << endl;
                
                
                Rcout << "mu_tmp_ystar = " << mu_tmp_ystar << endl;
                Rcout << "X_tmp*beta_samples.col(iter) = " << X_tmp*beta_samples.col(iter) << endl;
                Rcout << "Z_tmp*b_samples.slice(iter).col(i) = " << Z_tmp*b_samples.slice(iter).col(i) << endl;
                Rcout << "Ri_tmp = " << Ri_tmp << endl;
                Rcout << "Y_star_sample(span(0, tp-1), i)  = " << Y_star_sample(span(0, tp-1), i)  << endl;

                Rcout << "================After======================= " << endl;
                Rcout << "Sigma_tmp_b = " << Sigma_tmp_b << endl;
                Rcout << "res_b = " << res_b << endl;
                Rcout << "mu_tmp_b = " << mu_tmp_b << endl;
                Rcout << "b_samples.slice(iter+1).col(i) = " << b_samples.slice(iter+1).col(i) << endl;
            }
            */
            
            Sigma_tmp_beta += (X_tmp.t()*Ri_inv*X_tmp);
            
             //res = Y_star_sample((span(0, tp-1), i))- (Z.slice(i))*b_samples.slice(iter).col(i);
            
            res_beta = Y_star_sample(span(0, tp-1), i)- Z_tmp*b_samples.slice(iter+1).col(i);
            
            mu_tmp_beta += X_tmp.t()*Ri_inv*res_beta;
            
            b_vec = b_samples.slice(iter+1).col(i);
            
            //if(iter==121 && i==55){
            //    Rcout << i << "\t" << "b_vec = "<< b_vec << endl;
            //    Rcout << "nu = " << nu_samples(i, iter) << endl;
            //}
            //if(i==27){
            //    Rcout << i << "\t" << "nu = " <<  nu_samples(i, iter)<< endl;
            //    Rcout << i << "\t" << "b_vec" << b_vec<< endl;
            //}
            Sigma_tmp += nu_samples(i, iter)*(b_vec*b_vec.t());
            

            
          }
            
        //if(iter==121)
        //Rcout << "iter = " << iter << "\tBefore Sigma_tmp = " << Sigma_tmp << "\t" << Sigma_tmp.is_sympd() << endl;
        //if(iter<=1){
        //    Rcout << "mu_tmp_beta Before = " << mu_tmp_beta << endl;
        //}
        Sigma_tmp_beta.diag() += 1./sigma2_beta;
        Sigma_tmp_beta = inv_sympd(Sigma_tmp_beta);
        mu_tmp_beta = Sigma_tmp_beta * mu_tmp_beta;
        
        /*
        if(iter<=1){
            Rcout << "iter = Before" << iter << endl;
            Rcout << "Sigma_tmp_beta = " << Sigma_tmp_beta << endl;
            Rcout << "mu_tmp_beta = " << mu_tmp_beta << endl;
            Rcout << "beta_samples.col(iter+1) = " << beta_samples.col(iter+1) << endl;
        }
         */

        if(updatebeta)
            beta_samples.col(iter+1) = mvnrnd(mu_tmp_beta, Sigma_tmp_beta);
        else
            beta_samples.col(iter+1) = beta_samples.col(iter);
        
    
        //if(iter<=1){
        //    Rcout << "beta_samples.col(iter+1) = " << beta_samples.col(iter+1) << endl;
        //}


            //Rcout << "Sigma_tmp" << Sigma_tmp << endl;
        
         
        
        Sigma_tmp = (Sigma_tmp + Lambda);
        
        //Rcout << "Sigma_tmp = " << Sigma_tmp << "\t" << Sigma_tmp.is_sympd() << endl;

    //} while (!Sigma_tmp.is_sympd());
    

    //if(iter==121)
        //Rcout << "Sigma_tmp_inv = " << Sigma_tmp << "\t Num_of_obs + Vb =" << Num_of_obs + Vb << endl;
    if(!Sigma_tmp.is_symmetric() || Sigma_tmp.has_nan()){
        Sigma_tmp.eye();
        //Rcout << "Sigma_tmp = " << Sigma_tmp << endl;
    }
    //Rcout << "iter = " << iter << endl;
    //Rcout << "Sigma_tmp = " << Sigma_tmp << endl;
    //Rcout << "Sigma_tmp.is_symmetric() = " << Sigma_tmp.is_symmetric() << endl;
    if(updateSigma)
        Sigma_samples.slice(iter+1) = iwishrnd( Sigma_tmp, (Num_of_obs + Vb)); //1./randg( 1, distr_param(gamma_shape, gamma_scale) );// //riwish(df,
    else
        Sigma_samples.slice(iter+1) = Sigma_samples.slice(iter);
}

void ProbitMLModelSelection::Update_nu(int iter)
{
    //if(iter == 1)
     //   Rcout << "Update nu" << endl;
    double alpha_tmp, beta_tmp;
    vec b_vec;
    alpha_tmp = 0.5*(v_gamma + Num_of_RanEffs);

    for(int i=0; i<Num_of_obs; i++){
        b_vec = b_samples.slice(iter+1).col(i); //( span::all, span(i), span(iter+1));
        beta_tmp = 0.5*(as_scalar(b_vec.t()*Sigma_samples.slice(iter+1).i()*b_vec) + v_gamma);
        nu_samples(i, iter+1) = randg( 1, distr_param(alpha_tmp, 1./beta_tmp))(0);  //Rf_rgamma(alpha_tmp, 1./beta_tmp); //
    }
    
}


void ProbitMLModelSelection::Update_beta(int iter)
{
    if(iter % 100 == 0)
        Rcout << "Update beta" << endl;
    mat Sigma_tmp(Num_of_covariates, Num_of_covariates), R_inv;
    Sigma_tmp.zeros();
    mat X_tmp, Z_tmp;
    vec mu_tmp(Num_of_covariates);
    mu_tmp.zeros();
    vec res;
    int tp;
    for(int i=0; i<Num_of_obs; i++){
        //Rcout << "i = " << i << endl;
        tp = TimePointsAvailable(i);
        //Rcout << "Ri = " << endl <<  Ri(i, tp, delta_samples.col(iter)) << endl;
        R_inv = inv_sympd( Ri_Version2(i, tp, delta_samples.col(iter)) );
        
        
        X_tmp = X(span(0, tp-1), span(0, Num_of_covariates-1), span(i));
        Z_tmp = (Z.slice(i).rows(0, tp-1));
        
        Sigma_tmp += (X_tmp.t()*R_inv*X_tmp);
        
        //Rcout << "Sigma_tmp =" << Sigma_tmp << endl;
        //yi_star =  y_star_samples.slice(iter).col(i);
        //res = Y_star_sample((span(0, tp-1), i))- (Z.slice(i))*b_samples.slice(iter).col(i);
        
        res = Y_star_sample(span(0, tp-1), i)- Z_tmp*b_samples.slice(iter+1).col(i);
        
        mu_tmp += X_tmp.t()*R_inv*res;
    }
    
    Sigma_tmp.diag() += 1./sigma2_beta;
    Sigma_tmp = Sigma_tmp.i();
    mu_tmp = Sigma_tmp * mu_tmp;
    
    beta_samples.col(iter+1) = mvnrnd(mu_tmp, Sigma_tmp);
}



void ProbitMLModelSelection::Update_Sigma(int iter)
{
    if(iter % 100 == 0)
        Rcout << "Update Sigma" << endl;
    //Rcout << "Update Sigma" << endl;
    double df = Num_of_obs + Vb;
    mat Sigma_tmp(Num_of_RanEffs, Num_of_RanEffs);
    Sigma_tmp.zeros();
    vec b_vec;
    //double gamma_shape, gamma_scale;
    for(int i=0; i<Num_of_obs; i++){
        b_vec = b_samples.slice(iter+1).col(i);
        //Rcout << i << "\t" << "b_vec = "<< b_vec << endl;
        Sigma_tmp += nu_samples(i, iter)*(b_vec*b_vec.t());
    }
    //Rcout << "Sigma_tmp = " << endl << Sigma_tmp << endl;
    Sigma_tmp = (Sigma_tmp + Lambda);//Sigma_tmp + Lambda; //
    Rcout << "Sigma + Lambda = " << endl << Sigma_tmp << endl;
    //gamma_shape = df;
    //gamma_scale = 2./as_scalar(Sigma_tmp);
    Sigma_samples.slice(iter+1) =iwishrnd( Sigma_tmp, df); //1./randg( 1, distr_param(gamma_shape, gamma_scale) );// //riwish(df, Sigma_tmp);
    //Rcout << "End update Sigma" << endl;
}

void ProbitMLModelSelection::Update_delta(int iter)
{
    //if(iter == 0)
        //Rcout << "Update delta" << endl;
    double delta_den = 0., delta_num = 0.;
    vec delta_cand =  mvnrnd(delta_samples.col(iter), tuning_delta*Idelta_diag);
    
    //vectorise(rmvnorm(1, delta_samples.col(iter), tuning_delta*Idelta_diag));
    //Rcout << "delta_cand = " << delta_cand << endl;
    vec res;
    mat Ri_inv, X_tmp, Z_tmp;
    int tp;
    for(int i=0; i<Num_of_obs; i++){
        //Rcout << "i = " << i << endl;
        tp = TimePointsAvailable(i);
        X_tmp = X.slice(i).rows(0, tp-1);  //X(span(0, tp-1), span(0, Num_of_covariates-1), span(i));
        Z_tmp = Z.slice(i).rows(0, tp-1);

        res = Y_star_sample(span(0, tp-1), i) - X_tmp*beta_samples.col(iter+1)-Z_tmp*b_samples.slice(iter+1).col(i);

        Ri_inv = Ri_Version2(i, tp, delta_samples.col(iter));
        
        /*
        if(i==0 && iter == 10){
            Rcout << delta_samples.col(iter) << endl;
            Rcout <<Ri_Version2(i, tp, delta_samples.col(iter));
        }
        */
        
        if(!Ri_inv.is_sympd())
            Ri_inv.eye();
        else
            Ri_inv = inv_sympd(Ri_inv);
        //Ri_inv = inv_sympd( Ri(i, tp, delta_samples.col(iter)) );
        
        //Rcout << "Ri_inv_old = " << endl << Ri_inv << endl;
        
        delta_den += 0.5*log(det(Ri_inv)) - 0.5*as_scalar(res.t()* Ri_inv*res);
        
        //Ri_inv = inv_sympd( Ri(i, tp, delta_cand) );
        Ri_inv = Ri_Version2(i, tp, delta_cand);
        if(!Ri_inv.is_sympd())
            Ri_inv.eye();
        else
            Ri_inv = inv_sympd(Ri_inv);

        //Rcout << "Ri_inv_new = " << endl << Ri_inv << endl;
        
        delta_num += 0.5*log(det(Ri_inv)) - 0.5*as_scalar(res.t()* Ri_inv*res);
    }
    
    //Rcout << "delta_den = " << delta_den << "\t delta_num = " << delta_num << endl;
    
    delta_den = delta_den - 0.5*accu(square(delta_samples.col(iter)))/sigma2_delta;
    delta_num = delta_num - 0.5*accu(square(delta_cand))/sigma2_delta;
    
    //Rcout << "delta_den = " << delta_den << "\t delta_num = " << delta_num << "\t delta_num - delta_den =" << delta_num-delta_den << endl;
    if(log(Rf_runif(0., 1.)) < delta_num - delta_den ){
        delta_samples.col(iter+1) = delta_cand;
        //Rcout << "Accept" << endl;
        acc_rate_delta++;
    }
    else
        delta_samples.col(iter+1) = delta_samples.col(iter);
    
    if((iter+1)%500 == 0){
        //Rcout << "tuning_delta = " << tuning_delta << endl;
        //Rcout << "acc_rate_delta/iter = " << acc_rate_delta/iter << endl;
        if( acc_rate_delta/iter<0.25 )
            tuning_delta = tuning_delta/2.;
        if( (1.*acc_rate_delta)/iter>0.50 )
            tuning_delta = 2*tuning_delta;
        
        //Rcout << "tuning_delta = " << tuning_delta << endl;
        
    }
}


/*
void ProbitMLModelSelection::Update_delta_CompWise(int iter)
{
    if(iter % 100 == 0)
        Rcout << "Update delta" << endl;
    double delta_den = 0., delta_num = 0.;
    vec delta_cand = delta_samples.col(iter); //vectorise(rmvnorm(1, delta_samples.col(iter), tuning_delta*Idelta_diag));
    double delta_cand_double;
    //Rcout << "delta_cand = " << delta_cand << endl;
    vec res;
    mat Ri_inv, X_tmp, Z_tmp;
    int tp;
    for(int j=0; j<Num_of_deltas; j++){
        delta_cand(j) = delta_samples(j, iter) + randn()* sqrt(tuning_delta);
        delta_den = 0., delta_num = 0.;
        for(int i=0; i<Num_of_obs; i++){
            
            //Rcout << "i = " << i << endl;
            tp = TimePointsAvailable(i);
            X_tmp =  X(span(0, tp-1), span(0, Num_of_covariates-1), span(i));
            Z_tmp = (Z.slice(i).rows(0, tp-1));
            res = Y_star_sample(span(0, tp-1), i) - X_tmp*beta_samples.col(iter+1)-Z_tmp*b_samples.slice(iter+1).col(i);
        
            
            Ri_inv = inv_sympd( Ri(i, tp, delta_samples.col(iter)) );
        
            delta_den += 0.5*log(det(Ri_inv)) - 0.5*as_scalar(res.t()* Ri_inv*res);
            
            delta_den = delta_den - 0.5*square(delta_samples(j, iter))/sigma2_delta;
        
            delta_cand = delta_cand(j) + randn()* sqrt(tuning_delta);
        
            Ri_inv = inv_sympd( Ri(i, tp, delta_cand) );
        
        //Rcout << "Ri_inv_new = " << endl << Ri_inv << endl;
        
            delta_num += 0.5*log(det(Ri_inv)) - 0.5*as_scalar(res.t()* Ri_inv*res);
        
        
        //Rcout << "delta_den = " << delta_den << "\t delta_num = " << delta_num << endl;
        
            delta_den = delta_den - 0.5*accu(square(delta_samples.col(iter)))/sigma2_delta;
            delta_num = delta_num - 0.5*accu(square(delta_cand))/sigma2_delta;
        
            //Rcout << "delta_den = " << delta_den << "\t delta_num = " << delta_num << "\t delta_num - delta_den =" << delta_num-delta_den << endl;
            if(log(Rf_runif(0., 1.)) < delta_num - delta_den ){
                delta_samples.col(iter+1) = delta_cand;
                //Rcout << "Accept" << endl;
                acc_rate_delta++;
            }
            else
                delta_samples.col(iter+1) = delta_samples.col(iter);
        }
    
}
*/

/*
void ProbitMLModelSelection::Update_omega(int iter)
{
    if(iter == 100)
        Rcout << "Update omega" << endl;
    //Rcout << "Update Sigma" << endl;
    
    mat omega_old = omega_samples.slice(iter), omega_new;
    uvec omega_lower_indices = trimatl_ind( size(omega_old), -1);
    vec omega_old_vec = omega_old(omega_lower_indices)
    vec omega_new;
    omega_new_vec = omega_old + omega_new.randu(size(omega_old))*datum::pi;
    omega_new_vec = clamp(omega_new, 0, datum::pi);
    
    omega_new.zeros();
    omega_new(omega_lower_indices) = omega_old_vec;
    
    mat Ri_old, Ri_new;
    for(int i=0; i<Num_of_obs; i++){
        Ri_old = Ri_Unconstraint(i, TimePointsAvailable(i), omega_old);
        Ri_new = Ri_Unconstraint(i, TimePointsAvailable(i), omega_new);
        
    }
    
    mat Sigma_new_inv = Sigma_b(iter, omega_new).i();
    mat Sigma_old_inv = Sigma_b(iter, omega_old).i();
    //Sigma_tmp.zeros();
    vec b_vec;
    double exp_term_new = 0., exp_term_old =0.;
    for(int i=0; i<Num_of_obs; i++){
        b_vec = b_samples.slice(iter+1).col(i);
        exp_term_new += as_scalar(b_vec*Sigma_new_inv*b_vec);
        exp_term_old += as_scalar(b_vec*Sigma_old_inv*b_vec);
    }
    
    exp_term_old += 0.5*log(det(Sigma_old_inv));
    exp_term_new += 0.5*log(det(Sigma_new_inv));
    
    if(log(Rf_runif(0., 1.)) < exp_term_new-exp_term_old){
        omega_samples.slice(iter + 1) = omega_new;
        Sigma_samples.slice(iter+1) = Sigma_new_inv.i();
    }
    else{
        omega_samples.slice(iter + 1) = omega_old;
        Sigma_samples.slice(iter+1) = Sigma_old_inv.i();
    }
}

*/

void ProbitMLModelSelection::ParameterEstimation()
{
    //Rcout << "Time = " << TimePointsAvailable << endl;
     
    b_mean = mean(b_samples, 2);
    Sigma_mean = mean(Sigma_samples, 2);
    beta_mean = mean(beta_samples, 1);
    nu_mean = mean(nu_samples, 1);
    delta_mean = mean(delta_samples, 1);
    
    //Rcout <<"beta_mean =" << beta_mean << endl;
    
    //Rcout << "size(b_mean) = " << size(b_mean) << endl;
    rowvec X_tmp, Z_tmp, Ri_tmp;
    //vec mu_tmp;
    double pit, CPO_tmp, ESS=0, GP=0, ESS_GP_tmp, RJ1, RJ2;
    logL = 0.;
    mat Djt(Num_of_Timepoints, Num_of_covariates, fill::zeros);
    mat Omega_I(Num_of_covariates, Num_of_covariates, fill::zeros), M_LZ(Num_of_covariates, Num_of_covariates, fill::zeros);
    vec mu_it(Num_of_Timepoints), p_it(Num_of_Timepoints);
    mat A_sqrt, Cov_Y, V, V_inv, Omega_I_inv, V_LZ, Gamma_RJ, Djt_sub;
    
    mat CPO = zeros<mat>(Num_of_obs, TimePointsAvailable.max());
    int tp;
    
    CIC = 0.;
    RJ_R = 0.;
    ACC = 0.;

    //Rcout << "============== 0 ============"<<endl;
    for(int i=0; i<Num_of_obs; i++){
        //Rcout << "i = " << i << endl;
        tp = TimePointsAvailable(i);
        for(int t=0; t<tp; t++){
            
            
            
            
            X_tmp = X.slice(i).row(t);
            Z_tmp = Z.slice(i).row(t);
            //Rcout << "X_tmp = " << X_tmp << endl;
            //Rcout << "Z_tmp = " << Z_tmp << endl;
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
        if(updatedelta)
            V = A_sqrt*Ri_Version2(i, tp, delta_mean)*A_sqrt;
        else
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
    //Rcout << "============== 2 ============"<<endl;
    //Rcout << " Omega_I = " << endl << Omega_I << endl;
    Omega_I_inv = Omega_I.i();
    //Rcout << "============== 3 ============"<<endl;
    Gamma_RJ = Omega_I_inv*M_LZ;
    //Rcout << "============== 4 ============"<<endl;
    V_LZ = Gamma_RJ*Omega_I_inv;
    //Rcout << "============== 5 ============"<<endl;
    RJ1 = trace(Gamma_RJ)/Num_of_covariates;
    RJ2 = accu((diagvec(Gamma_RJ*Gamma_RJ)))/Num_of_covariates;

    RJ_R = sqrt((1-RJ1)*(1-RJ1)+(1-RJ2)*(1-RJ2));
    //SC = ESS/(N-P-a
    //GP = -0.5*GP

    CIC = trace(Omega_I*V_LZ);
    //Rcout << "RJ_R = " << RJ_R << "\tCIC = " << CIC << endl;
    //cat("RJ.R = ", RJ.R, "\t", "SC = ", SC, "\n")

    
    for(int i=0; i<Num_of_obs; i++){
        for(int t=0; t<TimePointsAvailable(i); t++){
            
            X_tmp = X.slice(i).row(t);
            Z_tmp = Z.slice(i).row(t);
                    
            
            for(int iter = Num_of_iterations/2; iter<Num_of_iterations; iter++){
                pit = normcdf(as_scalar( X_tmp*beta_samples.col(iter) + Z_tmp*b_samples.slice(iter).col(i) ) , 0., 1.);
                
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
    AIC = -2*logL + 2 * (Num_of_covariates+Num_of_obs*Num_of_RanEffs + Num_of_deltas);
    BIC = -2*logL + log(Num_of_obs) * (Num_of_covariates+Num_of_obs*Num_of_RanEffs+ Num_of_deltas);
    
}

SEXP ProbitMLModelSelection::MCMC_Procedure()
{
    Rcout << "Start running MCMC procedure:"<< endl;
    
    int percent = 0;
    //X1 <- rtmvnorm(n=20000, mean = c(-1, -2, 1, 3), sigma=C, lower=c(-Inf, -Inf, 0, 0), upper=c(0, 0, Inf, Inf), algorithm="gibbs", burn.in.samples=100, thinning=5)
    //cor(X1)
    
    //Y_star_sample(span(0, tp-1), i) = rtmvnorm_gibbs(1, mu_tmp, Ri_tmp, lower, upper, Y_star_sample.col(i) ).t();

    //Rcout << cor(rtmvnorm_gibbs(2000, vec , Ri_tmp, lower, upper, Y_star_sample.col(i) ))
    
    
    List PosteriorSamples;
    List PosteriorEstimates;
    List MH_AcceptanceRates;
    List Posterior;
    
    //time_t start = time(NULL);
    
    int iter = 0;

    while(iter < Num_of_iterations-1){
        Update_ystar_b_beta_Sigma(iter);
        
        if(updatenu)
            Update_nu(iter);
        else
            nu_samples.col(iter+1) = nu_samples.col(iter);
        

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
        if(updatedelta)
            Update_delta(iter);
        //else
            //delta_samples.col(iter+1) = delta_samples.col(iter);
        percent = (100 *iter) / (Num_of_iterations-2) ;
        iter++;
        /*
        if(percent%5==0){
            Rcout << "\r" << "[" << std::string(percent / 5, (char)61) << std::string(100 / 5 - percent / 5, ' ') << "]" << "\t" << percent << "%";
            //Rcout << percent << "%" << " [Iteration " << iter + 1 << " of " << Num_of_iterations << "]";
            std::cout.flush();
        }
        */
        
        //if(iter %100 == 0)
        //    Rcout << "iter = " << iter << endl;
            
        //Rcout << '\n' << std::flush;
        
        
        /*
        int step = 1;
        int displayNext = step;
        int percent = 0;

        cout << "Processing " << totalImagesCount << " images..." << endl;

        // loop through the image count
        for (size_t i = 0; i < totalImagesCount ; ++i)
        {
            // Individual image processing operations

                    // Formatted progress indicator
            percent = (100 * (i + 1)) / totalImagesCount ;
            if (percent >= displayNext)
            {
                cout << "\r" << "[" << std::string(percent / 5, (char)254u) << std::string(100 / 5 - percent / 5, ' ') << "]";
                cout << percent << "%" << " [Image " << i + 1 << " of " << totalImagesCount << "]";
                std::cout.flush();
                displayNext += step;
            }
        }

        */
        
        //printProgress( (double)percent);
        
        if(percent%2==0){
            Rcout << "\r" <<  "[" << std::string(percent / 2, (char)61) << std::string(100 / 2 - percent / 2, ' ') << "]" << "\t" << percent << "%";
            //Rcout << percent << "%" << " [Iteration " << iter + 1 << " of " << Num_of_iterations << "]";
            //std::cout.flush();
        }
        
        
        //if(iter %100 == 0)
            //Rcout << "iter = " << iter << endl;

    }
    Rcout << endl << "Finish MCMC Procedure." << endl;
  
    ParameterEstimation();
    
    //Rcout << endl << "============= MCMC: Done =============="<< endl;
    //time_t end = time(NULL);
    //Rcout << "Number of Iterations : " << Num_of_iterations << endl;
    //Rcout << "Execution Time: "<< (double)(end-start)<<" Seconds"<<std::endl;
    
    PosteriorSamples["ystar.samples"] = Y_star_sample;
    PosteriorSamples["b.samples"] = b_samples;
    PosteriorSamples["nu.samples"] = nu_samples;
    PosteriorSamples["beta.samples"] = beta_samples;
    PosteriorSamples["Sigma.samples"] = Sigma_samples;
    PosteriorSamples["delta.samples"] = delta_samples;
    if(updatedelta)
        PosteriorEstimates["delta.mean"] = delta_mean;

 
    PosteriorEstimates["beta.mean"] = beta_mean;
    PosteriorEstimates["nu.mean"] = nu_mean;
    
    
    PosteriorEstimates["b.mean"] = b_mean;
    PosteriorEstimates["Sigma.mean"] = Sigma_mean;


    PosteriorEstimates["AIC"] = AIC;
    PosteriorEstimates["BIC"] = BIC;
    PosteriorEstimates["CIC"] = CIC;
    PosteriorEstimates["logL"] = logL;

    PosteriorEstimates["DIC"] = DIC;
    PosteriorEstimates["RJR"] = RJ_R;

    PosteriorEstimates["MPL"] = MPL;
    PosteriorEstimates["ACC"] = ACC/accu(TimePointsAvailable);
    
//    PosteriorEstimates["pred.y"] = pred_y;
    MH_AcceptanceRates["Acceptance.rate.for.delta"] = acc_rate_delta/Num_of_iterations;
 
    Posterior["PosteriorEstimates"] = PosteriorEstimates;
    Posterior["PosteriorSamples"] = PosteriorSamples;
    Posterior["MH_AcceptanceRates"] = MH_AcceptanceRates;
    

    
    //Rcout << "Acceptance Rate in MH = " << acc_rate_delta/Num_of_iterations << endl;
    return (Posterior);
}

