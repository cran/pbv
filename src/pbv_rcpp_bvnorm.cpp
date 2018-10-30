//// File Name: pbv_rcpp_bvnorm.cpp
//// File Version: 0.541



// [[Rcpp::depends(RcppArmadillo)]]

// #include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::interfaces(r, cpp)]]



///********************************************************************
///** pbv_rcpp_pnorm0
// [[Rcpp::export]]
double pbv_rcpp_pnorm0( double z)
{
    double y = ::Rf_pnorm5(z, 0.0, 1.0, 1, 0);
    //--- OUTPUT
    return y;
}
///********************************************************************

///********************************************************************
///** pbv_rcpp_pnorm
// [[Rcpp::export]]
Rcpp::NumericVector pbv_rcpp_pnorm( Rcpp::NumericVector x)
{
    int N = x.size();
    Rcpp::NumericVector y(N);
    for (int nn=0; nn<N; nn++){
        y[nn] = pbv_rcpp_pnorm0( x[nn] );        
    }    
    //--- OUTPUT
    return y;
}
///********************************************************************


///********************************************************************
///** Drezner & Wesolowksy, 1990, JCSC
///** pbv_rcpp_pbvnorm0
// [[Rcpp::export]]
double pbv_rcpp_pbvnorm0( double h1, double hk, double r)
{
    double bv=0;
    int NX=5;
    Rcpp::NumericVector X(NX);
    Rcpp::NumericVector W(NX);
    // data
    // x = Array(0.04691008, 0.23076534, 0.5, 0.76923466, 0.95308992)
    // W = Array(0.018854042, 0.038088059, 0.0452707394, 0.038088059, 0.018854042)    
    X[0]=.04691008;
    X[1]=.23076534;
    X[2]=.5;
    X[3]=.76923466;
    X[4]=.95308992;
    W[0]=.018854042;
    W[1]=.038088059;
    W[2]=.0452707394;
    W[3]=.038088059;
    W[4]=.018854042;
    // declarations
    double r1=0;
    double r2=0;
    double rr=0;
    double rr2=0;
    double r3=0;
    double h3=0;
    double h5=0;
    double h6=0;
    double h7=0;    
    double cor_max = 0.7;
    
    // computation
    double h2 = hk;
    double h12 = (h1*h1+h2*h2)/2;        
    double r_abs = std::abs(r);    
    if (r_abs > cor_max){
        r2 = 1.0 - r*r;
        r3 = std::sqrt(r2);
        if (r<0){
            h2 = -h2;
        }
        h3 = h1*h2;
        h7 = std::exp( -h3 / 2.0);            
        if ( r_abs < 1){
            h6 = std::abs(h1-h2);
            h5 = h6*h6 / 2.0;
            h6 = h6 / r3;
            double aa = 0.5 - h3 / 8.0;
            double ab = 3.0 - 2.0 * aa * h5;        
            bv = 0.13298076*h6*ab*(1-pbv_rcpp_pnorm0(h6))-std::exp(-h5/r2)*(ab + aa*r2)*0.053051647;
            for (int ii=0; ii<NX; ii++){
                r1 = r3*X[ii];
                rr = r1*r1;
                r2 = std::sqrt( 1.0 - rr);
                bv += - W[ii]*std::exp(- h5/rr)*(std::exp(-h3/(1.0+r2))/r2/h7 - 1.0 - aa*rr);                
            }                    
        }        
        double h11 = std::min(h1,h2);
        bv = bv*r3*h7 + pbv_rcpp_pnorm0(h11);            
        if (r < 0){
            bv = pbv_rcpp_pnorm0(h1) - bv;
        }
        
    } else {        
        h3=h1*h2;
        for (int ii=0; ii<NX; ii++){
            r1 = r*X[ii];
            rr2 = 1.0 - r1*r1;
            bv += W[ii] * std::exp(( r1*h3 - h12)/rr2)/ std::sqrt(rr2);
        }
        bv = pbv_rcpp_pnorm0(h1)*pbv_rcpp_pnorm0(h2) + r*bv;
    }    
    //--- OUTPUT
    return bv;
}
///********************************************************************

///********************************************************************
///** pbv_rcpp_pbvnorm
// [[Rcpp::export]]
Rcpp::NumericVector pbv_rcpp_pbvnorm( Rcpp::NumericVector x, Rcpp::NumericVector y,
        Rcpp::NumericVector rho)
{
    int N = x.size();
    Rcpp::NumericVector res(N);
    for (int ii=0; ii<N; ii++){
        res[ii] = pbv_rcpp_pbvnorm0(x[ii], y[ii], rho[ii]);
    }    
    //--- OUTPUT
    return res;
}
///********************************************************************

const double pi = 3.1415926535897;

///********************************************************************
///** pbv_rcpp_dbvnorm0
// [[Rcpp::export]]
double pbv_rcpp_dbvnorm0( double x, double y, double rho, bool use_log)
{
    double pi2 = 2*pi;
    double r2 = 1-rho*rho;
    double r3 = std::sqrt(r2);
    double z = x*x - 2*rho*x*y + y*y;
    z = - z / r2 / 2.0;    
    if ( ! use_log ){
        z = std::exp(z) / pi2 / r3;        
    } else {
        z += - std::log(r3*pi2);
    }    
    //--- OUTPUT
    return z;
}
///********************************************************************

///********************************************************************
///** pbv_rcpp_dbvnorm
// [[Rcpp::export]]
Rcpp::NumericVector pbv_rcpp_dbvnorm( Rcpp::NumericVector x, Rcpp::NumericVector y,
        Rcpp::NumericVector rho, bool use_log)
{
    int N = x.size();
    Rcpp::NumericVector res(N);
    for (int ii=0; ii<N; ii++){
        res[ii] = pbv_rcpp_dbvnorm0(x[ii], y[ii], rho[ii], use_log);
    }    
    //--- OUTPUT
    return res;
}
///********************************************************************
