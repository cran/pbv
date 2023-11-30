//// File Name: RcppExports.cpp
//// File Version: 0.005047
// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <string>
#include <set>

using namespace Rcpp; using namespace arma;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// pbv_rcpp_pnorm0
double pbv_rcpp_pnorm0(double z);
static SEXP _pbv_pbv_rcpp_pnorm0_try(SEXP zSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< double >::type z(zSEXP);
    rcpp_result_gen = Rcpp::wrap(pbv_rcpp_pnorm0(z));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _pbv_pbv_rcpp_pnorm0(SEXP zSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_pbv_pbv_rcpp_pnorm0_try(zSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error("%s", CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// pbv_rcpp_pnorm
Rcpp::NumericVector pbv_rcpp_pnorm(Rcpp::NumericVector x);
static SEXP _pbv_pbv_rcpp_pnorm_try(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(pbv_rcpp_pnorm(x));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _pbv_pbv_rcpp_pnorm(SEXP xSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_pbv_pbv_rcpp_pnorm_try(xSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error("%s", CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// pbv_rcpp_pbvnorm0
double pbv_rcpp_pbvnorm0(double h1, double hk, double r);
static SEXP _pbv_pbv_rcpp_pbvnorm0_try(SEXP h1SEXP, SEXP hkSEXP, SEXP rSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< double >::type h1(h1SEXP);
    Rcpp::traits::input_parameter< double >::type hk(hkSEXP);
    Rcpp::traits::input_parameter< double >::type r(rSEXP);
    rcpp_result_gen = Rcpp::wrap(pbv_rcpp_pbvnorm0(h1, hk, r));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _pbv_pbv_rcpp_pbvnorm0(SEXP h1SEXP, SEXP hkSEXP, SEXP rSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_pbv_pbv_rcpp_pbvnorm0_try(h1SEXP, hkSEXP, rSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error("%s", CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// pbv_rcpp_pbvnorm
Rcpp::NumericVector pbv_rcpp_pbvnorm(Rcpp::NumericVector x, Rcpp::NumericVector y, Rcpp::NumericVector rho);
static SEXP _pbv_pbv_rcpp_pbvnorm_try(SEXP xSEXP, SEXP ySEXP, SEXP rhoSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type rho(rhoSEXP);
    rcpp_result_gen = Rcpp::wrap(pbv_rcpp_pbvnorm(x, y, rho));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _pbv_pbv_rcpp_pbvnorm(SEXP xSEXP, SEXP ySEXP, SEXP rhoSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_pbv_pbv_rcpp_pbvnorm_try(xSEXP, ySEXP, rhoSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error("%s", CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// pbv_rcpp_dbvnorm0
double pbv_rcpp_dbvnorm0(double x, double y, double rho, bool use_log);
static SEXP _pbv_pbv_rcpp_dbvnorm0_try(SEXP xSEXP, SEXP ySEXP, SEXP rhoSEXP, SEXP use_logSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type y(ySEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< bool >::type use_log(use_logSEXP);
    rcpp_result_gen = Rcpp::wrap(pbv_rcpp_dbvnorm0(x, y, rho, use_log));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _pbv_pbv_rcpp_dbvnorm0(SEXP xSEXP, SEXP ySEXP, SEXP rhoSEXP, SEXP use_logSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_pbv_pbv_rcpp_dbvnorm0_try(xSEXP, ySEXP, rhoSEXP, use_logSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error("%s", CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// pbv_rcpp_dbvnorm
Rcpp::NumericVector pbv_rcpp_dbvnorm(Rcpp::NumericVector x, Rcpp::NumericVector y, Rcpp::NumericVector rho, bool use_log);
static SEXP _pbv_pbv_rcpp_dbvnorm_try(SEXP xSEXP, SEXP ySEXP, SEXP rhoSEXP, SEXP use_logSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< bool >::type use_log(use_logSEXP);
    rcpp_result_gen = Rcpp::wrap(pbv_rcpp_dbvnorm(x, y, rho, use_log));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _pbv_pbv_rcpp_dbvnorm(SEXP xSEXP, SEXP ySEXP, SEXP rhoSEXP, SEXP use_logSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_pbv_pbv_rcpp_dbvnorm_try(xSEXP, ySEXP, rhoSEXP, use_logSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error("%s", CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}

// validate (ensure exported C++ functions exist before calling them)
static int _pbv_RcppExport_validate(const char* sig) { 
    static std::set<std::string> signatures;
    if (signatures.empty()) {
        signatures.insert("double(*pbv_rcpp_pnorm0)(double)");
        signatures.insert("Rcpp::NumericVector(*pbv_rcpp_pnorm)(Rcpp::NumericVector)");
        signatures.insert("double(*pbv_rcpp_pbvnorm0)(double,double,double)");
        signatures.insert("Rcpp::NumericVector(*pbv_rcpp_pbvnorm)(Rcpp::NumericVector,Rcpp::NumericVector,Rcpp::NumericVector)");
        signatures.insert("double(*pbv_rcpp_dbvnorm0)(double,double,double,bool)");
        signatures.insert("Rcpp::NumericVector(*pbv_rcpp_dbvnorm)(Rcpp::NumericVector,Rcpp::NumericVector,Rcpp::NumericVector,bool)");
    }
    return signatures.find(sig) != signatures.end();
}

// registerCCallable (register entry points for exported C++ functions)
RcppExport SEXP _pbv_RcppExport_registerCCallable() { 
    R_RegisterCCallable("pbv", "_pbv_pbv_rcpp_pnorm0", (DL_FUNC)_pbv_pbv_rcpp_pnorm0_try);
    R_RegisterCCallable("pbv", "_pbv_pbv_rcpp_pnorm", (DL_FUNC)_pbv_pbv_rcpp_pnorm_try);
    R_RegisterCCallable("pbv", "_pbv_pbv_rcpp_pbvnorm0", (DL_FUNC)_pbv_pbv_rcpp_pbvnorm0_try);
    R_RegisterCCallable("pbv", "_pbv_pbv_rcpp_pbvnorm", (DL_FUNC)_pbv_pbv_rcpp_pbvnorm_try);
    R_RegisterCCallable("pbv", "_pbv_pbv_rcpp_dbvnorm0", (DL_FUNC)_pbv_pbv_rcpp_dbvnorm0_try);
    R_RegisterCCallable("pbv", "_pbv_pbv_rcpp_dbvnorm", (DL_FUNC)_pbv_pbv_rcpp_dbvnorm_try);
    R_RegisterCCallable("pbv", "_pbv_RcppExport_validate", (DL_FUNC)_pbv_RcppExport_validate);
    return R_NilValue;
}

static const R_CallMethodDef CallEntries[] = {
    {"_pbv_pbv_rcpp_pnorm0", (DL_FUNC) &_pbv_pbv_rcpp_pnorm0, 1},
    {"_pbv_pbv_rcpp_pnorm", (DL_FUNC) &_pbv_pbv_rcpp_pnorm, 1},
    {"_pbv_pbv_rcpp_pbvnorm0", (DL_FUNC) &_pbv_pbv_rcpp_pbvnorm0, 3},
    {"_pbv_pbv_rcpp_pbvnorm", (DL_FUNC) &_pbv_pbv_rcpp_pbvnorm, 3},
    {"_pbv_pbv_rcpp_dbvnorm0", (DL_FUNC) &_pbv_pbv_rcpp_dbvnorm0, 4},
    {"_pbv_pbv_rcpp_dbvnorm", (DL_FUNC) &_pbv_pbv_rcpp_dbvnorm, 4},
    {"_pbv_RcppExport_registerCCallable", (DL_FUNC) &_pbv_RcppExport_registerCCallable, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_pbv(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
