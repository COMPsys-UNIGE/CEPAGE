
#ifndef ST_MODEL_H
#define ST_MODEL_H

#include "neuron_model.hpp"
#include <cmath>

class ST_model : public neuron_model {
    
private:
    double g_Na;
    double g_K;
    double g_L;
    double g_A;
    double g_T;
    double g_KCa;
    double g_HVA;
    double epsilon;
    double E_Na;
    double E_K;
    double E_L;
    double E_Ca;
    double C;
    double A;
    double y0;
    double V_c;
    double w;
    double k_Ca;
    double alpha;
    double k;
    double nu_m;
    double nu_h;
    double nu_n;
    double nu_nA;
    double nu_hA;
    double nu_mT;
    double nu_hT;
    double nu_mHVA;
    double s_m;
    double s_h;
    double s_n;
    double s_nA;
    double s_hA;
    double s_mT;
    double s_hT;
    double s_mHVA;
    double tau_nA;
    double tau_hA;
    double tau_hT;
    double tau_mHVA;
    double I_app;

    
public:
    ST_model();
    ST_model(const ST_model &h);
    ST_model(double g_Na,double g_K,double g_L,double g_A,double g_T,double g_KCa,double g_HVA,double epsilon,double E_Na,double E_K,double E_L,double E_Ca,double C,double A,double y0,double V_c,double w,double k_Ca,double alpha,double k,double nu_m,double nu_h,double nu_n,double nu_nA,double nu_hA,double nu_mT,double nu_hT,double nu_mHVA,double s_m,double s_h,double s_n,double s_nA,double s_hA,double s_mT,double s_hT,double s_mHVA,double tau_nA,double tau_hA,double tau_hT,double tau_mHVA,double I_app);
    
    virtual ST_model *clone() const;
    
    void getXdot(double t, double *x, double *xdot,double *Isyn);
    void getJacobian(double t, double *x, double **J);
    
    double dfdi(double t, double *x);
    
    virtual ~ST_model();

};

#endif