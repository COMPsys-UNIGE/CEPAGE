
#ifndef PLANT_MODEL_H
#define PLANT_MODEL_H

#include "neuron_model.hpp"
#include <cmath>

class Plant_model : public neuron_model {
    
private:
    double g_I;
    double g_K;
    double g_T;
    double g_L;
    double g_KCa;
    double V_I;
    double V_K;
    double V_L;
    double V_Ca;
    double K_c;
    double rho;
    double Iapp;
    
public:
    Plant_model();
    Plant_model(const Plant_model &h);
    Plant_model(double g_I, double g_K, double g_T, double g_L, double g_KCa, double V_I, double V_K, double V_L, double V_Ca, double K_c, double rho,double Iapp);
    
    virtual Plant_model *clone() const;
    
    void getXdot(double t, double *x, double *xdot,double *Isyn);
    void getJacobian(double t, double *x, double **J);
    
    double dfdi(double t, double *x);
    
    virtual ~Plant_model();

};

#endif /* Plant_model_H */