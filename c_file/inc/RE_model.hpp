
#ifndef RE_MODEL_H
#define RE_MODEL_H

#include "neuron_model.hpp"
#include <cmath>

class RE_model : public neuron_model {
    
private:
    double gna; 
    double ENa; 
    double gk;
    double Ek;
    double gca;
    double gl;
    double El; 
    double C; 
    double Iapp;
    double Ca0; 
    double d;
    double KT;
    double Kd;
    double gd;
    double D;
    double EsynEx;
    double tau;
    double tTm1;
    double tTm2;
    double tTh1;
    double tTh2;

    
public:
    RE_model();
    RE_model(const RE_model &h);
    RE_model(double gna, double ENa, double gk, double Ek, double gca, double gl, double El, double C, double Iapp, double Ca0,double d, double KT, double Kd, double gd, double D, double EsynEx, double tau,double tTm1,double tTm2,double tTh1,double tTh2);
    
    virtual RE_model *clone() const;
    
    void getXdot(double t, double *x, double *xdot,double *Isyn);
    void getJacobian(double t, double *x, double **J);
    
    double dfdi(double t, double *x);
    
    virtual ~RE_model();

};

#endif /* RE_MODEL_H */