
#ifndef A_MODEL_H
#define A_MODEL_H

#include "neuron_model.hpp"
#include <cmath>

class A_model : public neuron_model {
    
private:
    double Ca_shift; 
    double x_shift; 
    double gh;
    double Vhh;
    double Iext;

public:
    A_model();
    A_model(const A_model &h);
    A_model(double Ca_shift, double x_shift, double gh, double Vhh, double Iext);
    
    virtual A_model *clone() const;
    
    void getXdot(double t, double *x, double *xdot,double *Isyn);
    void getJacobian(double t, double *x, double **J);
    
    double dfdi(double t, double *x);
    
    virtual ~A_model();

};

#endif /* A_MODEL_H */