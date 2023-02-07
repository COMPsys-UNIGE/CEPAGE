/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   FTM_synapse_model.h
 * Author: teo
 *
 * Created on July 12, 2017, 10:21 AM
 */

#ifndef A_DELAYED_SYNAPSE_MODEL_H
#define A_DELAYED_SYNAPSE_MODEL_H

#include "synapse_model.hpp"
#include <cmath>

class A_delayed_synapse_model : public synapse_model {
    
private:
    double alpha;
    double beta;
    double nu;
    double theta;
    double g2;
    double g1;
public:
    A_delayed_synapse_model();
    A_delayed_synapse_model(const A_delayed_synapse_model &other);
    A_delayed_synapse_model(double alpha, double beta,double nu,double theta,double del,double g2);
    
    virtual A_delayed_synapse_model *clone() const;
    
    void getXdot(double t, double *x, double *xdot, double Vpre, double **xold, double *VpreOld);
    void getXdot(double t, double *x, double *xdot, double Vpre);
    void getJacobian(double t, double *x, double **J, double Vpre);
    void getJacobian(double t, double *x, double **J, double Vpre, double **xold, double *VpreOld);
    
    double getActivation(double *x, double Vpre);
    void getDA(double *x, double **DA,double Vpre);
    
    virtual ~A_delayed_synapse_model();

};


#endif /* FTM_SYNAPSE_MODEL_H */