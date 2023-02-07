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

#ifndef ALPHA_DELAYED_SYNAPSE_MODEL_H
#define ALPHA_DELAYED_SYNAPSE_MODEL_H

#include "synapse_model.hpp"
#include <cmath>

class alpha_delayed_synapse_model : public synapse_model {
    
private:
    double alpha;
    double beta;
    double nu;
    double theta;
    double g2;
    double g1;
    
public:
    alpha_delayed_synapse_model();
    alpha_delayed_synapse_model(const alpha_delayed_synapse_model &other);
    alpha_delayed_synapse_model(double alpha, double beta,double nu,double theta,double del,double g2);
    
    virtual alpha_delayed_synapse_model *clone() const;
    
    void getXdot(double t, double *x, double *xdot, double Vpre, double **xold, double *VpreOld);
    void getXdot(double t, double *x, double *xdot, double Vpre);
    void getJacobian(double t, double *x, double **J, double Vpre);
    void getJacobian(double t, double *x, double **J, double Vpre, double **xold, double *VpreOld);
    
    double getActivation(double *x, double Vpre);
    double getActivation(double *x, double Vpre, double *VpreOld);
    void getDA(double *x, double **DA,double Vpre);
    
    virtual ~alpha_delayed_synapse_model();

};


#endif /* ALPHA_DELAYED_SYNAPSE_MODEL_H */

