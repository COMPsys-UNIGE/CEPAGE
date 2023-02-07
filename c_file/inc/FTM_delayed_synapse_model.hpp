/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/*
 * File:   FTM_delayed_synapse_model.h
 * Author: teo
 *
 * Created on July 12, 2017, 10:21 AM
 */

#ifndef FTM_DELAYED_SYNAPSE_MODEL_H
#define FTM_DELAYED_SYNAPSE_MODEL_H

#include "synapse_model.hpp"
#include <cmath>


class FTM_delayed_synapse_model : public synapse_model {
    
private:
    double nu;
    double theta;
    double g1;
    double g2;
    
public:
    FTM_delayed_synapse_model();
    FTM_delayed_synapse_model(const FTM_delayed_synapse_model& other);
    FTM_delayed_synapse_model(double nu,double theta, double g1, double g2, double del);
    
    virtual FTM_delayed_synapse_model *clone() const;    
    
    void getXdot(double t, double *x, double *xdot,double Vpre);
    double getActivation(double *x, double Vpre);
    double getActivation(double *x, double Vpre, double *VpreOld);
    void getJacobian(double t, double *x, double **J,double Vpre);
    void getDA(double *x, double **DA,double Vpre);
    
    virtual ~FTM_delayed_synapse_model();
    
};


#endif /* FTM_SYNAPSE_MODEL_H */

