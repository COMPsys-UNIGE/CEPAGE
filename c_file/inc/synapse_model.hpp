/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/*
 * File:   neuron_model.h
 * Author: picio
 *
 * Created on 14 maggio 2017, 19.28
 */

#ifndef SYNAPSE_MODEL_H
#define SYNAPSE_MODEL_H


class synapse_model
{
    
protected:
    int nx;
    
    int Ndelay = 0;
    double *delays = {};
    
public:
    
    virtual synapse_model *clone() const = 0;
    
    virtual double getActivation(double *x,double Vpre, double *VpreOld);
    virtual double getActivation(double *x,double Vpre) = 0;
    
    /* Xold is a Ndelays x Nstates matrix */
    virtual void getXdot(double t, double *x, double *xdot, double Vpre, double **xold, double *VpreOld);
    virtual void getXdot(double t, double *x, double *xdot, double Vpre) = 0;
  
    virtual void getJacobian(double t, double *x, double **J, double Vpre, double **xold, double *VpreOld);
    virtual void getJacobian(double t, double *x, double **J, double Vpre) = 0;
    
    virtual void getDA(double *x, double **DA,double Vpre) = 0;
    //virtual bool getResetConditions(double *x, double Vpre, double *VpreOld);
    virtual bool getResetConditions(double *x, double Vpre);
    
    //virtual void resetStates(double *x, double Vpre, double *VpreOld);
    virtual void resetStates(double *x, double Vpre);
    
    //virtual void getJacobian(double t, double *x, double **J) = 0;
    //virtual void getDA(double t, double *x) = 0;
    
    int getnx();
    
    
    int getNdelay();
    void getDelays(double *del);
    
    virtual ~synapse_model()
    {}
    
};

#endif /* NEURON_MODEL_H */

