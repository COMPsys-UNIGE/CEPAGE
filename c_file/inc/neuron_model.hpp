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

#ifndef NEURON_MODEL_H
#define NEURON_MODEL_H

#include "dynSys.hpp"

class neuron_model : public dynSys{
    
    protected:
        int nx;

    
    
public:
    
    virtual neuron_model *clone() const = 0;
    
    virtual void getXdot(double t, double *x, double *xdot,double *Iext) = 0;
    
    /* Xprec is a Ndelays x Nstates matrix */
    virtual void getXdot(double t, double *x, double *xdot,double *Iext,double **Xprec);

    virtual void getJacobian(double t, double *x, double **J) = 0;
    virtual void getJacobian(double t, double *x, double **J,double **Xprec);
    
    virtual double dfdi(double t, double *x)  = 0;
    
    virtual bool getResetConditions(double *x);
    virtual void resetStates(double *x);
    
        
    int getnx();

    virtual ~neuron_model()
    {}
    
};

#endif /* NEURON_MODEL_H */

