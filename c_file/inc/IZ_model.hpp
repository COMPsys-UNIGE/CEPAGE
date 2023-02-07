/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   IZ_model.h
 * Author: teo
 *
 * Created on July 12, 2017, 10:06 AM
 */

#ifndef IZ_MODEL_H
#define IZ_MODEL_H


#include "neuron_model.hpp"

class IZ_model : public neuron_model {
    
private:
    double a;
    double b;
    double c;
    double d;
    double I;
    double gL;
    double El;
    
public:
    IZ_model();
    IZ_model(const IZ_model &iz);
    IZ_model(double a, double b, double c, double d, double I, double gL, double El);
    
    IZ_model *clone() const;
    
    void getXdot(double t, double *x, double *xdot,double *Iext);
    void getJacobian(double t, double *x, double **J);
    
    double dfdi(double t, double *x);
    
    bool getResetConditions(double *x);
    void resetStates(double *x);
    
    virtual ~IZ_model();

};

#endif /* IZ_MODEL_H */

