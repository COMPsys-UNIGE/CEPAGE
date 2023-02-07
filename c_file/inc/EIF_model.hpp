/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   EIF_model.h
 * Author: teo
 *
 * Created on July 12, 2017, 10:06 AM
 */

#ifndef EIF_MODEL_H
#define EIF_MODEL_H


#include "neuron_model.hpp"
#include <cmath>

class EIF_model : public neuron_model {
    
private:
    double a;
    double b;
    double C;
    double VT;
    double DT;
    double I;
    double gL;
    double EL;
    double Vr;
    double tw;   
    double gexp;
    double gd;
    double D;
    double EsynEx;
    
public:
    EIF_model();
    EIF_model(double C,double a,double b,double VT,double DT,double I,double gL,double EL,double Vr,double tw,double gexp,double gd,double D,double EsynEx);
    EIF_model(const EIF_model &e);
    
    virtual EIF_model *clone() const;
    
    void getXdot(double t, double *x, double *xdot,double *Iext);
    void getJacobian(double t, double *x, double **J);
    
    double dfdi(double t, double *x);
    
    bool getResetConditions(double *x);
    void resetStates(double *x);
    
    virtual ~EIF_model();

};

#endif /* EIF_MODEL_H */

