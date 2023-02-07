/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   FN_relaxation_model.h
 * Author: teo
 *
 * Created on July 12, 2017, 9:48 AM
 */

#ifndef FN_RELAXATION_MODEL_H
#define FN_RELAXATION_MODEL_H



#include "neuron_model.hpp"
#include <cmath>

class FN_relaxation_model : public neuron_model {
    
private:
    double eps;
    double I;
    
public:
    FN_relaxation_model();
    FN_relaxation_model(double I,double eps);
    FN_relaxation_model(const FN_relaxation_model &f);
    
    FN_relaxation_model *clone() const;
    
    void getXdot(double t, double *x, double *xdot,double *Iext);
    void getJacobian(double t, double *x, double **J);
    
    double dfdi(double t, double *x);
    
    virtual ~FN_relaxation_model();

};

#endif /* FN_RELAXATION_MODEL_H */

