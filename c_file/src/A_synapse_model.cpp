/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   FTM_synapse_model.cpp
 * Author: picio
 * 
 * Created on 14 maggio 2017, 19.31
 */

#include "../inc/A_synapse_model.hpp"
#include <mex.h>

A_synapse_model::A_synapse_model()
{
    this->nx = 1;
    this->alpha = 0;
    this->beta = 0;
    this->nu = 0;
    this->theta = 0;
}

A_synapse_model::A_synapse_model(const A_synapse_model &other)
{
    this->nx = other.nx;
    this->alpha = other.alpha;
    this->beta = other.beta;
    this->nu = other.nu;
    this->theta = other.theta;
}

A_synapse_model::A_synapse_model(double alpha, double beta, double nu, double theta)
{
    this->nx = 1;
    this->alpha = alpha;
    this->beta = beta;
    this->nu = nu;
    this->theta = theta;
}


A_synapse_model *A_synapse_model::clone() const
{
    return new A_synapse_model(*this);
}


double A_synapse_model::getActivation(double *x,double Vpre)
{
    return x[0];
}


void A_synapse_model::getXdot(double t, double *x, double *xdot,double Vpre)
{
    double nu = this->nu;
    double theta = this->theta;    
    double alpha = this->alpha;
    double beta = this->beta;
    
//     mexPrintf("alf =%f\n",alpha);
    
    xdot[0] = alpha*x[0]*(1.0-x[0])/(1.0+exp(-nu*(Vpre-theta)))-beta*x[0];
}

void A_synapse_model::getJacobian(double t, double *x, double **J,double Vpre)
{
    double nu = this->nu;
    double theta = this->theta;    
    double alpha = this->alpha;
    double beta = this->beta;
    
    J[0][0] = alpha*x[0]*(1.0-x[0])*((nu/((1.0+exp(-nu*(Vpre-theta)))*(1.0+exp(-nu*(Vpre-theta))))) - (nu/((1.0+exp(-nu*(Vpre-theta))))));
    J[0][1] = -alpha/(1.0+exp(-nu*(Vpre-theta)))*(1-2*x[0]) - beta;
}

void A_synapse_model::getDA(double *x, double **DA,double Vpre)
{
    DA[0][0] = 0;
    DA[0][1] = 1;
}

A_synapse_model::~A_synapse_model() {
}

