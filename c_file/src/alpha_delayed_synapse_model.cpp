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

#include "../inc/alpha_delayed_synapse_model.hpp"

alpha_delayed_synapse_model::alpha_delayed_synapse_model()
{
    this->nx = 2;
    this->alpha = 0;
    this->beta = 0;
    this->nu = 0;
    this->theta = 0;
    this->g2 = 0;
    this->g1 = 0;
    this->Ndelay = 1;
    
    this->delays = new double[1];
    
    this->delays[0] = 0;
}

alpha_delayed_synapse_model::alpha_delayed_synapse_model(const alpha_delayed_synapse_model &other)
{
    this->nx = other.nx;
    this->alpha = other.alpha;
    this->beta = other.beta;
    this->nu = other.nu;
    this->theta = other.theta;
    this->g2 = other.g2;
    this->g1 = other.g1;
    this->Ndelay = other.Ndelay;
    
    this->delays = new double[1];
    
    this->delays[0] = other.delays[0];
}

alpha_delayed_synapse_model::alpha_delayed_synapse_model(double alpha, double beta, double nu, double theta, double del,double g2)
{
    this->nx = 2;
    this->alpha = alpha;
    this->beta = beta;
    this->nu = nu;
    this->theta = theta;
    this->g2 = g2;
    this->g1 = 1;
    this->Ndelay = 1;
    this->delays = new double[1];
    
    this->delays[0] = del;
}


alpha_delayed_synapse_model *alpha_delayed_synapse_model::clone() const
{
    return new alpha_delayed_synapse_model(*this);
}


double alpha_delayed_synapse_model::getActivation(double *x,double Vpre)
{
    double g1 = this->g1;
    double g2 = this->g2;   
    double alpha = this->alpha;
    double beta = this->beta;
    return g1*((alpha+beta)/alpha)*x[0]+g2*((alpha+beta)/alpha)*x[1];
}

double alpha_delayed_synapse_model::getActivation(double *x,double Vpre, double *VpreOld)
{
    double nu = this->nu;
    double theta = this->theta;   
    double alpha = this->alpha;
    double beta = this->beta;
    double g1 = this->g1;
    double g2 = this->g2;  
    
    return g1*((alpha+beta)/alpha)*x[0] + g2 *  1.0/(1.0+exp(-nu*(VpreOld[0]-(-50))));
}



void alpha_delayed_synapse_model::getXdot(double t, double *x, double *xdot, double Vpre, double **xold, double *VpreOld)
{
    double nu = this->nu;
    double theta = this->theta;    
    double alpha = this->alpha;
    double beta = this->beta;
    
    xdot[0] = alpha*(1.0-x[0])/(1.0+exp(-nu*(Vpre-theta)))-beta*x[0];
    xdot[1] = alpha*(1.0-x[1])/(1.0+exp(-nu*(VpreOld[0]-theta)))-beta*x[1];
}
void alpha_delayed_synapse_model::getXdot(double t, double *x, double *xdot, double Vpre)
{
}

void alpha_delayed_synapse_model::getJacobian(double t, double *x, double **J, double Vpre, double **xold, double *VpreOld)
{   
    J[0][0] = 0;
    J[0][1] = 0;
    J[1][0] = 0;
    J[1][1] = 0;
}

void alpha_delayed_synapse_model::getJacobian(double t, double *x, double **J, double Vpre)
{   
}

void alpha_delayed_synapse_model::getDA(double *x, double **DA,double Vpre)
{
    double g1 = this->g1;
    double g2 = this->g2;   
    double alpha = this->alpha;
    double beta = this->beta;
    
    DA[0][0] = 0;
    DA[0][1] = g1*((alpha+beta)/alpha);
    DA[0][2] = g2*((alpha+beta)/alpha);
}

alpha_delayed_synapse_model::~alpha_delayed_synapse_model() 
{
    delete[](this->delays);
}

