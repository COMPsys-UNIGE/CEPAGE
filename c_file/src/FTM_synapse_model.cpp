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

#include "../inc/FTM_synapse_model.hpp"

FTM_synapse_model::FTM_synapse_model()
{
    this->nx = 0;
    this->nu = 0;
    this->theta = 0;
    this->Ndelay = 0;
}

FTM_synapse_model::FTM_synapse_model(double nu, double theta)
{
    this->nx = 0;
    this->nu = nu;
    this->theta = theta;
    this->Ndelay = 0;
}

FTM_synapse_model::FTM_synapse_model(const FTM_synapse_model &other)
{
    this->nx = other.nx;
    this->nu = other.nu;
    this->theta = other.theta;
    this->Ndelay = other.Ndelay;
}

FTM_synapse_model *FTM_synapse_model::clone() const
{
    return new FTM_synapse_model(*this);
}



double FTM_synapse_model::getActivation(double *x,double Vpre)
{
    double nu = this->nu;
    double theta = this->theta;    
    
    return 1.0/(1.0+exp(-nu*(Vpre-theta)));
}

void FTM_synapse_model::getDA(double *x, double **DA,double Vpre)
{
     double nu = this->nu;
    double theta = this->theta;    
    
    DA[0][0] = (nu/pow(1.0+exp(-nu*(Vpre-theta)),2)) - (nu/((1.0+exp(-nu*(Vpre-theta)))));
}

void FTM_synapse_model::getXdot(double t, double *x, double *xdot,double Vpre)
{
    
}

void FTM_synapse_model::getJacobian(double t, double *x, double **J,double Vpre)
{
    
}

FTM_synapse_model::~FTM_synapse_model() {
}

