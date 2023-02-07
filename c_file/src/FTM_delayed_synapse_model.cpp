/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   FTM_delayed_synapse_model.cpp
 * Author: picio
 * 
 * Created on 14 maggio 2017, 19.31
 */

#include "../inc/FTM_delayed_synapse_model.hpp"

FTM_delayed_synapse_model::FTM_delayed_synapse_model()
{
    this->nx = 0;
    this->nu = 0;
    this->theta = 0;
    this->g1 = 0;
    this->g2 = 0;
    this->Ndelay = 1;
    
    this->delays = new double[1];
    
    this->delays[0] = 0;
}

FTM_delayed_synapse_model::FTM_delayed_synapse_model(const FTM_delayed_synapse_model& other)
{

    this->nx = other.nx;
    this->nu = other.nu;
    this->theta = other.theta;
    this->g1 = other.g1;
    this->g2 = other.g2;
    this->Ndelay = other.Ndelay;
    
    this->delays = new double[1];
    
    this->delays[0] = other.delays[0];
}


FTM_delayed_synapse_model::FTM_delayed_synapse_model(double nu, double theta,double g1,double g2, double del)
{
    this->nx = 0;
    this->nu = nu;
    this->theta = theta;
    this->g1 = g1;
    this->g2 = g2;
    this->Ndelay = 1;
    this->delays = new double[1];
    
    this->delays[0] = del;
}


FTM_delayed_synapse_model *FTM_delayed_synapse_model::clone() const
{
    return new FTM_delayed_synapse_model(*this);
}


double FTM_delayed_synapse_model::getActivation(double *x,double Vpre)
{
    double nu = this->nu;
    double theta = this->theta;    
    
    return g1 * 1.0/(1.0+exp(-nu*(Vpre-theta)));
}

double FTM_delayed_synapse_model::getActivation(double *x,double Vpre, double *VpreOld)
{
    double nu = this->nu;
    double theta = this->theta;    
    
    return g1 * 1.0/(1.0+exp(-nu*(Vpre-theta))) + g2 *  1.0/(1.0+exp(-nu*(VpreOld[0]-theta)));
}

void FTM_delayed_synapse_model::getJacobian(double t, double *x, double **J,double Vpre)
{
    
}

void FTM_delayed_synapse_model::getDA(double *x, double **DA,double Vpre)
{

}


void FTM_delayed_synapse_model::getXdot(double t, double *x, double *xdot,double Vpre)
{

}

FTM_delayed_synapse_model::~FTM_delayed_synapse_model()
{
    delete[](this->delays);
}

