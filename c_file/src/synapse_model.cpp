/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   synapse_model.cpp
 * Author: picio
 * 
 * Created on 14 maggio 2017, 19.28
 */

#include "../inc/synapse_model.hpp"


int synapse_model::getnx() 
{
    return this->nx;
}


/*bool synapse_model::getResetConditions(double *x, double Vpre, double *VpreOld)
{
    this->getResetConditions(x,Vpre);
}*/

bool synapse_model::getResetConditions(double *x, double Vpre)
{
    return false;
}



/*void synapse_model::resetStates(double *x, double Vpre, double *VpreOld)
{
    this->resetStates(x,Vpre);
}*/

void synapse_model::resetStates(double *x, double Vpre)
{
    return;
}


void synapse_model::getXdot(double t, double *x, double *xdot, double Vpre, double **xold, double *VpreOld)
{
    this->getXdot(t,x,xdot,Vpre);
}

void synapse_model::getJacobian(double t, double *x, double **J, double Vpre, double **xold, double *VpreOld)
{
    this->getJacobian(t,x,J,Vpre);
}



double synapse_model::getActivation(double *x,double Vpre, double *VpreOld)
{
    return this->getActivation(x,Vpre);
}


int synapse_model::getNdelay()
{
    return Ndelay;
}

void synapse_model::getDelays(double *del)
{
    int i;
    for(i=0;i<Ndelay;i++)
        del[i] = this->delays[i];
}