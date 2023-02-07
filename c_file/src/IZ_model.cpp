/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/*
 * File:   IZ_model.cpp
 * Author: teo
 *
 * Created on July 12, 2017, 10:06 AM
 */

#include "../inc/IZ_model.hpp"

IZ_model::IZ_model()
{
    this->nx = 2;
    this->a = 0;
    this->b = 0;
    this->c = 0;
    this->d = 0;
    this->I = 0;
    this->gL = 0;
    this->El = 0;
}

IZ_model::IZ_model(const IZ_model &iz)
{
    this->nx = iz.nx;
    this->a = iz.a;
    this->b = iz.b;
    this->c = iz.c;
    this->d = iz.d;
    this->I = iz.I;
    this->gL = iz.gL;
    this->El = iz.El;
}


IZ_model::IZ_model(double a, double b, double c, double d, double I, double gL, double El)
{
    this->nx = 2;
    this->a = a;
    this->b = b;
    this->c = c;
    this->d = d;
    this->I = I;
    this->gL = gL;
    this->El = El;
}

IZ_model *IZ_model::clone() const
{
    return new IZ_model(*this);
}

void IZ_model::getXdot(double t, double *x, double *xdot,double *Iext)
{
    double a = this->a;
    double b = this->b;
    double I = this->I;
    double gL = this->gL;
    double El =this->El;
    
    double v = x[0];
    double u = x[1];
    
    
    xdot[0] = 0.04*v*v+5*v+140-u+I-gL*(v-El)+Iext[0];
    xdot[1] = a*(b*v-u);
    
}

void IZ_model::getJacobian(double t, double *x, double **J)
{
    double a = this->a;
    double b = this->b;
    double I = this->I;
    double gL = this->gL;
    double El =this->El;
    
    J[0][0] = 2.0*0.04*x[0]+5-gL;
    J[0][1] = -1.0;
    J[1][0] = a*b;
    J[1][1] = -a;
    
}

double IZ_model::dfdi(double t, double *x)
{
    return 1;
}


bool IZ_model::getResetConditions(double *x)
{
    if(x[0] > 30)
        return true;
    else
        return false;
}

void IZ_model::resetStates(double *x)
{
    if(this->getResetConditions(x))
    {
        x[0] = this->c;
        x[1] += this->d;
    }
    return;
}

IZ_model::~IZ_model()
{
}

