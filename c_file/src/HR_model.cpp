/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   HR_model.cpp
 * Author: picio
 * 
 * Created on 14 maggio 2017, 19.31
 */

#include "../inc/HR_model.hpp"

HR_model::HR_model()
{
    this->nx = 3;
    this->b = 0;
    this->I = 0;
    this->mu =  0;
    this->s = 0;
    this->x_rest = 0;
}

HR_model::HR_model(const HR_model &h) 
{
    this->nx = h.nx;
    this->b = h.b;
    this->I = h.I;
    this->mu =  h.mu;
    this->s = h.s;
    this->x_rest = h.x_rest;
}

HR_model::HR_model(double b, double I, double mu, double s, double x_rest) 
{
    this->nx = 3;
    this->b = b;
    this->I = I;
    this->mu =  mu;
    this->s = s;
    this->x_rest = x_rest;
}

HR_model *HR_model::clone() const
{
    return new HR_model(*this);
}


void HR_model::getXdot(double t, double *x, double *xdot,double *Iext)
{
    double b = this->b;
    double I = this->I;
    double mu = this->mu;
    double s = this->s;
    double x_rest = this->x_rest;
    
    double xx = x[0];
    double yy = x[1];
    double zz = x[2];
            
    double x2 = xx*xx;
    
    xdot[0] = yy - xx*x2 + b*x2 - zz + I +Iext[0]; 
    xdot[1] = 1-5*x2 - yy;
    xdot[2] = mu*(s*(xx - x_rest) - zz);
    
    
}


void HR_model::getJacobian(double t, double *x, double **J)
{
    double b = this->b;
    double mu = this->mu;
    double s = this->s;
    
J[0][0] = -3.0*x[0]*x[0]+2.0*b*x[0];
J[0][1] = 1;
J[0][2] = -1;
J[1][0] = -10.0*x[0];
J[1][1] = -1;
J[1][2] = 0;
J[2][0] = mu*s;
J[2][1] = 0;
J[2][2] = -mu;
}

double HR_model::dfdi(double t, double *x)
{
    return 1;
}

HR_model::~HR_model() {
}

