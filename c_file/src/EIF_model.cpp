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

#include "../inc/EIF_model.hpp"

EIF_model::EIF_model() 
{
    this->nx = 2;
    this->a = 0;
    this->b = 0;
    this->C = 0;
    this->VT = 0;
    this->DT = 0;
    this->I = 0;
    this->gL = 0;
    this->EL = 0;
    this->Vr = 0;
    this->tw = 0;
    this->gexp = 0;
    this->gd = 0;
    this->D = 0;
    this->EsynEx = 0;
}

EIF_model::EIF_model(double C,double a,double b,double VT,double DT,double I,double gL,double EL,double Vr,double tw, double gexp,double gd,double D,double EsynEx) 
{
    this->nx = 2;
    this->a = a;
    this->b = b;
    this->C = C;
    this->VT = VT;
    this->DT = DT;
    this->I = I;
    this->gL = gL;
    this->EL = EL;
    this->Vr = Vr;
    this->tw = tw;
    this->gexp = gexp;
    this->gd = gd;
    this->D = D;
    this->EsynEx = EsynEx;
}

EIF_model::EIF_model(const EIF_model &e)
{
this->nx = e.nx;
    this->a = e.a;
    this->b = e.b;
    this->C = e.C;
    this->VT = e.VT;
    this->DT = e.DT;
    this->I = e.I;
    this->gL = e.gL;
    this->EL = e.EL;
    this->Vr = e.Vr;
    this->tw = e.tw;
    this->gexp = e.gexp;
    this->gd = e.gd;
    this->D = e.D;
    this->EsynEx = e.EsynEx;
}

EIF_model *EIF_model::clone() const
{
    return new EIF_model(*this);
}

void EIF_model::getXdot(double t, double *x, double *xdot,double *Iext)
{
    double a = this->a;
    double b = this->b;
    double C = this->C;
    double VT = this->VT;
    double DT = this->DT;
    double I = this->I;
    double gL = this->gL;
    double EL = this->EL;
    double Vr = this->Vr;
    double tw = this->tw;
    double gexp = this->gexp;
    double gd = this->gd;
    double D = this->D;
    double EsynEx = this->EsynEx;
    
    double V = x[0];
    double u = x[1];
            	
	
	
    xdot[0] = 1/C*(-gL*(V-EL) + gexp * exp((V-VT)/DT) - u + I + Iext[0] + gd*D*(EsynEx-V));	
    xdot[1] = 1/tw*(a*(V-EL) - u);
    
}

void EIF_model::getJacobian(double t, double *x, double **J)
{
}

double EIF_model::dfdi(double t, double *x)
{
}


bool EIF_model::getResetConditions(double *x)
{
    if(x[0] > 20)
        return true;
    else
        return false;
}

void EIF_model::resetStates(double *x)
{
    if(this->getResetConditions(x))
    {
        x[0] = this->Vr;
        x[1] += this->b;
    }
    return;
}

EIF_model::~EIF_model() 
{
}

