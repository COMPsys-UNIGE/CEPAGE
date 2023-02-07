/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   RE_model.cpp
 * Author: teo
 * 
 * Created on July 12, 2017, 10:33 AM
 */

#include "../inc/RE_model.hpp"

RE_model::RE_model() 
{
    this->nx = 7;
    this->gna = 0; 
    this->ENa = 0; 
    this->gk = 0;
    this->Ek = 0;
    this->gca = 0;
    this->gl = 0;
    this->El = 0; 
    this->C = 0; 
    this->Iapp = 0;
    this->Ca0 = 0; 
    this->d = 0;
    this->KT = 0;
    this->Kd = 0;
    this->gd = 0;
    this->D = 0;
    this->EsynEx = 0;
    this->tau=0;
    this->tTm1 = 0;
    this->tTm2 = 0;
    this->tTh1 = 0;
    this->tTh2 = 0;
}

RE_model::RE_model(const RE_model &h)
{
    this->nx = h.nx;
    this->gna = h.gna; 
    this->ENa = h.ENa; 
    this->gk = h.gk;
    this->Ek = h.Ek;
    this->gca = h.gca;
    this->gl = h.gl;
    this->El = h.El; 
    this->C = h.C; 
    this->Iapp = h.Iapp;
    this->Ca0 = h.Ca0; 
    this->d = h.d;
    this->KT = h.KT;
    this->Kd = h.Kd;
    this->gd = h.gd;
    this->D = h.D;
    this->EsynEx = h.EsynEx;
    this->tau=h.tau;
    this->tTm1 = h.tTm1;
    this->tTm2 = h.tTm2;
    this->tTh1 = h.tTh1;
    this->tTh2 = h.tTh2;
}

RE_model::RE_model(double gna, double ENa, double gk, double Ek, double gca, double gl, double El, double C, double Iapp, double Ca0, double d, double KT, double Kd, double gd, double D, double EsynEx, double tau,double tTm1,double tTm2,double tTh1,double tTh2)
{
    this->nx = 7;
    this->gna = gna; 
    this->ENa = ENa; 
    this->gk = gk;
    this->Ek = Ek;
    this->gca = gca;
    this->gl = gl;
    this->El = El; 
    this->C = C; 
    this->Iapp = Iapp;
    this->Ca0 = Ca0; 
    this->d = d;
    this->KT = KT;
    this->Kd = Kd;
    this->gd = gd;
    this->D = D;
    this->EsynEx = EsynEx;
    this->tau = tau;
    this->tTm1 = tTm1;
    this->tTm2 = tTm2;
    this->tTh1 = tTh1;
    this->tTh2 = tTh2;
}

RE_model *RE_model::clone() const
{
    return new RE_model(*this);
}


void RE_model::getXdot(double t, double *x, double *xdot,double *Isyn)
{
double gna = this->gna; 
double ENa = this->ENa;
double gk = this->gk;
double Ek = this->Ek;
double gca = this->gca;
double gl = this->gl;
double El = this->El;
double C = this->C;
double Iapp = this->Iapp;
double Ca0 = this->Ca0;
double d = this->d;
double KT = this->KT;
double Kd = this->Kd;
double gd = this->gd;
double D = this->D;
double EsynEx = this->EsynEx;
double tau = this->tau;
double tTm1=this->tTm1;
double tTm2=this->tTm2;
double tTh1=this->tTh1;
double tTh2=this->tTh2;

double R = 8.31441; 
double T = 309.15;  
double F = 96469;  
double k0 = 1000;
double k = 0.1;

double mTinf = 1/(1 + exp(-(x[0]+52)/7.4));
double tTm = 0.44+0.15/(tTm1*exp((x[0]+27)/10)+tTm2*exp(-(x[0]+102)/15));
double hTinf = 1/(1+exp((x[0]+80)/5));
double tTh = 62.7+0.27/(tTh1*exp((x[0]+48)/4)+tTh2*exp(-(x[0]+407)/50));
double ECa = k0*(R*T/(2*F))*log(Ca0/x[6]);

double Il = gl*(x[0]-El);
double INa = gna*x[2]*x[2]*x[2]*x[1]*(x[0]-ENa);
double IK = gk*x[3]*x[3]*x[3]*x[3]*(x[0]-Ek);
double IT = gca*x[4]*x[4]*x[5]*(x[0]-ECa);

xdot[0] = ((-Iapp-IT-Il-INa-IK+Isyn[0]+gd*D*(EsynEx-x[0]))/C)/tau;
xdot[1] = (0.128*exp((17-x[0])/18)*(1-x[1])-(4/(exp(-0.2*(x[0]-40))+1))*x[1])/tau;
xdot[2] = (((0.32*(13-x[0]))/(exp(0.25*(13-x[0]))-1))*(1-x[2])-((0.28*(x[0]-40))/(exp(0.2*(x[0]-40))-1))*x[2])/tau;
xdot[3] = (((0.032*(15-x[0]))/(exp(0.2*(15-x[0]))-1))*(1-x[3])-0.5*exp((10-x[0])/40)*x[3])/tau;
xdot[4] = (-(x[4]-mTinf)/(tTm))/tau;
xdot[5] = (-(x[5]-hTinf)/(tTh))/tau;
xdot[6] = (-(k*IT)/(2*F*d)-(KT*x[6])/(x[6]+Kd))/tau;

}


void RE_model::getJacobian(double t, double *x, double **J)
{
double gna = this->gna; 
double ENa = this->ENa;
double gk = this->gk;
double Ek = this->Ek;
double gca = this->gca;
double gl = this->gl;
double El = this->El;
double C = this->C;
double Iapp = this->Iapp;
double Ca0 = this->Ca0;
double d = this->d;
double KT = this->KT;
double Kd = this->Kd;
double tau = this->tau;
double tTm1=this->tTm1;
double tTm2=this->tTm2;
double tTh1=this->tTh1;
double tTh2=this->tTh2;


double V = x[0];
double h = x[1];
double m = x[2];
double n = x[3];
double mT = x[4];
double hT = x[5];
double Ca = x[6];

J[0][0] = 0;
J[0][1] = 0;
J[0][2] = 0;
J[0][3] = 0;
J[0][4] = 0;
J[0][5] = 0;
J[0][6] = 0;
J[1][0] = 0;
J[1][1] = 0;
J[1][2] = 0;
J[1][3] = 0;
J[1][4] = 0;
J[1][5] = 0;
J[1][6] = 0;
J[2][0] = 0;
J[2][1] = 0;
J[2][2] = 0;
J[2][3] = 0;
J[2][4] = 0;
J[2][5] = 0;
J[2][6] = 0;
J[3][0] = 0;
J[3][1] = 0;
J[3][2] = 0;
J[3][3] = 0;
J[3][4] = 0;
J[3][5] = 0;
J[3][6] = 0;
J[4][0] = 0;
J[4][1] = 0;
J[4][2] = 0;
J[4][3] = 0;
J[4][4] = 0;
J[4][5] = 0;
J[4][6] = 0;
J[5][0] = 0;
J[5][1] = 0;
J[5][2] = 0;
J[5][3] = 0;
J[5][4] = 0;
J[5][5] = 0;
J[5][6] = 0;
J[5][0] = 0;
J[6][1] = 0;
J[6][2] = 0;
J[6][3] = 0;
J[6][4] = 0;
J[6][5] = 0;
J[6][6] = 0;
}

double RE_model::dfdi(double t, double *x)
{
    return this->C;
}

RE_model::~RE_model() {
}

