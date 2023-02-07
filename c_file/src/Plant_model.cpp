/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Plant_model.cpp
 * Author: teo
 * 
 * Created on July 12, 2017, 10:33 AM
 */

#include "../inc/Plant_model.hpp"

Plant_model::Plant_model() 
{       
    this->nx = 5;
    this->g_I = 0;
    this->g_K = 0;
    this->g_T = 0;
    this->g_L = 0;
    this->g_KCa = 0;
    this->V_I = 0;
    this->V_K = 0;
    this->V_L = 0;
    this->V_Ca = 0;
    this->K_c = 0;
    this->rho = 0;
    this->Iapp = 0;
}

Plant_model::Plant_model(const Plant_model &h)
{
    this->nx = h.nx;
    this->g_I = h.g_I;
    this->g_K = h.g_K;
    this->g_T = h.g_T;
    this->g_L = h.g_L;
    this->g_KCa = h.g_KCa;
    this->V_I = h.V_I;
    this->V_K = h.V_K;
    this->V_L = h.V_L;
    this->V_Ca = h.V_Ca;
    this->K_c = h.K_c;
    this->rho = h.rho;
    this->Iapp = h.Iapp;
}

Plant_model::Plant_model(double g_I, double g_K, double g_T, double g_L, double g_KCa, double V_I, double V_K, double V_L, double V_Ca, double K_c, double rho, double Iapp)
{
    this->nx = 5;
    this->g_I = g_I;
    this->g_K = g_K;
    this->g_T = g_T;
    this->g_L = g_L;
    this->g_KCa = g_KCa;
    this->V_I = V_I;
    this->V_K = V_K;
    this->V_L = V_L;
    this->V_Ca = V_Ca;
    this->K_c = K_c;
    this->rho = rho;
    this->Iapp = Iapp;
}

Plant_model *Plant_model::clone() const
{
    return new Plant_model(*this);
}


void Plant_model::getXdot(double t, double *x, double *xdot,double *Isyn)
{
    double g_I = this->g_I;
    double g_K = this->g_K;
    double g_T= this->g_T;
    double g_L= this->g_L;
    double g_KCa= this->g_KCa;
    double V_I= this->V_I;
    double V_K= this->V_K;
    double V_L= this->V_L;
    double V_Ca= this->V_Ca;
    double K_c= this->K_c;
    double rho= this->rho;
    double Iapp= this->Iapp;
            
    double V_s = (127*x[0]/105 + 8265/105);
    double a_m = 0.1*(50 - V_s)/(exp((50 - V_s)/10) - 1);
    double b_m = 4*exp((25 - V_s)/18);
    double a_h = 0.07*exp((25 - V_s)/20);
    double b_h = 1.0/(1 + exp((55 - V_s)/10));
    double a_n = 0.01*(55 - V_s)/(exp((55 - V_s)/10) - 1);
    double b_n = 0.125*exp((45 - V_s)/80); 
    double m_inf = a_m/(a_m + b_m);
    double h_inf = a_h/(a_h + b_h);
    double n_inf = a_n/(a_n + b_n);
    double chi_inf = 1/(exp(0.15*(-x[0]-50)) + 1);
    double tau_h = 12.5/(a_h + b_h);
    double tau_n = 12.5/(a_n + b_n);
    double tau_chi = 235; 

    double I_Na = g_I*m_inf*m_inf*m_inf*x[2]*(V_I - x[0]);
    double I_K = g_K*x[3]*x[3]*x[3]*x[3]*(V_K - x[0]);
    double I_T = g_T*x[4]*(V_I-x[0]);
    double I_KCa = g_KCa*x[1]/(.5 + x[1])*(V_K - x[0]);
    double I_L = g_L*(V_L - x[0]);

    xdot[0] = I_Na+I_K+I_T+I_KCa+I_L+Iapp+Isyn[0];
    xdot[1] = rho*(K_c*x[4]*(V_Ca - x[0]) - x[1]);
    xdot[2] = (h_inf - x[2])/tau_h;
    xdot[3] = (n_inf - x[3])/tau_n;
    xdot[4] = (chi_inf - x[4])/tau_chi;

}


void Plant_model::getJacobian(double t, double *x, double **J)
{
    double g_I = this->g_I;
    double g_K = this->g_K;
    double g_T= this->g_T;
    double g_L= this->g_L;
    double g_KCa= this->g_KCa;
    double V_I= this->V_I;
    double V_K= this->V_K;
    double V_L= this->V_L;
    double V_Ca= this->V_Ca;
    double K_c= this->K_c;
    double rho= this->rho;
    double Iapp= this->Iapp;


double V = x[0];
double Ca = x[1];
double h = x[2];
double n = x[3];
double chi = x[4];


J[0][0] = 0;
J[0][1] = 0;
J[0][2] = 0;
J[0][3] = 0;
J[0][4] = 0;
J[1][1] = 0;
J[1][2] = 0;
J[1][3] = 0;
J[1][4] = 0;
J[2][1] = 0;
J[2][2] = 0;
J[2][3] = 0;
J[2][4] = 0;
J[3][1] = 0;
J[3][2] = 0;
J[3][3] = 0;
J[3][4] = 0;
J[4][1] = 0;
J[4][2] = 0;
J[4][3] = 0;
J[4][4] = 0;
}

double Plant_model::dfdi(double t, double *x)
{
    return 1;
}

Plant_model::~Plant_model() {
}

