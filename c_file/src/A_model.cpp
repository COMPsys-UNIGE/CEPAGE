

#include "../inc/A_model.hpp"

A_model::A_model() 
{
    this->nx = 6;
    this->Ca_shift = 0; 
    this->x_shift = 0; 
    this->gh = 0;
    this->Vhh = 0;
    this->Iext = 0;
}

A_model::A_model(const A_model &h)
{
    this->nx = h.nx;
    this->Ca_shift = h.Ca_shift; 
    this->x_shift = h.x_shift; 
    this->gh = h.gh;
    this->Vhh = h.Vhh;
    this->Iext = h.Iext;
}

A_model::A_model(double Ca_shift, double x_shift, double gh, double Vhh, double Iext)
{
    this->nx = 6;
    this->Ca_shift = Ca_shift; 
    this->x_shift = x_shift; 
    this->gh = gh;
    this->Vhh = Vhh;
    this->Iext = Iext;
}

A_model *A_model::clone() const
{
    return new A_model(*this);
}


void A_model::getXdot(double t, double *x, double *xdot,double *Isyn)
{
double Ca_shift = this->Ca_shift; 
double x_shift = this->x_shift;
double gh = this->gh;
double Vhh = this->Vhh;
double Iext = this->Iext;

xdot[0] = 8*((0.1*(50-(127*x[0]/105+8265/105))/(exp((50 - (127*x[0]/105+8265/105))/10) - 1))/((0.1*(50 - (127*x[0]/105 + 8265/105))/(exp((50 - (127*x[0]/105 + 8265/105))/10) - 1))+(4*exp((25 - (127*x[0]/105 + 8265/105))/18))))*((0.1*(50-(127*x[0]/105+8265/105))/(exp((50 - (127*x[0]/105+8265/105))/10) - 1))/((0.1*(50 - (127*x[0]/105 + 8265/105))/(exp((50 - (127*x[0]/105 + 8265/105))/10) - 1))+(4*exp((25 - (127*x[0]/105 + 8265/105))/18))))*((0.1*(50-(127*x[0]/105+8265/105))/(exp((50 - (127*x[0]/105+8265/105))/10) - 1))/((0.1*(50 - (127*x[0]/105 + 8265/105))/(exp((50 - (127*x[0]/105 + 8265/105))/10) - 1))+(4*exp((25 - (127*x[0]/105 + 8265/105))/18))))*x[1]*(30 - x[0]) + 1.3*x[2]*x[2]*x[2]*x[2]*(-75-x[0])+0.01*x[3]*(30-x[0])+0.03*x[5]/(0.5 + x[5])*(-75 - x[0])+0.003*(-40 - x[0])+gh*((1/(1+exp(-(x[0]+63)/7.8)))*(1/(1+exp(-(x[0]+63)/7.8)))*(1/(1+exp(-(x[0]+63)/7.8))))*x[4]*(120-x[0])+Isyn[0]+Iext;

xdot[1] = ((1-x[1])*(0.07*exp((25 - (127*x[0]/105 + 8265/105))/20))-x[1]*(1.0/(1 + exp((55 - (127*x[0]/105 + 8265/105))/10))))/12.5;

xdot[2] = ((1-x[2])*(0.01*(55 - (127*x[0]/105 + 8265/105))/(exp((55 - (127*x[0]/105 + 8265/105))/10) - 1))-x[2]*(0.125*exp((45 - (127*x[0]/105 + 8265/105))/80)))/12.5;

xdot[3] = ((1/(exp(0.15*(-x[0]-50+x_shift))+1))-x[3])/235;

xdot[4] = 0.5*((1/(1+exp(10*(x[0]-Vhh))))-x[4])/(7.1+10.4/(1+exp((x[0]+68)/2.2)));

xdot[5] = 0.0001*(0.0085*x[3]*(140-x[0]+Ca_shift)-x[5]);
}


void A_model::getJacobian(double t, double *x, double **J)
{
double Ca_shift = this->Ca_shift; 
double x_shift = this->x_shift;
double gh = this->gh;
double Vhh = this->Vhh;
double Iext = this->Iext;


double V = x[0];
double h = x[1];
double n = x[2];
double chi = x[3];
double ksi = x[4];
double Ca = x[5];

J[0][0] = 0;
J[0][1] = 0;
J[0][2] = 0;
J[0][3] = 0;
J[0][4] = 0;
J[0][5] = 0;
J[1][0] = 0;
J[1][1] = 0;
J[1][2] = 0;
J[1][3] = 0;
J[1][4] = 0;
J[1][5] = 0;
J[2][0] = 0;
J[2][1] = 0;
J[2][2] = 0;
J[2][3] = 0;
J[2][4] = 0;
J[2][5] = 0;
J[3][0] = 0;
J[3][1] = 0;
J[3][2] = 0;
J[3][3] = 0;
J[3][4] = 0;
J[3][5] = 0;
J[4][0] = 0;
J[4][1] = 0;
J[4][2] = 0;
J[4][3] = 0;
J[4][4] = 0;
J[4][5] = 0;
J[5][0] = 0;
J[5][1] = 0;
J[5][2] = 0;
J[5][3] = 0;
J[5][4] = 0;
J[5][5] = 0;
}

double A_model::dfdi(double t, double *x)
{
    return 0;
}

A_model::~A_model() {
}

