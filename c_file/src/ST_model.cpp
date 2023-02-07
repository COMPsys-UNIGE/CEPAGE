

#include "../inc/ST_model.hpp"

ST_model::ST_model() 
{
    this->nx = 8;
    this->g_Na = 0;
    this->g_K = 0;
    this->g_L = 0;
    this->g_A = 0;
    this->g_T = 0;
    this->g_KCa = 0;
    this->g_HVA = 0;
    this->epsilon = 0;
    this->E_Na = 0;
    this->E_K = 0;
    this->E_L = 0;
    this->E_Ca = 0;
    this->C = 0;
    this->A = 0;
    this->y0 = 0;
    this->V_c = 0;
    this->w = 0;
    this->k_Ca = 0;
    this->alpha = 0;
    this->k = 0;
    this->nu_m = 0;
    this->nu_h = 0;
    this->nu_n = 0;
    this->nu_nA = 0;
    this->nu_hA = 0;
    this->nu_mT = 0;
    this->nu_hT = 0;
    this->nu_mHVA = 0;
    this->s_m = 0;
    this->s_h = 0;
    this->s_n = 0;
    this->s_nA = 0;
    this->s_hA = 0;
    this->s_mT = 0;
    this->s_hT = 0;
    this->s_mHVA = 0;
    this->tau_nA = 0;
    this->tau_hA = 0;
    this->tau_hT = 0;
    this->tau_mHVA = 0;
    this->I_app = 0;
}

ST_model::ST_model(const ST_model &h)
{
    this->nx = h.nx;
    this->g_Na = h.g_Na;
    this->g_K = h.g_K;
    this->g_L = h.g_L;
    this->g_A = h.g_A;
    this->g_T = h.g_T;
    this->g_KCa = h.g_KCa;
    this->g_HVA = h.g_HVA;
    this->epsilon = h.epsilon;
    this->E_Na = h.E_Na;
    this->E_K = h.E_K;
    this->E_L = h.E_L;
    this->E_Ca = h.E_Ca;
    this->C = h.C;
    this->A = h.A;
    this->y0 = h.y0;
    this->V_c = h.V_c;
    this->w = h.w;
    this->k_Ca = h.k_Ca;
    this->alpha = h.alpha;
    this->k = h.k;
    this->nu_m = h.nu_m;
    this->nu_h = h.nu_h;
    this->nu_n = h.nu_n;
    this->nu_nA = h.nu_nA;
    this->nu_hA = h.nu_hA;
    this->nu_mT = h.nu_mT;
    this->nu_hT = h.nu_hT;
    this->nu_mHVA = h.nu_mHVA;
    this->s_m = h.s_m;
    this->s_h = h.s_h;
    this->s_n = h.s_n;
    this->s_nA = h.s_nA;
    this->s_hA = h.s_hA;
    this->s_mT = h.s_mT;
    this->s_hT = h.s_hT;
    this->s_mHVA = h.s_mHVA;
    this->tau_nA = h.tau_nA;
    this->tau_hA = h.tau_hA;
    this->tau_hT = h.tau_hT;
    this->tau_mHVA = h.tau_mHVA;
    this->I_app = h.I_app;
}

ST_model::ST_model(double g_Na,double g_K,double g_L,double g_A,double g_T,double g_KCa,double g_HVA,double epsilon,double E_Na,double E_K,double E_L,double E_Ca,double C,double A,double y0,double V_c,double w,double k_Ca,double alpha,double k,double nu_m,double nu_h,double nu_n,double nu_nA,double nu_hA,double nu_mT,double nu_hT,double nu_mHVA,double s_m,double s_h,double s_n,double s_nA,double s_hA,double s_mT,double s_hT,double s_mHVA,double tau_nA,double tau_hA,double tau_hT,double tau_mHVA,double I_app)
{
    this->nx = 8;
    this->g_Na = g_Na;
    this->g_K = g_K;
    this->g_L = g_L;
    this->g_A = g_A;
    this->g_T = g_T;
    this->g_KCa = g_KCa;
    this->g_HVA = g_HVA;
    this->epsilon = epsilon;
    this->E_Na = E_Na;
    this->E_K = E_K;
    this->E_L = E_L;
    this->E_Ca = E_Ca;
    this->C = C;
    this->A = A;
    this->y0 = y0;
    this->V_c = V_c;
    this->w = w;
    this->k_Ca = k_Ca;
    this->alpha = alpha;
    this->k = k;
    this->nu_m = nu_m;
    this->nu_h = nu_h;
    this->nu_n = nu_n;
    this->nu_nA = nu_nA;
    this->nu_hA = nu_hA;
    this->nu_mT = nu_mT;
    this->nu_hT = nu_hT;
    this->nu_mHVA = nu_mHVA;
    this->s_m = s_m;
    this->s_h = s_h;
    this->s_n = s_n;
    this->s_nA = s_nA;
    this->s_hA = s_hA;
    this->s_mT = s_mT;
    this->s_hT = s_hT;
    this->s_mHVA = s_mHVA;
    this->tau_nA = tau_nA;
    this->tau_hA = tau_hA;
    this->tau_hT = tau_hT;
    this->tau_mHVA = tau_mHVA;
    this->I_app = I_app;
}

ST_model *ST_model::clone() const
{
    return new ST_model(*this);
}


void ST_model::getXdot(double t, double *x, double *xdot,double *Isyn)
{
    double g_Na = this->g_Na;
    double g_K = this->g_K;
    double g_L = this->g_L;
    double g_A = this->g_A;
    double g_T = this->g_T;
    double g_KCa = this->g_KCa;
    double g_HVA = this->g_HVA;
    double epsilon = this->epsilon;
    double E_Na = this->E_Na;
    double E_K = this->E_K;
    double E_L = this->E_L;
    double E_Ca = this->E_Ca;
    double C = this->C;
    double A = this->A;
    double y0 = this->y0;
    double V_c = this->V_c;
    double w = this->w;
    double k_Ca = this->k_Ca;
    double alpha = this->alpha;
    double k = this->k;
    double nu_m = this->nu_m;
    double nu_h = this->nu_h;
    double nu_n = this->nu_n;
    double nu_nA = this->nu_nA;
    double nu_hA = this->nu_hA;
    double nu_mT = this->nu_mT;
    double nu_hT = this->nu_hT;
    double nu_mHVA = this->nu_mHVA;
    double s_m = this->s_m;
    double s_h = this->s_h;
    double s_n = this->s_n;
    double s_nA = this->s_nA;
    double s_hA = this->s_hA;
    double s_mT = this->s_mT;
    double s_hT = this->s_hT;
    double s_mHVA = this->s_mHVA;
    double tau_nA = this->tau_nA;
    double tau_hA = this->tau_hA;
    double tau_hT = this->tau_hT;
    double tau_mHVA = this->tau_mHVA;
    double I_app = this->I_app;
    


double m_inf = 1/(1+exp(-(x[0]-nu_m)/s_m));
double h_inf = 1/(1+exp(-(x[0]-nu_h)/s_h));
double n_inf = 1/(1+exp(-(x[0]-nu_n)/s_n));
double n_Ainf = 1/(1+exp(-(x[0]-nu_nA)/s_nA));
double h_Ainf = 1/(1+exp(-(x[0]-nu_hA)/s_hA));
double m_Tinf = 1/(1+exp(-(x[0]-nu_mT)/s_mT));
double h_Tinf = 1/(1+exp(-(x[0]-nu_hT)/s_hT));
double m_HVAinf = 1/(1+exp(-(x[0]-nu_mHVA)/s_mHVA));

double tau_h = y0+2*A*w/(4*3.1415*(x[0]-V_c)*(x[0]-V_c)+w*w);
double tau_n = 6/(1+exp((x[0]+23)/15));

double I_Na = g_Na*m_inf*m_inf*m_inf*x[1]*(x[0]-E_Na);
double I_K = g_K*x[2]*x[2]*x[2]*x[2]*(x[0]-E_K);
double I_L = g_L*(x[0]-E_L);
double I_A = g_A*x[3]*x[4]*(x[0]-E_K);
double I_T = g_T*m_Tinf*x[5]*(x[0]-E_Ca);
double I_HVA = g_HVA*x[6]*(x[0]-E_Ca);
double I_KCa = g_KCa*(x[7]*x[7]*x[7]*x[7]*x[7]/(k_Ca*k_Ca*k_Ca*k_Ca*k_Ca+x[7]*x[7]*x[7]*x[7]*x[7]))*(x[0]-E_K);

xdot[0] = (1/C)*(I_app-I_Na-I_K-I_L-I_A-I_T-I_KCa-I_HVA+Isyn[0]);
xdot[1] = (h_inf-x[1])/tau_h;
xdot[2] = (n_inf-x[2])/tau_n;
xdot[3] = (n_Ainf-x[3])/tau_nA;
xdot[4] = (h_Ainf-x[4])/tau_hA;
xdot[5] = (h_Tinf-x[5])/tau_hT;
xdot[6] = (m_HVAinf-x[6])/tau_mHVA;
xdot[7] = -epsilon*(k*alpha*(I_T+I_HVA)+k_Ca);

}


void ST_model::getJacobian(double t, double *x, double **J)
{
    double g_Na = this->g_Na;
    double g_K = this->g_K;
    double g_L = this->g_L;
    double g_A = this->g_A;
    double g_T = this->g_T;
    double g_KCa = this->g_KCa;
    double g_HVA = this->g_HVA;
    double epsilon = this->epsilon;
    double E_Na = this->E_Na;
    double E_K = this->E_K;
    double E_L = this->E_L;
    double E_Ca = this->E_Ca;
    double C = this->C;
    double A = this->A;
    double y0 = this->y0;
    double V_c = this->V_c;
    double w = this->w;
    double k_Ca = this->k_Ca;
    double alpha = this->alpha;
    double k = this->k;
    double nu_m = this->nu_m;
    double nu_h = this->nu_h;
    double nu_n = this->nu_n;
    double nu_nA = this->nu_nA;
    double nu_hA = this->nu_hA;
    double nu_mT = this->nu_mT;
    double nu_hT = this->nu_hT;
    double nu_mHVA = this->nu_mHVA;
    double s_m = this->s_m;
    double s_h = this->s_h;
    double s_n = this->s_n;
    double s_nA = this->s_nA;
    double s_hA = this->s_hA;
    double s_mT = this->s_mT;
    double s_hT = this->s_hT;
    double s_mHVA = this->s_mHVA;
    double tau_nA = this->tau_nA;
    double tau_hA = this->tau_hA;
    double tau_hT = this->tau_hT;
    double tau_mHVA = this->tau_mHVA;
    double I_app = this->I_app;


double V = x[0];
double h = x[1];
double n = x[2];
double n_A = x[3];
double h_A = x[4];
double h_T = x[5];
double m_HVA = x[6];
double Ca = x[7];

J[0][0] = 0;
J[0][1] = 0;
J[0][2] = 0;
J[0][3] = 0;
J[0][4] = 0;
J[0][5] = 0;
J[0][6] = 0;
J[0][7] = 0;
J[1][0] = 0;
J[1][1] = 0;
J[1][2] = 0;
J[1][3] = 0;
J[1][4] = 0;
J[1][5] = 0;
J[1][6] = 0;
J[1][7] = 0;
J[2][0] = 0;
J[2][1] = 0;
J[2][2] = 0;
J[2][3] = 0;
J[2][4] = 0;
J[2][5] = 0;
J[2][6] = 0;
J[2][7] = 0;
J[3][0] = 0;
J[3][1] = 0;
J[3][2] = 0;
J[3][3] = 0;
J[3][4] = 0;
J[3][5] = 0;
J[3][6] = 0;
J[3][7] = 0;
J[4][0] = 0;
J[4][1] = 0;
J[4][2] = 0;
J[4][3] = 0;
J[4][4] = 0;
J[4][5] = 0;
J[4][6] = 0;
J[4][7] = 0;
J[5][0] = 0;
J[5][1] = 0;
J[5][2] = 0;
J[5][3] = 0;
J[5][4] = 0;
J[5][5] = 0;
J[5][6] = 0;
J[5][0] = 0;
J[5][7] = 0;
J[6][1] = 0;
J[6][2] = 0;
J[6][3] = 0;
J[6][4] = 0;
J[6][5] = 0;
J[6][6] = 0;
J[6][7] = 0;
J[7][1] = 0;
J[7][2] = 0;
J[7][3] = 0;
J[7][4] = 0;
J[7][5] = 0;
J[7][6] = 0;
J[7][7] = 0;
}

double ST_model::dfdi(double t, double *x)
{
    return this->C;
}

ST_model::~ST_model() {
}

