function x_dot = getXdot(object,t,x,varargin)


% synapsis input
Isyn= 0;
if nargin >= 4
    if isscalar(varargin{end})
        Isyn = varargin{end};
    end
end


x_dot = zeros(8,1);

% Define param
g_Na = object.g_Na;
g_K = object.g_K;
g_L = object.g_L;
g_A = object.g_A;
g_T = object.g_T;
g_KCa = object.g_KCa;
g_HVA = object.g_HVA;
epsilon = object.epsilon;
E_Na = object.E_Na;
E_K = object.E_K;
E_L = object.E_L;
E_Ca = object.E_Ca;
C = object.C;
A = object.A;
y0 = object.y0;
V_c = object.V_c;
w = object.w;
k_Ca = object.k_Ca;
alpha = object.alpha;
k = object.k;
nu_m = object.nu_m;
nu_h = object.nu_h;
nu_n = object.nu_n;
nu_nA = object.nu_nA;
nu_hA = object.nu_hA;
nu_mT = object.nu_mT;
nu_hT = object.nu_hT;
nu_mHVA = object.nu_mHVA;
s_m = object.s_m;
s_h = object.s_h;
s_n = object.s_n;
s_nA = object.s_nA;
s_hA = object.s_hA;
s_mT = object.s_mT;
s_hT = object.s_hT;
s_mHVA = object.s_mHVA;
tau_nA = object.tau_nA;
tau_hA = object.tau_hA;
tau_hT = object.tau_hT;
tau_mHVA = object.tau_mHVA;
I_app = object.I_app;



% equation
V = x(1);
h = x(2);
n = x(3);
n_A = x(4);
h_A = x(5);
h_T = x(6);
m_HVA = x(7);
Ca = x(8);

% compute help variable 

m_inf = 1/(1+exp(-(V-nu_m)/s_m));
h_inf = 1/(1+exp(-(V-nu_h)/s_h));
n_inf = 1/(1+exp(-(V-nu_n)/s_n));
n_Ainf = 1/(1+exp(-(V-nu_nA)/s_nA));
h_Ainf = 1/(1+exp(-(V-nu_hA)/s_hA));
m_Tinf = 1/(1+exp(-(V-nu_mT)/s_mT));
h_Tinf = 1/(1+exp(-(V-nu_hT)/s_hT));
m_HVAinf = 1/(1+exp(-(V-nu_mHVA)/s_mHVA));

tau_h = y0+2*A*w/(4*3.1415*(V-V_c)*(V-V_c)+w*w);
tau_n = 6/(1+exp((V+23)/15));

I_Na = g_Na*m_inf*m_inf*m_inf*h*(V-E_Na);
I_K = g_K*n*n*n*n*(V-E_K);
I_L = g_L*(V-E_L);
I_A = g_A*n_A*h_A*(V-E_K);
I_T = g_T*m_Tinf*h_T*(V-E_Ca);
I_HVA = g_HVA*m_HVA*(V-E_Ca);
I_KCa = g_KCa*(Ca*Ca*Ca*Ca*Ca/(k_Ca*k_Ca*k_Ca*k_Ca*k_Ca+Ca*Ca*Ca*Ca*Ca))*(V-E_K);


V_dot = (1/C)*(I_app-I_Na-I_K-I_L-I_A-I_T-I_KCa-I_HVA+Isyn);
h_dot = (h_inf-h)/tau_h;
n_dot = (n_inf-n)/tau_n;
n_A_dot = (n_Ainf-n_A)/tau_nA;
h_A_dot = (h_Ainf-h_A)/tau_hA;
h_T_dot = (h_Tinf-h_T)/tau_hT;
m_HVA_dot = (m_HVAinf-m_HVA)/tau_mHVA;
Ca_dot = -epsilon*(k*alpha*(I_T+I_HVA)+k_Ca);


x_dot = [V_dot; h_dot; n_dot; n_A_dot; h_A_dot; h_T_dot; m_HVA_dot; Ca_dot];
end