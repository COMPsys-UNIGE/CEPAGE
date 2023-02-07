function x_dot = getXdot(object,t,x,varargin)
% getXdot    Computes the derivative of the state 
%
%  x_dot = getXdot(object,t,x)
%   compute the time derivative of the model at time instant t, in state x
%   and with I_{ext} = 0;
%
%  x_dot = getXdot(object,t,x,I_{ext})
%   compute the time derivative of the model at time instant t, in state x
%   and with I_{ext} specified by the user;


% synapsis input
Isyn= 0;
if nargin >= 4
    if isscalar(varargin{end})
        Isyn = varargin{end};
    end
end


x_dot = zeros(5,1);

% Define param
g_I = object.g_I;
g_K = object.g_K;
g_T = object.g_T;
g_L = object.g_L;
g_KCa = object.g_KCa;
V_I = object.V_I;
V_K = object.V_K;
V_L = object.V_L;
V_Ca = object.V_Ca;
K_c = object.K_c;
rho = object.rho;
Iapp = object.Iapp;

% equation
V = x(1);
Ca = x(2);
h = x(3);
n = x(4);
chi = x(5);


% compute help variable 
V_s = (127*V/105 + 8265/105);
a_m = 0.1*(50 - V_s)/(exp((50 - V_s)/10) - 1);
b_m = 4*exp((25 - V_s)/18);
a_h = 0.07*exp((25 - V_s)/20);
b_h = 1.0/(1 + exp((55 - V_s)/10));
a_n = 0.01*(55 - V_s)/(exp((55 - V_s)/10) - 1);
b_n = 0.125*exp((45 - V_s)/80); 
m_inf = a_m/(a_m + b_m);
h_inf = a_h/(a_h + b_h);
n_inf = a_n/(a_n + b_n);
chi_inf = 1/(exp(0.15*(-V-50)) + 1);
tau_h = 12.5/(a_h + b_h);
tau_n = 12.5/(a_n + b_n);
tau_chi = 235; 

I_Na = g_I*m_inf*m_inf*m_inf*h*(V_I - V);
I_K = g_K*n*n*n*n*(V_K - V);
I_T = g_T*chi*(V_I-V);
I_KCa = g_KCa*Ca/(.5 + Ca)*(V_K - V);
I_L = g_L*(V_L - V);

V_dot = I_Na+I_K+I_T+I_KCa+I_L+Iapp+Isyn;
Ca_dot = rho*(K_c*chi*(V_Ca - V) - Ca);
h_dot = (h_inf - h)/tau_h;
n_dot = (n_inf - n)/tau_n;
chi_dot = (chi_inf - chi)/tau_chi;


x_dot = [V_dot; Ca_dot; h_dot; n_dot, chi_dot];
end