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


x_dot = zeros(7,1);

% Define param
gna = object.gna;
ENa = object.ENa;
gk = object.gk;
Ek = object.Ek; 
gca = object.gca;
gl = object.gl;
El = object.El;
C = object.C;
Iapp = object.Iapp;
Ca0 = object.Ca0;
d = object.d;
KT = object.KT;
Kd = object.Kd;
gd = object.gd;
D = object.D;
EsynEx = object.EsynEx;
tau = object.tau;
tTm1 = object.tTm1;
tTm2 = object.tTm2;
tTh1 = object.tTh1;
tTh2 = object.tTh2;

% equation
V = x(1);
h = x(2);
m = x(3);
n = x(4);
mT = x(5);
hT = x(6);
Ca = x(7);

% constants
R = 8.31441;  %[J/K*mol]
T = 309.15;  %[K]
F = 96469;  %[C/mol]
k0 = 1000; %(to have E_Ca in mV)
k = 0.1;

% compute help variable 
mTinf = 1/(1 + exp(-(V+52)/7.4));
tTm = 0.44+0.15/(tTm1*exp((V+27)/10)+tTm2*exp(-(V+102)/15));
hTinf = 1/(1+exp((V+80)/5));
tTh = 62.7+0.27/(tTh1*exp((V+48)/4)+tTh2*exp(-(V+407)/50));
ECa = k0*(R*T/(2*F))*log(Ca0/Ca);

Il = gl*(V-El);
INa = gna*m*m*m*h*(V-ENa);
IK = gk*n*n*n*n*(V-Ek);
IT = gca*mT*mT*hT*(V-ECa);


V_dot = ((-Iapp-IT-Il-INa-IK+Isyn+gd*D*(EsynEx-V))/C)/tau;
h_dot = (0.128*exp((17-V)/18)*(1-h)-(4/(exp(-0.2*(V-40))+1))*h)/tau;
m_dot = (((0.32*(13-V))/(exp(0.25*(13-V))-1))*(1-m)-(0.28*(V-40)/(exp(0.2*(V-40))-1))*m)/tau;
n_dot = (((0.032*(15-V))/(exp(0.2*(15-V))-1))*(1-n)-0.5*exp((10-V)/40)*n)/tau;
mT_dot = (-(mT-mTinf)/(tTm))/tau;
hT_dot = (-(hT-hTinf)/(tTh))/tau;
Ca_dot = (-(k*IT)/(2*F*d)-(KT*Ca)/(Ca+Kd))/tau;


x_dot = [V_dot; h_dot; m_dot; n_dot; mT_dot; hT_dot; Ca_dot];
end