function J = getJacobian(object,t,x)
nx = object.nx;

if size(x,1) ~= nx
    error('input vector must be a nx x npoints vector');
end

np = size(x,2);

J = cell(np,1);

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

for i=1:np

xx = x(:,i);

V = xx(1);
h = xx(2);
n = xx(3);
n_A = xx(4);
h_A = xx(5);
h_T = xx(6);
m_HVA = xx(7);
Ca = xx(8);


J{i} = [];

end

if size(x,2) == 1
    JJ = J{1};
    J = [];
    J = JJ;
end