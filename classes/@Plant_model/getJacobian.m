function J = getJacobian(object,t,x)
nx = object.nx;

if size(x,1) ~= nx
    error('input vector must be a nx x npoints vector');
end

np = size(x,2);

J = cell(np,1);

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

for i=1:np

xx = x(:,i);

V = xx(1);
Ca = xx(2);
h = xx(3);
n = xx(4);
chi = xx(5);


J{i} = [];

end

if size(x,2) == 1
    JJ = J{1};
    J = [];
    J = JJ;
end