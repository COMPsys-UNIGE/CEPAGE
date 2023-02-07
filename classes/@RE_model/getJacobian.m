function J = getJacobian(object,t,x)
nx = object.nx;

if size(x,1) ~= nx
    error('input vector must be a nx x npoints vector');
end

np = size(x,2);

J = cell(np,1);

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

for i=1:np

xx = x(:,i);

V = xx(1);
h = xx(2);
m = xx(3);
n = xx(4);
mT = xx(5);
hT = xx(6);
Ca = xx(7);


J{i} = [];

end

if size(x,2) == 1
    JJ = J{1};
    J = [];
    J = JJ;
end