function J = getJacobian(object,t,x,Vpre)

nx = object.nx;

if size(x,1) ~= nx
    error('input vector must be a nx x npoints vector');
end

np = size(x,2);


J = cell(np,1);

for i=1:np
    J{i} = [];
end

if size(x,2) == 1
    JJ = J{1};
    J = [];
    J = JJ;
end