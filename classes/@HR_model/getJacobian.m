function J = getJacobian(object,t,x)

% getJacobian    Computes the Jacobian of the vector field in point x and
% time t
%
%  J = getJacobian(object,t,x,Isyn)
%   compute theJacobian of the model in state x ant time t
%  J = getJacobian(object,t,x,Isyn)
%   compute theJacobian of the model in state x,time t and with synaptic 
%   current Isyn
%
% Contributors:
%
% Matteo Lodi (matteo.lodi@edu.unige.it)
%
% Copyright (C) 2016 University of Genoa, Italy.

% Legal note:
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public
% License as published by the Free Software Foundation; either
% version 2.1 of the License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
%
% You should have received a copy of the GNU General Public
% License along with this library; if not, write to the
% Free Software Foundation, Inc.,
% 59 Temple Place, Suite 330,
% Boston, MA  02111-1307  USA

Isyn = 0;
if nargin == 4
    Isyn = varargin{1};
end

nx = object.nx;

if size(x,1) ~= nx
    error('input vector must be a nx x npoints vector');
end

np = size(x,2);

J = cell(np,1);

mu = object.mu;
s = object.s;
b = object.b;

for i=1:np

xx = x(:,i);
    
J1 = -3*xx(1)^2+2*b*xx(1) + Isyn;
J2 = 1;
J3 = -1;
J4 = -10*xx(1);
J5 = -1;
J6 = 0;
J7 = mu*s;
J8 = 0;
J9 = -mu;
J{i} = [J1 J2 J3; J4 J5 J6; J7 J8 J9];
end

if size(x,2) == 1
    JJ = J{1};
    J = [];
    J = JJ;
end
