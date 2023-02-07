function J = getJacobian(object,t,x)

% getJacobian    Computes the Jacobian of the vector field in point x and
% time t
%
%  J = getJacobian(object,x)
%   compute theJacobian of the model in state x
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


nx = object.nx;

if size(x,1) ~= nx
    error('input vector must be a nx x npoints vector');
end

np = size(x,2);

J = cell(np,1);

gK2 = object.gk2;
gL  = object.gl;
gNa  = object.gna;
C  = object.C;
EK2  = object.Ek;
ENa  = object.ENa;
EL  = object.El;
Iapp  = object.Iapp;
tauNa  = object.tNa;
tauK2  = object.tk2;
Vshift = object.VshiftK2;

for i=1:np

xx = x(:,i);

V = x(1);
mK2 = x(2);
hNa = x(3);


J{i} = [ -(gL + gK2*mK2^2 + (gNa*hNa)/(exp(- 150*V - 183/40) + 1)^3 - (450*gNa*hNa*exp(- 150*V - 183/40)*(ENa - V))/(exp(- 150*V - 183/40) + 1)^4)/C, (2*gK2*mK2*(EK2 - V))/C, (gNa*(ENa - V))/(C*(exp(- 150*V - 183/40) + 1)^3);
                                                    (83*exp(- 83*V - 83*Vshift - 747/500))/(tauK2*(exp(- 83*V - 83*Vshift - 747/500) + 1)^2),                -1/tauK2,                                                 0;
                                                                                  -(500*exp(500*V + 65/4))/(tauNa*(exp(500*V + 65/4) + 1)^2),                       0,                                          -1/tauNa];

end

if size(x,2) == 1
    JJ = J{1};
    J = [];
    J = JJ;
end

