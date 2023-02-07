function J = getJacobian(object,t,x,Vpre)

% getJacobian    Computes the Jacobian of the vector field in point x and
% time t and presynaptic voltage Vpre with respect to [Vpre x]'
%
%  J = getJacobian(object,t,x)
%   compute theJacobian of the model in state x and time t and with
%   presynaptic potential 0
%
%  J = getJacobian(object,t,x,Vpre)
%   compute theJacobian of the model in state x, time t and with
%   presynaptic potential Vpre
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

for i=1:np
    J{i} = [];
end

if size(x,2) == 1
    JJ = J{1};
    J = [];
    J = JJ;
end
