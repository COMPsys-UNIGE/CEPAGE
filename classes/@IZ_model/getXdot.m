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

% synapsis input
Isyn = 0;
if nargin >= 4
    if isscalar(varargin{end})
        Isyn = varargin{end};
    end
end



a = object.a;
b = object.b;

I = object.I;

gL = object.gL;
EL = object.EL;

v = x(1);
u = x(2);

x_dot(1,1) = 0.04*v^2+5*v+140-u+I-gL*(v-EL)+Isyn;
x_dot(2,1) = a*(b*v-u);

    
end