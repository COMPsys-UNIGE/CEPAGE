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


x_dot = zeros(3,1);

% Define param
gna = object.gna;
ENa = object.ENa;
gk2 = object.gk2;
Ek = object.Ek; 
gl = object.gl;
El = object.El;
tNa = object.tNa;
tk2 = object.tk2;
C = object.C;
Iapp = object.Iapp;
VshiftK2 = object.VshiftK2;

% equation
V = x(1);
h = x(2);
m = x(3);

% compute help variable 
hInf = sigmf(V,[-500, -0.0333]);
nInf = sigmf(V,[150, -0.0305]);
mInf = sigmf(V,[83, -(0.018+VshiftK2)]);

n = nInf;
INa = gna*n*n*n*h*(V-ENa);
Ik2 = gk2*m*m*(V-Ek);
Il = gl*(V-El);

% tna0 = 0.0405;
% tauMax = 0.0705;
% VhalfTau = -0.03;
% kTau = 15e-3;
% tNa = tna0;% +((tauMax-tna0)./cosh((V-VhalfTau)/kTau));


V_dot = (-INa-Ik2-Il-Iapp+Isyn)/C;
h_dot = (hInf- h)/tNa;
m_dot = (mInf- m)/tk2;

x_dot = [V_dot; h_dot; m_dot];
end