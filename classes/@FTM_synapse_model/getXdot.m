function x_dot = getXdot(object,t,x,varargin)
% getXdot    Computes the derivative of the state 
%
%
%  x_dot = getXdot(object,t,x)
%   compute the time derivative of the model 
%
%  x_dot = getXdot(object,t,x,Vpre)
%   compute the time derivative of the model 
%
%  x_dot = getXdot(object,t,x,Vpre,Xold,Vpreold)
%   compute the time derivative of the model 
%
% t is the current time , x is the current state, Vpre is the pre-synaptic
% neuron membrane potential, Xold is the delayed version of x and Vpreold
% is the delayed version of Vpre
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

x_dot = zeros(0);