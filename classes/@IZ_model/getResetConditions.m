function [position,isterminal,direction] = getResetConditions(object,t,y,varargin)

% getResetConditions    get the reset condition for neuron model
%
% Compute the reset condition V = 30 for neuron model
%
% [position,isterminal,direction] = getResetConditions(object,t,x)
% where t is current time and x is current state 
%
% position, isterminal and direction have the same role as the ode event
% function. See https://it.mathworks.com/help/matlab/math/ode-event-location.html
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



position = [y(1)-30;-1]; % The value that we want to be zero
isterminal = [1;0];  % Halt integration 
direction = [1;0];   % The zero can be approached from either direction

