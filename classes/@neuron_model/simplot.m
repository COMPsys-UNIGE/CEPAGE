function [T,X] = simplot(object,Tspan,x0,varargin)
% simplot    Simulates the neuron and plots results
%
% Performs a simulation by means of method neuron_model/sim and plots the
% results. The syntax is the same of method neuron_model/sim.
% [T,X] = simplot(OBJ,T,X0)
% [T,X] = simplot(OBJ,T,X0,OPTS)
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
if nargin == 4
    [TT,XX] = object.sim(Tspan,x0,varargin{1});
else
    [TT,XX] = object.sim(Tspan,x0);
end
names = object.getStateNames;
nx = object.getnx;
figure

for i=1:nx
    subplot(nx,1,i);
    plot(TT,XX(:,i));
    ylabel(names{i});
end
xlabel('T');

subplot(nx,1,1);

title([object.modelName,' model']);


T = TT;
X = XX;