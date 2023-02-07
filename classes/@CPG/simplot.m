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

if size(x0,1) > 1
    error(['x0 must have one row; in order to simulate network behaviour',...
        'starting from different initial condition you must use sim method']);
end

if nargin == 4
    [TT,XX] = object.sim(Tspan,x0,varargin{1});
else
    [TT,XX] = object.sim(Tspan,x0);
end

nx = object.neurons{1}.getnx;
N = object.N;
figure
for j = 1:N
    names = object.neurons{j}.getStateNames;
    for i=1:nx
        subplot(nx,1,i);
        hold on
        plot(TT,XX(:,i+(j-1)*nx));
        ylabel(names{i});
    end
end
xlabel('T');

subplot(nx,1,1);

title('Neurons network model');


T = TT;
X = XX;