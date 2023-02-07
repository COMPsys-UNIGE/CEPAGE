function CI = getCIfromPhi(object,deltaPhi,varargin)
% getCIfromPhi   Get the initial condition of a CPG in which
% the phase difference between neuron 1 and neuron i+1 is deltaPhi(i);
%
% CI = getCIfromPhi(object,deltaPhi)
% CI is a cell array with N-1 elements describing the evolution of the
% phase difference between neuron 1 and neuron i+1
% deltaPhi is a Mx(N-1) matrix; each row represent a different initial
% condition of the network.
% 
% CI = getCIfromPhi(object,deltaPhi,options)
% options is a struct with fields:
% - Vth - threshold that define events, it must be between max and minimum
%         value of the membrane potential (default 0)
% - Ttrans - Transitory time of the neuron (default 100)
% - x0 - Initial condition of the neuoron (default 0)
% - integrator: string indicating the solver used to integrate the
%                differential equations. It can be either 'ode45','ode23',
%                'ode113','ode15s','ode23s','ode23t' , 'ode23tb', 'odeint' or 'eulero.
%                 If you choose eulero or odeint the simulation is made in C throught mex
%                 file and a supported mex compiler is required.
%                Default: 'ode45'.
% - integratorOptions: options to provide to the ODE solver. Type doc odeset to
%                 get help. If you choose eulero the variable must have a
%                 field dt that describe the integration step size.
%                 Default: odeset.
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
N = object.N;
neur = object.neurons{1};
x0 = zeros(neur.getnx,1);
nx = object.neurons{1}.getnx;

Ttrans = 100;
Vth = 0;

integrator = 'ode45';
integratorOptions = odeset;

if nargin == 3
    if isstruct(varargin{1})
        if isfield(varargin{1},'Vth')
            Vth = varargin{1}.Vth;
        end
        if isfield(varargin{1},'Ttrans')
            Ttrans = varargin{1}.Ttrans;
        end
        if isfield(varargin{1},'x0')
            if numel(x0) == neur.getnx
            x0 = varargin{1}.x0;
            end
        end
        if isfield(varargin{1},'integrator')
            integrator = varargin{1}.integrator;
        end
        if isfield(varargin{1},'integratorOptions')
            integratorOptions = varargin{1}.integratorOptions;
        end
    end
end



if ((numel(deltaPhi(:,1)) == N-1) || (numel(deltaPhi(1,:)) == N-1))
    %Se è messo male lo giro
    if numel(deltaPhi(:,1)) == N-1
        deltaPhi = deltaPhi';
    end
else
    error 'input array must be Mx(N-1)';
end

if ~(all(deltaPhi(:) <= 1) && all(deltaPhi(:) >= 0))
    error 'phi must be between 0 and 1';
end

[M, ~] = size(deltaPhi);

CI = zeros(M,object.totState);


% Calcolo il periodo T e trovo le CI1
opt.integratorOptions = integratorOptions;
opt.integrator = integrator;


[T, X] = neur.sim([0 Ttrans],x0,opt);

up = find(diff(X(:,1) > Vth) == 1);  

Te = T(up);

% Prendo le ultime così sono più sicuro che il transitorio sia esaurito
% e calcolo periodo e CI
Tau = Te(end)-Te(end-1);
CI1 = X(up(end) , 1:neur.getnx);
CI(:,1:neur.getnx) = repmat(CI1,M,1);
% Quindi calcolo le altre

x0 = CI1;


[T, X] = neur.sim([0 Tau*1.5],x0,opt);

[T,ii,c] = unique(T);
X = X(ii,:);

TT = Tau*(1-deltaPhi);
TT(TT > Tau) = Tau;
for i=1:N-1
    XX = interp1(T,X,TT(:,i));
    CI(:,(i*nx+1):(i*nx+nx)) = XX;
end


end