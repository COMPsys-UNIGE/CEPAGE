function [PRC,limitCycle] = computePRC(object,nPoints,Ttrans,varargin)

% computePRC   Compute the Phase resetting curve for the neuron
%
% [PRC,limitCycle] = computePRC(object,nPoints,Ttrans)
% nPoints is the number of points in which the function eval the PRC, 
% Ttrans is the transitory time used before finding a limit cycle and 
% object is the neuron_model object.
%
% [PRC,limitCycle] = computePRC(object,nPoints,Ttrans,opt)
% A structure opt can be provided with % the following fields:
%   - x0: integrator starting condition
%   - Vth: voltage threshold used to select a single period
% To compute PRC the toolbox use the method explained in 
% NoviÄ?enko, Viktor, and Kestutis Pyragas. "Computation of phase response 
% curves via a direct method adapted to infinitesimal perturbations." 
% Nonlinear Dynamics 67.1 (2012): 517-526.
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

if ~object.isContinuous
    error('PRC can be computed only for continuous systems.');
end
nx = object.getnx;

opt = [];

Tlimit = 1e7;

Vth = 0;
x0 = zeros(nx,1);


if nargin == 4
    if isstruct(varargin{1})
        opt = varargin{1};
    else
        error('option must be a struct');
    end
end

if isfield(opt,'Vth')
    Vth = opt.Vth;
elseif isfield(opt,'x0')
    x0 = opt.x0;
elseif isfield(opt,'Tlimit')
    Tlimit = opt.Tlimit;
end

% Build limit cycle

eventFun = @(T, Y) eventFun_(T, Y, Vth);

[~,X] = ode45(@object.getXdot,[0 Ttrans],x0,odeset);
[~,X] = ode45(@object.getXdot,[0 Tlimit],X(end,:),odeset('events',eventFun));
CI = X(end,:);

[Torbit, Xorbit] = ode45(@object.getXdot,[0 Tlimit],CI,odeset('events',eventFun,'RelTol',1e-10,'AbsTol',1e-25));

period = Torbit(end);

[~,ind] = max(Xorbit(:,1));
% [~,ind] = min(Xorbit(:,1));
% tmp1 = [0;diff(Xorbit(:,1)>=0)];
% ind = find(tmp1==1);
rightCI = Xorbit(ind,:);

[Torbit, Xorbit] = ode45(@object.getXdot,[0 period],rightCI,odeset('RelTol',1e-15,'AbsTol',1e-20));

lc.T = Torbit;
lc.X = Xorbit;

limitCycle = lc;
% Compute fundamenta matrix

dphi = period/nPoints;

fundMatr = object.getFundMatrix(nPoints,lc);

% Extract eigenvector and get PRC


PRC_ = zeros(nPoints,nx);
for i=1:nPoints
    X = interp1(Torbit,Xorbit,(i-1)*dphi);
    V = object.getXdot(0,X);
    fMatr_i = fundMatr{i};
    
    L1 = zeros(nx,1);
    L1(1) = 1;
    
    L1tr = zeros(1,nx);
    L1tr(1,1) = 1;
    
    A = fMatr_i(2:end,2:end)-eye(1);
    B = fMatr_i(1,2:end);
    
    L1tr(1,2:end) = -B*inv(A);
    
    L1 = L1tr';
    
    tmp = L1 * inv(L1tr * V);
    PRC_(i,:) = tmp;
end

PRC.phi = linspace(0,1,nPoints);
PRCdotdot=[0;0;diff(diff(PRC_(:,1)))];
[~,locs] = findpeaks(PRCdotdot,'MinPeakHeight',0.3);
PRC_=[PRC_;PRC_];
locs(1)=locs(1)+nPoints;
PRC_(locs(1)-10:locs(1)+10,1)=linspace(PRC_(locs(1)-10,1),PRC_(locs(1)+10,1),21);
PRC_(locs(2)-10:locs(2)+10,1)=linspace(PRC_(locs(2)-10,1),PRC_(locs(2)+10,1),21);
PRC_(1:20,1)=PRC_(nPoints+1:nPoints+20,1);

% mm = (PRC_(end)+PRC_(1))/2;
% PRC_(1) = mm;
% PRC_(end) = mm;


PRC.PRC = PRC_(1:nPoints,:);

    function  [value, isterminal, direction] = eventFun_(t,y,Vth_)
        value = y(1)-Vth_;
        isterminal = 1;
        direction = 1;
    end

end

