function [w_in,w_ex,w_el,M_in,M_ex,M_el] = computeMw(object,PRC,limitCycle,varargin)
% computeMw   Compute M and w functions from Phase Resetting Curve
% using standard phase reduction method. For mathematical explanation see
% the pdf toolbox manual
%
% [w_in,w_ex,w_el,M_in,M_ex,M_el] = computeMw(object,PRC,limitCycle)
% w_in, w_ex and w_el are struct with fields deltaPhi and value that 
% represent the x and the value of the function; the value can be computed
% through interp1 function.
% M_in, M_ex and M_el are struct with fields deltaPhi1i, deltaPhi1j and 
% value that represent the x, the y and the value of the function; the 
% value can be computed through interp2 function.
% PRC is a struct with fields phi and PRC that describe the phase resetting
% curve of the neuron model; orbit is a struct with fields T and X that
% describe the time evolution of a signle orbit of the neuron. To compute
% PRC and orbit you can use method computePRC of neuron_model object.
%
% [w_in,w_ex,w_el,M_in,M_ex,M_el] = computeMw(object,PRC,limitCycle,opts)
% opts is a struct with two fields:
% - nDeltaPhi: the number of suddivision over each of the N-1 dimension of 
% the phase difference domain (default: 10)
% - nStepIntegral: the step over which the integral is computed (default: 20)
% 
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


if ~object.is_continuous
    error('matrices M and w can be computed only for continuous CPGs');
end

warndlg('Warning! This method works only for CPG with equal synapses model.')

nStepIntegral = 20;
nDeltaPhi = 10;

if nargin == 4
    opt = varargin{1};
    if ~isstruct(opt)
        error('opt must be a struct');
    else
        if isfield(opt,'nStepIntegral')
            nStepIntegral = opt.nStepIntegral;
        end
        if isfield(opt,'nDeltaPhi')
            nDeltaPhi = opt.nDeltaPhi;
        end
    end
end

EsynEx = object.EsynEx;
EsynIn = object.EsynIn;

inhActivation = object.inhActivation;
excActivation = object.excActivation;

% Compute functions w
phi1 = linspace(0,1,nStepIntegral);
deltaPhi1jw = linspace(0,1,nDeltaPhi);

w_ex_ = zeros(nDeltaPhi,1);
w_in_ = zeros(nDeltaPhi,1);
w_el_ = zeros(nDeltaPhi,1);

dphi = diff(phi1);
dphi = dphi(1);

T = limitCycle.T(end);

disp('Compute w');

parfor i=1:nDeltaPhi
    for k=1:nStepIntegral
        phii = phi1(k)-deltaPhi1jw(i);
        if(phii < 0)
            phii = phii+1;
        end
        Q = interp1(PRC.phi,PRC.PRC(:,1),phi1(k),'spline');
        V1 = interp1(limitCycle.T,limitCycle.X(:,1),phi1(k)*T,'spline');
        Vj = interp1(limitCycle.T,limitCycle.X(:,1),phii*T,'spline');
        
        w_ex_(i) = w_ex_(i)+Q*(-V1+EsynEx)*excActivation{1,1}.getActivation(Vj);
        
        w_in_(i) = w_in_(i)+ Q*(-V1+EsynIn)*inhActivation{1,1}.getActivation(Vj);
        
        w_el_(i) = w_el_(i)+ Q*(Vj-V1);
    end
end

w_in_ = w_in_*dphi;
w_ex_ = w_ex_*dphi;
w_el_ = w_el_*dphi;

w_in.value = w_in_;
w_in.deltaPhi = deltaPhi1jw;

w_ex.value = w_ex_;
w_ex.deltaPhi = deltaPhi1jw;

w_el.value = w_el_;
w_el.deltaPhi = deltaPhi1jw;

% Compute functions M


phi1 = linspace(0,1,nStepIntegral);

phiVec = linspace(0,1,nDeltaPhi);

[deltaPhi1i,deltaPhi1j] = meshgrid(phiVec,phiVec);




M_ex_ = zeros(nDeltaPhi,nDeltaPhi);
M_in_ = zeros(nDeltaPhi,nDeltaPhi);
M_el_ = zeros(nDeltaPhi,nDeltaPhi);

T = limitCycle.T(end);

disp('Compute M');

for i=1:nDeltaPhi
    parfor j=1:nDeltaPhi
        for k=1:nStepIntegral
            phii = phi1(k)-deltaPhi1i(i,j);
            if(phii < 0)
                phii = phii+1;
            end
            phij = phi1(k)-deltaPhi1j(i,j);
            if(phij < 0)
                phij = phij+1;
            end
            Q = interp1(PRC.phi,PRC.PRC(:,1),phii,'spline');
            Vi = interp1(limitCycle.T,limitCycle.X(:,1),phii*T,'spline');
            Vj = interp1(limitCycle.T,limitCycle.X(:,1),phij*T,'spline');
            
            M_ex_(i,j) = M_ex_(i,j)+Q*(-Vi+EsynEx)*excActivation{1,1}.getActivation(Vj);
            
            M_in_(i,j) = M_in_(i,j)+ Q*(-Vi+EsynIn)*inhActivation{1,1}.getActivation(Vj);
            
            M_el_(i,j) = M_el_(i,j)+ Q*(-Vi+Vj);
        end
    end
    waitbar(i/nDeltaPhi);
end

M_in_ = M_in_*dphi;
M_ex_ = M_ex_*dphi;
M_el_ = M_el_*dphi;

M_in.value = M_in_;
M_in.deltaPhi1i = deltaPhi1i;
M_in.deltaPhi1j = deltaPhi1j;

M_ex.value = M_ex_;
M_ex.deltaPhi1i = deltaPhi1i;
M_ex.deltaPhi1j = deltaPhi1j;

M_el.value = M_el_;
M_el.deltaPhi1i = deltaPhi1i;
M_el.deltaPhi1j = deltaPhi1j;
