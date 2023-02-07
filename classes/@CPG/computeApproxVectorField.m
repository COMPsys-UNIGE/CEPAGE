function [deltaPhiMatrix, deltaPhiDot] = computeApproxVectorField(object,PRC,orbit,varargin)
% computeApproxVectorialField   Get the phase difference vectorial field 
% using standard phase reduction method
%
% [deltaPhiMatrix, deltaPhiDot] = computeApproxVectorialField(object,PRC,orbit)
% deltaPhiMatrix is a nStepPhase^N x N-1 matrix describe the phase
% different point in which the phase evolution are computed
% deltaPhiDot is a nStepPhase^N x N-1 matrix that describe the evolution of
% the phase difference 
% PRC is a struct with fields phi and PRC that describe the phase resetting
% curve of the neuron model; orbit is a struct with fields T and X that
% describe the time evolution of a signle orbit of the neuron. To compute
% PRC and orbit you can use method computePRC of neuron_model object.
%
% [deltaPhiMatrix, deltaPhiDot] = computeApproxVectorialField(object,PRC,orbit,optd)
% opts is a struct with two fields:
% - nStepPhase: the number of suddivision over each of the N-1 dimension of 
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
N = object.N;

nStepIntegrale = 20;
nStepFase = 10;

if nargin == 4
    if isfield(varargin{1},'nStepPhase')
        nStepFase = varargin{1}.nStepPhase;
    end
     if isfield(varargin{1},'nStepIntegral')
        nStepIntegrale = varargin{1}.nStepIntegral;
    end
end

deltaPhiSingle = linspace(0,1,nStepFase);
phaseVec1 = linspace(0,1,nStepIntegrale);

deltaPhiCellTmp = cell(N-1,1);
deltaPhiCell = cell(N-1,1);
for i=1:N-1
    deltaPhiCellTmp{i} = deltaPhiSingle;
end

[deltaPhiCell{1:N-1}] = ndgrid(deltaPhiCellTmp{:});

Qfase = PRC.phi;
Qval = PRC.PRC(:,1);



Torbit = orbit.T;
Xorbit = orbit.X;
phiOrbit = linspace(0,1,numel(Torbit));

np = numel(deltaPhiCell{1});

deltaPhiMatrix = zeros(numel(deltaPhiCell{1}),N-1);

for i=1:N-1
    deltaPhiMatrix(:,i) = deltaPhiCell{i}(:);
end

deltaPhiDot = zeros(size(deltaPhiMatrix));

% loop on evry neuron
for i=1:np
    
    deltaPhi = deltaPhiMatrix(i,:);
    
    % loop on evry delta
    integrale = zeros(N-1,1);
    h = waitbar(i/np);
    for j=1:nStepIntegrale
        phi1 = phaseVec1(j);
        phij = phi1-deltaPhi;
        
        phi =[phi1;phij(:)];
        phi(phi<0) = phi(phi<0)+1;
        Qtmp = interp1(Qfase,Qval,phi);
        Vtmp = interp1(Torbit,Xorbit(:,1),phi*Torbit(end));
        
        integrale = integrale + f_function(object,Qtmp,Vtmp(:,1));
    end
    deltaPhiDot(i,:) = -integrale(:)';
end
close(h)



end

