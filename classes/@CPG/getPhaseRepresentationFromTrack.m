function [Tphi,phi] = getPhaseRepresentationFromTrack(object,T,X,Vth )
% getPhaseRepresentationFromTrack   Gets the phase evolution of the network
%
% [Tphi,phi] = getPhaseRepresentationFromTrack(object,T,X,Vth )
% phi is a cell array with N-1 elements describing the evolution of the
% phase difference between neuron 1 and neuron i+1
% Tphi is a cell array with N-1 elements describing the time
% T and X are the time and state vectors that describe the time evolution of
% the network, Vth is the value of the membrane potential that correspond to 
% an event (it must be between max and minimum value of the membrane 
% potential), object is the CPG object.
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


totState = object.totState;
incrementalIndexState = object.incrementalIndexState;


nc = size(X,2);

if nc ~= totState
    error(['X must have ',num2str(N*nx),'columns']);
end

ii = incrementalIndexState(1:N);

X = X(:,ii);

Te = cell(N,1);

low = X < Vth;
up = X >= Vth;

len = Inf;

for i=1:N
    try
    evDif = low(1:end-1,i)-up(2:end,i);
    evSum = low(1:end-1,i)+up(2:end,i);
    Te{i} = T(evDif == 0 & evSum == 2);
    len = min(len,numel(Te{i}));
    catch
    end
end

T1 = Te{1}(end)-Te{1}(end-1);

phi = zeros(N-1,len);
Tphi = zeros(len,1);
for i=2:N
    phi(i-1,1:len) = mod((Te{i}(1:len)-Te{1}(1:len))/T1,1);
    Tphi(1:len) = Te{1}(1:len);
end

phi = phi';
