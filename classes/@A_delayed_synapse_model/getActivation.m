function act = getActivation(object,Vpre,varargin)
% getActivation    Computes the activation function of the synapse
%
%  act = getActivation(object,Vpre)
%   compute the time activation function of the synapse; Vpre is the
%   pre-synaptic neuron membrane potential
%
%  act = getActivation(object,Vpre,x)
%   compute the time activation function of the synapse; Vpre is the
%   pre-synaptic neuron membrane potential and x are the synapse states
%
%  act = getActivation(object,Vpre,x,VpreOld)
%   compute the time activation function of the synapse; Vpre is the
%   pre-synaptic neuron membrane potential, x are the synapse states and
%   VpreOld is the delayed pre-synaptic neuron membrane potential

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

if nargin == 3
    x1 = varargin{1}(1);
    x2 = varargin{1}(2);
else
	x1 = 0;
    x2 = 0;
end

act = object.g1*x1+object.g2*x2;
end
