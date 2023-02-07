function neur = getNeuron(object,varargin)
% getNeuron   Gets the neurons in the network
%
% neur = getNeuron(object)
% neur is a cell array containing the neurons in the network, object is 
% the CPG object.
%
% neur = getNeuron(object,ii)
% neur is the ii-th neurons in the network, object is 
% the CPG object.
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

if nargin == 1
    neur = object.neurons;
    
elseif nargin == 2
    neur = object.neurons(varargin{1});
else
    error('Wrong number of input arguments');
end
end