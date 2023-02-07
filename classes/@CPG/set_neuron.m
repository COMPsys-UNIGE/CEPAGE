function object = set_neuron(object,newNeur,ii)
% set_neuron   Set the neuron ii-th of the network
%
% object = set_neuron(object,newNeur,ii)
% object is a neuron_network object, newNeur is neuron_model object, ii is
% the index of the nueron to change
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
if ~isa(newNeur,'neuron_model')
    error('newNeur must be a neruon_model object ')
end
object.neurons{ii} = newNeur;