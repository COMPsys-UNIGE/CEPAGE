function str = getCbuilder(object)

% getCbuilder   return the neuron model C_builder
%
% str = getCbuilder(object)
% object is the neruon object and str is a string containing the C code to
% create a neuron model object.
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

    str = sprintf('IZ_model(%e,%e,%e,%e,%e,%e,%e)',...
    object.a,object.b,object.c,object.d,object.I,object.gL,object.EL);
end