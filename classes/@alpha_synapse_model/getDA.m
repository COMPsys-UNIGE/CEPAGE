function J = getDA(object,x,Vpre,varargin)

% getDA    Computes the differential of the acrivation function with
%           respect to presynaptic potential Vpre and states x
%
%  J = getDA(object,Vpre)
%   Computes the differential of the acrivation function with respect to 
%   presynaptic potential Vpre and with state x=0;
%
%  J = getDA(object,Vpre,x)
%   Computes the differential of the acrivation function with respect to 
%   presynaptic potential Vpre and with state x
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


J = [0 1];