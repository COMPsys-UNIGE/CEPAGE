function disp(object)
% disp   Displays some information about the FTM_synapse_model object
%
% disp(OBJ)
% OBJ is the FTM_synapse_model object.

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
disp('Alpha synapse model.');
disp('The synapse activation function is');
disp('   ');
disp('x');
disp('   ');
disp('The state equations is:');
disp(' alpha*(1-s)*1/(1+exp(-nu*(Vpre-theta))-beta*s');

disp(' ');
disp('Parameters value:');
disp(['     alpha = ',num2str(object.al)]);
disp(['      beta = ',num2str(object.beta)]);
disp(['     theta = ',num2str(object.theta)]);
disp(['     nu    = ',num2str(object.nu)]);

end