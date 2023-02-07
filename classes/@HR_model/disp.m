function disp(object)
% disp   Displays some information about the HR_model object
%
% disp(OBJ)
% OBJ is the HR_model object.

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

    disp('Hindmarsh Rose based neuron model.');
    disp('The Neuron is modeld as follow');
    disp('   _');
    disp('  |');
    disp('  | dx/dt = y − x^3 + bx^2 − z + I + I_{ext}');
    disp('  | ');
    disp(' <  dẏ/dt = 1 − 5x^2 − y');
    disp('  |');
    disp('  | dz/dt = \mu (s (x-x_{rest}) - z)');
    disp('  |_');
    disp(' ');
    disp('Parameters value:');
    disp(['     b = ',num2str(object.b)]);
    disp(['     I = ',num2str(object.I)]);
    disp(['     s = ',num2str(object.s)]);
    disp(['    mu = ',num2str(object.mu)]);
    disp(['x_rest = ',num2str(object.x_rest)]);
end