function disp(object)
% TO DO

% disp   Displays some information about the IZ_model object
%
% disp(OBJ)
% OBJ is the IZ_model object.

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
disp('Izhikevich neuron model.');
disp('The Neuron is modeld as follow');
disp('   _');
disp('  |');
disp('  | dv/dt = 0.04*v^2+5*v+140-u+I+gL(v-EL)-I_{ext}');
disp(' <');
disp('  | du/dt = a*(b*v-u);');
disp('  |_-');
disp('   ');
disp(' with the reset condition');
disp('                 _');
disp('                |  v = c;');
disp(' if v > 30     <');
disp('                |_ u = u+d;');
disp('   ');
disp('Parameters value:');
disp(['     a = ',num2str(object.a)]);
disp(['     b = ',num2str(object.b)]);
disp(['     c = ',num2str(object.c)]);
disp(['     d = ',num2str(object.d)]);
disp(['     I = ',num2str(object.I)]);
disp(['     gL = ',num2str(object.gL)]);
disp(['     EL = ',num2str(object.EL)]);
end