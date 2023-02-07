function disp(object)
% disp   Displays some information about the HH_model object
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
    disp('Hodgkin-Huxley based neuron model.');
    disp('The Neuron is modeld as follow');
    disp('   _');
    disp('  |');
    disp('  | dV/dt = (-I_{Na}-I_{k2}-I_l-I_{app}-I_{ext})/C; ');
    disp('  | ');
    disp('  <  dh/dt = (h_{Inf}- h)/t_{Na}');
    disp('  |');
    disp('  | dm/dt = (m_{Inf}- m)/t_{k2}');
    disp('  |_');
    disp('  ');
    disp('  h_{Inf} = 1/(1 + exp(500*(X+0.333)))');
    disp('  n_{Inf} = 1/(1 + exp(-150*(X+0.305)))');
    disp('  m_{Inf} = 1/(1 + exp(-83*(X+(0.018+VshiftK2))))');
    disp('  ');
    disp('  n = n_{Inf};');
    disp('  I_{Na} = gna*n^3*h*(V-E_{Na})');
    disp('  I_{k2} = gk2*m^2*(V-E_k)');
    disp('  I_l = g_l*(V-E_l)');
    disp(' ');
    disp('Parameters value:');

    disp(['     gna = ',num2str(object.gna)]);
    disp(['     ENa = ',num2str(object.ENa)]);
    disp(['     gk2 = ',num2str(object.gk2)]);
    disp(['      Ek = ',num2str(object.Ek)]);
    disp(['      gl = ',num2str(object.gl)]);
    disp(['      El = ',num2str(object.El)]);
    disp(['     tNa = ',num2str(object.tNa)]);
    disp(['     tk2 = ',num2str(object.tk2)]);
    disp(['       C = ',num2str(object.C)]);
    disp(['    Iapp = ',num2str(object.Iapp)]);
    disp(['VshiftK2 = ',num2str(object.VshiftK2)]);

end