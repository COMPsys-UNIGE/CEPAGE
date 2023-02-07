function generateC(object,fileName,varargin)
% generateC Generates C files for the computation of model vector field
%
% generateC(OBJ,filename)
% Generates the C filename.c and filename.h files for the computation of
% model vector field
%
% generateC(OBJ,filename,folder)
% Generates the C filename.c and filename.h files for the computation of
% model vector field in the directory folder
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
if nargin >= 3
    folder = [varargin{1},'/'];
    if(exist(varargin{1},'dir') ~= 7)
        mkdir(folder);
    end
else
    folder = '';
end
gpustring = '';
if nargin == 4
    if(varargin{2})
        gpustring = ' __device__ ';
    end
end
fout = fopen([folder,fileName,'.h'],'w+');
fprintf(fout,['#ifndef ',fileName,'_H\n']);
fprintf(fout,['#define ',fileName,'_H\n\n']);
fprintf(fout,'#include <math.h>\n');
fprintf(fout,['void ',gpustring,fileName,'(double *x,double *xdot,double Iext);\n\n']);
fprintf(fout,['void ',gpustring,fileName,'_jac(double *x,double **jac,double Iext);\n\n']);
fprintf(fout,'#endif');
fclose(fout);

fout = fopen([folder,fileName,'.c'],'w+');
fprintf(fout,['#include "',[fileName,'.h'],'"\n\n']);
fprintf(fout,['void ',gpustring,fileName,'(double *x,double *xdot,double Iext)\n']);
fprintf(fout,'{\n');
fprintf(fout,'\tconst double a = %f;\n',object.a);
fprintf(fout,'\tconst double b = %f;\n',object.b);
fprintf(fout,'\tconst double c = %f;\n',object.c);
fprintf(fout,'\tconst double d = %f;\n',object.d);
fprintf(fout,'\tconst double I = %f;\n',object.I);

fprintf(fout,'\tconst double gL = %f;\n',object.gL);
fprintf(fout,'\tconst double El = %f;\n',object.El);


fprintf(fout,'\n');
fprintf(fout,'xdot[0] = 0.04*v*v+5*v+140-u+I+gL*(v-El)+Iext;;\n');
fprintf(fout,'xdot[1] = a*(b*v-u);\n');
    fprintf(fout,'}\n\n');

fprintf(fout,['void ',gpustring,fileName,'_jac(double *x,double **jac,double Iext)\n']);
fprintf(fout,'{\n');
fprintf(fout,'}');

fclose(fout);


end