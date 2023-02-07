function writeApproxVectorField(object,w_in,w_ex,w_el,M_in,M_ex,M_el,varargin)

% writeApproxVectorField Generates files for the computation of 
%   approximate CPG vector field
%
% writeApproxVectorField(OBJ,w_in,w_ex,w_el,M_in,M_ex,M_el)
% Generates the  filename.m file for the computation of the approximate 
% vector field
%
% writeApproxVectorField(OBJ,w_in,w_ex,w_el,M_in,M_ex,M_el,opt )
% Generates files for the computation of the approximate vector field; opt
% is a struct with the following fields:
% - code: a string describing the code language to be generated. It can be
% 'c' (c code), 'm' (matlab code), 'auto' (c+auto07 code), 'matcont' 
%  ('m'+matcont code).
% - filename : name of the function and of the file that describe the
% vectorial field
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

if ~object.is_continuous
    error('approximate CPG vector field can be computed only for continuous CPGs');
end


filename = 'vectorialField';
code = 'm';
folder = [];

if nargin == 8
    opt = varargin{1};
    if ~isstruct(opt)
        error('opt must be a struct');
    else
        if isfield(opt,'code')
            code = opt.code;
        end
        if isfield(opt,'filename')
            filename = opt.filename;
        end
        if isfield(opt,'folder')
            folder = opt.folder;
        end
    end
end

if ~isempty(folder)
    if ~strcmp(folder(end),'\') && ~strcmp(folder(end),'/')
        folder = [folder,'/'];
    end
    
    if strcmp(folder(end),'\')
        folder(end) = '/';
    end
    
    mkdir(folder(1:end-1));
end

%%%%%%% MATLAB MATCONT CODE %%%%%%%%
if strcmp(code,'m') || strcmp(code,'matcont')
    
    save([folder,'Mw_cont.mat'],'w_in','w_ex','w_el','M_in','M_ex','M_el');
    
    if strcmp(code,'matcont')
    
    fin = fopen([folder,'matcont_bifurcation_format.m'],'w');
    
    fprintf(fin,'function out = matcont_bifurcation_format(w_in,w_ex,w_el,M_in,M_ex,M_el)\n');
    fprintf(fin,'out{1} = @init;\n');
    fprintf(fin,'out{2} = @(t,x,p)fun_eval(t,x,p,w_in,w_ex,w_el,M_in,M_ex,M_el);\n');
    fprintf(fin,'out{3} =[];\n');
    fprintf(fin,'out{4} = [];\n');
    fprintf(fin,'out{5} = [];\n');
    fprintf(fin,'out{6} = [];\n');
    fprintf(fin,'out{7} = [];\n');
    fprintf(fin,'out{8} = [];\n');
    fprintf(fin,'out{9} = [];\n\n');
    fprintf(fin,'%% --------------------------------------------------------------------------\n');
    fprintf(fin,'function dydt = fun_eval(time,kmrgd,bif_param,w_in,w_ex,w_el,M_in,M_ex,M_el)\n');
    
    fprintf(fin,'%%---------------- WRITE HERE NETWORK CONFIGURATION AND BIFURCATION PROBLEM ----------------%%\n');
    fprintf(fin,'g_in =  [];\n');
    fprintf(fin,'g_ex =  [];\n');
    fprintf(fin,'g_el =  [];\n\n');
    
    
    fprintf(fin,'N = (numel(kmrgd)/2)+1;\n');
    fprintf(fin,'bigPhi = mod(kmrgd,1);\n\n');
    
    fprintf(fin,'\tdydt = %s(time,bigPhi,g_in,g_ex,g_el,w_in,w_ex,w_el,M_in,M_ex,M_el);\n\n',filename);
    
    
    fprintf(fin,'%% --------------------------------------------------------------------------\n');
    fprintf(fin,'function [tspan,y0,options] = init\n');
    fprintf(fin,'%%---------------- WRITE HERE INITIAL CONDITION ----------------%%\n');
    fprintf(fin,'y0 = [];\n');
    fprintf(fin,'options = odeset;\n');
    fprintf(fin,'handles = feval(odefile);\n');
    fprintf(fin,'tspan = [0 10];\n');
        
    fclose(fin);
    
    end
    %%%%%%% VECTORIAL FIELD CODE %%%%%%%%
    
    fin = fopen([folder,filename,'.m'],'w');
    fprintf(fin,'function xdot = %s(t,x,g_in,g_ex,g_el,w_in,w_ex,w_el,M_in,M_ex,M_el)\n',filename);
    fprintf(fin,'N = %d;\n',object.getN);
    fprintf(fin,'xdot = zeros(%d,1);\n',object.N-1);
    fprintf(fin,'bigPhi = [0;x];\n');
    
    fprintf(fin,'constantContrib = 0;\n');
    fprintf(fin,'for j=1:N\n');
    fprintf(fin,'\tphij = bigPhi(j);\n');
    fprintf(fin,'\tphij = mod(phij,1);\n');
    fprintf(fin,'\tconstantContrib = constantContrib+g_in(1,j)*interp1(w_in.deltaPhi,w_in.value,phij);\n');
    fprintf(fin,'\tconstantContrib = constantContrib+g_ex(1,j)*interp1(w_ex.deltaPhi,w_ex.value,phij);\n');
    fprintf(fin,'\tconstantContrib = constantContrib+g_el(1,j)*interp1(w_el.deltaPhi,w_el.value,phij);\n');
    fprintf(fin,'end\n\n');
    
    fprintf(fin,'for i=2:N\n');
    fprintf(fin,'\ttmp = 0;\n');
    fprintf(fin,'\tfor j=1:N\n');
    fprintf(fin,'\t\tphi1i = bigPhi(i);\n');
    fprintf(fin,'\t\tphi1j = bigPhi(j);\n');
    fprintf(fin,'\t\tphi1i = mod(phi1i,1);\n');
    fprintf(fin,'\t\tphi1j = mod(phi1j,1);\n');
    fprintf(fin,'\t\ttmp = tmp+g_in(i,j)*interp2(M_in.deltaPhi1i,M_in.deltaPhi1j,M_in.value,phi1i,phi1j);\n');
    fprintf(fin,'\t\ttmp = tmp+g_ex(i,j)*interp2(M_ex.deltaPhi1i,M_ex.deltaPhi1j,M_ex.value,phi1i,phi1j);\n');
    fprintf(fin,'\t\ttmp = tmp+g_el(i,j)*interp2(M_el.deltaPhi1i,M_el.deltaPhi1j,M_el.value,phi1i,phi1j);\n\n');
%     fprintf(fin,'\t\ttmp = tmp+g_in(i,j)*lininterp2(M_in.deltaPhi1j,M_in.deltaPhi1i,M_in.value,phi1j,phi1i);\n');
%     fprintf(fin,'\t\ttmp = tmp+g_ex(i,j)*lininterp2(M_ex.deltaPhi1j,M_ex.deltaPhi1i,M_ex.value,phi1j,phi1i);\n');
%     fprintf(fin,'\t\ttmp = tmp+g_el(i,j)*lininterp2(M_el.deltaPhi1j,M_el.deltaPhi1i,M_el.value,phi1j,phi1i);\n\n');

    
    fprintf(fin,'\tend\n');
    fprintf(fin,'\txdot(i-1) = constantContrib-tmp;\n');
    fprintf(fin,'end\n');
    
    fclose(fin);
    
    
    %%%%%%% C CODE %%%%%%%%
elseif strcmp(code,'c')||strcmp(code,'auto')
    
    %%%%%%% HEADER FILE %%%%%%%%
    fin = fopen([folder,filename,'.h'],'w');
    nphi = numel(w_in.deltaPhi);
    pp = spline(w_in.deltaPhi,w_in.value);
    
    delta = diff(w_in.deltaPhi);
    delta = delta(1);
    fprintf(fin,'#define delta %ff\n',delta);
    fprintf(fin,'#define N %d\n\n',object.getN);

    
    fprintf(fin,'double w_in(double phi);\n');
    fprintf(fin,'double w_ex(double phi);\n');
    fprintf(fin,'double w_el(double phi);\n\n');
    fprintf(fin,'double BilinearInterpolation(double q11, double q12, double q21, double q22, double x1, double x2, double y1, double y2, double x, double y);\n\n');
    fprintf(fin,'double M_in(double phii,double phij);\n');
    fprintf(fin,'double M_ex(double phii,double phij);\n');
    fprintf(fin,'double M_el(double phii,double phij);\n\n');

    fprintf(fin,'void %s(double t,double *x,double *xdot,double g_in[N][N],double g_ex[N][N],double g_el[N][N]);\n\n',filename);
   
    
    fclose(fin);
    
    
    %%%%%%% VECTORIAL FIELD FILE %%%%%%%%
    
    fin = fopen([folder,filename,'.c'],'w');
    fprintf(fin,['#include "',filename,'.h"\n']);
    fprintf(fin,['#include <stdlib.h>\n']);
	if strcmp(code,'auto')
    		fprintf(fin,'#include "auto_f2c.h"\n\n');
	else
		fprintf(fin,'#include <stdio.h>\n\n');
    end
    
         fprintf(fin,'const double phi_vect[%d] = {',nphi);
    for i=1:nphi-1
        fprintf(fin,'%.64e,',w_in.deltaPhi(i));
    end
    fprintf(fin,'%.64e};\n\n',w_in.deltaPhi(end));
    
    npp = size(pp.coefs,1);
    fprintf(fin,'const double w_in_vect[%d][4] = {{',npp);
    for i=1:npp-1
        for j=1:3
            fprintf(fin,'%.64e,',pp.coefs(i,j));
        end
        fprintf(fin,'%.64e},\n\t\t{',pp.coefs(i,end));
    end
    for j=1:3
        fprintf(fin,'%.64e,',pp.coefs(end,j));
    end
    fprintf(fin,'%.64e}};\n\n',pp.coefs(end,end));
    
    
    nphi = numel(w_ex.deltaPhi);
    pp = spline(w_ex.deltaPhi,w_ex.value);
    
    npp = size(pp.coefs,1);
    fprintf(fin,'const double w_ex_vect[%d][4] = {{',npp);
    for i=1:npp-1
        for j=1:3
            fprintf(fin,'%.64e,',pp.coefs(i,j));
        end
        fprintf(fin,'%.64e},\n\t\t{',pp.coefs(i,end));
    end
    for j=1:3
        fprintf(fin,'%.64e,',pp.coefs(end,j));
    end
    fprintf(fin,'%.64e}};\n\n',pp.coefs(end,end));
    
    nphi = numel(w_el.deltaPhi);
    pp = spline(w_el.deltaPhi,w_el.value);
    
    npp = size(pp.coefs,1);
    fprintf(fin,'const double w_el_vect[%d][4] = {{',npp);
    for i=1:npp-1
        for j=1:3
            fprintf(fin,'%.64e,',pp.coefs(i,j));
        end
        fprintf(fin,'%.64e},\n\t\t{',pp.coefs(i,end));
    end
    for j=1:3
        fprintf(fin,'%.64e,',pp.coefs(end,j));
    end
    fprintf(fin,'%.64e}};\n\n',pp.coefs(end,end));
    
    nM = size(M_in.value);
    nPhii = linspace(0,1,nM(1));
    nPhij = linspace(0,1,nM(2));
    nphi = numel(w_el.deltaPhi);
    
    fprintf(fin,'const double phii_vect[%d] = {',nM(1));
    for i=1:nM(1)-1
        fprintf(fin,'%.64e,',nPhii(i));
    end
    fprintf(fin,'%.64e};\n',nPhii(end));
    
    fprintf(fin,'const double phij_vect[%d] = {',nM(2));
    for i=1:nM(2)-1
        fprintf(fin,'%.64e,',nPhij(i));
    end
    fprintf(fin,'%.64e};\n',nPhij(end));
    %
    %
    %     writeBicubicInterp(M_in,fin,'M_in_vect');
    %     writeBicubicInterp(M_ex,fin,'M_ex_vect');
    %     writeBicubicInterp(M_el,fin,'M_el_vect');
    
    
    M_in.value = M_in.value';
    fprintf(fin,'const double M_in_vect[%d][%d] = {{',nM(1),nM(2));
    for i=1:nM(1)-1
        for j=1:nM(2)-1
            fprintf(fin,'%.64e,',M_in.value(i,j));
        end
        fprintf(fin,'%.64e},\n\t\t{',M_in.value(i,end));
    end
    for j=1:nM(1)-1
        fprintf(fin,'%.64e,',M_in.value(end,j));
    end
    fprintf(fin,'%.64e}};\n\n',M_in.value(end,end));
    
    nM = size(M_ex.value);
    nPhii = linspace(0,1,nM(1));
    nPhij = linspace(0,1,nM(2));
    nphi = numel(w_el.deltaPhi);
    
    M_ex.value = M_ex.value';
    fprintf(fin,'const double M_ex_vect[%d][%d] = {{',nM(1),nM(2));
    for i=1:nM(1)-1
        for j=1:nM(2)-1
            fprintf(fin,'%.64e,',M_ex.value(i,j));
        end
        fprintf(fin,'%.64e},\n\t\t{',M_ex.value(i,end));
    end
    for j=1:nM(1)-1
        fprintf(fin,'%.64e,',M_ex.value(end,j));
    end
    fprintf(fin,'%.64e}};\n\n',M_ex.value(end,end));
    
    nM = size(M_el.value);
    nPhii = linspace(0,1,nM(1));
    nPhij = linspace(0,1,nM(2));
    nphi = numel(w_el.deltaPhi);
    
    
    M_el.value = M_el.value';
    fprintf(fin,'const double M_el_vect[%d][%d] = {{',nM(1),nM(2));
    for i=1:nM(1)-1
        for j=1:nM(2)-1
            fprintf(fin,'%.64e,',M_el.value(i,j));
        end
        fprintf(fin,'%.64e},\n\t\t{',M_el.value(i,end));
    end
    for j=1:nM(1)-1
        fprintf(fin,'%.64e,',M_el.value(end,j));
    end
    fprintf(fin,'%.64e}};\n\n',M_el.value(end,end));
    
    % function w_in
    fprintf(fin,'\tdouble w_in(double phi)\n\t{\n');
    
    
    fprintf(fin,'\t\tint i,ii;\n');
    fprintf(fin,'\t\tdouble tmp,mul;\n\n');
    
    fprintf(fin,'\t\tfor(i=0;i<%d;i++)\n',nphi-1);
    fprintf(fin,'\t\t{\n');
    fprintf(fin,'\t\t\tif(phi < phi_vect[i])\n');
    fprintf(fin,'\t\t\t\tbreak;\n');
    fprintf(fin,'\t\t}\n');
    fprintf(fin,'\t\tii = i-1;\n');
    fprintf(fin,'\t\tphi = phi-phi_vect[ii];\n');
    fprintf(fin,'\t\ttmp = w_in_vect[ii][0];\n');
    fprintf(fin,'\t\tfor(i=1;i<4;i++)\n');
    fprintf(fin,'\t\t{\n');
    fprintf(fin,'\t\t\ttmp = phi*tmp+w_in_vect[ii][i];\n');
    fprintf(fin,'\t\t}\n');
    fprintf(fin,'\t\treturn tmp;\n\t}\n\n');
    
    % function w_ex
    fprintf(fin,'\tdouble w_ex(double phi)\n\t{\n');
    fprintf(fin,'\t\tint i,ii;\n');
    fprintf(fin,'\t\tdouble tmp,mul;\n\n');
    
    
    fprintf(fin,'\t\tfor(i=0;i<%d;i++)\n',nphi-1);
    fprintf(fin,'\t\t{\n');
    fprintf(fin,'\t\t\tif(phi < phi_vect[i])\n');
    fprintf(fin,'\t\t\t\tbreak;\n');
    fprintf(fin,'\t\t}\n');
    fprintf(fin,'\t\tii = i-1;\n');
    fprintf(fin,'\t\tphi = phi-phi_vect[ii];\n');
    fprintf(fin,'\t\ttmp = w_ex_vect[ii][0];\n');
    fprintf(fin,'\t\tfor(i=1;i<4;i++)\n');
    fprintf(fin,'\t\t{\n');
    fprintf(fin,'\t\t\ttmp = phi*tmp+w_ex_vect[ii][i];\n');
    fprintf(fin,'\t\t}\n');
    fprintf(fin,'\t\treturn tmp;\n\t}\n\n');
    
    % function w_el
    fprintf(fin,'\tdouble w_el(double phi)\n\t{\n');
    fprintf(fin,'\t\tint i,ii;\n');
    fprintf(fin,'\t\tdouble tmp,mul;\n\n');
    
    fprintf(fin,'\t\tfor(i=0;i<%d;i++)\n',nphi-1);
    fprintf(fin,'\t\t{\n');
    fprintf(fin,'\t\t\tif(phi < phi_vect[i])\n');
    fprintf(fin,'\t\t\t\tbreak;\n');
    fprintf(fin,'\t\t}\n');
    fprintf(fin,'\t\tii = i-1;\n');
    fprintf(fin,'\t\tphi = phi-phi_vect[ii];\n');
    fprintf(fin,'\t\ttmp = w_el_vect[ii][0];\n');
    fprintf(fin,'\t\tfor(i=1;i<4;i++)\n');
    fprintf(fin,'\t\t{\n');
    fprintf(fin,'\t\t\ttmp = phi*tmp+w_el_vect[ii][i];\n');
    fprintf(fin,'\t\t}\n');
    fprintf(fin,'\t\treturn tmp;\n\t}\n\n');
    
    
    % Bilinear interpolation
    fprintf(fin,'\tdouble BilinearInterpolation(double q11, double q12, double q21, double q22, double x1, double x2, double y1, double y2, double x, double y)\n');
    fprintf(fin,'\t{\n');
    fprintf(fin,'\t\tdouble x2x1, y2y1, x2x, y2y, yy1, xx1;\n');
    fprintf(fin,'\t\tx2x1 = x2 - x1;\n');
    fprintf(fin,'\t\ty2y1 = y2 - y1;\n');
    fprintf(fin,'\t\tx2x = x2 - x;\n');
    fprintf(fin,'\t\ty2y = y2 - y;\n');
    fprintf(fin,'\t\tyy1 = y - y1;\n');
    fprintf(fin,'\t\txx1 = x - x1;\n\n');
    
    fprintf(fin,'\t\treturn 1.0 / (x2x1 * y2y1) * (q11 * x2x * y2y + q21 * xx1 * y2y + q12 * x2x * yy1 +q22 * xx1 * yy1);\n');
    fprintf(fin,'\t}\n\n');
    
    % function M_in
    fprintf(fin,'\tdouble M_in(double phii,double phij)\n\t{\n');
    
    fprintf(fin,'\t\tint i,j,jj,ii;\n\t\tdouble tmp;\n');
    
    
    fprintf(fin,'\t\tfor(i=0;i<%d;i++)\n\t\t{\n\t\t\tif(phii < phii_vect[i])\n\t\t\t\tbreak;\n\t\t}\n\t\tii = i-1;\n',nM(1)-1);
    fprintf(fin,'\t\tfor(i=0;i<%d;i++)\n\t\t{\n\t\t\tif(phij < phij_vect[i])\n\t\t\t\tbreak;\n\t\t}\n\t\tjj = i-1;\n',nM(1)-1);
%     fprintf(fin,'\t\tphii = (phii-phii_vect[ii])/delta;\n\t\tphij = (phij-phij_vect[jj])/delta;\n');
    %     fprintf(fin,'\t\ttmp = 0;\n\t\tfor(i=0;i<4;i++)\n\t\t\tfor(j=0;j<4;j++)\n\t\t\t\ttmp += (M_in_vect[ii][jj][i][j])*pow(phii,i)*pow(phij,j);\n\t\treturn tmp;\n\t}\n');
    
    fprintf(fin,'\t\treturn BilinearInterpolation(M_in_vect[ii][jj],M_in_vect[ii][jj+1],M_in_vect[ii+1][jj],M_in_vect[ii+1][jj+1],phii_vect[ii],phii_vect[ii+1],phij_vect[jj],phij_vect[jj+1],phii,phij);\n\t}\n\n');
    
    % function M_ex
    fprintf(fin,'\tdouble M_ex(double phii,double phij)\n\t{\n');
    fprintf(fin,'\t\tint i,j,jj,ii;\n\t\tdouble tmp;\n');
    
    
    fprintf(fin,'\t\tfor(i=0;i<%d;i++)\n\t\t{\n\t\t\tif(phii < phii_vect[i])\n\t\t\t\tbreak;\n\t\t}\n\t\tii = i-1;\n',nM(1)-1);
    fprintf(fin,'\t\tfor(i=0;i<%d;i++)\n\t\t{\n\t\t\tif(phij < phij_vect[i])\n\t\t\t\tbreak;\n\t\t}\n\t\tjj = i-1;\n',nM(1)-1);
%     fprintf(fin,'\t\tphii = (phii-phii_vect[ii])/delta;\n\t\tphij = (phij-phij_vect[jj])/delta;\n');
    %     fprintf(fin,'\t\ttmp = 0;\n\t\tfor(i=0;i<4;i++)\n\t\t\tfor(j=0;j<4;j++)\n\t\t\t\ttmp += (M_ex_vect[ii][jj][i][j])*pow(phii,i)*pow(phij,j);\n\t\treturn tmp;\n\t}\n');
    
    fprintf(fin,'\t\treturn BilinearInterpolation(M_ex_vect[ii][jj],M_ex_vect[ii][jj+1],M_ex_vect[ii+1][jj],M_ex_vect[ii+1][jj+1],phii_vect[ii],phii_vect[ii+1],phij_vect[jj],phij_vect[jj+1],phii,phij);\n\t}\n\n');
    
    % function M_el
    fprintf(fin,'\tdouble M_el(double phii,double phij)\n\t{\n');
    fprintf(fin,'\t\tint i,j,jj,ii;\n\t\tdouble tmp;\n');
    
    
    fprintf(fin,'\t\tfor(i=0;i<%d;i++)\n\t\t{\n\t\t\tif(phii < phii_vect[i])\n\t\t\t\tbreak;\n\t\t}\n\t\tii = i-1;\n',nM(1)-1);
    fprintf(fin,'\t\tfor(i=0;i<%d;i++)\n\t\t{\n\t\t\tif(phij < phij_vect[i])\n\t\t\t\tbreak;\n\t\t}\n\t\tjj = i-1;\n',nM(1)-1);
%     fprintf(fin,'\t\tphii = (phii-phii_vect[ii])/delta;\n\t\tphij = (phij-phij_vect[jj])/delta;\n');
    %     fprintf(fin,'\t\ttmp = 0;\n\t\tfor(i=0;i<4;i++)\n\t\t\tfor(j=0;j<4;j++)\n\t\t\t\ttmp += (M_el_vect[ii][jj][i][j])*pow(phii,i)*pow(phij,j);\n\t\treturn tmp;\n\t}\n');
    
    fprintf(fin,'\t\treturn BilinearInterpolation(M_el_vect[ii][jj],M_el_vect[ii][jj+1],M_el_vect[ii+1][jj],M_el_vect[ii+1][jj+1],phii_vect[ii],phii_vect[ii+1],phij_vect[jj],phij_vect[jj+1],phii,phij);\n\t}\n\n');
    
     fprintf(fin,'void %s(double t,double *x,double *xdot,double g_in[N][N],double g_ex[N][N],double g_el[N][N])\n{\n',filename);
   
    
    
    fprintf(fin,'\tint i,j;\n');
    fprintf(fin,'\tdouble constantContrib,tmp;\n');
    fprintf(fin,'\tdouble *bigDeltaPhi;\n\n');
    fprintf(fin,'\tbigDeltaPhi = (double *)malloc(sizeof(double)*(N));\n\n');
    fprintf(fin,'\tbigDeltaPhi[0] = 0;\n');
    fprintf(fin,'\tfor(i=0;i<N-1;i++)\n');
    fprintf(fin,'\t{\n');
    fprintf(fin,'\t\twhile(x[i] >= 1)\n');
    fprintf(fin,'\t\t\tx[i]-= 1;\n');
    fprintf(fin,'\t\twhile(x[i] < 0)\n');
    fprintf(fin,'\t\t\tx[i]+= 1;\n');
    fprintf(fin,'\t\tbigDeltaPhi[i+1] = x[i];\n');
    fprintf(fin,'\t}\n\n');
    
    fprintf(fin,'\tconstantContrib = 0;\n');
    fprintf(fin,'\tfor(j=0;j<N;j++)\n');
    fprintf(fin,'\t{\n');
    fprintf(fin,'\t\tconstantContrib = constantContrib+g_in[0][j]*w_in(bigDeltaPhi[j]);\n');
    fprintf(fin,'\t\tconstantContrib = constantContrib+g_ex[0][j]*w_ex(bigDeltaPhi[j]);\n');
    fprintf(fin,'\t\tconstantContrib = constantContrib+g_el[0][j]*w_el(bigDeltaPhi[j]);\n');
    fprintf(fin,'\t}\n\n');
    
    
    fprintf(fin,'\tfor(i=0;i<N-1;i++)\n');
    fprintf(fin,'\t{\n');
    fprintf(fin,'\t\t\ttmp=0;\n');
    fprintf(fin,'\t\tfor(j=0;j<N;j++)\n');
    fprintf(fin,'\t\t{\n');
    fprintf(fin,'\t\t\ttmp = tmp+g_in[i+1][j]*M_in(bigDeltaPhi[i+1],bigDeltaPhi[j]);\n');
    fprintf(fin,'\t\t\ttmp = tmp+g_ex[i+1][j]*M_ex(bigDeltaPhi[i+1],bigDeltaPhi[j]);\n');
    fprintf(fin,'\t\t\ttmp = tmp+g_el[i+1][j]*M_el(bigDeltaPhi[i+1],bigDeltaPhi[j]);\n');
    fprintf(fin,'\t\t}\n');
    fprintf(fin,'\txdot[i] = constantContrib-tmp;\n');
    fprintf(fin,'\t}\n');
    
    fprintf(fin,'}');
    
    %%%%%%% AUTO INTERFACE FILE (if needed) %%%%%%%%
    
    if strcmp(code,'auto')
        fprintf(fin,'\n\n\n');
        N = object.getN;
        g_in = object.get_g_in;
        g_ex = object.get_g_ex;
        g_el = object.get_g_el;
        fprintf(fin,'int func(integer ndim, const doublereal *u, const integer *icp,\n');
        fprintf(fin,'\tconst doublereal *par, integer ijac, doublereal *f, doublereal *dfdu, \n');
        fprintf(fin,'\tdoublereal *dfdp)\n');
        fprintf(fin,'{\n');
        fprintf(fin,'\tint j;\n');
        fprintf(fin,'\tconst double g_in[N][N] = {');
        for i=1:N
            fprintf(fin,'\t\t{');
            for j=1:N-1
                fprintf(fin,'%e,',g_in(i,j));
            end
            fprintf(fin,'%e}',g_in(i,end));
            if i ~= N
                fprintf(fin,',\n');
            end
        end
        fprintf(fin,'};\n');
        
        fprintf(fin,'\tconst double g_ex[N][N] = {');
        for i=1:N
            fprintf(fin,'\t\t{');
            for j=1:N-1
                fprintf(fin,'%e,',g_ex(i,j));
            end
            fprintf(fin,'%e}',g_ex(i,end));
            if i ~= N
                fprintf(fin,',\n');
            end
        end
        fprintf(fin,'};\n');
        
        fprintf(fin,'\tconst double g_el[N][N] = {');
        for i=1:N
            fprintf(fin,'\t\t{');
            for j=1:N-1
                fprintf(fin,'%e,',g_el(i,j));
            end
            fprintf(fin,'%e}',g_el(i,end));
            if i ~= N
                fprintf(fin,',\n');
            end
        end
        fprintf(fin,'};\n\tdouble x[N-1],xdot[N-1];\n\n');
        fprintf(fin,'\t/**** Write here how bifurcation parameters change synamptic matrices ****/\n\n\n');
        
        fprintf(fin,'\tfor(j=0;j<N-1;j++)\n\t\tx[j] = u[j];\n\n');
        
        fprintf(fin,['\t',filename,'(0,x,xdot,g_in,g_ex,g_el);\n\n']);
        
        fprintf(fin,'\tfor(j=0;j<N-1;j++)\n\t\tf[j] = xdot[j];\n\n');
        
        fprintf(fin,'\treturn 0;\n}\n\n');
        
        
        fprintf(fin,'int stpnt(integer ndim, doublereal t, doublereal *u, doublereal *par)\n');
        fprintf(fin,'{\n');
        
        fprintf(fin,'\tint j;\n\t/*******Write here initial condition *****/\ndouble deltaPhi0[] = {};\n\n');
        fprintf(fin,'\tfor(j=0;j<N-1;j++)\n\t\tu[j] = deltaPhi0[j];\n\n');
        fprintf(fin,'return 0;\n}\n');
        
        
        fprintf(fin,'int pvls(integer ndim, const doublereal *u, doublereal *par)\n');
        fprintf(fin,'{\n');
        fprintf(fin,'\treturn 0;\n');
        fprintf(fin,'}\n\n');
        
        
        fprintf(fin,'int bcnd(integer ndim, const doublereal *par, const integer *icp,\n');
        fprintf(fin,'\tinteger nbc, const doublereal *u0, const doublereal *u1, integer ijac,\n');
        fprintf(fin,'\tdoublereal *fb, doublereal *dbc)\n');
        fprintf(fin,'{\n');
        fprintf(fin,'\treturn 0;\n');
        fprintf(fin,'}\n\n');
        
        
        fprintf(fin,'int icnd(integer ndim, const doublereal *par, const integer *icp,\n');
        fprintf(fin,'\tinteger nint, const doublereal *u, const doublereal *uold,\n');
        fprintf(fin,'\tconst doublereal *udot, const doublereal *upold, integer ijac,\n');
        fprintf(fin,'\tdoublereal *fi, doublereal *dint)\n');
        fprintf(fin,'{\n');
        fprintf(fin,'\treturn 0;\n');
        fprintf(fin,'}\n\n\n');
        
        
        fprintf(fin,'int fopt(integer ndim, const doublereal *u, const integer *icp,\n');
        fprintf(fin,'\tconst doublereal *par, integer ijac,\n');
        fprintf(fin,'\tdoublereal *fs, doublereal *dfdu, doublereal *dfdp)\n');
        fprintf(fin,'{\n');
        fprintf(fin,'\treturn 0;\n');
        fprintf(fin,'}\n');
        
        
        fin2 = fopen([folder,'c.',filename],'w');
        fprintf(fin2,'unames = {');
        for i=1:N-2
            fprintf(fin2,'%d: ''delta1%d'',',i,i+1);
        end
        fprintf(fin2,'%d: ''delta1%d''}\n',N-1,N);
        fprintf(fin2,'parnames = {}\n');
        fprintf(fin2,'NDIM=   %d, IPS =  1, IRS =   0, ILP =   1\n',N-1);
        fprintf(fin2,'ICP =  [1]\n');
        fprintf(fin2,'NTST=  20, NCOL=   4, IAD =   1, ISP =   0, ISW = 1, IPLT= 0, NBC= 0, NINT=0\n');
        fprintf(fin2,'NMX=  200, NPR=  100, MXBF=  10, IID =   2, ITMX=   10, ITNW= 10, NWTN= 3, JAC= 0\n');
        fprintf(fin2,'EPSL= 1e-07, EPSU = 1e-07, EPSS = 1e-05\n');
        fprintf(fin2,'DS  =  0.01, DSMIN= 0.00001, DSMAX=   0.5, IADS=   1\n');
        fprintf(fin2,'NPAR=  1, THL =  {}, THU =  {}\n');
        fclose(fin2);
        
        fin2 = fopen([folder,filename,'.auto'],'w');
        fprintf(fin2,'import numpy as np\n');
        fprintf(fin2,'import scipy.io as sio\n');
        
        fprintf(fin2,'#***Compute a solution family***\n');
        fprintf(fin2,'equilCont=run(e=''vectorialField'',c=''vectorialField'')\n\n');
        fprintf(fin2,'#and save\nsave(equilCont,''equilCont'')\n\n');
        fprintf(fin2,'#save in matlab format\n');
        fprintf(fin2,'sio.savemat(''equil_cont.mat'', {''gbif'':equilCont[0][''gbif''],\\\n');
        for i=1:N-2
            fprintf(fin2,'''delta1%d'':equilCont[0][''delta1%d''],\\\n',i+1,i+1);
        end
        fprintf(fin2,'''delta1%d'':equilCont[0][''delta1%d'']})\n',N,N);
        
        fclose(fin2);
       
        msgbox(['Warning! In order to compute continuation you must set continuation problem and ',...
                'initial condition in file ',filename,'.c and set the continuation constant in ',...
                'c.',filename,' file.']);
        
    end
    
    fclose(fin);
    % TO DO
    
    %
    % xx = mod(xx,1);
    % xdot = zeros(N-1,1);
    % bigPhi = [0 ;xx(:)];
    %
    % constantContrib = 0;
    % for j=1:N
    %     phij = bigPhi(j);
    %     constantContrib = constantContrib+g_in(1,j)*interp1(deltaPhi1jw,w_in,phij);
    %     constantContrib = constantContrib+g_ex(1,j)*interp1(deltaPhi1jw,w_ex,phij);
    %     constantContrib = constantContrib+g_el(1,j)*interp1(deltaPhi1jw,w_el,phij);
    % end
    %
    % for i=2:N
    %     tmp = 0;
    %     for j=1:N
    %         phi1i = bigPhi(i);
    %         phi1j = bigPhi(j);
    %         tmp = tmp+g_in(i,j)*interp2(deltaPhi1i,deltaPhi1j,M_in,phi1i,phi1j);
    %         tmp = tmp+g_ex(i,j)*interp2(deltaPhi1i,deltaPhi1j,M_ex,phi1i,phi1j);
    %         tmp = tmp+g_el(i,j)*interp2(deltaPhi1i,deltaPhi1j,M_el,phi1i,phi1j);
    %         if isnan(tmp)
    %             disp('');
    %         end
    %     end
    %     xdot(i-1) = -constantContrib+tmp;
    % end
    
else
    error('code must be c, auto or m');
end
end


