function [T,X] = sim(object,Tspan,x0,varargin)
% sim    Simulates the neuron
%
%
% [T,X] = sim(OBJ,T,X0)
% Simulates the neuron model OBJ, starting from initial condition X0 between
% Tspan(1) and Tspan(2). Output T is time and output X is the state evolution
%
% [T,X] = sim(OBJ,T,X0,OPTS)
% A structure OPTS can be provided with % the following fields:
% - integrator: string indicating the solver used to integrate the
%              differential equations. It can be either 'ode45','ode23',
%              'ode113','ode15s','ode23s','ode23t' , 'ode23tb','dde23','odeint',
%               'implicit_eulero' or 'explicit_eulero'.
%               If you choose explicit_eulero the simulation is made in C throught mex
%               file and a supported mex compiler is required. See the
%              documentation of ode45 to get help (doc ode45).
%              Default: 'ode45' for non delayed systems and 'dde23' for delayed systems.
% - integrProp: options to provide to the ODE solver. Type doc odeset to
%               get help. 
%               If you choose explicit_eulero the variable must have a
%               field dt that describe the integration step size.
%               If you choose implicit_eulero the variable must have a
%               field dt that describe the integration step size.
%               If you choose odeint integrator the variable must be an
%               odeset structure.
%               Default: odeset.
% - x0_delayed: delayed initial condition from lower to higher delays. x0_delayed must be a Ndelay x
%               Nstates matrix. default : 0.
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

if ~object.is_delayed
    possibleSolver = {'ode45','ode23','ode113','ode15s','ode23s',...
        'ode23t','ode23tb','explicit_eulero','implicit_eulero','odeint'};
    
    nx = object.getnx;
    
    if numel(x0) ~= nx
        error(['x0 must be a vector with ',num2str(nx),' elements']);
    end
    
    if numel(Tspan) == 1
        Tspan(2) = Tspan(1);
        Tspan(1) = 0;
    elseif numel(Tspan) > 2
        error('Tspan must be a vector with two elements [Tstart, Tend]');
    end
    
    integrator = 'ode45';
    integratorOptions = odeset;
    
    if nargin == 4
        if isfield(varargin{1},'integratorOptions')
            integratorOptions = varargin{1}.integratorOptions;
        end
        if isfield(varargin{1},'integrator')
            integrator = varargin{1}.integrator;
        end
    end
    
    intOk = -1;
    for i=1:numel(possibleSolver)
        if strcmp(integrator,possibleSolver{i})
            intOk = 1;
            break;
        end
    end
    
    if intOk == -1
        error('Choosen integrator not allowed');
    end
    
    
    if strcmp(integrator,'explicit_eulero')
        if ~isfield(integratorOptions,'dt')
            error('Integrator step dt must be provided when using explicit_eulero integrator');
        end
        
        oldFolder = cd;
        
        ii = 0;
        nameFolder = 'tmp';
        
        while(exist(nameFolder,'dir') == 7)
            ii = ii+1;
            nameFolder = ['tmp',num2str(ii)];
        end
        
        mkdir(nameFolder);
        
        
        dt = integratorOptions.dt;
        cd(nameFolder);
        %     object.generateC('vectorField');
        
        % create vectorField.cpp
        fout = fopen('vectorField.cpp','w');
        fprintf(fout,'#include "vectorField.hpp"\n');
        fprintf(fout,'void initVectorField(dynSys **vf)\n{\n');
        fprintf(fout,['\t*vf = new ',object.getCbuilder,';\n}']);
        fclose(fout);
        
        cpth = getcpath();
        
        if ispc
            copyfile([cpth,'explicit_eulero.obj']);
            eval(['mex -silent explicit_eulero.obj vectorField.cpp "-I',cpth,'" "-L',cpth,'" -lCEPAGE -output explicit_eulero']);
        elseif isunix
            copyfile([cpth,'explicit_eulero.o']);
            eval(['mex -silent explicit_eulero.o vectorField.cpp "-I',cpth,'" "-L',cpth,'" -lCEPAGE -output explicit_eulero']);
        else
            error('Unsopported operative system');
        end
        
        
        nStep = (Tspan(2)-Tspan(1))/dt;
        X = explicit_eulero(nx,nStep,dt,x0);
        
        clear explicit_eulero;
        
        cd(oldFolder);
        
        a = 0;
        tic
        while (a == 0) && (toc < 5)
            a = rmdir(nameFolder,'s');
        end
        
        if a ~= 1
            warning(['Cannot delete ',nameFolder,' folder, remove manually!']);
        end
        
        X = X';
        T = linspace(Tspan(1),Tspan(2),nStep);
    elseif strcmp(integrator,'implicit_eulero')
        if ~isfield(integratorOptions,'dt')
            error('Integrator step dt must be provided when using implicit_eulero integrator');
        end
        
        oldFolder = cd;
        
        ii = 0;
        nameFolder = 'tmp';
        
        while(exist(nameFolder,'dir') == 7)
            ii = ii+1;
            nameFolder = ['tmp',num2str(ii)];
        end
        
        mkdir(nameFolder);
        
        
        dt = integratorOptions.dt;
        cd(nameFolder);
        %     object.generateC('vectorField');
        
        % create vectorField.cpp
        fout = fopen('vectorField.cpp','w');
        fprintf(fout,'#include "vectorField.hpp"\n');
        fprintf(fout,'void initVectorField(dynSys **vf)\n{\n');
        fprintf(fout,['\t*vf = new ',object.getCbuilder,';\n}']);
        fclose(fout);
        
        cpth = getcpath();
        
        if ispc
            copyfile([cpth,'implicit_eulero.obj']);
            eval(['mex -silent implicit_eulero.obj vectorField.cpp "-I',cpth,'" "-L',cpth,'" -lCEPAGE -output implicit_eulero']);
        elseif isunix
            copyfile([cpth,'implicit_eulero.o']);
            eval(['mex -silent implicit_eulero.o vectorField.cpp "-I',cpth,'" "-L',cpth,'" -lCEPAGE -output implicit_eulero']);
        else
            error('Unsopported operative system');
        end
        
        
        nStep = (Tspan(2)-Tspan(1))/dt;
        X = implicit_eulero(nx,nStep,dt,x0,1e-9);
        
        clear implicit_eulero;
        
        cd(oldFolder);
        
        a = 0;
        tic
        while (a == 0) && (toc < 5)
            a = rmdir(nameFolder,'s');
        end
        
        if a ~= 1
            warning(['Cannot delete ',nameFolder,' folder, remove manually!']);
        end
        
        X = X';
        T = linspace(Tspan(1),Tspan(2),nStep);
    elseif strcmp(integrator,'odeint')
        
        useBoost = getCEPAGEPar();
        useBoost = useBoost.useBoost;
        
        if ~useBoost
            error('Boost c++ integrator is not installed');
        end
        
        % default values
        dt = 1e-3;
        relTol = 1e-6;
        absTol = 1e-10;
        
        if isfield(integratorOptions,'MaxStep')
            if ~isempty(integratorOptions.MaxStep)
                dt = integratorOptions.MaxStep;
            end
        end
        
        if isfield(integratorOptions,'RelTol')
            if ~isempty(integratorOptions.RelTol)
                relTol = integratorOptions.RelTol;
            end
        end
        
        
        if isfield(integratorOptions,'AbsTol')
            if ~isempty(integratorOptions.AbsTol)
                absTol = integratorOptions.AbsTol;
            end
        end
        
        oldFolder = cd;
        
        ii = 0;
        nameFolder = 'tmp';
        
        while(exist(nameFolder,'dir') == 7)
            ii = ii+1;
            nameFolder = ['tmp',num2str(ii)];
        end
        
        mkdir(nameFolder);
        cd(nameFolder);
        %     object.generateC('vectorField');
        %     movefile('vectorField.c','vectorField.cpp');
        
        % create vectorField.cpp
        fout = fopen('vectorField.cpp','w');
        fprintf(fout,'#include "vectorField.hpp"\n');
        fprintf(fout,'void initVectorField(dynSys **vf)\n{\n');
        fprintf(fout,['\t*vf = new ',object.getCbuilder,';\n}']);
        fclose(fout);
        
        cpth = getcpath();
        
        boostDir = getCEPAGEPar();
        boostDir = boostDir.boostDir;
        
        if ispc
            copyfile([cpth,'odeint.obj']);
            eval(['mex -silent odeint.obj vectorField.cpp "-I',cpth,'" -L"',cpth,'" -lCEPAGE -output odeint_v2']);
        elseif isunix
            copyfile([cpth,'odeint.o']);
            eval(['mex -silent odeint.o vectorField.cpp "-I',cpth,'" -L"',cpth,'" -lCEPAGE -output odeint']);
        else
            error('Unsopported operative system');
        end
        
        
        [T,X] = odeint_v2(nx,Tspan(2)-Tspan(1),dt,x0,absTol,relTol);
        T = T+Tspan(1);
        clear odeint_v2;
        cd(oldFolder);
        
        a = 0;
        tic
        while (a == 0) && (toc < 5)
            a = rmdir(nameFolder,'s');
        end
        
        if a ~= 1
            warning(['Cannot delete ',nameFolder,' folder, remove manually!']);
        end
    else
        %     integratorOptions.Jacobian = @object.getJacobian;
        integratorOptions.Events = @object.getResetConditions;
        
        currentT = Tspan(1);
        
        T = [0];
        X = [zeros(1,object.nx)];
        
        while(currentT < Tspan(2))
            
            command = [integrator,'(@object.getXdot,[currentT Tspan(2)],x0,integratorOptions);'];
            [Ttmp,Xtmp] = eval(command);
            [x0,object] = object.resetStates(Ttmp(end),Xtmp(end,:));
            
            T(end) = [];
            X(end,:) = [];
            
            T = [T;Ttmp];
            X = [X;Xtmp];
            
            currentT = T(end);
            
        end
        
    end
    
else
    
    %% TIME DELAYED
    
    possibleSolver = {'dde23','explicit_eulero','odeint'};
    
    nx = object.getnx;
    
    if numel(x0) ~= nx
        error(['x0 must be a vector with ',num2str(nx),' elements']);
    end
    
    if numel(Tspan) == 1
        Tspan(2) = Tspan(1);
        Tspan(1) = 0;
    elseif numel(Tspan) > 2
        error('Tspan must be a vector with two elements [Tstart, Tend]');
    end
    
    integrator = 'dde23';
    integratorOptions = ddeset;
    
    if nargin == 4
        if isfield(varargin{1},'integratorOptions')
            integratorOptions = varargin{1}.integratorOptions;
        end
        if isfield(varargin{1},'integrator')
            integrator = varargin{1}.integrator;
        end
    end
    
    intOk = -1;
    for i=1:numel(possibleSolver)
        if strcmp(integrator,possibleSolver{i})
            intOk = 1;
            break;
        end
    end
    
    if intOk == -1
        error('Choosen integrator not allowed');
    end
    
    
    nDelay = numel(object.delays);
    if nargin == 4
        if isfield(varargin{1},'x0_delayed')
            x0_del = varargin{1}.x0_delayed;
            
            if(any(size(x0_del) ~= [nDelay,nx]))
                if (all(size(x0_del') == [nDelay,nx]))
                    x0_del = x0_del';
                else
                    x0_del = zeros(nx,nDelay);
                end
            end
            
        else
            x0_del = zeros(nx,nDelay);
        end
        
    else
        x0_del = zeros(nx,nDelay);
    end
    
    
    if strcmp(integrator,'explicit_eulero')
        if ~isfield(integratorOptions,'dt')
            error('Integrator step dt must be provided when using explicit_eulero integrator');
        end
        
       oldFolder = cd;
        
        ii = 0;
        nameFolder = 'tmp';
        
        while(exist(nameFolder,'dir') == 7)
            ii = ii+1;
            nameFolder = ['tmp',num2str(ii)];
        end
        
        mkdir(nameFolder);
        
        
        dt = integratorOptions.dt;
        cd(nameFolder);
        %     object.generateC('vectorField');
        
        % create vectorField.cpp
        fout = fopen('vectorField.cpp','w');
        fprintf(fout,'#include "vectorField.hpp"\n');
        fprintf(fout,'void initVectorField(dynSys **vf)\n{\n');
        fprintf(fout,['\t*vf = new ',object.getCbuilder,';\n}']);
        fclose(fout);
        
        cpth = getcpath();
        
        
        if ispc
            copyfile([cpth,'eulero_delayed.obj']);
            eval(['mex -silent eulero_delayed.obj vectorField.cpp "-I',cpth,'" "-L',cpth,'" -lCEPAGE -output eulero_delayed']);
        elseif isunix
            copyfile([cpth,'eulero_delayed.o']);
            eval(['mex -silent eulero_delayed.o vectorField.cpp "-I',cpth,'" "-L',cpth,'" -lCEPAGE -output eulero_delayed']);
        else
            error('Unsopported operative system');
        end
        
        nStep = (Tspan(2)-Tspan(1))/dt;
        X = eulero_delayed(nx,nStep,dt,x0,x0_del');
        cd(oldFolder);
        
        clear eulero_delayed;
        
        rmdir(nameFolder,'s');
        X = X';
        T = linspace(Tspan(1),Tspan(2),nStep);
    elseif strcmp(integrator,'implicit_eulero')
            error('Implicit eulero integrator not implemented for CPGs with delay');
    elseif strcmp(integrator,'odeint')
             error('odeint integrator not implemented for CPGs with delay');
        useBoost = getCEPAGEPar();
        useBoost = useBoost.useBoost;
        
        if ~useBoost
            error('Boost c++ integrator is not installed');
        end
        
        % default values
        dt = 1e-3;
        relTol = 1e-6;
        absTol = 1e-10;
        
        if isfield(integratorOptions,'MaxStep')
            if ~isempty(integratorOptions.MaxStep)
                dt = integratorOptions.MaxStep;
            end
        end
        
        if isfield(integratorOptions,'RelTol')
            if ~isempty(integratorOptions.RelTol)
                relTol = integratorOptions.RelTol;
            end
        end
        
        
        if isfield(integratorOptions,'AbsTol')
            if ~isempty(integratorOptions.AbsTol)
                absTol = integratorOptions.AbsTol;
            end
        end
        
        
        oldFolder = cd;
        
        ii = 0;
        nameFolder = 'tmp';
        
        while(exist(nameFolder,'dir') == 7)
            ii = ii+1;
            nameFolder = ['tmp',num2str(ii)];
        end
        
        mkdir(nameFolder);
        
        cd(nameFolder);
        
        %     object.generateC('vectorField');
        %     movefile('vectorField.c','vectorField.cpp');
        
        % create vectorField.cpp
        fout = fopen('vectorField.cpp','w');
        fprintf(fout,'#include "vectorField.hpp"\n');
        fprintf(fout,'void initVectorField(dynSys **vf)\n{\n');
        fprintf(fout,['\t*vf = new ',object.getCbuilder,';\n}']);
        fclose(fout);
        
        cpth = getcpath();
        
        
        boostDir = getCEPAGEPar();
        boostDir = boostDir.boostDir;
        
        if ispc
            copyfile([cpth,'odeint_delayed.obj']);
            eval(['mex -silent odeint_delayed.obj vectorField.cpp "-I',cpth,'" "-I',boostDir,'/include" -L"',cpth,'" -lCEPAGE -output odeint_delayed']);
            
        elseif isunix
            copyfile([cpth,'odeint_delayed.o']);
            eval(['mex -silent odeint_delayed.o vectorField.cpp "-I',cpth,'" "-I',boostDir,'/include" -L"',cpth,'" -lCEPAGE -output odeint_delayed']);
        else
            error('Unsopported operative system');
        end
        
        [T,X] = odeint_delayed(nx,Tspan(2)-Tspan(1),dt,x0,x0_del,absTol,relTol);
        T = T+Tspan(1);
        
        clear odeint_delayed;
        
        cd(oldFolder);
        rmdir(nameFolder,'s');
        
    else
        %     integratorOptions.Jacobian = @object.getJacobian;
        integratorOptions.Events = @object.getResetConditions;
        
        currentT = Tspan(1);
        
        delays = object.delays(:);
        
        
        T = [-delays(end:-1:1);0]';
        X = [x0_del(:,end:-1:1),x0(:)]';
        
        
        jumps = [];
        inity0 = [];
        
        while(currentT < Tspan(2))
            
            hyst = @(t)object.hystoryFun(t,T,X);
            integratorOptions.Jumps = jumps;
            integratorOptions.InitialY = inity0;
            command = [integrator,'(@object.getXdot,[delays],hyst, [currentT Tspan(2)],integratorOptions);'];
            sol = eval(command);
            
            Ttmp = sol.x';
            Xtmp = sol.y';
            
            T(end) = [];
            X(end,:) = [];
            
            T = [T(:);Ttmp(:)];
            X = [X;Xtmp];
            
            Zpast = object.hystoryFun(T(end)-delays,T,X);
            
            [x0,object] = object.resetStates(Ttmp(end),Xtmp(end,:),Zpast);
            
            
            
            currentT = T(end);
            
            jumps = [jumps;currentT];
            inity0 = x0;
        end
        
    end
    
    
end

end

