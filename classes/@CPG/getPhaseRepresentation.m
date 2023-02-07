function [Tphi,phi,Xfinale] = getPhaseRepresentation(object,Tspan,CI,varargin )
% getPhaseRepresentation   Gets the phase evolution of the network
%
% [Tphi,phi] = getPhaseRepresentation(object,T,CI )
% Tphi is a cell array with N-1 elements describing the time
% phi is a cell array with N-1 elements describing the evolution of the
% phase difference between neuron 1 and neuron i+1
% The network is simulated for T seconds.
% CI is the initial conditions of the networks and must be a 1xN array. If
% CI is a Ns x N matrix the function compute the phase evolutions of the
% networks starting from the Ns starting conditions.
% If only 1 starting condition is considered Tphi and phi are vectors.
%
% phi = getPhaseRepresentation(object,T,CI,OPTS)
% A structure OPTS can be provided with % the following fields:
% - integrator: string indicating the solver used to integrate the
%              differential equations. It can be either 'ode45','ode23',
%              'ode113','ode15s','ode23s','ode23t' , 'ode23tb', 'odeint',
%               'implicit_eulero' or 'explicit_eulero.
%               If you choose explicit_eulero or odeint the simulation is made in C throught mex
%               file and a supported mex compiler is required.
%               If you choose implicit_eulero the variable must have a
%               field dt that describe the integration step size.
%              Default: 'ode45'.
% - integratorOptions: options to provide to the ODE solver. Type doc odeset to
%               get help. If you choose explicit_eulero or implicit_eulero the variable must have a
%               field dt that describe the integration step size.
%               Default: odeset.
% - Vth :  is the value of the membrane potential that correspond to an
%        event (it must be between max and minimum value of the membrane
%        potential)
%        Default : 0
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

Xfinale = [];

if ~object.is_delayed
    
    possibleSolver = {'ode45','ode23','ode113','ode15s','ode23s','ode23t',...
        'ode23tb','explicit_eulero','odeint','implicit_eulero'};
    
    % nx = object.neurons{1}.getnx;
    N  = object.N;
    Nstati = object.totState;
    
    if size(CI,2) ~= Nstati
        error(['initial conditions must be a vector with ',num2str(Nstati),' columns']);
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
    
    
    Vth = 0;
    stopThreshold = eps;
    if nargin == 4
        if isfield(varargin{1},'Vth')
            Vth = varargin{1}.Vth;
        end
        if isfield(varargin{1},'stopThreshold')
            stopThreshold = varargin{1}.stopThreshold;
        end
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
        
        fout = fopen('vectorField.cpp','w');
        fprintf(fout,'#include "vectorField.hpp"\n');
        fprintf(fout,'void initVectorField(dynSys **vf)\n{\n');
        fprintf(fout,[object.getCbuilder,';\n}']);
        fclose(fout);
        
        cpth = getcpath();
        
        if ispc
            copyfile([cpth,'explicit_euleroEvents.obj']);
            eval(['mex -silent explicit_euleroEvents.obj vectorField.cpp "-I',cpth,'" "-L',cpth,'" -lCEPAGE -output explicit_euleroEvents']);
        elseif isunix
            copyfile([cpth,'explicit_euleroEvents.o']);
            eval(['mex -silent explicit_euleroEvents.o vectorField.cpp "-I',cpth,'" "-L',cpth,'" -lCEPAGE -output explicit_euleroEvents']);
        else
            error('Unsopported operative system');
        end
        
        
        nStep = (Tspan(2)-Tspan(1))/dt;
        
        if size(CI,1) > 1
            Tphi = cell(size(CI,1),1);
            phi = cell(size(CI,1),1);
            parfor kk = 1:size(CI,1)
                [Tphi{kk},phi{kk}] = explicit_euleroEvents(Nstati,nStep,dt,CI(kk,:),Vth,stopThreshold,N);
                phi{kk} = mod(phi{kk},1);
            end
        else
            [Tphi,phi] = explicit_euleroEvents(Nstati,nStep,dt,CI,Vth,stopThreshold,N);
            phi = mod(phi,1);
        end
        
        clear explicit_euleroEvents;
        
        cd(oldFolder);
        a = 0;
        tic
        while (a == 0) && (toc < 5)
            a = rmdir(nameFolder,'s');
        end
        
        if a ~= 1
            warning(['Cannot delete ',nameFolder,' folder, remove manually!']);
        end
        
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
        
        fout = fopen('vectorField.cpp','w');
        fprintf(fout,'#include "vectorField.hpp"\n');
        fprintf(fout,'void initVectorField(dynSys **vf)\n{\n');
        fprintf(fout,[object.getCbuilder,';\n}']);
        fclose(fout);
        
        cpth = getcpath();
        
        if ispc
            copyfile([cpth,'implicit_euleroEvents.obj']);
            eval(['mex -silent implicit_euleroEvents.obj vectorField.cpp "-I',cpth,'" "-L',cpth,'" -lCEPAGE -output implicit_euleroEvents']);
        elseif isunix
            copyfile([cpth,'implicit_euleroEvents.o']);
            eval(['mex -silent implicit_euleroEvents.o vectorField.cpp "-I',cpth,'" "-L',cpth,'" -lCEPAGE -output implicit_euleroEvents']);
        else
            error('Unsopported operative system');
        end
        
        
        nStep = (Tspan(2)-Tspan(1))/dt;
        
        if size(CI,1) > 1
            phi = cell(size(CI,1),1);
            Tphi = cell(size(CI,1),1);
            parfor kk = 1:size(CI,1)
                [Tphi{kk},phi{kk}] = implicit_euleroEvents(Nstati,nStep,dt,CI(kk,:),Vth,stopThreshold,N,1e-9);
                phi{kk} = mod(phi{kk},1);
            end
        else
            [Tphi,phi] = implicit_euleroEvents(Nstati,nStep,dt,CI,Vth,stopThreshold,N,1e-9);
            phi = mod(phi,1);
        end
        
        clear implicit_euleroEvents;
        
        cd(oldFolder);
        a = 0;
        tic
        while (a == 0) && (toc < 5)
            a = rmdir(nameFolder,'s');
        end
        
        if a ~= 1
            warning(['Cannot delete ',nameFolder,' folder, remove manually!']);
        end
        
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
        
        fout = fopen('vectorField.cpp','w');
        fprintf(fout,'#include "vectorField.hpp"\n');
        fprintf(fout,'void initVectorField(dynSys **vf)\n{\n');
        fprintf(fout,[object.getCbuilder,';\n}']);
        fclose(fout);
        
        boostDir = getCEPAGEPar();
        boostDir = boostDir.boostDir;
        
        cpth = getcpath();
        
        if ispc
            copyfile([cpth,'odeintEvents.obj']);
            eval(['mex -silent odeintEvents.obj vectorField.cpp "-I',cpth,'" "-I',boostDir,'/include" -L"',cpth,'" -lCEPAGE -output odeintEvents']);
        elseif isunix
            copyfile([cpth,'odeintEvents.o']);
            eval(['mex -silent odeintEvents.o vectorField.cpp "-I',cpth,'" "-I',boostDir,'/include" -L"',cpth,'" -lCEPAGE -output odeintEvents']);
        else
            error('Unsopported operative system');
        end
        
        if size(CI,1) > 1
            phi = cell(size(CI,1),1);
            Tphi = cell(size(CI,1),1);
            parfor kk = 1:size(CI,1)
                [Tphi{kk},phi{kk}] = odeintEvents(Nstati,Tspan(2),dt,CI(kk,:),Vth,N,absTol,relTol);
                phi{kk} = mod(phi{kk},1);
            end
        else
            [Tphi,phi] = odeintEvents(Nstati,Tspan(2),dt,CI,Vth,N,absTol,relTol);
            phi = mod(phi,1);
        end
        
        clear odeintEvents;
        
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
        
        if ~object.is_continuous
            error(['Pahse representation employing MATLAB integrator '...
                'can be used only for continuous system']);
        end
        
        Te = cell(size(CI,1),1);
        ie = cell(size(CI,1),1);
        eventFun = @(T, Y) object.eventsTh(T, Y, Vth);
        integratorOptions.Events = eventFun;
        
        
        if size(CI,1) > 1
            Tphi = cell(size(CI,1),1);
            phi = cell(size(CI,1),1);
            
            parfor i=1:size(CI,1)
                [Tphi{i},phi{i}] = object.getPhaseRepresentation(Tspan,CI(i,:),struct('Vth',Vth,'integrator',integrator,'integratorOptions',integratorOptions));
            end
            
        else
            x0  = CI;
            
            command = [integrator,'(@object.getXdot,[Tspan(1) Tspan(2)],x0,integratorOptions);'];
            [~,~,Te,~,ie] = eval(command);
            
            
            len = Inf;
            
            for n=1:N
                len = min([len,sum(ie == n)]);
            end
            
            phiTmp = zeros(len-1,N-1);
            
            T1 = Te(ie == 1);
            
            for n=2:N
                Tn = Te(ie == n);
                phiTmp(:,n-1) = mod((Tn(2:len)-T1(2:len))./(T1(2:len)-T1(1:len-1)),1);
            end
            
            Tphi = T1(1:len-1);
            phi = phiTmp;
        end
        
        
        
    end
    
    
else
    
    
    possibleSolver = {'dde23','explicit_eulero','odeint'};
    
    % nx = object.neurons{1}.getnx;
    N  = object.N;
    Nstati = object.totState;
    
    if size(CI,2) ~= Nstati
        error(['initial conditions must be a vector with ',num2str(Nstati),' columns']);
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
    
    
    Vth = 0;
    stopThreshold = eps;
    if nargin == 4
        if isfield(varargin{1},'Vth')
            Vth = varargin{1}.Vth;
        end
        if isfield(varargin{1},'stopThreshold')
            stopThreshold = varargin{1}.stopThreshold;
        end
    end
    
    nDelay = numel(object.delays);
    if nargin == 4
        if isfield(varargin{1},'x0_delayed')
            x0_del = varargin{1}.x0_delayed;
            
            if(any(size(x0_del) ~= [nDelay,object.totState]))
                if (all(size(x0_del') == [nDelay,object.totState]))
                    x0_del = x0_del';
                else
                    x0_del = zeros(object.totState,nDelay);
                end
            end
            
        else
            x0_del = zeros(object.totState,nDelay);
        end
        
    else
        x0_del = zeros(object.totState,nDelay);
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
        
        fout = fopen('vectorField.cpp','w');
        fprintf(fout,'#include "vectorField.hpp"\n');
        fprintf(fout,'void initVectorField(dynSys **vf)\n{\n');
        fprintf(fout,[object.getCbuilder,';\n}']);
        fclose(fout);
        
        cpth = getcpath();
        
        if ispc
            copyfile([cpth,'explicit_euleroEvents_delayed.obj']);
            eval(['mex -silent explicit_euleroEvents_delayed.obj vectorField.cpp "-I',cpth,'" "-L',cpth,'" -lCEPAGE -output explicit_euleroEvents_delayed']);
        elseif isunix
            copyfile([cpth,'explicit_euleroEvents_delayed.o']);
            eval(['mex -silent explicit_euleroEvents_delayed.o vectorField.cpp "-I',cpth,'" "-L',cpth,'" -lCEPAGE -output explicit_euleroEvents_delayed']);
        else
            error('Unsopported operative system');
        end
        
        nStep = (Tspan(2)-Tspan(1))/dt;
        
        if size(CI,1) > 1
            Tphi = cell(size(CI,1),1);
            phi = cell(size(CI,1),1);
            Xfinale = cell(size(CI,1),1);
            parfor kk = 1:size(CI,1)
                [Tphi{kk},phi{kk},Xfinale{kk}] = explicit_euleroEvents_delayed(Nstati,nStep,dt,CI(kk,:),Vth,stopThreshold,N,x0_del);
                phi{kk} = mod(phi{kk},1);
            end
        else
            [Tphi,phi,Xfinale] = explicit_euleroEvents_delayed(Nstati,nStep,dt,CI,Vth,stopThreshold,N,x0_del);
            phi = mod(phi,1);
        end
        
        clear explicit_euleroEvents_delayed;
        
        cd(oldFolder);
        a = 0;
        tic
        while (a == 0) && (toc < 5)
            a = rmdir(nameFolder,'s');
        end
        
        if a ~= 1
            warning(['Cannot delete ',nameFolder,' folder, remove manually!']);
        end
        
    elseif strcmp(integrator,'implicit_eulero')
        error('Implicit eulero integrator not implemented for CPGs with delay');
    elseif strcmp(integrator,'odeint')
        %error('odeint integrator not implemented for CPGs with delay');
        
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
        
        fout = fopen('vectorField.cpp','w');
        fprintf(fout,'#include "vectorField.hpp"\n');
        fprintf(fout,'void initVectorField(dynSys **vf)\n{\n');
        fprintf(fout,[object.getCbuilder,';\n}']);
        fclose(fout);
        
        boostDir = getCEPAGEPar();
        boostDir = boostDir.boostDir;
        
        cpth = getcpath();
        
        if ispc
            copyfile([cpth,'odeint_delayedEvents.obj']);
            eval(['mex -silent odeint_delayedEvents.obj vectorField.cpp "-I',cpth,'" "-I',boostDir,'/include" -L"',cpth,'" -lCEPAGE -output odeint_delayedEvents']);
        elseif isunix
            copyfile([cpth,'odeint_delayedEvents.o']);
            eval(['mex -silent odeint_delayedEvents.o vectorField.cpp "-I',cpth,'" "-I',boostDir,'/include" -L"',cpth,'" -lCEPAGE -output odeint_delayedEvents']);
        else
            error('Unsopported operative system');
        end
        
        if size(CI,1) > 1
            Tphi = cell(size(CI,1),1);
            phi = cell(size(CI,1),1);
            Xfinale = cell(size(CI,1),1);
            parfor kk = 1:size(CI,1)
                [Tphi{kk},phi{kk},Xfinale{kk}] = odeint_delayedEvents(Nstati,Tspan(2),dt,CI(kk,:),Vth,N,x0_del,absTol,relTol);
                phi{kk} = mod(phi{kk},1);
            end
        else
            [Tphi,phi,Xfinale] = odeint_delayedEvents(Nstati,Tspan(2),dt,CI,Vth,N,x0_del,absTol,relTol);
            phi = mod(phi,1);
        end
        
        clear odeint_delayedEvents;
        
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
        
        if ~object.is_continuous
            error(['Pahse representation employing MATLAB integrator '...
                'can be used only for continuous system']);
        end
        
        Te = cell(size(CI,1),1);
        ie = cell(size(CI,1),1);
        eventFun = @(T, Y) object.eventsTh(T, Y, Vth);
        integratorOptions.Events = eventFun;
        
        
        if size(CI,1) > 1
            Tphi = cell(size(CI,1),1);
            phi = cell(size(CI,1),1);
            
            parfor i=1:size(CI,1)
                [Tphi{i},phi{i}] = object.getPhaseRepresentation(Tspan,CI(i,:),struct('x0_delayed',x0_del,'Vth',Vth,'integrator',integrator,'integratorOptions',integratorOptions));
            end
            
        else

            x0  = CI;
            
            command = [integrator,'(@object.getXdot,[Tspan(1) Tspan(2)],x0,integratorOptions);'];
            [~,~,Te,~,ie] = eval(command);
            
            
            len = Inf;
            
            for n=1:N
                len = min([len,sum(ie == n)]);
            end
            
            phiTmp = zeros(len-1,N-1);
            
            T1 = Te(ie == 1);
            
            for n=2:N
                Tn = Te(ie == n);
                phiTmp(:,n-1) = mod((Tn(2:len)-T1(2:len))./(T1(2:len)-T1(1:len-1)),1);
            end
            
            
            phi = phiTmp;
            Tphi = T1(1:len-1);
        end
        
        
        
    end
    
    
end

