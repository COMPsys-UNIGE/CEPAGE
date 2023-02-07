classdef CPG
    
    % CPG   Generic neuron network model
    %
    % This object represents a generic Central Pattern Generator model.
    % The network is characterized by a set of states (x)
    %
    % Each neuron i-th of the network is influenced by neuron j-th by mean
    % the synapsis current:
    %
    %   $$ Isyn = \sum_j g_{in}(i,j)(V_i-EsynIn) inhActivation(V_j) +
    %   g_{ex}(i,j)(V_i-EsynEx)  inhActivation(V_j)+
    %   g_{el}(V_j-V_i) $$
    %
    % OBJ = CPG(N)
    % Builds a CPG composed of N cells object OBJ with all parameters equal to 0.
    %
    %  OBJ = CPG(N,neurons,g_in,g_ex,g_el,Esyn_In,Esyn_Ex,inhActivation,excActivation)
    %  Builds a CPG object OBJ with user assigned parameters value.
    %   - neurons could is a cell array containing N neuron_model objects, each one
    %   representing the cells of the CPG. If neurons is a scalar, all the 
    %   cells within the CPG are described by the same model.
    %   - g_in,g_ex and g_el are N x N matrices containing synapses weights
    %   - Esyn_In (Esyn_Ex) is the chemical inhibitory (excitatory) 
    %   synapses reverse potential.
    %   - inhActivation (excActivation) is the chemical inhibitory (excitatory) 
    %   activation model. It could be a N x N cell array of FTM_synapses
    %   objects or a single FTM_synapses model.
    %
    %   CPG methods:
    %   getPhaseRepresentation - Gets the phase evolution of the CPG%
    %   getCIfromPhi - Get the initial condition of the CPG
    %   getApproxPhaseRepresentation - Get the phase difference evolution
    %                                  of the network using standard
    %                                  phase reduction method
    %   getXdot - computes the derivative of the state
    %   sim - simulates the neuron.
    %   simplot - simulates the neuron system and plots time evolution of
    %             states.
    %   disp - displays some information about the HR_model object
    
    %   getN - Gets the number of neurons in the CPG
    
    %
    % See also HH_model, neuron_model
    %
    % Contributors:
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
    
    %Properties
    
    properties (Access = protected)
        N = 0; % number of neurons in the network
        neurons = []; % Array containing neuron object
        g_in = []; % matrix describing chimical inibhitory synapsis link
        g_ex = [];  % matrix describing chimical excitatory synapsis link
        g_el = [];  % matrix describing electrical synapsis link
        EsynIn = 0; % Chimical inhibitory synapsis reverse potential
        EsynEx = 0;  % Chimical excitatory synapsis reverse potential
        inhActivation = {};
        excActivation = {};
        
        
        delays = [];
        
        delayIndexNeur = {};
        delayInhSyn = {};
        delayExcSyn = {};
        
        totState = 0;
        incrementalIndexState = [];
        
    end
    
    %methods
    methods
        
        %Constructor
        function object = CPG(N,varargin)
            
            object.N = N;
            object.neurons = cell(N,1);
            object.g_in = zeros(N,N);
            object.g_ex = zeros(N,N);
            object.g_el = zeros(N,N);
            object.inhActivation = cell(N,N);
            object.excActivation = cell(N,N);
            
            object.incrementalIndexState = zeros(N+2*N*N+1,1);
            object.incrementalIndexState(1) = 1;
            if nargin >= 2
                %set neurons model
                nModel = varargin{1};
                if isscalar(nModel) && isa(nModel,'neuron_model')
                    
                    for i=1:N
                        object.neurons{i} = nModel;
                        object.incrementalIndexState(i+1) = object.incrementalIndexState(i)+nModel.getnx;
                    end
                elseif numel(nModel) == N && iscell(nModel)
                    for i=1:N
                        if ~isa(nModel{i},'neuron_model')
                            error('Error in object input');
                        end
                        object.neurons{i} = nModel{i};
                        object.incrementalIndexState(i+1) = object.incrementalIndexState(i)+nModel{i}.getnx;
                    end
                else
                    error('Error in object input');
                end
            end
            
            if nargin == 9
                %set synaps strenght
                object.g_in = varargin{2};
                object.g_ex = varargin{3};
                object.g_el = varargin{4};
                object.EsynIn = varargin{5};
                object.EsynEx = varargin{6};
                
                act = varargin{7};
                
                ii = 2;
                
                if isscalar(act) && isa(act,'synapse_model')
                    for i=1:N
                        for j=1:N
                            object.inhActivation{i,j} = act;
                            object.incrementalIndexState(N+ii) = object.incrementalIndexState(N+ii-1)+act.getnx;
                            ii = ii+1;
                        end
                    end
                elseif all(size(act) == [N,N]) && iscell(act)
                    for i=1:N
                        for j=1:N
                            if ~isa(act{i,j},'synapse_model')
                                error('Error in object input');
                            end
                            
                            object.inhActivation{i,j} = act{i,j};
                            object.incrementalIndexState(N+ii) = object.incrementalIndexState(N+ii-1)+act{i,j}.getnx;
                            ii = ii+1;
                        end
                    end
                else
                    error('Error in object input');
                end
                
                
                act = varargin{8};
                
                if isscalar(act) && isa(act,'synapse_model')
                    for i=1:N
                        for j=1:N
                            object.excActivation{i,j} = act;
                            object.incrementalIndexState(N+ii) = object.incrementalIndexState(N+ii-1)+act.getnx;
                            ii = ii+1;
                        end
                    end
                elseif all(size(act) == [N,N]) && iscell(act)
                    for i=1:N
                        for j=1:N
                            if ~isa(act{i,j},'synapse_model')
                                error('Error in object input');
                            end
                            
                            object.excActivation{i,j} = act{i,j};
                            object.incrementalIndexState(N+ii) = object.incrementalIndexState(N+ii-1)+act{i,j}.getnx;
                            ii = ii+1;
                        end
                    end
                else
                    error('Error in object input');
                end
                
                object.totState = object.incrementalIndexState(end)-1;
                
                % Set up delays
                
                delays_  = [];
                
                for i=1:N
                    delays_ = [delays_;object.neurons{i}.getDelays()];
                end
                
                for i=1:N
                    for j=1:N
                        delays_ = [delays_;object.inhActivation{i,j}.getDelays()];
                    end
                end
                
                for i=1:N
                    for j=1:N
                        delays_ = [delays_;object.excActivation{i,j}.getDelays()];
                    end
                end
                
                
                object.delays = sort(unique(delays_));
                
                
                object.delayIndexNeur = cell(N,1);
                object.delayInhSyn = cell(N,N);
                object.delayExcSyn = cell(N,N);
                
                for i=1:N
                    object.delayIndexNeur{i} = find(ismember(object.delays, object.neurons{i}.getDelays()));
                end
                
                for i=1:N
                    for j=1:N
                        object.delayInhSyn{i,j} = find(ismember(object.delays, object.inhActivation{i,j}.getDelays()));
                    end
                end
                
                for i=1:N
                    for j=1:N
                        object.delayExcSyn{i,j} = find(ismember(object.delays, object.excActivation{i,j}.getDelays()));
                    end
                end
                
                
                
                
            end
        end
        
        x_dot = getXdot(object,t,x,varargin);
        J = getJacobian(object,t,x);
        str = getCbuilder(object);
        plot(object);
        disp(object);
        CI = getCIfromPhi(object,deltaPhi,Ttrans);
        [Tphi,phi] = getPhaseRepresentationFromTrack(object,T,X,Vth)
        [Tphi,phi,Xfinale] = getPhaseRepresentation(object,Tspan,CI,varargin )
        [TphiApprox,phiApprox] = getApproxPhaseRepresentation(object,PRC,limitCycle,g_in,g_ex,g_el,Tspan,phi0,varargin)
        [deltaPhiMatrix, deltaPhiDot] = computeApproxVectorField(object,PRC,orbit,varargin)
        writeApproxVectorField(w_in,w_ex,w_el,M_in,M_ex,M_el,varargin);
        % get methods
        N = getN(object);
        
        neur = getNeuron(object,varargin)
        
        inhAct = get_inhAct(object);
        excAct = get_excAct(object);
        
        g_in = get_g_in(object);
        g_ex = get_g_ex(object);
        g_el = get_g_el(object);
        % set methods
        object = set_g_in(object,g_in);
        object = set_g_ex(object,g_ex);
        object = set_g_el(object,g_el);
        [w_in,w_ex,w_el,M_in,M_ex,M_el] = computeMw(object,PRC,limitCycle,varargin);
        
        [position,isterminal,direction] = getResetConditions(object,t,y,varargin);
        [xreset,object] = resetStates(object,t,x,ie,varargin);
        cont = is_continuous(object)
    end
    
    
    
    %     methods (Access = private)
    %         [value,isterminal,direction] = eventsTh(object,t,y,varargin);
    %     end
end
