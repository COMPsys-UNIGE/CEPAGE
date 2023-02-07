classdef synapse_model
    
    % synapse_model Generic synapse model
    %
    % This object represents a generic synapse model.
    % The synapse is characterized by an activation function f(V_j,x) and
    % possibly by a state vector x. V_j is the pre-synaptic neuron membrane
    % potential
    %
    % This is an abstract class which cannot be istantiated.
    %
    %
    % synapse_model methods:
    %   getStateNames - gets the mnemonical name of the synapse states
    %   getnx - Gets the number of state varible nx
    %
    %   getActivation - Compute the activation function
    %   getXdot       - Compute the synapse states derivative
    %
    %   isContinuous - Report if the model is time continuous
    %   is_delayed - Report if the model is time delayed
    %   getDelays - Gets the delays
    %
    %   resetStates - reset the state variables
    %   getResetConditions - get the reset condition
    %
    %   getCbuilder - Generates a string of C code for the computation 
    %                 of model vector field
    %
    % See also FTM_synapse_model
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
        nx = 0 % number of state variable
        xnames = cell(0); % States name
        modelName = '';
        isContinuous = false;
        delays = [];
    end
    
    methods (Abstract)
        % Other methods
        act = getActivation(object,Vpre,varargin);  
        JA = getDA(object,x,Vpre,varargin);    
        J = getJacobian(object,t,x,Vpre);
        x_dot = getXdot(object,t,x,varargin);
        str = getCbuilder(object);
    end
    

    
    methods
        %constructor
        function object = synapse_model(varargin)
            % abstract class empty cconstructor
        end
        
        
        nx = getnx(object);
        names = getStateNames(object);
        [position,isterminal,direction] = getResetConditions(object,t,x,Vpre)
        [xreset,object] = resetStates(object,t,x,Vpre)
        cont = is_continuous(object);
        del = is_delayed(object);
        del = getDelays(object)
    end
    
end
