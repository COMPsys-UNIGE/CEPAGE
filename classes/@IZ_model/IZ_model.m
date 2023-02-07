classdef IZ_model < neuron_model
    
    % IZ_model   Neuron represented by the Izhikevich model
    %
    % The Neuron is modeld as follow:
    %   _
    %  |
    %  | dv/dt = 0.04*v^2+5*v+140-u+I-gL(v-EL)-I_{ext}
    % <
    %  | du/dt = a*(b*v-u);
    %  |_-
    %
    % with the reset condition
    %                 _
    %                |  v = c;
    % if v > 30     <
    %                |_ u = u+d;
    %
    % Ther model has two state variable; v represent the membrane potential.
    % a,b,c,d,I,gL,EL are parameters that allows to change neuron behaviours.
    %
    % OBJ = IZ_model()
    % Builds an IZ_model object OBJ with all parameters equal to 0.
    %
    % OBJ = IZ_model(a,b,c,d,I,gL,EL)
    % Builds an IZ_model object OBJ with user assigned parameters value
    %
    %
    %   IZ_model methods:
    %   getXdot - omputes the derivative of the state
    %   sim - simulates the neuron.
    %   simplot - simulates the neuron system and plots time evolution of
    %             states.
    %   disp - displays some information about the FN_relaxation_model object
    %   get_I - gets the value of parameter I
    %   set_I - sets the value of parameter I
    %   getJacobian -Computes the Jacobian of the vector field in point x
    %   generateC - Generates C files for the computation of model vector field
    %
    % The IZ_model object is derived from neuron_model and
    % inherits all its methods.
    %
    % See also FN_relaxation_model, neuron_model
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
    
    %properties
    properties(Access = protected)
        a = 0;
        b = 0;
        c = 0;
        d = 0;
        I = 0;
        gL = 0;
        EL = 0;
    end
    
    
    methods
        function obj = IZ_model(varargin)
            obj.nx = 2;
            obj.xnames = {'v','u'};
            obj.modelName = 'Morris-Lecar';
            obj.isContinuous = false;
            if nargin == 7
                obj.a = varargin{1};
                obj.b = varargin{2};
                obj.c = varargin{3};
                obj.d = varargin{4};
                obj.I = varargin{5};
                obj.gL = varargin{6};
                obj.EL = varargin{7};
            end
        end
        
        % Other methods
        x_dot = getXdot(object,t,x,varargin);
        disp(object);
        str = getCbuilder(object);
        dfdi = getdfdi(object,x); % compute the differential of f with respect to Isyn
        
        [position,isterminal,direction] = getResetConditions(object,t,y,varargin);
        [xreset,object] = resetStates(object,t,x,varargin);
        
        % get methods
        gL = get_gL(object);
        EL = get_EL(object);
        I = get_I(object);
        
        % set methods
        object = set_gL(object,gL);
        object = set_EL(object,EL);
        object = set_I(object,I);
    end
    
end
