classdef FN_relaxation_model < neuron_model
    
    % FN_relaxation_model   Neuron represented by the FitzHugh-Nagumo relaxation model
    %
    % The Neuron is modeld as follow:
    %   _
    %  |
    %  | dV/dt = V − V^3 − x + I - I_{ext}
    % <
    %  | dx/dt = eps ((1/1+e^-10V)-X)
    %  |_-
    %
    % Ther model has two state variable; V represent the membrane potential.
    % eps and I are parameters that allows to change neuron behaviours.
    %
    % OBJ = FN_relaxation_model()
    % Builds an FN_relaxation_model object OBJ with all parameters equal to 0.
    %
    % OBJ = FN_relaxation_model(I,eps)
    % Builds an FN_relaxation_model object OBJ with user assigned parameters value
    %
    %
    %   FN_relaxation_model methods:
    %   getXdot - omputes the derivative of the state
    %   sim - simulates the neuron.
    %   simplot - simulates the neuron system and plots time evolution of
    %             states.
    %   disp - displays some information about the FN_relaxation_model object
    %   get_eps - gets the value of parameter eps
    %   get_I - gets the value of parameter I
    %   set_eps - sets the value of parameter eps
    %   set_I - sets the value of parameter I
    %   getJacobian -Computes the Jacobian of the vector field in point x
    %   generateC - Generates C files for the computation of model vector field
    %
    % The FN_relaxation_model object is derived from neuron_model and
    % inherits all its methods.
    %
    % See also HR_model, neuron_model
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
        I = 0;
        eps = 0;
    end
    
    methods
        function obj = FN_relaxation_model(varargin)
            obj.nx = 2;
            obj.xnames = {'V','x'};
            obj.modelName = 'FitzHugh-Nagumo relaxation';
            obj.isContinuous = true;
            if nargin == 2
                obj.I = varargin{1};
                obj.eps = varargin{2};
            end
        end
        
        % Other methods
        x_dot = getXdot(object,t,x,varargin);
        disp(object);
        str = getCbuilder(object);
        dfdi = getdfdi(object,x); % compute the differential of f with respect to Isyn
        
        % get methods
        I = get_I(object);
        eps = get_eps(object);
        
        % set methods
        object = set_I(object,I);
        object = set_eps(object,eps);
    end
    
end