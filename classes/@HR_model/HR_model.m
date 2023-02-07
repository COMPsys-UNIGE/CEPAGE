classdef HR_model < neuron_model
    
    % HR_model   Neuron represented by the Hindmarsh Rose model
    %
    % The Neuron is modeld as follow:
    %   _
    %  |
    %  | dx/dt = y - x^3 + bx^2 - z + I - I_{ext}
    %  |
    % <  dy/dt = 1-5x^2 - y
    %  |
    %  | dz/dt = \mu (s (x-x_{rest}) - z)
    %  |_
    %
    % Ther model has three state variable; x represent the membrane potential.
    % b,I,\mu,s and x_{rest} are parameters that allows to change neuron
    % behaviours.
    %
    % OBJ = HR_model()
    % Builds an HR_model object OBJ with all parameters equal to 0.
    %
    % OBJ = HR_model(b,mu,s,I,x_{rest})
    % Builds an HR_model object OBJ with user assigned parameters value
    %
    %
    %   HR_model methods:
    %   getXdot - computes the derivative of the state
    %   sim - simulates the neuron.
    %   simplot - simulates the neuron system and plots time evolution of
    %             states.
    %   disp - displays some information about the HR_model object
    %   getnx - gets the number of state variables
    %   get_b - gets the value of parameter b
    %   get_s - gets the value of parameter s
    %   get_I - gets the value of parameter I
    %   get_x_rest - gets the value of parameter x_{rest}
    %   get_mu - gets the value of parameter \mu
    %   set_b - sets the value of parameter b
    %   set_s - sets the value of parameter s
    %   set_I - sets the value of parameter I
    %   set_x_rest - sets the value of parameter x_{rest}
    %   set_mu - sets the value of parameter \mu
    %   getJacobian -Computes the Jacobian of the vector field in point x
    %   getCbuilder - Generates a string of C code for the computation
    %                 of model vector field
    %
    % The HR_model object is derived from neuron_model and inherits all its methods.
    %
    % See also HH_model, neuron_moedl
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
        b = 0;  % parameter b
        mu = 0; % parameter mu
        s = 0;  % parameter s
        I = 0;  % parameter I
        x_rest = 0; % parameter x_rest
    end
    
    methods
        % Constructor
        function obj = HR_model(varargin)
            obj.nx = 3;
            obj.xnames = {'x','y','z'};
            obj.modelName = 'Hindmarsh-Rose';
            obj.isContinuous = true;
            if nargin == 5
                obj.b = varargin{1};
                obj.mu = varargin{2};
                obj.s = varargin{3};
                obj.I = varargin{4};
                obj.x_rest = varargin{5};
            end
        end
        
        % Other methods
        str = getCbuilder(object);
        x_dot = getXdot(object,t,x,varargin);
        disp(object);
        dfdi = getdfdi(object,x); % compute the differential of f with respect to Isyn
        
        
        % Get methods
        b = get_b(object);
        s = get_s(object);
        I = get_I(object);
        x_rest = get_x_rest(object);
        mu = get_mu(object);
        
        % Set methods
        object = set_b(object,b);
        object = set_s(object,s);
        object = set_I(object,I);
        object = set_x_rest(object,x_rest);
        object = set_mu(object,mu);
    end
    
    
end