classdef alpha_synapse_model < synapse_model
    
    % Alphan synapse model
    %
    % This object represents an alpha synapse
    % The synapse is characterized by an activation function
    %
    %   f(t,x) = x
    %
    % The states equations are:
    %  dx/dt = alpha * (1-x) * 1/(1+exp(-nu*(V-theta))-beta*x
    %
    %
    % OBJ = alpha_synapse_model()
    % Builds an alpha_synapse_model object OBJ with all parameters equal to 0.
    %
    % OBJ = alpha_synapse_model(alpha,beta,nu,theta)
    % Builds an alpha_synapse_model object OBJ with user assigned parameters value
    %
    % alpha_synapse_model methods:
    %   getActivation - Compute the activation function
    %
    % See also synapse_model
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
        nu = 0;
        theta = 0;
        al = 0;
        beta = 0;
    end
    
    
    
    methods
        %constructor
        function object = alpha_synapse_model(varargin)
            if nargin == 0
                object.nu = 0;
                object.theta = 0;
                object.al = 0;
                object.beta = 0;
                object.nx = 1;
                object.modelName = 'alpha_synapse';
                object.isContinuous = true;
            elseif nargin == 4
                object.al = varargin{1};
                object.beta = varargin{2};
                object.nu = varargin{3};
                object.theta = varargin{4};
                object.nx = 1;
                object.modelName = 'alpha_synapse';
                object.isContinuous = true;
            else
                error('Wrong number of input arguments');
            end
            
            
        end
        
        act = getActivation(object,Vpre,varargin);
        str = getCbuilder(object);
        disp(object);
        x_dot = getXdot(object,t,x,varargin);
        
    end
    
end
