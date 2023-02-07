classdef FTM_delayed_synapse_model < synapse_model
    
    % Fast threshold modulation delayedsynapse model
    %
    % This object represents a synapse described by the Fast threshold 
    % modulation model.
    % The synapse is characterized by an activation function 
    %
    %   f(x) = g1*1/(1+exp(-nu*(V(t)-theta))+g2*1/(1+exp(-nu*(V(t-t0)-theta))
    %
    %
    % OBJ = FTM_synapse_model()
    % Builds an FTM_synapse_model object OBJ with all parameters equal to 0.
    %
    % OBJ = FTM_synapse_model(nu,theta,t0,g1,g2)
    % Builds an FTM_synapse_model object OBJ with user assigned parameters value
    %
    % FTM_synapse_model methods:
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
        g2 = 0;
        g1 = 0;
    end
       

    
    methods
        %constructor
        function object = FTM_delayed_synapse_model(varargin)
            if nargin == 0
                object.nu = 0;
                object.theta = 0;
                object.nx = 0;
                object.modelName = 'FTM_delayed_synapse';
                object.isContinuous = true;
            elseif nargin == 4
                object.nu = varargin{1};
                object.theta = varargin{2};
                object.delays = varargin{3};
%                 object.g1 = varargin{4};
%                 object.g2 = varargin{5};
                object.g2 = varargin{4};
                object.g1 = 1;
                object.nx = 0;
                object.modelName = 'FTM_delayed_synapse';
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
