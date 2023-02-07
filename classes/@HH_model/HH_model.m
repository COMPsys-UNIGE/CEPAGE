% HH_model   Neuron represented by the Hodgkin-Huxley model
    %
    % The Neuron is modeld as follow:
    %   _
    %  |
    %  | dV/dt = (-I_{Na}-I_{k2}-I_l-I_{app}-I_{ext})/C; 
    %  | 
    % <  dh/dt = (h_{Inf}- h)/t_{Na}
    %  |
    %  | dm/dt = (m_{Inf}- m)/t_{k2}
    %  |_
    %
    %   h_{Inf} = 1/(1 + exp(500*(X+0.333)))
    %   n_{Inf} = 1/(1 + exp(-150*(X+0.305)))
    %   m_{Inf} = 1/(1 + exp(-83*(X+(0.018+VshiftK2))))
    %
    %   n = n_{Inf};
    %   I_{Na} = gna*n^3*h*(V-E_{Na})
    %   I_{k2} = gk2*m^2*(V-E_k)
    %   I_l = g_l*(V-E_l)
    %
    % Ther model has three state variable; V represent the membrane potential.
    % gna,ENa,gk2,Ek,gl,El,tNa,tk2,C,Iapp,VshiftK2 are parameters that 
    % allows to change neuron behaviours.
    %
    % OBJ = HH_model()
    % Builds an HH_model object OBJ with all parameters equal to 0.
    %
    % OBJ = HH_model(gna,ENa,gk2,Ek,gl,El,tNa,tk2,C,Iapp,VshiftK2)
    % Builds an HH_model object OBJ with user assigned parameters value
    %
    %
    %   HH_model methods:
    %   getXdot - omputes the derivative of the state
    %   sim - simulates the neuron.
    %   simplot - simulates the neuron system and plots time evolution of
    %             states.
    %   disp - displays some information about the HH_model object
    %   get_gna - gets the value of parameter gna
    %   get_ENa - gets the value of parameter ENa
    %   get_gk2 - gets the value of parameter gk2
    %   get_Ek - gets the value of parameter Ek
    %   get_gl - gets the value of parameter gl
    %   get_El - gets the value of parameter El
    %   get_tNa - gets the value of parameter tNa
    %   get_tk2 - gets the value of parameter tk2
    %   get_C - gets the value of parameter C
    %   get_Iapp - gets the value of parameter Iapp
    %   get_VshiftK2 - gets the value of parameter VshiftK2
    %   set_gna - sets the value of parameter gna
    %   set_ENa - sets the value of parameter ENa
    %   set_gk2 - sets the value of parameter gk2
    %   set_Ek - sets the value of parameter Ek
    %   set_gl - sets the value of parameter gl
    %   set_El - sets the value of parameter El
    %   set_tNa - sets the value of parameter tNa
    %   set_tk2 - sets the value of parameter tk2
    %   set_C - sets the value of parameter C
    %   set_Iapp - sets the value of parameter Iapp
    %   set_VshiftK2 - sets the value of parameter VshiftK2
    %   getJacobian -Computes the Jacobian of the vector field in point x 
    %   generateC - Generates C files for the computation of model vector field
    %
    % The HH_model object is derived from neuron_model and inherits all its methods.
    %
    % See also ML_model, neuron_model
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
    
classdef HH_model < neuron_model
    
    %properties
    properties(Access = protected)
        gna = 0; % parameter gna
        ENa = 0; % parameter ENa
        gk2 = 0; % parameter gk2
        Ek = 0; % parameter Ek
        gl = 0; % parameter gl
        El = 0; % parameter El
        tNa = 0; % parameter tNa
        tk2 = 0; % parameter tk2
        C = 0; % parameter C
        Iapp = 0; % parameter Iapp
        VshiftK2 = 0; % parameter VshiftK2
    end
    
    methods
        function obj = HH_model(varargin)
            obj.nx = 3;
            obj.xnames = {'V','h','m'};
            obj.modelName = 'Hodgkin-Huxley';
            obj.isContinuous = true;
            if nargin == 11
                obj.gna = varargin{1};
                obj.ENa = varargin{2};
                obj.gk2 = varargin{3};
                obj.Ek = varargin{4};
                obj.gl = varargin{5};
                obj.El = varargin{6};
                obj.tNa = varargin{7};
                obj.tk2 = varargin{8};
                obj.C = varargin{9};
                obj.Iapp = varargin{10};
                obj.VshiftK2 = varargin{11};
            end
        end
        
        x_dot = getXdot(object,t,x,varargin);
        disp(object);
        str = getCbuilder(object);
        dfdi = getdfdi(object,x); % compute the differential of f with respect to Isyn
        
        % get methods
        gna = get_gna(object);
        ENa = get_ENa(object);
        gk2 = get_gk2(object);
        Ek = get_Ek(object);
        gl = get_gl(object);
        El = get_El(object);
        tNa = get_tNa(object);
        tk2 = get_tk2(object);
        C = get_C(object);
        Iapp = get_Iapp(object);
        VshiftK2 = get_VshiftK2(object);

        
        %set methods
        object = set_gna(object,gna);
        object = set_EnN(object,ENa);
        object = set_gk2(object,gk2);
        object = set_Ek(object,Ek);
        object = set_gl(object,gl);
        object = set_El(object,El);
        object = set_tNa(object,tNa);
        object = set_tk2(object,tk2);
        object = set_C(object,C);
        object = set_Iapp(object,Iapp);
        object = set_VshiftK2(object,VshiftK2);
    end
    
end