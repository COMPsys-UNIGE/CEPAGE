% A_model  
classdef A_model < neuron_model
    
    %properties
    properties(Access = protected)
        Ca_shift = 0;
        x_shift = 0;
        gh = 0;
        Vhh = 0;
        Iext = 0;
    end
    
    methods
        function obj = A_model(varargin)
            obj.nx = 6;
            obj.xnames = {'V','h','n','chi','ksi','Ca'};
            obj.modelName = 'Andrey';
            obj.isContinuous = true;
            if nargin == 5
                obj.Ca_shift = varargin{1};
                obj.x_shift = varargin{2};
                obj.gh = varargin{3};
                obj.Vhh = varargin{4};
                obj.Iext = varargin{5};
            end
        end
        
        x_dot = getXdot(object,t,x,varargin);
        disp(object);
        str = getCbuilder(object);
        dfdi = getdfdi(object,x); % compute the differential of f with respect to Isyn
        
        % get methods
        Ca_shift = get_Ca_shift(object);
        x_shift = get_x_shift(object);
        gh = get_gh(object);
        Vhh = get_Vhh(object);
        Iext = get_Iext(object);
        
        %set methods
        object = set_Ca_shift(object,Ca_shift);
        object = set_x_shift(object,x_shift);
        object = set_gh(object,gh);
        object = set_Vhh(object,Vhh);
        object = set_Iext(object,Iext)
    end
    
end