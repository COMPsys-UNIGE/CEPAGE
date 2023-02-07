% Plant_model  
   
    
classdef Plant_model < neuron_model
    
    %properties
    properties(Access = protected)
        g_I = 0;
        g_K = 0;
        g_T = 0;
        g_L = 0;
        g_KCa = 0;
        V_I = 0;
        V_K = 0;
        V_L = 0;
        V_Ca = 0;
        K_c = 0;
        rho = 0;
        Iapp = 0;
    end
    
    methods
        function obj = Plant_model(varargin)
            obj.nx = 5;
            obj.xnames = {'V','Ca','h','n','chi'};
            obj.modelName = 'Plant';
            obj.isContinuous = true;
            if nargin == 12
                obj.g_I= varargin{1};
                obj.g_K= varargin{2};
                obj.g_T= varargin{3};
                obj.g_L= varargin{4};
                obj.g_KCa= varargin{5};
                obj.V_I= varargin{6};
                obj.V_K= varargin{7};
                obj.V_L= varargin{8};
                obj.V_Ca= varargin{9};
                obj.K_c= varargin{10};
                obj.rho= varargin{11};
                obj.Iapp= varargin{12};

            end
        end
        
        x_dot = getXdot(object,t,x,varargin);
        disp(object);
        str = getCbuilder(object);
        dfdi = getdfdi(object,x); % compute the differential of f with respect to Isyn
        
        % get methods
        g_I = get_g_I(object);
        g_K = get_g_K(object);
        g_T = get_g_T(object);
        g_L = get_g_L(object);
        g_KCa = get_g_KCa(object);
        V_I = get_V_I(object);
        V_K = get_V_K(object);
        V_L = get_V_L(object);
        V_Ca = get_V_Ca(object);
        K_c = get_K_c(object);
        rho = get_rho(object);
        Iapp = get_Iapp(object);
     
        
        %set methods
        object = set_g_I(object,g_I);
        object = set_g_K(object,g_K);
        object = set_g_T(object,g_T);
        object = set_g_L(object,g_L);
        object = set_g_KCa(object,g_KCa);
        object = set_V_I(object,V_I);
        object = set_V_K(object,V_K);
        object = set_V_L(object,V_L);
        object = set_V_Ca(object,V_Ca);
        object = set_K_c(object,K_c);
        object = set_rho(object,rho);
        object = set_Iapp(object,Iapp);
    end
    
end