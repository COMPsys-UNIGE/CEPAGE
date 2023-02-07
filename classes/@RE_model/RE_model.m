    
classdef RE_model < neuron_model
    
    %properties
    properties(Access = protected)
        gna = 0; % parameter gna
        ENa = 0; % parameter ENa
        gk = 0; % parameter gk
        Ek = 0; % parameter Ek
        gca = 0; % parameter gca
        gl = 0; % parameter gl
        El = 0; % parameter El
        C = 0; % parameter C
        Iapp = 0; % parameter Iapp
        Ca0 = 0; % parameter Ca0
        d = 0;
        KT = 0;
        Kd = 0;
        gd = 0;
        D = 0;
        EsynEx = 0;
        tau = 0;
        tTm1 = 0;
        tTm2 = 0;
        tTh1 = 0;
        tTh2 = 0;
    end
    
    methods
        function obj = RE_model(varargin)
            obj.nx = 7;
            obj.xnames = {'V','h','m','n','mT','hT','Ca'};
            obj.modelName = 'Reticular';
            obj.isContinuous = true;
            if nargin == 21
                obj.gna = varargin{1};
                obj.ENa = varargin{2};
                obj.gk = varargin{3};
                obj.Ek = varargin{4};
                obj.gca = varargin{5};
                obj.gl = varargin{6};
                obj.El = varargin{7};
                obj.C = varargin{8};
                obj.Iapp = varargin{9};
                obj.Ca0 = varargin{10};
                obj.d = varargin{11};
                obj.KT = varargin{12};
                obj.Kd = varargin{13};
                obj.gd = varargin{14};
                obj.D = varargin{15};
                obj.EsynEx = varargin{16};
                obj.tau = varargin{17};
                obj.tTm1 = varargin{18};
                obj.tTm2 = varargin{19};
                obj.tTh1 = varargin{20};
                obj.tTh2 = varargin{21};
            end
        end
        
        x_dot = getXdot(object,t,x,varargin);
        disp(object);
        str = getCbuilder(object);
        dfdi = getdfdi(object,x); % compute the differential of f with respect to Isyn
        
        % get methods
        gna = get_gna(object);
        ENa = get_ENa(object);
        gk = get_gk(object);
        Ek = get_Ek(object);
        gca = get_gca(object);
        gl = get_gl(object);
        El = get_El(object);
        C = get_C(object);
        Iapp = get_Iapp(object);
        Ca0 = get_Ca0(object);
        d = get_d(object);
        KT = get_KT(object);
        Kd = get_Kd(object);
        gd = get_gd(object);
        D = get_D(object); 
        EsynEx = get_EsynEx(object);
        tau = get_tau(object);
        tTm1 = get_tTm1(object);
        tTm2 = get_tTm2(object);
        tTh1 = get_tTh1(object);
        tTh2 = get_tTh2(object);
        

        
        %set methods
        object = set_gna(object,gna);
        object = set_ENa(object,ENa);
        object = set_gk(object,gk);
        object = set_Ek(object,Ek);
        object = set_gca(object,gca);
        object = set_gl(object,gl);
        object = set_El(object,El);
        object = set_C(object,C);
        object = set_Iapp(object,Iapp);
        object = set_Ca0(object,Ca0);
        object = set_d(object,d);
        object = set_KT(object,KT);
        object = set_Kd(object,Kd);
        object = set_gd(object,gd);
        object = set_D(object,D); 
        object = set_EsynEx(object,EsynEx);
        object = set_tau(object,tau);
        object = set_tTm1(object,tTm1);
        object = set_tTm2(object,tTm2);
        object = set_tTh1(object,tTh1);
        object = set_tTh2(object,tTh2);
    end
    
end