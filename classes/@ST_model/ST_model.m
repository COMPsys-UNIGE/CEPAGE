
    
classdef ST_model < neuron_model
    
    %properties
    properties(Access = protected)
g_Na = 0;
g_K = 0;
g_L = 0;
g_A = 0;
g_T = 0;
g_KCa = 0;
g_HVA = 0;
epsilon = 0;
E_Na = 0;
E_K = 0;
E_L = 0;
E_Ca = 0;
C = 0;
A = 0;
y0 = 0;
V_c = 0;
w = 0;
k_Ca = 0;
alpha = 0;
k = 0;
nu_m = 0;
nu_h = 0;
nu_n = 0;
nu_nA = 0;
nu_hA = 0;
nu_mT = 0;
nu_hT = 0;
nu_mHVA = 0;
s_m = 0;
s_h = 0;
s_n = 0;
s_nA = 0;
s_hA = 0;
s_mT = 0;
s_hT = 0;
s_mHVA = 0;
tau_nA = 0;
tau_hA = 0;
tau_hT = 0;
tau_mHVA = 0;
I_app = 0;
    end
    
    methods
        function obj = ST_model(varargin)
            obj.nx = 8;
            obj.xnames = {'V_dot', 'h_dot', 'n_dot', 'n_A_dot', 'h_A_dot', 'h_T_dot', 'm_HVA_dot', 'Ca_dot'};
            obj.modelName = 'ST';
            obj.isContinuous = true;
            if nargin == 41
                obj.g_Na = varargin{1};
                obj.g_K = varargin{2};
                obj.g_L = varargin{3};
                obj.g_A = varargin{4};
                obj.g_T = varargin{5};
                obj.g_KCa = varargin{6};
                obj.g_HVA = varargin{7};
                obj.epsilon = varargin{8};
                obj.E_Na = varargin{9};
                obj.E_K = varargin{10};
                obj.E_L = varargin{11};
                obj.E_Ca = varargin{12};
                obj.C = varargin{13};
                obj.A = varargin{14};
                obj.y0 = varargin{15};
                obj.V_c = varargin{16};
                obj.w = varargin{17};
                obj.k_Ca = varargin{18};
                obj.alpha = varargin{19};
                obj.k = varargin{20};
                obj.nu_m = varargin{21};
                obj.nu_h = varargin{22};
                obj.nu_n = varargin{23};
                obj.nu_nA = varargin{24};
                obj.nu_hA = varargin{25};
                obj.nu_mT = varargin{26};
                obj.nu_hT = varargin{27};
                obj.nu_mHVA = varargin{28};
                obj.s_m = varargin{29};
                obj.s_h = varargin{30};
                obj.s_n = varargin{31};
                obj.s_nA = varargin{32};
                obj.s_hA = varargin{33};
                obj.s_mT = varargin{34};
                obj.s_hT = varargin{35};
                obj.s_mHVA = varargin{36};
                obj.tau_nA = varargin{37};
                obj.tau_hA = varargin{38};
                obj.tau_hT = varargin{39};
                obj.tau_mHVA = varargin{40};
                obj.I_app = varargin{41};
                
            end
        end
        
        x_dot = getXdot(object,t,x,varargin);
        disp(object);
        str = getCbuilder(object);
        dfdi = getdfdi(object,x); % compute the differential of f with respect to Isyn
        
        % get methods
        g_Na = get_g_Na(object);
        g_K = get_g_K(object);
        g_L = get_g_L(object);
        g_A = get_g_A(object);
        g_T = get_g_T(object);
        g_KCa = get_g_KCa(object);
        g_HVA = get_g_HVA(object);
        epsilon = get_epsilon(object);
        E_Na = get_E_Na(object);
        E_K = get_E_K(object);
        E_L = get_E_L(object);
        E_Ca = get_E_Ca(object);
        C = get_C(object);
        A = get_A(object);
        y0 = get_y0(object);
        V_c = get_V_c(object);
        w = get_w(object);
        k_Ca = get_k_Ca(object);
        alpha = get_alpha(object);
        k = get_k(object);
        nu_m = get_nu_m(object);
        nu_h = get_nu_h(object);
        nu_n = get_nu_n(object);
        nu_nA = get_nu_nA(object);
        nu_hA = get_nu_hA(object);
        nu_mT = get_nu_mT(object);
        nu_hT = get_nu_hT(object);
        nu_mHVA = get_nu_mHVA(object);
        s_m = get_s_m(object);
        s_h = get_s_h(object);
        s_n = get_s_n(object);
        s_nA = get_s_nA(object);
        s_hA = get_s_hA(object);
        s_mT = get_s_mT(object);
        s_hT = get_s_hT(object);
        s_mHVA = get_s_mHVA(object);
        tau_nA = get_tau_nA(object);
        tau_hA = get_tau_hA(object);
        tau_hT = get_tau_hT(object);
        tau_mHVA = get_tau_mHVA(object);       
        I_app = get_I_app(object);       

        
        %set methods
        object = set_g_Na(object,g_Na);
        object = set_g_K(object,g_K);
        object = set_g_L(object,g_L);
        object = set_g_A(object,g_A);
        object = set_g_T(object,g_T);
        object = set_g_KCa(object,g_KCa);
        object = set_g_HVA(object,g_HVA);
        object = set_epsilon(object,epsilon);
        object = set_E_Na(object,E_Na);
        object = set_E_K(object,E_K);
        object = set_E_L(object,E_L);
        object = set_E_Ca(object,E_Ca);
        object = set_C(object,C);
        object = set_A(object,A);
        object = set_y0(object,y0);
        object = set_V_c(object,V_c);
        object = set_w(object,w);
        object = set_k_Ca(object,k_Ca);
        object = set_alpha(object,alpha);
        object = set_k(object,k);
        object = set_nu_m(object,nu_m);
        object = set_nu_h(object,nu_h);
        object = set_nu_n(object,nu_n);
        object = set_nu_nA(object,nu_nA);
        object = set_nu_hA(object,nu_hA);
        object = set_nu_mT(object,nu_mT);
        object = set_nu_hT(object,nu_hT);
        object = set_nu_mHVA(object,nu_mHVA);
        object = set_s_m(object,s_m);
        object = set_s_h(object,s_h);
        object = set_s_n(objects,s_n);
        object = set_s_nA(object,s_nA);
        object = set_s_hA(object,s_hA);
        object = set_s_mT(object,s_mT);
        object = set_s_hT(object,s_hT);
        object = set_s_mHVA(object,s_mHVA);
        object = set_tau_nA(object,tau_nA);
        object = set_tau_hA(object,tau_hA);
        object = set_tau_hT(object,tau_hT);
        object = set_tau_mHVA(object,tau_mHVA); 
        object = set_I_app(object,I_app); 
    end
    
end