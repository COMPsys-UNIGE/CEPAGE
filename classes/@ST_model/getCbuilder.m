function str = getCbuilder(object)

%   getCbuilder Generates a string of C code for the computation 
%                 of model vector field 
%
    str = sprintf('ST_model(%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e)',...
    object.g_Na,object.g_K,object.g_L,object.g_A,object.g_T,...
    object.g_KCa,object.g_HVA,object.epsilon,object.E_Na,object.E_K,object.E_L,object.E_Ca,...
    object.C,object.A,object.y0,object.V_c,object.w,object.k_Ca,object.alpha,object.k,...
    object.nu_m,object.nu_h,object.nu_n,object.nu_nA,object.nu_hA,object.nu_mT,object.nu_hT,...
    object.nu_mHVA,object.s_m,object.s_h,object.s_n,object.s_nA,object.s_hA,object.s_mT,...
    object.s_hT,object.s_mHVA,object.tau_nA,object.tau_hA,object.tau_hT,object.tau_mHVA,object.I_app);
end