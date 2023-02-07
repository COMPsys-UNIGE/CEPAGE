function str = getCbuilder(object)

%   getCbuilder Generates a string of C code for the computation 
%                 of model vector field 
%
    str = sprintf('Plant_model(%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e)',...
    object.g_I,object.g_K,object.g_T,object.g_L,object.g_KCa,object.V_I,...
    object.V_K,object.V_L,object.V_Ca,object.K_c,object.rho,object.Iapp);
end