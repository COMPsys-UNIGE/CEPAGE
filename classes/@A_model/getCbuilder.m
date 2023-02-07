function str = getCbuilder(object)

%   getCbuilder Generates a string of C code for the computation 
%                 of model vector field 
%
    str = sprintf('A_model(%e,%e,%e,%e,%e)',...
    object.Ca_shift,object.x_shift,object.gh,object.Vhh,object.Iext);
end
