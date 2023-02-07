function str = getCbuilder(object)
    str = sprintf('EIF_model(%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e)',...
    object.C,object.a,object.b,object.VT,object.DT,object.I,object.gL,object.EL,object.Vr,object.tw,object.gexp,object.gd,object.D,object.EsynEx);
end