function [xreset,object] = resetStates(object,t,x,varargin)

xreset(1) = object.Vr;
xreset(2) = x(2)+object.b;

end

