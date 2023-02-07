function J = getDA(object,Vpre,varargin)



J = object.nu*exp(-object.nu*(Vpre-object.theta))./((1+exp(-object.nu*(Vpre-object.theta))).^2);