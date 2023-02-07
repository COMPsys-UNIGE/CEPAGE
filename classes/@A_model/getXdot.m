function x_dot = getXdot(object,t,x,varargin)
% getXdot    Computes the derivative of the state 
%
%  x_dot = getXdot(object,t,x)
%   compute the time derivative of the model at time instant t, in state x
%   and with I_{ext} = 0;
%
%  x_dot = getXdot(object,t,x,I_{ext})
%   compute the time derivative of the model at time instant t, in state x
%   and with I_{ext} specified by the user;


% synapsis input
Isyn= 0;
if nargin >= 5
    if isscalar(varargin{end})
        Isyn = varargin{end};
    end
end


x_dot = zeros(6,1);

% Define param
Ca_shift = object.Ca_shift;
x_shift = object.x_shift;
gh = object.gh;
Vhh = object.Vhh;
Iext = object.Iext;

% equation
V = x(1);
h = x(2);
n = x(3);
chi = x(4);
ksi = x(5);
Ca = x(6);


% compute help variable


V_dot = 8*((0.1*(50-(127*V/105+8265/105))/(exp((50 - (127*V/105 ...
    +8265/105))/10) - 1))/((0.1*(50 - (127*V/105 + 8265/105))/(exp((50 - (127*V/105 + 8265/105))/10) - 1))+...
    (4*exp((25 - (127*V/105 + 8265/105))/18))))^3*h*(30 - V) + 1.3*n^4*(-75-V)+0.01*chi*(30-V) ...
    +0.03*Ca/(.5 + Ca)*(-75 - V)+0.003*(-40 - V)+gh*((1/(1+exp(-(V+63)/7.8)))^3)*ksi*(120-V)+Isyn+Iext;

h_dot = ((1-h)*(0.07*exp((25 - (127*V/105 + 8265/105))/20))-h*(1.0/(1 + exp((55 - (127*V/105 + 8265/105))/10))))/12.5;

n_dot = ((1-n)*(0.01*(55 - (127*V/105 + 8265/105))/(exp((55 - (127*V/105 + 8265/105))/10) - 1))-n*(0.125*exp((45 - (127*V/105 + 8265/105))/80)))/12.5;

chi_dot = ((1/(exp(0.15*(-V-50+x_shift))+1))-chi)/235;

ksi_dot = 0.5*((1/(1+exp(10*(V-Vhh))))-ksi)/(7.1+10.4/(1+exp((V+68)/2.2)));

Ca_dot = 0.0001*(0.0085*chi*(140-V+Ca_shift)-Ca);


x_dot = [V_dot; h_dot; n_dot; chi_dot; ksi_dot; Ca_dot];
end