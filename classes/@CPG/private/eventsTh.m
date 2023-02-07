%% EVENTSTH
% Funzione per intercettare l'evento V_membrana = Vth in una rete di neuroni con modello LH


function [value,isterminal,direction] = eventsTh(object,t,y,varargin)

Vth = 0;
if nargin == 4
    if isscalar(varargin{1})
        Vth = varargin{1};
    else
        error('Threshold must be a scalar');
    end
end

N = object.N;

startIndex = object.incrementalIndexState;

value = zeros(N,1);
for i = 1:N
    index = startIndex(i);
    value(i,1) = y(index) - Vth;
end
isterminal = zeros(N,1);
direction = ones(N,1);

end