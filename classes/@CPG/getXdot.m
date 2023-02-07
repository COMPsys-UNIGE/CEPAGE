
function x_dot = getXdot(object,t,x,varargin)
% getXdot    Computes the derivative of the state
%
%  x_dot = getXdot(object,t,x)
%   compute the time derivative of the model at time instant t, in state x
%   and with I_{ext} = 0;
%
%  x_dot = getXdot(object,t,x,x_delayed)
%   compute the time derivative of the model at time instant t, in state x
%   and with I_{ext} = 0; x_delayed are the state values at t-delays
%
%  x_dot = getXdot(object,t,x,x_delayed,I_{ext})
%   compute the time derivative of the model at time instant t and in state
%   x. x_delayed are the state values at t-delays. I_{ext} is a vector
%   containing the external currents that act on each neuron.
%
% Contributors:
%
% Matteo Lodi (matteo.lodi@edu.unige.it)
%
% Copyright (C) 2016 University of Genoa, Italy.

% Legal note:
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public
% License as published by the Free Software Foundation; either
% version 2.1 of the License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
%
% You should have received a copy of the GNU General Public
% License along with this library; if not, write to the
% Free Software Foundation, Inc.,
% 59 Temple Place, Suite 330,
% Boston, MA  02111-1307  USA



delays = object.delays;

if nargin == 4
    Iext = zeros(object.N,1);
    Z = varargin{1};
elseif nargin == 5
    Iext = varargin(2);
    Z = varargin{1};
else
    Iext = zeros(object.N,1);
    Z = zeros(object.totState,delays);
end

N = object.N;
g_in = object.g_in;
g_ex = object.g_ex;
g_el = object.g_el;

inhActivation = object.inhActivation;
excActivation = object.excActivation;

neurons = object.neurons;

%The neurons must be all the same name of state variable
% nx = neurons{1}.getnx();

totStates = object.totState;
incrementalIndexState = object.incrementalIndexState;

x_dot = zeros(totStates,1);

%Valori sinapsi
EsynIn = object.EsynIn;
EsynEx = object.EsynEx;



if numel(delays) == 0
    
    
    % Update neurons states
    for i=1:N
        Isyn = 0;
        beginI = incrementalIndexState(i);
        endI = incrementalIndexState(i+1)-1;
        %Compute exc syn
        for j=1:N
            if g_ex(i,j) ~= 0
                
                synBeginI = incrementalIndexState(N+N*N+(i-1)*N+j);
                synEndI = incrementalIndexState(N+N*N+(i-1)*N+j+1)-1;
                
                otherI = incrementalIndexState(j);
                Isyn = Isyn + g_ex(i,j)*(-x(beginI)+EsynEx)*excActivation{i,j}.getActivation(x(otherI),x(synBeginI:synEndI));
                x_dot(synBeginI:synEndI) = excActivation{i,j}.getXdot(t,x(synBeginI:synEndI),x(otherI));
            end
        end
        %compute inh syn
        for j=1:N
            if g_in(i,j) ~= 0
                synBeginI = incrementalIndexState(N+(i-1)*N+j);
                synEndI = incrementalIndexState(N+(i-1)*N+j+1)-1;
                otherI = incrementalIndexState(j);
                Isyn = Isyn + g_in(i,j)*(-x(beginI)+EsynIn)*inhActivation{i,j}.getActivation(x(otherI),x(synBeginI:synEndI));
                x_dot(synBeginI:synEndI) = inhActivation{i,j}.getXdot(t,x(synBeginI:synEndI),x(otherI));
                
            end
        end
        %compute el syn
        for j=1:N
            if g_el(i,j) ~= 0
                otherI = incrementalIndexState(j);
                Isyn = Isyn + g_el(i,j)*(+x(otherI)-x(beginI));
            end
        end
        
        Isyn = Isyn+Iext(i);
        
        x_dot(beginI:endI) = neurons{i}.getXdot(t,x(beginI:endI),Isyn);
    end
    
    
else
    
    
    delayIndexNeur = object.delayIndexNeur;
    delayInhSyn = object.delayInhSyn;
    delayExcSyn = object.delayExcSyn;
    
    
    % Update neurons states
    for i=1:N
        Isyn = 0;
        beginI = incrementalIndexState(i);
        endI = incrementalIndexState(i+1)-1;
        %Compute exc syn
        for j=1:N
            if g_ex(i,j) ~= 0
                
                synBeginI = incrementalIndexState(N+N*N+(i-1)*N+j);
                synEndI = incrementalIndexState(N+N*N+(i-1)*N+j+1)-1;
                
                otherI = incrementalIndexState(j);
                Ztmp = Z(otherI,delayExcSyn{i,j});
                Isyn = Isyn + g_ex(i,j)*(-x(beginI)+EsynEx)*excActivation{i,j}.getActivation(x(otherI),x(synBeginI:synEndI),Ztmp);
                
                Ztmp2 = Z(:,delayExcSyn{i,j});
                x_dot(synBeginI:synEndI) = excActivation{i,j}.getXdot(t,x(synBeginI:synEndI),x(otherI),Ztmp2,Ztmp);
                
            end
        end
        %compute inh syn
        for j=1:N
            if g_in(i,j) ~= 0
                
                synBeginI = incrementalIndexState(N+(i-1)*N+j);
                synEndI = incrementalIndexState(N+(i-1)*N+j+1)-1;
                
                otherI = incrementalIndexState(j);
                Ztmp = Z(otherI,delayInhSyn{i,j});
                Isyn = Isyn + g_in(i,j)*(-x(beginI)+EsynIn)*inhActivation{i,j}.getActivation(x(otherI),x(synBeginI:synEndI),Ztmp);
                Ztmp2 = Z(:,delayInhSyn{i,j});
                x_dot(synBeginI:synEndI) = inhActivation{i,j}.getXdot(t,x(synBeginI:synEndI),x(otherI),Ztmp2,Ztmp);
            end
        end
        %compute el syn
        for j=1:N
            if g_el(i,j) ~= 0
                otherI = incrementalIndexState(j);
                Isyn = Isyn + g_el(i,j)*(+x(otherI)-x(beginI));
            end
        end
        
        Isyn = Isyn+Iext(i);
        
        Ztmp = Z(:,delayIndexNeur{i});
        
        x_dot(beginI:endI) = neurons{i}.getXdot(t,x(beginI:endI),Ztmp,Isyn);
    end
    
end


end
