function [position,isterminal,direction] = getResetConditions(object,t,y,varargin)

% getResetConditions    get the reset condition for neuron model
%
% Compute the reset condition x = xbar for neuron model
%
% [position,isterminal,direction] = getResetConditions(object,t,x)
% where t is current time and x is current state 
%
% position, isterminal and direction have the same role as the ode event
% function. See https://it.mathworks.com/help/matlab/math/ode-event-location.html
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

N = object.N;

totState = object.totState;
incrementalIndexState = object.incrementalIndexState;

position = zeros(totState,1); % The value that we want to be zero
isterminal = zeros(totState,1);  % Halt integration
direction = zeros(totState,1);   % The zero can be approached from either direction


inhActivation = object.inhActivation;
excActivation = object.excActivation;

neurons = object.neurons;


if ~object.is_delayed
    
    for i=1:N
        beginI = incrementalIndexState(i);
        endI = incrementalIndexState(i+1)-1;
        

        
        [position(beginI:endI),isterminal(beginI:endI),direction(beginI:endI)] = neurons{i}.getResetConditions(t,y(beginI:endI));
    end
    
    
    % synapses reset
    ii = 1;
    for i=1:N
        for j=1:N
            beginI = incrementalIndexState(N+ii);
            endI = incrementalIndexState(N+ii+1)-1;
            
                    otherI = incrementalIndexState(j);
            
            [position(beginI:endI),isterminal(beginI:endI),direction(beginI:endI)] = inhActivation{i,j}.getResetConditions(t,y(beginI:endI),y(otherI));
            ii = ii+1;
        end
    end
    
    for i=1:N
        for j=1:N
            beginI = incrementalIndexState(N+ii);
            endI = incrementalIndexState(N+ii+1)-1;
            
            otherI = incrementalIndexState(j);
            
            [position(beginI:endI),isterminal(beginI:endI),direction(beginI:endI)] = excActivation{i,j}.getResetConditions(t,y(beginI:endI),y(otherI));
            ii = ii+1;
        end
    end
    
    
else
    
    
    if nargin == 4
        Z = varargin{1};
    else
        Z = zeros(object.totState,numel(object.delays));
    end
    
    
    delayIndexNeur = object.delayIndexNeur;
    delayInhSyn = object.delayInhSyn;
    delayExcSyn = object.delayExcSyn;
    
    
    for i=1:N
        beginI = incrementalIndexState(i);
        endI = incrementalIndexState(i+1)-1;
        
        Ztmp = Z(:,delayIndexNeur{i});
        
        
        [position(beginI:endI),isterminal(beginI:endI),direction(beginI:endI)] = neurons{i}.getResetConditions(t,y(beginI:endI),Ztmp);
    end
    
    
    % synapses reset
    ii = 1;
    for i=1:N
        for j=1:N
            beginI = incrementalIndexState(N+ii);
            endI = incrementalIndexState(N+ii+1)-1;
            
            otherI = incrementalIndexState(j);
            
            Ztmp = Z(otherI,delayInhSyn{i,j});
            
            
            
%             [position(beginI:endI),isterminal(beginI:endI),direction(beginI:endI)] = inhActivation{i,j}.getResetConditions(t,y(beginI:endI),x(otherI),Ztmp);
            ii = ii+1;
        end
    end
    
    for i=1:N
        for j=1:N
            beginI = incrementalIndexState(N+ii);
            endI = incrementalIndexState(N+ii+1)-1;
            
            otherI = incrementalIndexState(j);
            
            Ztmp = Z(:,delayExcSyn{i,j});
            
%             [position(beginI:endI),isterminal(beginI:endI),direction(beginI:endI)] = excActivation{i,j}.getResetConditions(t,y(beginI:endI),x(otherI),Ztmp);
            ii = ii+1;
        end
    end
    
end

