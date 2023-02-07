function [xreset,object] = resetStates(object,t,x,ie,varargin)

% resetStates    reset the states 
%
% Compute the resetted state x = x_r
%
% [xreset,object] = resetStates(object,t,x)
% where t is current time and x is current state 
%
% xreset is the resetted state
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

xreset = zeros(size(x));

N = object.N;

incrementalIndexState = object.incrementalIndexState;

inhActivation = object.inhActivation;
excActivation = object.excActivation;

neurons = object.neurons;



if numel(object.delays) == 0
    
    for i=1:N
        beginI = incrementalIndexState(i);
        endI = incrementalIndexState(i+1)-1;
        
        if any(ie >= beginI & ie <= endI)
            xreset(beginI:endI) = neurons{i}.resetStates(t,x(beginI:endI));
        else
            xreset(beginI:endI) = x(beginI:endI);
        end
    end
    
    
    % synapses reset
    ii = 1;
    for i=1:N
        for j=1:N
            beginI = incrementalIndexState(N+ii);
            endI = incrementalIndexState(N+ii+1)-1;
            
            otherI = incrementalIndexState(j);
            
            if any(ie >= beginI & ie <= endI)
                xreset(beginI:endI) = inhActivation{i,j}.resetStates(t,x(beginI:endI),x(otherI));
            else
                xreset(beginI:endI) = x(beginI:endI);
            end
            ii = ii+1;
        end
    end
    
    for i=1:N
        for j=1:N
            beginI = incrementalIndexState(N+ii);
            endI = incrementalIndexState(N+ii+1)-1;
            
            otherI = incrementalIndexState(j);
            
            if any(ie >= beginI & ie <= endI)
                xreset(beginI:endI) = excActivation{i,j}.resetStates(t,x(beginI:endI),x(otherI));
            else
                xreset(beginI:endI) = x(beginI:endI);
            end
            ii = ii+1;
        end
    end
    
    
else
    
    
    if nargin == 5
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
        
        
        if any(ie >= beginI & ie <= endI)
            Ztmp = Z(:,delayIndexNeur{i});
            xreset(beginI:endI) = neurons{i}.resetStates(t,x(beginI:endI),Ztmp);
        else
            xreset(beginI:endI) = x(beginI:endI);
        end
    end
    
    
    % synapses reset
    ii = 1;
    for i=1:N
        for j=1:N
            beginI = incrementalIndexState(N+ii);
            endI = incrementalIndexState(N+ii+1)-1;
            
            otherI = incrementalIndexState(j);
            
            if any(ie >= beginI & ie <= endI)
                Ztmp = Z(otherI,delayInhSyn{i,j});
                xreset(beginI:endI) = inhActivation{i,j}.resetStates(t,x(beginI:endI),x(otherI),Ztmp);
            else
                xreset(beginI:endI) = x(beginI:endI);
            end
            ii = ii+1;
        end
    end
    
    for i=1:N
        for j=1:N
            beginI = incrementalIndexState(N+ii);
            endI = incrementalIndexState(N+ii+1)-1;
            
            otherI = incrementalIndexState(j);
            
            if any(ie >= beginI & ie <= endI)
                Ztmp = Z(otherI,delayExcSyn{i,j});
                xreset(beginI:endI) = excActivation{i,j}.resetStates(t,x(beginI:endI),x(otherI),Ztmp);
            else
                xreset(beginI:endI) = x(beginI:endI);
            end
            ii = ii+1;
        end
    end
    
    
    
    
    
end

