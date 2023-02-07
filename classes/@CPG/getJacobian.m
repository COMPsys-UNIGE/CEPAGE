function J = getJacobian(object,t,x)

% getJacobian    Computes the Jacobian of the vector field in point x
%
%  J = getJacobian(object,x)
%   compute theJacobian of the model in state x
%
% Contributor:
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

N = object.getN;

g_in = object.get_g_in;
g_ex = object.get_g_ex;
g_el = object.get_g_el;

inhActivation = object.inhActivation;
excActivation = object.excActivation;

Esyn_Ex = object.EsynEx;
Esyn_In = object.EsynIn;


ndim_vector = zeros(N,1);

for i=1:N
    ndim_vector(i) = object.neurons{i}.getnx;
end

incrementalIndexState = object.incrementalIndexState;

ndim = incrementalIndexState(end)-1;

if size(x,1) ~= ndim
    error('input vector must be a nx x npoints vector');
end

np = size(x,2);

J = cell(np,1);


for i=1:np
    tmpJ = zeros(ndim,ndim);
    xx = x(:,i);
    
    currentPointer = 1;
    
    
    
    for k=1:N
        neurNx = ndim_vector(k);
        xsingolo = xx(currentPointer:currentPointer+neurNx-1);
        neuronJac = object.neurons{k}.getJacobian(0,xsingolo);
        
        %uncoupled neuron
        tmpJ(currentPointer:currentPointer+neurNx-1,currentPointer:currentPointer+neurNx-1) = neuronJac;
        
        % Add synaptic contribute
        % Add dVidot/dVi
        otherI = 1;
        for j=1:N
            if g_in(k,j) ~= 0
                
                inhInd = k*N+j;
               indN = incrementalIndexState(inhInd):(incrementalIndexState(inhInd+1)-1);     
                
                DA = inhActivation{k,j}.getDA(xx(indN),xx(otherI));
                
                tmpJ(currentPointer,currentPointer) =  tmpJ(currentPointer,currentPointer) - ...
                    (g_in(k,j)*inhActivation{k,j}.getActivation(xx(otherI),xx(indN)))*object.neurons{k}.getdfdi(xsingolo);
                
                tmpJ(currentPointer,otherI) =  tmpJ(currentPointer,otherI) + ...
                    ((g_in(k,j)*(Esyn_In-xx(currentPointer))*DA(1)))*object.neurons{k}.getdfdi(xsingolo);
                
                if inhActivation{k,j}.getnx > 0
                    inhInd = k*N+j;%(i-1)*N+j+N;
                    indN = incrementalIndexState(inhInd):(incrementalIndexState(inhInd+1)-1);
                    
                    tmpJ(currentPointer,indN) = ((g_in(k,j)*(Esyn_In-xx(currentPointer))*DA(2:end)'))*object.neurons{k}.getdfdi(xsingolo);
                end
                
            end
            if g_ex(k,j) ~= 0
                
                inhExc = k*N+j;
               indN = incrementalIndexState(inhExc):(incrementalIndexState(inhExc+1)-1);
               
                
                DA = excActivation{k,j}.getDA(xx(indN),xx(otherI));
                
                tmpJ(currentPointer,currentPointer) =  tmpJ(currentPointer,currentPointer) - ...
                    (g_ex(k,j)*excActivation{k,j}.getActivation(xx(otherI),xx(indN)))*object.neurons{k}.getdfdi(xsingolo);
                
                tmpJ(currentPointer,otherI) =  tmpJ(currentPointer,otherI) + ...
                    ((g_ex(k,j)*(Esyn_Ex-xx(currentPointer))*DA(1)))*object.neurons{k}.getdfdi(xsingolo);
                
                if excActivation{k,j}.getnx > 0
                    excExc = k*N+j+N*N;
                    indN = incrementalIndexState(excExc):(incrementalIndexState(excExc+1)-1);
                    
                    tmpJ(currentPointer,indN) = ((g_ex(k,j)*(Esyn_Ex-xx(currentPointer))*DA(2:end)'))*object.neurons{k}.getdfdi(xsingolo);
                end
                
            end
            if g_el(k,j) ~= 0
                tmpJ(currentPointer,currentPointer) =  tmpJ(currentPointer,currentPointer) - ...
                    g_el(k,j)*object.neurons{k}.getdfdi(xsingolo);
            end
            %                 tmpJ(currentPointer,currentPointer) =  tmpJ(currentPointer,currentPointer) - ...
            %                     ((g_in(k,j)*inhActivation{k,j}.getActivation(xx(otherI)) - ...
            %                     g_ex(k,j)*excActivation{k,j}.getActivation(xx(otherI)) -...
            %                     g_el(k,j)))*object.neurons{k}.getdfdi(xsingolo);
            otherI = otherI + ndim_vector(j);
        end
        
        %             % Add dVidot/dVj
        %             otherI = 1;
        %             for j=1:N
        %                 if g_in(k,j) ~= 0
        %                      tmpJ(currentPointer,otherI) =  tmpJ(currentPointer,otherI) + ...
        %                     ((g_in(k,j)*(Esyn_In-xx(currentPointer))*inhActivation{k,j}.getDA(xx(otherI))))*object.neurons{k}.getdfdi(xsingolo);
        %                 end
        %                 if g_ex(k,j) ~= 0
        %                 tmpJ(currentPointer,otherI) =  tmpJ(currentPointer,otherI) + ...
        %                     ((g_ex(k,j)*(Esyn_Ex-xx(currentPointer))*excActivation{k,j}.getDA(xx(otherI))))*object.neurons{k}.getdfdi(xsingolo);
        %                 end
        %                 if g_el(k,j) ~= 0
        %                     tmpJ(currentPointer,otherI) =  tmpJ(currentPointer,otherI) + g_el(k,j)*object.neurons{k}.getdfdi(xsingolo);
        %                 end
        %
        % %                 tmpJ(currentPointer,otherI) =  tmpJ(currentPointer,otherI) + ...
        % %                     ((g_in(k,j)*(Esyn_In-xx(currentPointer))*inhActivation{k,j}.getDA(xx(otherI)) + ...
        % %                     g_ex(k,j)*(Esyn_Ex-xx(currentPointer))*excActivation{k,j}.getDA(xx(otherI)) + ...
        % %                     g_el(k,j)))*object.neurons{k}.getdfdi(xsingolo);
        %                 otherI = otherI + ndim_vector(j);
        %             end
        
        
        currentPointer = currentPointer+ndim_vector(k);
        
    end
    
    % Add synapses states Jac
    for k=1:N
       for j=1:N
           if inhActivation{k,j}.getnx > 0 && g_in(k,j) ~= 0
               
               inhInd = k*N+j;
               indN = incrementalIndexState(inhInd):(incrementalIndexState(inhInd+1)-1);
               Jn = inhActivation{k,j}.getJacobian(t,xx(indN),xx(incrementalIndexState(j)));
               tmpJ(indN,incrementalIndexState(j)) = Jn(:,1);
               tmpJ(indN,indN) = Jn(:,2:end);
           end
           
           if excActivation{k,j}.getnx > 0 && g_ex(k,j) ~= 0
               
               inhExc = k*N+j+N*N;
               indN = incrementalIndexState(inhExc):(incrementalIndexState(inhExc+1)-1);
               Jn = excActivation{k,j}.getJacobian(t,xx(indN),xx(incrementalIndexState(j)));
               tmpJ(indN,incrementalIndexState(j)) = Jn(:,1);
               tmpJ(indN,indN) = Jn(:,2:end);
           end
           
       end
    end
    
    
    
    J{i} = tmpJ;
    
end

if size(x,2) == 1
    JJ = J{1};
    J = [];
    J = JJ;
end
