function h = plotPhaseSpace(Tphi,phi)
% plotPhaseSpace Plot the time evolution of the phase difference of the
%                neurons in the network described by variable phi. If
%                network is composed by 3 or 4 neurons the function plot
%                also the phase space evolution of the network.
%
% h = simplot(object,phi)
% h is the is the figure(s) handle. phi must be a cell array; each element
% of phi must be a matrix describing the phase different evolution of the
% network
%
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

if ~iscell(Tphi) || ~iscell(phi)
    error('Tphi and phi must be cell arrays');
end

nCI = numel(phi);

soglia = 0.8; % threshold for turn over torus 

N = size(phi{1},2)+1;

possibleFP = zeros(nCI,N-1);

eps = 5e-2;

for i=1:nCI
    possibleFP(i,:) = phi{i}(end,:);
end

roundedFP = exp(1i*2*pi*possibleFP);
roundedFP = round(roundedFP/eps)*eps;
[~,IA,~] = unique(roundedFP,'rows');
stableFP = possibleFP(IA,:);



Nstab = size(stableFP,1);
cmap = hsv(N-1);

colorBasin = hsv(Nstab);


if N == 3
    h = cell(2,1);
    h{2} = figure;
    hold on
    for i=1:nCI
        
        [~,indexColor] = min(sum(abs(exp(1i*2*pi*repmat(phi{i}(end,:),Nstab,1)) - exp(1i*2*pi*stableFP)),2));
        
        differenza = diff(phi{i});
        ii = [1; find(any(abs(differenza) > soglia,2));size(phi{i},1)];
        
        for k=1:numel(ii)-1
            plot(phi{i}(ii(k)+1:ii(k+1),1),phi{i}(ii(k)+1:ii(k+1),2),'color',colorBasin(indexColor,:))
        end
        
    end
    xlabel('\Delta \phi_{12}');
    ylabel('\Delta \phi_{13}');
    axis([0 1 0 1]);
    
%     h{3} = figure;
%     hold on
%     
%     R = 11;
%     r = 1;
%     
%     u=linspace(0,2*pi,100);
%     v=linspace(0,2*pi,100);
%     
%     [u,v]=meshgrid(u,v);
%     
%     x=(R+r.*cos(v)).*cos(u);
%     y=(R+r.*cos(v)).*sin(u);
%     z=r.*sin(v);
%     
%     m = surf(x,y,z);
%     
%     set(m,'FaceColor',[0.7 0.7 0.7],'FaceAlpha',0.7,'edgecolor','none');
%     %     alpha(m,0)
%     
%     for i=1:nCI
%         u = 2*pi*phi{i}(:,1);
%         v = 2*pi*phi{i}(:,2);
%         
%         x=(R+r.*cos(v)).*cos(u);
%         y=(R+r.*cos(v)).*sin(u);
%         z=r.*sin(v);
%         [~,indexColor] = min(sum(abs(repmat(phi{i}(end,:),Nstab,1) - stableFP)));
%         
%         l = plot3(x,y,z,'r','linewidth',2,'color',colorBasin(indexColor,:));
%     end
%     
%     % plot axes
%     zer = zeros(1,100);
%     ls = linspace(0,1,100);
%     
%     u = 2*pi*zer;
%     v = 2*pi*ls;
%     
%     x=(R+r.*cos(v)).*cos(u);
%     y=(R+r.*cos(v)).*sin(u);
%     z=r.*sin(v);
%     
%     
%     plot3(x,y,z,'k','linewidth',2);
%     
%     u = 2*pi*ls;
%     v = 2*pi*zer;
%     
%     x=(R+r.*cos(v)).*cos(u);
%     y=(R+r.*cos(v)).*sin(u);
%     z=r.*sin(v);
%     
%     
%     plot3(x,y,z,'k','linewidth',2);
    
elseif N == 4
    h = cell(2,1);
    h{2} = figure;
    hold on
    for i=1:nCI
        [~,indexColor] = min(sum(abs(repmat(phi{i}(end,:),Nstab,1) - stableFP)));
        differenza = diff(phi{i});
        ii = [1; find(any(abs(differenza) > soglia,2));size(phi{i},1)];
        
        for k=1:numel(ii)-1
            plot3(phi{i}(ii(k)+1:ii(k+1),1),phi{i}(ii(k)+1:ii(k+1),2),phi{i}(ii(k)+1:ii(k+1),3),'color',colorBasin(indexColor,:))
        end
        
    end
    xlabel('\Delta \phi_{12}');
    ylabel('\Delta \phi_{13}');
    zlabel('\Delta \phi_{14}');
    axis([0 1 0 1 0 1]);
    view(-30,35)
else
    h = cell(1,1);
end

h{1} = figure;
hold on
for j=1:N-1
    plot(0,0,'color',cmap(j,:))
end

for i=1:nCI
    for j=1:N-1
        differenza = diff(phi{i}(:,j));
        ii = [1; find(abs(differenza) > soglia);size(phi{i},1)];
        
        for k=1:numel(ii)-1
            plot(Tphi{i}(ii(k)+1:(ii(k+1))),phi{i}(ii(k)+1:ii(k+1),j),'color',cmap(j,:));
        end
    end
end

str = '';
for i=2:N-1
    str = ([str,'''neuron ',num2str(i),''',']);
end
str = ([str,'''neuron ',num2str(N),'''']);
eval(['legend(',str,')']);

end

