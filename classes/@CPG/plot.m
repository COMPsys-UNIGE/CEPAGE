function plot(object)
% plot   Plot the neurons network topology
%
% disp(OBJ)
% plot(object) is the CPG object.

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
g_in = object.g_in;
g_ex = object.g_ex;
g_el = object.g_el;

[tmp,tmp1,inhEl] = find(unique(g_in));
[tmp,tmp1,exEl] = find(unique(g_ex));
[tmp,tmp1,elEl] = find(unique(g_el));

inhEl = sort(inhEl);
exEl = sort(exEl);
elEl = sort(elEl);

nIn = numel(inhEl);
nEx = numel(exEl);
nEl = numel(elEl);

lineWidthIn = linspace(0.5,0.5*2*nIn,nIn);
lineWidthEx = linspace(0.5,0.5*2*nEx,nEx);
lineWidthEl = linspace(0.5,0.5*2*nEl,nEl);

elColor = 'b';
inColor = 'r';
exColor = 'g';

neurCoord = cell(N,1);


R = 4;
ro = 2*pi/N;
markSize = R/8;

figure
hold on

plot(0,0 ,elColor);
plot(0,0,inColor);
plot(0,0,exColor);

for i=1:N
    x = R*cos(ro*i);
    y = R*sin(ro*i);
    neurCoord{i} = struct('x',x,'y',y);
    t = linspace(0,2*pi,1000);
    xx = markSize*cos(t)+x;
    yy = markSize*sin(t)+y;
    %        plot(x,y,'o','markersize',markSize);
    plot(xx,yy,'k');
    text(x,y,num2str(i))
end



%plot el syn
for i=1:N-1
    for j=i:N
        if g_el(i,j) ~= 0
            ii = find(elEl == g_el(i,j));
            X = [neurCoord{i}.x neurCoord{j}.x];
            Y = [neurCoord{i}.y neurCoord{j}.y];
            
            alpha = atan2(diff(Y),diff(X));
            
            X = [neurCoord{i}.x+markSize*cos(alpha),neurCoord{j}.x+markSize*cos(alpha-pi)];
            Y = [neurCoord{i}.y+markSize*sin(alpha),neurCoord{j}.y+markSize*sin(alpha-pi)];
            
            Xtext = mean(X);%neurCoord{i}.x+markSize*cos(pi/2);
            Ytext = mean(Y);%neurCoord{i}.y+markSize*sin(pi/2);
            
            text(Xtext,Ytext,num2str(g_el(i,j)));
            
            h1 = plot(X,Y,elColor,'linewidth',lineWidthEl(ii));
        end
    end
end





%plot inh syn
for i=1:N
    for j=1:N
        if g_in(i,j) ~= 0
            
            ii = find(inhEl == g_in(i,j));
            p1 = neurCoord{i};
            p2 = neurCoord{j};
            
            a = sqrt((p1.x-p2.x)^2+(p1.y-p2.y)^2)/2;
            b = markSize*2/3;
            
            tmp1 = (1-b^2/a^2);
            tmp2 = -2*a;
            tmp3 = b^2-markSize^2+a^2;
            
            xtmp = roots([tmp1,tmp2,tmp3]);
            xtmp = min(xtmp);
            ytmp = sqrt(markSize^2-(xtmp-a)^2);
            
            deltaT = atan2(ytmp*a,xtmp*b);
            
            x0 = (p1.x+p2.x)/2;
            y0 = (p1.y+p2.y)/2;
            
            alpha = atan2(p1.y-p2.y,p1.x-p2.x);
            
            t = linspace(pi+deltaT,2*pi-deltaT,1000);
            
            xx = a*cos(t)*cos(alpha)-b*sin(t)*sin(alpha)+x0;
            yy = a*cos(t)*sin(alpha)+b*sin(t)*cos(alpha)+y0;
            
            Xtext = a*cos(pi/2)*cos(alpha)-b*sin(pi/2)*sin(alpha)+x0;
            Ytext = a*cos(pi/2)*sin(alpha)+b*sin(pi/2)*cos(alpha)+y0;
            
            
            if i > j
                text(Xtext,Ytext,num2str(g_in(i,j)),'HorizontalAlignment','left')
            else
                text(Xtext,Ytext,num2str(g_in(i,j)),'HorizontalAlignment','right')
            end
            h2 = plot(xx,yy,inColor,'linewidth',lineWidthIn(ii));
            plot(xx(end),yy(end),'o','MarkerFaceColor',inColor,'MarkerEdgeColor',inColor);
        end
    end
end


%plot ex syn
for i=1:N
    for j=1:N
        if g_ex(i,j) ~= 0
            
            ii = find(exEl == g_ex(i,j));
            p1 = neurCoord{i};
            p2 = neurCoord{j};
            
            a = sqrt((p1.x-p2.x)^2+(p1.y-p2.y)^2)/2;
            b = markSize*2/3;
            
            tmp1 = (1-b^2/a^2);
            tmp2 = -2*a;
            tmp3 = b^2-markSize^2+a^2;
            
            xtmp = roots([tmp1,tmp2,tmp3]);
            xtmp = min(xtmp);
            ytmp = sqrt(markSize^2-(xtmp-a)^2);
            
            deltaT = atan2(ytmp*a,xtmp*b);
            
            x0 = (p1.x+p2.x)/2;
            y0 = (p1.y+p2.y)/2;
            
            alpha = atan2(p1.y-p2.y,p1.x-p2.x);
            
            a = sqrt((p1.x-p2.x)^2+(p1.y-p2.y)^2)/2;
            b = markSize*2/3;
            
            t = linspace(pi+deltaT,2*pi-deltaT,1000);
            
            xx = a*cos(t)*cos(alpha)-b*sin(t)*sin(alpha)+x0;
            yy = a*cos(t)*sin(alpha)+b*sin(t)*cos(alpha)+y0;
            
            Xtext = a*cos(pi/2)*cos(alpha)-b*sin(pi/2)*sin(alpha)+x0;
            Ytext = a*cos(pi/2)*sin(alpha)+b*sin(pi/2)*cos(alpha)+y0;
            
            
            if i > j
                text(Xtext,Ytext,num2str(g_ex(i,j)),'HorizontalAlignment','left')
            else
                text(Xtext,Ytext,num2str(g_ex(i,j)),'HorizontalAlignment','right')
            end
            
            
            h3 = plot(xx,yy,exColor,'linewidth',lineWidthEx(ii));
            plot(xx(end),yy(end),'s','MarkerFaceColor',exColor,'MarkerEdgeColor',exColor);
        end
    end
end


axis([-5 5 -5 5])

legend('electric','inhibitory','excitatory');

% hh = [];
% stringMsg = {};
% if exist('h1')
%     hh = [hh, h1];
%     stringMsg{numel(stringMsg)+1} = 'electric';
% end
% if exist('h2')
%     hh = [hh, h2];
%     stringMsg{numel(stringMsg)+1} = 'inhibitory';
% end
% if exist('h3')
%     hh = [hh, h3];
%     stringMsg{numel(stringMsg)+1} = 'excitatory';
% end
% legend(hh,stringMsg);
end