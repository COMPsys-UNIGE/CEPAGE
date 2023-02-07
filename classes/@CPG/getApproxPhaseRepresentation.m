function [TphiApprox,phiApprox] = getApproxPhaseRepresentation(object,PRC,limitCycle,g_in,g_ex,g_el,Tspan,phi0,varargin)

EsynEx = object.EsynEx;
EsynIn = object.EsynIn;

inhActivation = object.inhActivation;
excActivation = object.excActivation;

N = object.N;
Q = PRC.PRC(:,1);
T = limitCycle.T(end);

[TphiApprox,phiApprox] = ode45(@vectorialField, Tspan, phi0);

    function xdot = vectorialField(t,x)   
        xdot = zeros(N-1,1);
        phi_current = [0;x];
        phi_current = mod(phi_current,1);
        Q_current = interp1(PRC.phi,Q,phi_current,'spline');
        V_current = interp1(limitCycle.T,limitCycle.X(:,1),phi_current*T,'spline');
        contribIn = 0;
        contribEx = 0;
        contribEl = 0;
        for j=1:N
            contribIn = contribIn + g_in(1,j)*(-V_current(1)+EsynIn)*inhActivation{1,1}.getActivation(V_current(j));
            contribEx = contribEx + g_ex(1,j)*(-V_current(1)+EsynEx)*excActivation{1,1}.getActivation(V_current(j));
            contribEl = contribEl + g_el(1,j)*(V_current(j)-V_current(1)); 
        end
        constantContrib = contribIn+contribEx+contribEl;

        for i=2:N
            tmpIn = 0;
            tmpEx = 0;
            tmpEl = 0;
            for j=1:N
                tmpIn = tmpIn + g_in(i,j)*(-V_current(i)+EsynIn)*inhActivation{1,1}.getActivation(V_current(j));
                tmpEx = tmpEx + g_ex(i,j)*(-V_current(i)+EsynEx)*excActivation{1,1}.getActivation(V_current(j));
                tmpEl = tmpEl + contribEl + g_el(i,j)*(V_current(j)-V_current(i));             
            end
            tmp = tmpIn+tmpEx+tmpEl;
            xdot(i-1) = -Q_current(1)*constantContrib+Q_current(i)*tmp;
        end
    end
end