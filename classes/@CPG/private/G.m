function res = G(object,V,i,j)

N = object.N;
g_in = object.g_in;
g_ex = object.g_ex;
g_el = object.g_el;

inhActivation = object.inhActivation;
excActivation = object.excActivation;


%Valori sinapsi
EsynIn = object.EsynIn;
EsynEx = object.EsynEx;

res = 0;

res = res + g_ex(i,j)*(-V(i)+EsynEx)*excActivation{i,j}.getActivation(V(j));

res = res + g_in(i,j)*(-V(i)+EsynIn)*inhActivation{i,j}.getActivation(V(j));

res = res + g_el(i,j)*(-V(j)+V(i));

