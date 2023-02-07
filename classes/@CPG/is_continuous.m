function cont = is_continuous(object)
% isContinuous   Report if the model is time continuous
%
% cont = isContinuous(object)
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

neurons = object.neurons;

cont = true;

for i=1:N
    if ~neurons{i}.is_continuous()
        cont = false;
        break;
    end
end


for i=1:N
    for j=1:N
        if ~object.inhActivation{i,j}.is_continuous()
            cont = false;
            break;
        end
    end
end

for i=1:N
    for j=1:N
        if ~object.excActivation{i,j}.is_continuous()
            cont = false;
            break;
        end
    end
end



end