function [position,isterminal,direction] = getResetConditions(object,t,y,varargin)
position = [y(1)-20;-1]; % The value that we want to be zero
isterminal = [1;0];  % Halt integration 
direction = [1;0];   % The zero can be approached from either direction

