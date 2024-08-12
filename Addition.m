function [output] = Addition(s,C)
% Calculates the value of a state
%e.g. C=[1 1.1 1.3];
%e.g. s=[1 0 1];
output=sum(s(1,1:end).*C);
