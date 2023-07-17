function [g_v] = RC_g(v)
% function for nonlinear resistors
g_v = exp(40*v) + v - 1;

end