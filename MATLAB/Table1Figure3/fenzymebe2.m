function ft=fenzymebe2(u)
% FTIME 
% Nonlinear residual for time-dependent problem in Chapter 2.
% This code has the zero boundary conditions built in.
% The time step and solution are passed as globals.
%
global w_oldd k M1

 d2u = M1*u;
 fun = -u./(1+u);
 
% Nonlinear residual for implicit Euler discretization.

ft= d2u - w_oldd - k*fun;

