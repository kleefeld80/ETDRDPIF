function ft=fenzymebe(u)
% FTIME 
% Nonlinear residual for time-dependent problem in Chapter 2.
% This code has the zero boundary conditions built in.
% The time step and solution are passed as globals.
%
global uold dt M2

 d2u = M2*u;
 fun = -u./(1+u);
 
% Nonlinear residual for implicit Euler discretization.

ft= d2u - uold - dt*fun;

