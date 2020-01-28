function ft=fenzymebdf2(u)
% fallenbdf2 
% Nonlinear residual for Allen Cahn Equation
% This code has the zero boundary conditions built in
% The time step and solution are passed as globals.
% BDF2 is used in the nonlinear function formulation
%
global uold dt M1 uold1

d2u = M1*u; 
funu = -u./(1+u);

% Nonlinear residual for BDF2 discretization.
ft=d2u - (4/3)*uold + (1/3)*uold1 - (2*dt/3)*funu;