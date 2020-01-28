function ft=fbrusselatorbdf2(u)
% fbrusselatorbdf2 
% Nonlinear residual for Brusselator
% This code has the zero Neumann boundary conditions built in
% The time step and solution are passed as globals.
% BDF2 is used in the nonlinear function formulation
%
global W_old dt M2 W_old1 A B nb

 d2u = M2*u; 
 funu = zeros(nb,1);
 u_1 = u(1:2:nb-1); u_2 = u(2:2:nb);
 f1 = A+u_1.^2.*u_2 -(B+1)*u_1;
 f2 = B*u_1-u_1.^2.*u_2;
 funu(1:2:nb-1) = f1; funu(2:2:nb) = f2;

% Nonlinear residual for BDF2 discretization.
ft=d2u-(4/3)*W_old + (1/3)*W_old1-(2*dt/3)*funu;