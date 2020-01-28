function ft=fbrusselatorbe(u)
% fbrusselatorbe
% Nonlinear residual for time-dependent brusselator.
% This code has the neumann boundary condition built in.
% The time step and solution are passed as globals.
%
global W_old dt M1 A B nb

 d2u = M1*u;
 funu = zeros(nb,1);
 u_1 = u(1:2:nb-1); u_2 = u(2:2:nb);
 f1 = A+u_1.^2.*u_2 -(B+1)*u_1;
 f2 = B*u_1-u_1.^2.*u_2;
 funu(1:2:nb-1) = f1; funu(2:2:nb) = f2;
 
% Nonlinear residual for implicit Euler discretization.
ft= d2u - W_old - dt*funu;

