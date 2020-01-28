%Emmanuel Asante-Asamani
% 02-26-15

function [cputime,u_soln] = enzymekinetics2Dbdf2(step,n,tol)
% enzymekinetics2Dbdf2(n,step,tol) solves the Enzyme Kinetics Equation
% u_t = u_xx -u/(1+u); u(0,t) =0; u(1,t) = 1; 
% u(x,0) = 1; t>0.
%
% with the BDF2 discretization. Newton's method is used
% for the nonlinear solver. The Jacobian is tridiagonal, so we
% use the banded differencing function.
%
% The value of u at the current time and the time step are passed
% to the nonlinear residual as MATLAB global variables.
%
% INPUTS: step: time step
%         n: number of internal nodes
%         tol: tolerance for nonlinear residual
% OUTPUTS: cputime: duration of simulation
%          error: Linfinity error
%***********************************************************************
close all;clc

global uold dt M1 M2 uold1 

dt = step;

% diffusion coefficient
d1 = 0.2; d2 = 0.2; 

% create nodes
nnodes = n^2;
x = linspace(0,1,n+2); y=x;
xint = x(2:n+1); yint=xint;
h = abs(x(1)-x(2)); 
nodes = zeros(nnodes,2);
j = 1;
for k = 1 : n
        for i = 1:n
               nodes(j,:) = [xint(i) yint(k)];
            j = j+1;
        end
end

% discretize time interval
t = 0:dt:1; nt = length(t);

% initial condition for w_old
%smooth
uold = zeros(nnodes,1);
for i = 1: nnodes
 pn = nodes(i,:); % extract point
 xn = pn(1);
 yn = pn(2);
 uold(i) = sin(pi*xn)*sin(pi*yn);
end

% discontinuity at boundary
uold = ones(nnodes,1);

%% Assemble diffusion matrix
e = ones(n,1); II = eye(nnodes);I = eye(n);
r1=d1/(h^2); r2 = d2/(h^2);
As1 = spdiags([-r1*e 2*r1*e -r1*e], -1:1, n, n);
As2 = spdiags([-r2*e 2*r2*e -r2*e], -1:1, n, n);
A1 = kron(I,As1); A2 = kron(As2,I); A = A1+A2;
M1 = sparse(II+(2/3)*dt*A); M2 = sparse(II+dt*A);

%% Use tight tolerances, Newton's method, and a tridiagonal Jacobian.
tol=[tol,tol];
%parms=[40, 1, 0, 1, 1, 1]; % Uses Newtons method
parms=[40, 1000, 0.5, 1]; % uses Modified Newtons method

it_hist = zeros(nt-1,1); stoptol = it_hist;
tic
% backward euler step for u_1
     [unew,ithist,~,stol] = nsold(uold,'fenzymebe',tol,parms);
      it_hist(1) = ithist;
      stoptol(1) = stol;
      uold1 = uold;
      uold = unew;
      
% bdf2 step from here on      
for it=2:nt-1
      [unew, ithist,~,stol] = nsold(uold,'fenzymebdf2',tol,parms);
       it_hist(it) = ithist;
       stoptol(it) = stol;
       uold1 = uold;
       uold=unew; 
end
cputime = toc;
tolvec = stoptol;
u_soln= unew;

%% Plot the results.
% U = zeros(n+2);
% U(2:n+1,2:n+1) = reshape(u_soln,n,n);
% 
% figure
% surf(x,y,U')
% xlabel('x')
% ylabel('y')
% zlabel('U')
% title('\bf\fontsize{20} U solution ')
% 
%  figure
%  plot(it_hist,'b')
%  hold on
%  plot(tolvec,'r')
%  grid on
%  legend('fnorm','tolerance')
