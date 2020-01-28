%Emmanuel Asante-Asamani
% 02-26-15

function [cputime,u_soln] = brusselator2Dbdf2(step,n,tol)
% brusselator2Dbdf2(step,n,tol) solves the Brusselator Equation
% u_t = eps1u_xx + A + u^2v-(B+1)u; 
% v_t = eps2v_xx + Bu - u^2v
% u(,t) =0; u(1,t) = 1; 
% u(x,y0) = 1/2 +y; t>0.
% v(x,y,0) = 1 + 5x;
% dudn = 0, dvdn = 0;
%
% with the BDF2 discretization. Newton's method is used
% for the nonlinear solver. The Jacobian is tridiagonal, so we
% use the banded differencing function.
%
% The value of u at the current time and the time step are passed
% to the nonlinear residual as MATLAB global variables.
%
% INPUTS: step: time step
%         n: number of spatial points in each coordinate direction
%         tol: tolerance for nonlinear residual
% OUTPUTS: cputime: duration of simulation
%          error: Linfinity error
%***********************************************************************
close all;clc

global W_old dt M1 M2 W_old1 A B nb
dt = step;
A = 1; B = 3.4;

% diffusion coefficient
epsln = 2.e-3; 

% create nodes
x = linspace(0,1,n); h = abs(x(1)-x(2)); 
y = x;
nnodes = n^2;
nodes = zeros(nnodes,2);
j = 1;
for k = 1 : n
        for i = 1:n
               nodes(j,:) = [x(i) y(k)];
            j = j+1;
        end
end
nb = 2*nnodes; % becuase we are solving a system of 2 RDE

% discretize time interval
t = 0:dt:2; nt = length(t);

% initial condition for u
u_old = 0.5 + nodes(:,2); 

% initial condition for v
v_old = 1 + 5*nodes(:,1); 

% Stacking nodes for evolution
W_old = zeros(nb,1);
W_old(1:2:nb-1) = u_old; W_old(2:2:nb) = v_old; 

%% Assemble diffusion matrix
ns = 2*n;
Z = zeros(2);
I = eye(2);
C = zeros(2);
C(1,1) = (epsln*dt)/h^2;
C(2,2) = (epsln*dt)/h^2;
I1 = blktridiag(-2*C,Z,Z,n);
I2 = blktridiag(-C,Z,Z,n);

Aq1 = (I+4*C);
Q1 = blktridiag(Aq1,-C,-C,n);
Q1(1:2,3:4) = -2*C; Q1(ns-1:ns,ns-3:ns-2) = -2*C;
M1 = blktridiag(Q1,I2,I2,n);
M1(1:ns,ns+1:2*ns)=I1; 
M1(nb-(ns-1):nb,nb-(2*ns-1):nb-ns)=I1;

Aq2 = (I+(8/3)*C);
Q2 = blktridiag(Aq2,-(2/3)*C,-(2/3)*C,n);
Q2(1:2,3:4) = -(4/3)*C; Q2(ns-1:ns,ns-3:ns-2) = -(4/3)*C;
M2 = blktridiag(Q2,(2/3)*I2,(2/3)*I2,n);
M2(1:ns,ns+1:2*ns)=(2/3)*I1; 
M2(nb-(ns-1):nb,nb-(2*ns-1):nb-ns)=(2/3)*I1;


%% Use tight tolerances, Newton's method, and a tridiagonal Jacobian.
tol=[tol,tol];
%parms=[40, 1, 0, 1, 1, 1]; % Uses Newtons method
parms=[40, 1000, 0.5, 1]; % uses Modified Newtons method

it_hist = zeros(nt-1,1); stoptol = it_hist;

tic
% backward euler step for u_1
     [unew,ithist,~,stol] = nsold(W_old,'fbrusselatorbe',tol,parms);
      it_hist(1) = ithist;
      stoptol(1) = stol;
      W_old1 = W_old;
      W_old = unew;
      
% bdf2 step from here on      
for it=2:nt-1
      [unew, ithist,~,stol] = nsold(W_old,'fbrusselatorbdf2',tol,parms);
       it_hist(it) = ithist;
       stoptol(it) = stol;
       W_old1 = W_old;
       W_old=unew; 
end
cputime = toc;
tolvec = stoptol;
u_soln = W_old(1:2:nb-1); 
v_soln = W_old(2:2:nb);
U = reshape(u_soln,n,n); V = reshape(v_soln,n,n);

%% Plot the results.
% contourf(x,y,U')
%  %title('\bf\fontsize{14} U solution ')
% 
% figure
% contourf(x,y,V')
%  %title('\bf\fontsize{14} V solution ')
% 
%  figure
%  plot(it_hist,'b')
%  hold on
%  plot(tolvec,'r')
%  grid on
%  legend('fnorm','tolerance')
