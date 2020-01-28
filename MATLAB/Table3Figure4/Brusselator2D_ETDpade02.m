%% Solution to 2D Brusselator Using ETD Pade 02 scheme
% E.O Asante-Asamani
% 05/08/2014

function [runtime,u_soln] = Brusselator2D_ETDpade02(dt,steps)
 clc;
% dt: time step (0.01)
% steps: number of spatial points in each coordinate direction (51)

%% Model Paramters and initial conditions
A = 1; B = 3.4;

% diffusion coefficient
epsln = 2.e-3; 

% create nodes
x = linspace(0,1,steps); h = abs(x(1)-x(2)); 
y = linspace(0,1,steps);
nnodes = steps^2;
nodes = zeros(nnodes,2);
j = 1;
for k = 1 : steps
        for i = 1:steps
               nodes(j,:) = [x(i) y(k)];
            j = j+1;
        end
end
nb = 2*nnodes;

% discretize time interval
t = 0:dt:2; tlen = length(t);

% initial condition for u
u_old = 0.5 + nodes(:,2); 

% initial condition for v
v_old = 1 + 5*nodes(:,1); 

% Stacking nodes for evolution
W_old = zeros(nb,1);
W_old(1:2:nb-1) = u_old; W_old(2:2:nb) = v_old; 

%% Block matrix Assembly
w1 = 0.5; w2=complex(0.5,-0.5); cm = complex(-1,1); w=complex(0,-1);
ns = 2*steps;
Z = zeros(2);
I = eye(2);
C = zeros(2);
C(1,1) = (epsln*dt)/h^2;
C(2,2) = (epsln*dt)/h^2;
I1 = blktridiag(-2*C,Z,Z,steps);
I2 = blktridiag(-C,Z,Z,steps);
Aq = (4*C-cm*I);
Q = blktridiag(Aq,-C,-C,steps);
Q(1:2,3:4) = -2*C; Q(ns-1:ns,ns-3:ns-2) = -2*C;
M = blktridiag(Q,I2,I2,steps);
M(1:ns,ns+1:2*ns)=I1; 
M(nb-(ns-1):nb,nb-(2*ns-1):nb-ns)=I1;

%% Time Evolution 
[L,U]=lu(M);

%hw = waitbar(0,'Simulating...');
for i = 2:tlen
     %waitbar(i/(tlen),hw) 
     Na = U\(L\(w*W_old + w2*dt*F(W_old)));
     b = 2*real(Na);
     Nb = U\(L\(w1*dt*(F(b)-F(W_old))));
     W_old = b + 2*real(Nb);
end
% uncomment this section to display solution
%****************************************************************

 u_soln = W_old(1:2:nb-1); 
%  v_soln = W_old(2:2:nb);
%  U = reshape(u_soln,steps,steps); V = reshape(v_soln,steps,steps);

%*******************************************************************
runtime = toc;

%% Plots
% uncomment this section to display solution
%*****************************************************************
% contourf(x,y,U')
%  title('\bf\fontsize{20} U solution ')
% 
% figure
% contourf(x,y,V')
%  title('\bf\fontsize{20} V solution ')

%******************************************************************






%****************function calls**************************************
function Fr = F(U)
 Fr = zeros(nb,1);
 u = U(1:2:nb-1); v = U(2:2:nb);
 f1 = A+u.^2.*v -(B+1)*u;
 f2 = B*u-u.^2.*v;
 Fr(1:2:nb-1) = f1; Fr(2:2:nb) = f2;
end



end