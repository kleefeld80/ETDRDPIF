%% Solution to 2D Brusselator Using ETD RRP with ETD1 predictor scheme
% Ref:Zegeling et al (2004)(see attached paper for details on initial and 
%     boundary conditions
% E.O Asante-Asamani
% 05/08/2014

function [runtime,u_soln] = Brusselator2D_IMEX_TR(dt,steps)
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
ns = 2*steps;
Z = zeros(2);
I = eye(2);
C = zeros(2);D = zeros(2);
C(1,1) = (epsln*dt)/h^2;D(1,1) = (epsln)/h^2;
C(2,2) = (epsln*dt)/h^2;D(2,2) = (epsln)/h^2;
I1 = blktridiag(-2*C,Z,Z,steps);I1d = blktridiag(-2*D,Z,Z,steps);
I2 = blktridiag(-C,Z,Z,steps);I2d = blktridiag(-D,Z,Z,steps);

% matrix 1
A_m = -4*D;
Q = blktridiag(A_m,D,D,steps);
Q(1:2,3:4) = 2*D; Q(ns-1:ns,ns-3:ns-2) = 2*D;
M = blktridiag(Q,-I2d,-I2d,steps);
M(1:ns,ns+1:2*ns)=-I1d; 
M(nb-(ns-1):nb,nb-(2*ns-1):nb-ns)=-I1d;


% matrix 2
A_m = (I+2*C);
Q1 = blktridiag(A_m,-0.5*C,-0.5*C,steps);
Q1(1:2,3:4) = -C; Q1(ns-1:ns,ns-3:ns-2) = -C;
M1 = blktridiag(Q1,0.5*I2,0.5*I2,steps);
M1(1:ns,ns+1:2*ns) = 0.5*I1; 
M1(nb-(ns-1):nb,nb-(2*ns-1):nb-ns) = 0.5*I1;

%% Time Evolution 
% LU decomposition to speed up computation
[L1,U1]=lu(M1);

%hw = waitbar(0,'Simulating...');
tic
 
for i = 1:tlen-1
     F_old = F(W_old);
     M_old = M*W_old;
     W_star = U1\(L1\(W_old + dt*F_old +0.5*dt*M_old));
     F_star = F(W_star);
     M_star = M*W_star;
     % Main Step
     W_old = W_old + 0.5*dt*(F_old + M_old) + 0.5*dt*(F_star + M_star);  
end

 u_soln = W_old(1:2:nb-1);
 v_soln = W_old(2:2:nb);
 U = reshape(u_soln,steps,steps); V = reshape(v_soln,steps,steps);
 
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