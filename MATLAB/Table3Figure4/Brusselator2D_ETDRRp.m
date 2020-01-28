%% Solution to 2D Brusselator Using ETD RRP scheme
% E.O Asante-Asamani
% 05/08/2014

function runtime = Brusselator2D_ETDRRp(dt,steps)
 clc;
% dt: time step (0.01)
% steps: number of spatial points in each coordinate direction (51)

%% Model Paramters and initial conditions
A = 1; B = 3.4;

% diffusion coefficient
epsln = 2.e-3; 

% create nodes
x = linspace(0,1,steps); h = abs(x(1)-x(2)); 
y = x;
nnodes = steps^2;
nodes = zeros(nnodes,2);
j = 1;
for k = 1 : steps
        for i = 1:steps
               nodes(j,:) = [x(i) y(k)];
            j = j+1;
        end
end
nb = 2*nnodes; % becuase we are solving a system of 2 RDE

% discretize time interval
t = 0:dt:11; tlen = length(t);

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
C = zeros(2);
C(1,1) = (epsln*dt)/h^2;
C(2,2) = (epsln*dt)/h^2;
I1 = blktridiag(-2*C,Z,Z,steps);
I2 = blktridiag(-C,Z,Z,steps);
% matrix 1
Aq1 = (I + (4/3)*C);
Q1 = blktridiag(Aq1,-(1/3)*C,-(1/3)*C,steps);
Q1(1:2,3:4) = -(2/3)*C; Q1(ns-1:ns,ns-3:ns-2) = -(2/3)*C;
M1 = blktridiag(Q1,(1/3)*I2,(1/3)*I2,steps);
M1(1:ns,ns+1:2*ns)=(1/3)*I1; 
M1(nb-(ns-1):nb,nb-(2*ns-1):nb-ns)=(1/3)*I1;
% matrix 2
Aq2 = (I + C);
Q2 = blktridiag(Aq2,-(1/4)*C,-(1/4)*C,steps);
Q2(1:2,3:4) = -(1/2)*C; Q2(ns-1:ns,ns-3:ns-2) = -(1/2)*C;
M2 = blktridiag(Q2,(1/4)*I2,(1/4)*I2,steps);
M2(1:ns,ns+1:2*ns)=(1/4)*I1; 
M2(nb-(ns-1):nb,nb-(2*ns-1):nb-ns)=(1/4)*I1;

%% Time Evolution 
[L1,U1] =lu(M1); [L2,U2]=lu(M2);
%hw = waitbar(0,'Simulating...');
tic
for i = 2:tlen
    % waitbar(i/(tlen),hw) 
     % RRp  predictor
     a = U1\(L1\(9*W_old + 3*dt*F(W_old)));
     b = U2\(L2\(-8*W_old - 2*dt*F(W_old)));
     c = a+b;
     % Correction
     a = U1\(L1\(dt*(F(c)-F(W_old))));
     b = U2\(L2\((dt/2)*(F(c)-F(W_old))));
     W_old = c+a-b; 
end

% uncomment this section to display solution
%****************************************************************

% u_soln = W_old(1:2:nb-1); v_soln = W_old(2:2:nb);
% U = reshape(u_soln,steps,steps); V = reshape(v_soln,steps,steps);

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