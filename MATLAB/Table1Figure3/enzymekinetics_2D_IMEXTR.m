%% Solution to 2D enzyme kinetics using IMEX-BDF2
% Ref:Bhatt, Khaliq 2014
% E.O Asante-Asamani
% 10/22/2014

function [runtime,wsoln] = enzymekinetics_2D_IMEXTR(dt,steps)
 % INPUTS: n is the number of intervals
%         k is the time step 
%         T is the final evolution time
% OUTPUT:  ftime cpu time for the evolution
%         u_soln is the final solution

clc;close all

%% Model Paramters and initial conditions

% diffusion coefficient
d1 = 0.2; d2 = 0.2; 
%eps1 = 2; eps2=2; % for advection

% create nodes
nnodes = steps^2;
x = linspace(0,1,steps+2); y=x;
xint = x(2:steps+1); yint=xint;

h = abs(x(1)-x(2)); 
nodes = zeros(nnodes,2);

j = 1;
for c = 1 : steps
        for i = 1:steps
               nodes(j,:) = [xint(i) yint(c)];
            j = j+1;
        end
end

% discretize time interval
t = 0:dt:1; tlen = length(t);

% initial condition for w_old(continuous)
w_old = zeros(nnodes,1);
for i = 1: nnodes
 pn = nodes(i,:); % extract point
 xn = pn(1);
 yn = pn(2);
 w_old(i) = sin(pi*xn)*sin(pi*yn);
end

w_old = ones(nnodes,1);

% w_old = zeros(nnodes,1);
% for i = 1:nnodes
%     if nodes(i,2) < 0.5
%         w_old(i)=-1;
%     else
%         w_old(i)=1;
%     end
% end

%%  matrix Assembly
e = ones(steps,1); Id = sparse(eye(nnodes));I = eye(steps);
As1 = (d1/h^2)*spdiags([e -2*e e], -1:1, steps, steps); % finite difference diffusion operator
As2 = (d2/h^2)*spdiags([e -2*e e], -1:1, steps, steps); % finite difference diffusion operator

A1 = kron(I,As1); A2 = kron(As2,I); A = A1+A2;

% for advection
%   e = ones(steps,1);
%  a1 = (eps2/h + d/h^2); a2 = (eps1/h+d/h^2); a3 = (eps1/h-d/h^2); a4= (eps2/h-d/h^2); a5 = d/h^2;
%  I1 = eye(steps)*a1; I2 = eye(steps)*a4;
%  As = spdiags([-e*a2 4*a5*e e*a3], -1:1, steps, steps); % finite difference diffusion operator
%  A = blktridiag(As,-I1,I2,steps); Id = sparse(eye(nnodes));
 %M1 = sparse(Id-gamma*dt*A);
 %M2 = sparse(Id + ddt*A);
% M3 = sparse(Id + 0.5*dt*A); 
A1 = sparse(Id - 0.5*dt*A);
%M1 = sparse(Id-dt*A);
%M2 = sparse(Id + dt*A);
%M3 = sparse(Id + 0.5*dt*A); 
%M4 = sparse(Id - 0.5*dt*A);

% M3 = spdiags([0.5*r*e (1-r)*e 0.5*r*e], -1:1, nx, nx);
% M4 = spdiags([-0.5*r*e (1+r)*e -0.5*r*e], -1:1, nx, nx);
 
%% Time Evolution
%hw = waitbar(0,'Simulating...');
%tol = dt^3;
%tol=[tol,tol];
%parms=[40, 1, 0, 1, 1, 1]; % Uses Newtons method
%parms=[40, 1000, 0.5, 1, 1, 1]; % uses Modified Newtons methodEvolution


tic
[L1,U1] = lu(A1);

% Fully implicit Backward Euler step   
  %[w_old,~,~,~] = nsold(w_old1,'fenzymebe2',tol,parms);
  
% Semi implicit Backward Euler  
  %w_old = M1\(w_old1 + dt*F(w_old1));
  
% % IMEX-theta (0.5)
%     w_old1 = M4\(M3*w_old + dt*F(w_old)); 
 
for i = 1:tlen-1
    F_old = F(w_old);
    A_old = A*w_old;
    w_star = U1\(L1\(w_old + dt*F_old + 0.5*dt*A_old));
    F_star = F(w_star);
    A_star = A*w_star;
    % Full Scheme
    w_old = w_old + 0.5*dt*(F_old+F_star) + 0.5*dt*(A_old +A_star);
end  
runtime = toc;
wsoln = w_old;
%% Plots
% U = zeros(steps+2);
% U(2:steps+1,2:steps+1) = reshape(w_old,steps,steps);
% figure
% surf(x,y,U')
% xlabel('x')
% ylabel('y')
% zlabel('U')
% title('\bf\fontsize{20} U solution ')
%  colorbar
% 
% figure
% contourf(x,y,V')
%  title('\bf\fontsize{20} V solution ')
% colorbar


%****************function calls**************************************
function Fr = F(u)
 Fr = -u./(1+u);
end


end