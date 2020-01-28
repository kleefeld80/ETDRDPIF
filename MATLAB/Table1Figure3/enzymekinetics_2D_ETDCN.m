%% Solution to 2D enzyme kinetics using ETDCN
% Ref:Bhatt, Khaliq 2014
% E.O Asante-Asamani
% 09/16/2014

function [runtime,w_old] = enzymekinetics_2D_ETDCN(dt,steps)
 
clc; close all
% dt: time step
% steps: number of interior spatial points in each coordinate direction

%% Model Paramters and initial conditions

% diffusion coefficient
d1 = 0.2; d2 = 0.2; 
%eps1 = 2; eps2=2;

% create nodes
nnodes = steps^2;
x = linspace(0,1,steps+2); y=x;
xint = x(2:steps+1); yint=xint;

h = abs(x(1)-x(2)); 
nodes = zeros(nnodes,2);

j = 1;
for k = 1 : steps
        for i = 1:steps
               nodes(j,:) = [xint(i) yint(k)];
            j = j+1;
        end
end

% discretize time interval
t = 0:dt:1; tlen = length(t);

% initial condition for w_old

% smooth
% w_old = zeros(nnodes,1);
% for i = 1: nnodes
%  pn = nodes(i,:); % extract point
%  xn = pn(1);
%  yn = pn(2);
%  w_old(i) = sin(pi*xn)*sin(pi*yn);
% end

% discontinuity at boundary
w_old = ones(nnodes,1);

% discontinuity at boundary and mid section of grid
% w_old = zeros(nnodes,1);
% for i = 1:nnodes
%     if nodes(i,2) < 0.5
%         w_old(i)=-1;
%     else
%         w_old(i)=1;
%     end
% end
%plot(w_old);

%%  matrix Assembly
e = ones(steps,1); II = eye(nnodes);I = eye(steps);
As1 = (d1/h^2)*spdiags([-e 2*e -e], -1:1, steps, steps); % finite difference diffusion operator
As2 = (d2/h^2)*spdiags([-e 2*e -e], -1:1, steps, steps); % finite difference diffusion operator
A1 = kron(I,As1); A2 = kron(As2,I); A = A1+A2;
M1 = sparse(II+0.5*dt*A); M2 = sparse(II-0.5*dt*A); M3 = sparse(II+dt*A);

% % for advective system
%   e = ones(steps,1);
%  a1 = (eps2/h + d/h^2); a2 = (eps1/h+d/h^2); a3 = (eps1/h-d/h^2); a4= (eps2/h-d/h^2); a5 = d/h^2;
%  I1 = eye(steps)*a1; I2 = eye(steps)*a4;
%  As = spdiags([-e*a2 4*a5*e e*a3], -1:1, steps, steps); % finite difference diffusion operator
%  A = blktridiag(As,-I1,I2,steps); Id = sparse(eye(nnodes));
%  M1 = sparse(Id+0.5*dt*A); M2 = sparse(Id-0.5*dt*A); %M3 = sparse(Id+dt*A);

%% Time Evolution
%hw = waitbar(0,'Simulating...');
tic
for i = 2:tlen
     %waitbar(i/(tlen),hw) 
      F_old = F(w_old);
%      if i<6
%          w_old = M3\(w_old + dt*F_old);
%      else
%          w_b= M1\(M2*w_old + dt*F_old);
%          w_a = M1\(0.5*dt*(F(w_b)-F_old));
%          w_old = w_a+w_b;  
%      end
        w_b= M1\(M2*w_old + dt*F_old);
        w_a = M1\(0.5*dt*(F(w_b)-F_old));
        w_old = w_a+w_b; 
end
runtime = toc;
%% Plots
U = zeros(steps+2);
U(2:steps+1,2:steps+1) = reshape(w_old,steps,steps);

figure
surf(x,y,U')
xlabel('x')
ylabel('y')
zlabel('U')
title('\bf\fontsize{20} U solution ')
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