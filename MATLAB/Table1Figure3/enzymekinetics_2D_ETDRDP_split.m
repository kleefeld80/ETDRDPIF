%% Solution to 2D enzyme kinetics using the new LOD-RRP by Prof Wade
% this algorithm uses RDP
% E.O Asante-Asamani
% 10/21/2014

function [runtime,w_old] = enzymekinetics_2D_ETDRDP_split(dt,steps)
 
clc; close all;
% dt: time step
% steps: number of interior spatial points in each coordinate direction

%% Model Paramters and initial conditions

% diffusion coefficient
d = 0.2; 

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

% initial condition for w_old(continuous)
w_old = zeros(nnodes,1);
for i = 1: nnodes
 pn = nodes(i,:); % extract point
 xn = pn(1);
 yn = pn(2);
 w_old(i) = sin(pi*xn)*sin(pi*yn);
end
w_old = ones(nnodes,1);
%%  matrix Assembly
e = ones(steps,1); II = eye(nnodes);I = speye(steps);
A = (d/h^2)*spdiags([-e 2*e -e], -1:1, steps, steps); % finite difference diffusion operator
A1 = kron(I,A);  % for x direction
A2 = kron(A,I);  % for y direction
M1 = sparse(II+dt*A1); M2 = sparse(II+(dt/3)*A1); M3 = sparse(II+(dt/4)*A1);
M11 = sparse(II+dt*A2); M22 = sparse(II+(dt/3)*A2); M33 = sparse(II+(dt/4)*A2);

%% Time Evolution
[L1,U1]= lu(M1); [L2,U2]=lu(M2); [L3,U3]=lu(M3);
[L11,U11]=lu(M11); [L22,U22]=lu(M22); [L33,U33]=lu(M33);
tic
for i = 2:tlen
   %  waitbar(i/(tlen),hw)     
     F_old = F(w_old);
     w_star = U11\(L11\(w_old + dt*F_old));
     w_star = U1\(L1\w_star);
     F_star = F(w_star);
     a_1 = U2\(L2\w_old);
     b_1 = U3\(L3\w_old);
     c_1 = 9*a_1 - 8*b_1;
     a_2 = U2\(L2\F_old);
     b_2 = U3\(L3\F_old);
     c_2 = 9*a_2-8*b_2;
     d_1 = U22\(L22\(9*c_1+2*dt*c_2 + dt*F_star));
     d_2 = U33\(L33\(-8*c_1-(3/2)*dt*c_2-(dt/2)*F_star));
     w_old = d_1+d_2;

% without lu
%      F_old = F(w_old);
%      w_star = M11\(w_old + dt*F_old);
%      w_star = M1\w_star;
%      F_star = F(w_star);
%      a_1 = M2\w_old;
%      b_1 = M3\w_old;
%      c_1 = 9*a_1 - 8*b_1;
%      a_2 = M2\F_old;
%      b_2 = M3\F_old;
%      c_2 = 9*a_2-8*b_2;
%      d_1 = M22\(9*c_1+2*dt*c_2 + dt*F_star);
%      d_2 = M33\(-8*c_1-(3/2)*dt*c_2-(dt/2)*F_star);
%      w_old = d_1+d_2; 

end
runtime = toc;
%% Plots
% U = zeros(steps+2);
% U(2:steps+1,2:steps+1) = reshape(w_old,steps,steps);
% figure
% surf(x,y,U')
% xlabel('x')
% ylabel('y')
% zlabel('U')

% figure
% contourf(x,y,U')
%  title('\bf\fontsize{20} U solution ')
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