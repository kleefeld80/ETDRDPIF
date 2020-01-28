function [runtime,zsoln] = Shrodinger_2DN_IFETDRDP(dt,steps)
% dt: time step
% steps: number of interior spatial points in each coordinate direction

%% Model Paramters and initial conditions 
% create nodes
nnodes = steps^2;
x = linspace(0,1,steps); y=x;
xint = x; yint=y;

h = abs(x(1)-x(2))
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

% initial condition for z_old
vold = cos(pi*nodes(:,1)).*cos(pi*nodes(:,2));
zold = zeros(2*nnodes,1);
zold(1:2:2*nnodes-1)=vold;

%%  matrix Assembly
% 1D regular diffusion matrix with Neumann boundary
e = ones(steps,1);
A = spdiags([e -2*e e], -1:1, steps, steps);
A(1,2) = 2;
A(steps,steps-1) = 2;

% 2x2 matrix for reversing solution ordering
B = zeros(2);
B(1,2)=1; B(2,1)=-1;

% 2D stacked diffusion matrix for each coordinate direction
I1 = sparse(eye(steps));
r=1/h^2;
Ex = r*kron(I1,kron(A,B));
Ey = r*kron(A,kron(I1,B));

% Final matrices for U 
I2 = sparse(eye(2*nnodes));
r1 = 1/3; r2 = 1/4;
Ex1 = (I2 + dt*Ex);
Ex3 = (I2 + r1*dt*Ex);
Ex4 = (I2 + r2*dt*Ex);
Ey1 = (I2 + dt*Ey);
Ey3 = (I2 + r1*dt*Ey);
Ey4 = (I2 + r2*dt*Ey);

%% Time Evolution
tic
[L1,U1]= lu(Ex1); [L2,U2]=lu(Ex3); [L3,U3]=lu(Ex4);
[L11,U11]=lu(Ey1); [L22,U22]=lu(Ey3); [L33,U33]=lu(Ey4);
for i = 2:tlen
     Fold = F(zold);
     zstar = U11\(L11\(zold + dt*Fold));
     zstar = U1\(L1\zstar);
     Fstar = F(zstar);
     a_1 = U2\(L2\zold);
     b_1 = U3\(L3\zold);
     c_1 = 9*a_1 - 8*b_1;
     a_2 = U2\(L2\Fold);
     b_2 = U3\(L3\Fold);
     c_2 = 9*a_2-8*b_2;
     d_1 = U22\(L22\(9*c_1+2*dt*c_2 + dt*Fstar));
     d_2 = U33\(L33\(-8*c_1-(3/2)*dt*c_2-(dt/2)*Fstar));
     zold = d_1+d_2;
         
     
end
runtime = toc;
%% Plots
vsoln  = zold(1:2:2*nnodes-1);
wsoln = zold(2:2:2*nnodes);
zsoln=vsoln+1i*wsoln;

zexact=exp(-1i*t(end))*cos(pi*nodes(:,1)).*cos(pi*nodes(:,2));
norm(zsoln-zexact,2)

zsoln = sqrt(vsoln.^2+wsoln.^2);
zexact = sqrt(cos(t(tlen)).^2.*cos(pi*nodes(:,1)).^2.*cos(pi*nodes(:,2)).^2 ...
+ sin(t(tlen)).^2.*(cos(pi*nodes(:,1)).^2.*cos(pi*nodes(:,2)).^2));
norm(zsoln-zexact,2)
zmat = reshape(zsoln,steps,steps);
zexact = reshape(zexact,steps,steps);

figure
subplot(2,2,1)
surf(x,y,zmat')
xlabel('x')
ylabel('y')
zlabel('|\Psi|')

subplot(2,2,2)
surf(x,y,zexact')
xlabel('x')
ylabel('y')
zlabel('|\Psi|')

subplot(2,2,3)
contourf(x,y,zmat')
 title('\bf\fontsize{20} Numerical solution ')
 colorbar

subplot(2,2,4)
contourf(x,y,zexact')
 title('\bf\fontsize{20} Exact solution ')
colorbar


%****************function calls**************************************
function Fr = F(soln)
 Fr = zeros(2*nnodes,1);
 v = soln(1:2:2*nnodes-1);
 w = soln(2:2:2*nnodes);
 fdiag = -(1-2*pi^2)*(1-cos(pi*nodes(:,1)).^2.*cos(pi*nodes(:,2)).^2+(v.^2+w.^2));
 Fr(1:2:2*nnodes-1) = -fdiag.*w;
 Fr(2:2:2*nnodes) = fdiag.*v;
end

end