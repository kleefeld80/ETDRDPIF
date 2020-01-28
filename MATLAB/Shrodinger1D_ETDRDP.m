function [runtime,zapprox] = Shrodinger1D_ETDRDP(dt,h)

alpha = 0.01;   % 0.01
c = 0.1;     % 0.1
q = 1; 
a = -80; b = 100;  % -80, 100

% create nodes
x =a:h:b; 
xlen = length(x);

% discretize time interval
t = 0:dt:108; tlen = length(t); %108

% initial condition 
v = sqrt(2*alpha/q)*sech(sqrt(alpha)*x).*cos(0.5*c*x);
w = sqrt(2*alpha/q)*sech(sqrt(alpha)*x).*sin(0.5*c*x);

% Stacking nodes for evolution
zmass = zeros(tlen,1);
zenergy = zeros(tlen,1);
solz = zeros(2*xlen,1);
solz(1:2:2*xlen-1) = v;
solz(2:2:2*xlen) = w;
zmod = v.^2 + w.^2;

% Initial mass
mass = 2*sum(zmod(2:xlen-1));
mass = mass + (zmod(1)+zmod(xlen));
zmass(1)=0.5*h*mass;

% Initial energy
zen = zmod.^2;
zmid = ((v(3:xlen)-v(1:xlen-2)).^2 + (w(3:xlen)-w(1:xlen-2)).^2)/(2*h)^2 - zen(2:xlen-1)/2;
zenergy(1) = 0.5*h*(-(zen(1)+zen(xlen))/2+2*sum(zmid));

% 1D regular diffusion matrix with Neumann boundary
e = ones(xlen,1);
A = spdiags([e -2*e e], -1:1, xlen, xlen);
A(1,2) = 2;
A(xlen,xlen-1) = 2;

% 2x2 matrix for reversing solution ordering
B = zeros(2);
B(1,2)=1; B(2,1)=-1;

% 1D stacked diffusion matrix
r=1/h^2;
E = r*kron(A,B);

% Final matrices for U 
I = sparse(eye(2*xlen));
r1 = 1/3; r2 = 1/4;
E1 = (I + dt*E);
E3 = (I + r1*dt*E);
E4 = (I + r2*dt*E);

%% Time Evolution 
tic
for i = 2:tlen
   % RDP implementation
     F_old = F(solz);
     zstar = E1\(solz + dt*F_old); 
     % Main RDP code  
     F_star = F(zstar);
     a_3 = E3\(9*solz+2*dt*F_old+dt*F_star);
     b_3 = E4\(-8*solz-(3/2)*dt*F_old-0.5*dt*F_star);
     solz = a_3 + b_3; 
     
     % Calculate Mass
     v = solz(1:2:2*xlen-1);
     w = solz(2:2:2*xlen);
     zmod = v.^2 + w.^2;
     mass = 2*sum(zmod(2:xlen-1));
     mass = mass + (zmod(1)+zmod(xlen));
     zmass(i)=0.5*h*mass;
     
     % Calculate Energy
     zen = zmod.^2;
     zmid = ((v(3:xlen)-v(1:xlen-2)).^2 + (w(3:xlen)-w(1:xlen-2)).^2)/(2*h)^2 - zen(2:xlen-1)/2;
     zenergy(i) = 0.5*h*(-(zen(1)+zen(xlen))/2+2*sum(zmid));
     
end
runtime = toc;
zapprox = sqrt(solz(1:2:2*xlen-1).^2 + solz(2:2:2*xlen).^2);
    

% Exact solution
vexact = (sqrt(2*alpha/q)*sech(sqrt(alpha)*(x-c*t(tlen))).*cos(0.5*c*x-(0.25*c^2-alpha)*t(tlen)))';
wexact = (sqrt(2*alpha/q)*sech(sqrt(alpha)*(x-c*t(tlen))).*sin(0.5*c*x-(0.25*c^2-alpha)*t(tlen)))';
zexact = sqrt(vexact.^2 + wexact.^2);

figure(1)
hold on
plot(x,solz(1:2:2*xlen-1),'r')
plot(x,solz(2:2:2*xlen),'b')
plot(x,vexact,'r--')
plot(x,wexact,'b--')
hold off
%norm([vexact-solz(1:2:2*xlen-1);wexact-solz(2:2:2*xlen)],2)
%norm(zexact-zapprox,2)
%error(solz(1:2:2*xlen-1),vexact)
%pause
zmass(1)
zmass(end)
zenergy(1)
zenergy(end)
%% Plots
% Organize solution 
figure
subplot(1,3,1)
plot(x,zapprox,'r','LineWidth',2)
hold on
plot(x,zexact,'b','LineWidth',2)
legend('ETDRDP','Exact')
hold off
xlabel('x')
ylabel('|u|')

subplot(1,3,2)
plot(t,zmass,'k','LineWidth',2)
xlabel('t')
ylabel('M(t)')
title('Mass')

subplot(1,3,3)
plot(t,zenergy,'k','LineWidth',2)
xlabel('t')
ylabel('E(t)')
title('Energy')


%****************function calls**************************************
function Fr = F(soln)
 Fr = zeros(2*xlen,1);
 v = soln(1:2:2*xlen-1);
 w = soln(2:2:2*xlen);
 fdiag = q*(v.^2+w.^2);
 Fr(1:2:2*xlen-1) = -fdiag.*w;
 Fr(2:2:2*xlen) = fdiag.*v;
end

end