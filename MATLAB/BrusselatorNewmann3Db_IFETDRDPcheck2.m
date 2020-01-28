function [runtime,u_soln] = BrusselatorNewmann3Db_IFETDRDPcheck2(dt,steps)

% dt: time step. Default is 0.001
% steps: number of spatial points in each coordinate direction. Default is 11

%% Model Paramters and initial conditions
a = 3; b = 1.0;
% diffusion coefficient
mu = 0.02; 

% create nodes
x = linspace(0,1,steps); h = abs(x(1)-x(2)); 
y = x;
z = y;
nnodes = steps^3;
nodes = zeros(nnodes,3);
m = 1;
for k = 1:steps
    for j = 1 : steps
            for i = 1:steps
                   nodes(m,:) = [x(i) y(j) z(k)];
                m = m+1;
            end
     end
end

% discretize time interval
t = 0:dt:40; 
tlen = length(t);

% initial condition for u
u_old = 1 + sin(2*pi*nodes(:,1)).*sin(2*pi*nodes(:,2)).*sin(2*pi*nodes(:,3)); 

% initial condition for v
v_old = 3*ones(nnodes,1); 

%% Block matrix Assembly
% 1D  matrix
I = speye(steps);
e = ones(steps,1);r=1/h^2;
B = spdiags([-r*e 2*r*e -r*e], -1:1, steps, steps);
B(1,2) = -2*r;
B(steps,steps-1) = -2*r;

% 3D  matrix for space
Ax = kron(I,kron(I,B));
Ay = kron(I,kron(B,I));
Az = kron(B,kron(I,I));

% System matrices
r1 = mu/3; r2 = mu/4; r3=mu;
Id = kron(I,kron(I,I));
A1x = (Id + r1*dt*Ax);
A2x = (Id + r2*dt*Ax);
A3x = (Id + r3*dt*Ax);
A1y = (Id + r1*dt*Ay);
A2y = (Id + r2*dt*Ay);
A3y = (Id + r3*dt*Ay);
A1z = (Id + r1*dt*Az);
A2z = (Id + r2*dt*Az);
A3z = (Id + r3*dt*Az);

clear Ax Ay Az I Id B

[L3x,U3x]=lu(A3x);
[L3y,U3y]=lu(A3y);
[L3z,U3z]=lu(A3z);

[L2x,U2x]=lu(A2x);
[L2y,U2y]=lu(A2y);
[L2z,U2z]=lu(A2z);

[L1x,U1x]=lu(A1x);
[L1y,U1y]=lu(A1y);
[L1z,U1z]=lu(A1z);

%%
 Uplot = zeros(tlen,1);
 Vplot = Uplot;
 Usoln = reshape(u_old,steps,steps,steps); 
 Vsoln = reshape(v_old,steps,steps,steps);
 Uplot(1,:) = Usoln(4,4,4);
 Vplot(1,:) = Vsoln(4,4,4);
 
tic
for i = 2:tlen
    F_old = F(u_old,v_old);
    % For u
    p1 = U3x\(L3x\F_old(:,1));
    p2 = U3y\(L3y\p1);
    p3u = U3z\(L3z\p2);
    % For v
    p1 = U3x\(L3x\F_old(:,2));
    p2 = U3y\(L3y\p1);
    p3v = U3z\(L3z\p2);
    
    % For u
    d1 = U3x\(L3x\u_old);
    d2 = U3y\(L3y\d1);
    d3u = U3z\(L3z\d2);
    u_star = d3u + dt*p3u;
    % For v
    d1 = U3x\(L3x\v_old);
    d2 = U3y\(L3y\d1);
    d3v = U3z\(L3z\d2);
    v_star = d3v + dt*p3v;
    F_star = F(u_star,v_star);
       
    % For u
    b1 = U1x\(L1x\F_old(:,1));
    b2 = U2x\(L2x\F_old(:,1));
    c2 = 9*b1-8*b2;
    b3 = U1y\(L1y\c2);
    b4 = U2y\(L2y\c2);
    c4u = 9*b3-8*b4;
    % For v
    b1 = U1x\(L1x\F_old(:,2));
    b2 = U2x\(L2x\F_old(:,2));
    c2 = 9*b1-8*b2;
    b3 = U1y\(L1y\c2);
    b4 = U2y\(L2y\c2);
    c4v = 9*b3-8*b4;
    
    % For u
    a1 = U1x\(L1x\u_old);
    a2 = U2x\(L2x\u_old);
    c1 = 9*a1-8*a2;
    a3 = U1y\(L1y\c1);
    a4 = U2y\(L2y\c1);
    c3u = 9*a3-8*a4;
    s1u = U1z\(L1z\(9*c3u+2*dt*c4u+dt*F_star(:,1)));
    s2u = U2z\(L2z\(8*c3u+(3/2)*dt*c4u+0.5*dt*F_star(:,1)));
    u_old = s1u-s2u;
    % For v
    a1 = U1x\(L1x\v_old);
    a2 = U2x\(L2x\v_old);
    c1 = 9*a1-8*a2;
    a3 = U1y\(L1y\c1);
    a4 = U2y\(L2y\c1);
    c3v = 9*a3-8*a4;
    s1v = U1z\(L1z\(9*c3v+2*dt*c4v+dt*F_star(:,2)));
    s2v = U2z\(L2z\(8*c3v+(3/2)*dt*c4v+0.5*dt*F_star(:,2)));
    v_old = s1v-s2v;
    
    % postprocessing for plot
    Usoln = reshape(u_old,steps,steps,steps); 
    Vsoln = reshape(v_old,steps,steps,steps);
    Uplot(i,:) = Usoln(4,4,4);
    Vplot(i,:) = Vsoln(4,4,4);
    
end
  runtime = toc;
  u_soln = u_old;

    figure(18)
    hold on
    plot(t,Uplot,'LineWidth',1)
    plot(t,Vplot,'LineWidth',1)
    xlabel('t')
    ylabel('concentration')
    grid on
    box on
    set(gca,'LineWidth', 1);
    set(gca,'FontSize',10);
    set(gca,'FontWeight','bold');
    legend('u_1','u_2')
    hold off
    axis([0,40,0,5])
    print -depsc2 profiles2.eps
    
function Fr = F(u,v)
 f1 = b+u.^2.*v -(a+1)*u;
 f2 = a*u-u.^2.*v;
 Fr = [f1 f2];
end

end