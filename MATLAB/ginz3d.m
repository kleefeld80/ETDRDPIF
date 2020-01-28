function [runtime] = ginz3d(dt,steps)

% dt: time step
% steps: number of spatial points in each coordinate direction

% diffusion coefficient
epsln1 = 1;

% create nodes
nnodes = steps^3;
x = linspace(0,200,steps+1); y=x; z=x;
xint = x(1:steps);
yint=xint;
zint=xint;

h = abs(x(1)-x(2))
nodes = zeros(nnodes,3);

j = 1;
for k = 1 : steps
        for i = 1:steps
            for ii=1:steps
               nodes(j,:) = [xint(ii) yint(i) zint(k)];
               j = j+1;
            end
        end
end

% discretize time interval
t = 0:dt:100; tlen = length(t)

% Stacking nodes for evolution
%u_old = reshape(0.1*randn(steps,steps,steps),steps*steps*steps,1);
u_old=exp(-((nodes(:,1)-  50).^2+(nodes(:,2)-  50).^2+(nodes(:,3)-  50).^2)/1000)...
     -exp(-((nodes(:,1)- 100).^2+(nodes(:,2)- 100).^2+(nodes(:,3)- 100).^2)/1000);

%% Block matrix Assembly
w = (epsln1*dt)/h^2;
Q = blktridiag(2,-1,-1,steps);
Q(1,end) = -1; Q(end,1) = -1; I = speye(steps);

A1 = w*kron(I,kron(I,Q));
A2 = w*kron(I,kron(Q,I));
A3 = w*kron(Q,kron(I,I));
clear Q

% System matrices
r1 = 1/3; r2 = 1/4;
Id = kron(I,kron(I,I));

A1x = (Id + r1*A1);
A2x = (Id + r2*A1);
A3x = (Id + A1);
clear A1

A1y = (Id + r1*A2);
A2y = (Id + r2*A2);
A3y = (Id + A2);
clear A2

A1z = (Id + r1*A3);
A2z = (Id + r2*A3);
A3z = (Id + A3);
clear A3
clear I
clear Id

[L1,U1] =lu(A3x);
clear A3x
[L2,U2] =lu(A3y);
clear A3y
[L3,U3] =lu(A3z);
clear A3z

 [L11,U11] =lu(A1x); 
 clear A1x
 [L22,U22] =lu(A1y);
 clear A1y
 [L33,U33] =lu(A1z);
 clear A1z
 
 [L111,U111] =lu(A2x); 
 clear A2x
 [L222,U222] =lu(A2y);
 clear A2y
 [L333,U333] =lu(A2z);
 clear A2z

tic
for i = 2:tlen
    
    % Processor 4 (could potentially be split between 2 processors)   
    F_old = F(u_old);
    % For u
    p1 = U1\(L1\(F_old));
    p2 = U2\(L2\(p1));
    p3u =U3\(L3\(p2));
    
    % Processor 3 (could potentially be split between 2 processors)
    % For u
    d1 = U1\(L1\(u_old));
    d2 = U2\(L2\(d1));
    d3u =U3\(L3\(d2));
    u_star = d3u + dt*p3u;
    F_star = F(u_star);
       
    % Processor 2
    % For u
    b1 = U11\(L11\(F_old));
    b2 = U111\(L111\(F_old));
    c2 = 9*b1-8*b2;
    b3 = U22\(L22\(c2));
    b4 = U222\(L222\(c2));
    c4u = 9*b3-8*b4;
   
    
    % Processor 1
    % For u
    a1 = U11\(L11\(u_old));
    a2 = U111\(L111\(u_old));
    c1 = 9*a1-8*a2;
    a3 = U22\(L22\(c1));
    a4 = U222\(L222\(c1));
    c3u = 9*a3-8*a4;
    s1u = U33\(L33\(9*c3u+2*dt*c4u+dt*F_star));
    s2u = U333\(L333\(8*c3u+(3/2)*dt*c4u+0.5*dt*F_star));
    u_old = s1u-s2u;
   
    
     %waitbar(i/tlen,hw)
end
%u_old
runtime = toc;
    v=real(reshape(u_old,steps,steps,steps));
    figure(1)
    [x,y,z]=meshgrid(xint,yint,zint);
    xmin = min(x(:)); 
    ymin = min(y(:)); 
    zmin = min(z(:));

    xmax = max(x(:)); 
    ymax = max(y(:)); 
    zmax = max(z(:));

    hold on
    hx = slice(x,y,z,v,xmin,[],[]);
    hx.FaceColor = 'interp';
    hx.EdgeColor = 'none';

    hy = slice(x,y,z,v,[],ymin,[]);
    hy.FaceColor = 'interp';
    hy.EdgeColor = 'none';

    hz = slice(x,y,z,v,[],[],zmax);
    hz.FaceColor = 'interp';
    hz.EdgeColor = 'none';
    daspect([1,1,1])
    axis tight
    view(-38.5,16)
    set(gca,'TickLength',[0 0]);
    set(gca,'LineWidth', 1);
    set(gca,'FontSize',10);
    set(gca,'FontWeight','bold');
    ax=gca;
    ax.BoxStyle = 'full';
    ax.CLim=[-1,1];
    xlabel('x')
    ylabel('y')
    zlabel('z')
    set(gca,'XTick',[0 50 100 150 200]);
    set(gca,'YTick',[0 50 100 150 200]);
    set(gca,'ZTick',[0 50 100 150 200]);
    colorbar
    box on
    hold off
    print -dpng ginzburg3real.png
    
    v=imag(reshape(u_old,steps,steps,steps));
    figure(2)
    [x,y,z]=meshgrid(xint,yint,zint);
    xmin = min(x(:)); 
    ymin = min(y(:)); 
    zmin = min(z(:));

    xmax = max(x(:)); 
    ymax = max(y(:)); 
    zmax = max(z(:));

    hold on
    hx = slice(x,y,z,v,xmin,[],[]);
    hx.FaceColor = 'interp';
    hx.EdgeColor = 'none';

    hy = slice(x,y,z,v,[],ymin,[]);
    hy.FaceColor = 'interp';
    hy.EdgeColor = 'none';

    hz = slice(x,y,z,v,[],[],zmax);
    hz.FaceColor = 'interp';
    hz.EdgeColor = 'none';
    daspect([1,1,1])
    axis tight
    view(-38.5,16)
    set(gca,'TickLength',[0 0]);
    set(gca,'LineWidth', 1);
    set(gca,'FontSize',10);
    set(gca,'FontWeight','bold');
    ax=gca;
    ax.BoxStyle = 'full';
    ax.CLim=[-1,1];
    xlabel('x')
    ylabel('y')
    zlabel('z')
    set(gca,'XTick',[0 50 100 150 200]);
    set(gca,'YTick',[0 50 100 150 200]);
    set(gca,'ZTick',[0 50 100 150 200]);
    colorbar
    box on
    hold off
    print -dpng ginzburg3imag.png

% %****************function calls**************************************
function Fr = F(u)
 Fr=u-(1+1i*1.3).*u.*abs(u).*abs(u);
end

end
