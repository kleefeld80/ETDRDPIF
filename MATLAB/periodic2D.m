function [runtime] = periodic2D(dt,steps)

% dt: time step
% steps: number of spatial points in each coordinate direction

% diffusion coefficient
epsln1 = 1;

% create nodes
nnodes = steps^2;
x = linspace(0,200,steps+1); y=x;
xint = x(1:steps);
yint=xint;

h = abs(x(1)-x(2));
fprintf('dt=%f h=%f\n',dt,h)
nodes = zeros(nnodes,2);


j = 1;
for k = 1 : steps
        for i = 1:steps
               nodes(j,:) = [xint(i) yint(k)];
            j = j+1;
        end
end

% discretize time interval
t = 0:dt:100; tlen = length(t);

% Stacking nodes for evolution
w_old = reshape(0.1*randn(steps,steps),steps*steps,1); 

%% Block matrix Assembly
w = (epsln1*dt)/h^2;
Q = blktridiag(2,-1,-1,steps);
Q(1,end) = -1; Q(end,1) = -1; I = speye(steps);
A1 = w*kron(I,Q);A2 = w*kron(Q,I); 
I = speye(nnodes);
M1 = sparse(I+A1); M2 = sparse(I+(1/3)*A1); M3 = sparse(I+(1/4)*A1);
M11 = sparse(I+A2); M22 = sparse(I+(1/3)*A2); M33 = sparse(I+(1/4)*A2);

%% Time Evolution 

[L1,U1] =lu(M1);
[L2,U2] =lu(M2);
[L3,U3] =lu(M3);
 [L11,U11] =lu(M11); 
 [L22,U22] =lu(M22);
 [L33,U33] =lu(M33);

tic
for i = 2:tlen
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
end
  U = reshape(w_old,steps,steps);
  right=U(:,end);
  U=[right U];
  upper=U(1,:);
  U=[U;upper];
runtime = toc;

 %% Plots
figure(19)
imagesc(x,y,real(U),[-1,1])
set(gca, 'YDir', 'normal');
set(gca,'TickLength',[0 0]);
set(gca,'LineWidth', 1);
set(gca,'FontSize',10);
pbaspect(gca,[1 1 1])
set(gca,'XTick',[0 20 40 60 80 100 120 140 160 180 200]);
set(gca,'YTick',[0 20 40 60 80 100 120 140 160 180 200]);
set(gca,'FontWeight','bold');
colorbar
print -deps2c ginzburgreal.eps

figure(20)
imagesc(x,y,imag(U),[-1,1])
set(gca, 'YDir', 'normal');
set(gca,'TickLength',[0 0]);
set(gca,'LineWidth', 1);
set(gca,'FontSize',10);
pbaspect(gca,[1 1 1])
set(gca,'XTick',[0 20 40 60 80 100 120 140 160 180 200]);
set(gca,'YTick',[0 20 40 60 80 100 120 140 160 180 200]);
set(gca,'FontWeight','bold');
colorbar
print -deps2c ginzburgimag.eps

% %****************function calls**************************************
function Fr = F(u)
 Fr=u-(1+1i*1.3).*u.*abs(u).*abs(u);
end

end
