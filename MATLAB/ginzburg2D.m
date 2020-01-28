% Matlab code template
% 2D {3D} complex Gizburg-Landau equation on a periodic domain.
% u_t = u -(1+iA)u*abs(u)^2 + D(u_xx + u_yy)
% For a 3D code, replace relevant sections in code with code in { }.
% For a different equation, change the diffusion paramter, D, and
% nonlinear function, g. In Nv, Na, Nb, Nc, g must be called properly.
%=========================  CUSTOM SET UP =========================
A=1.3; lam=200; N=400; D=1; h=1/20; %domain/grid size, diff par, timestep
g=inline('u-(1+1i*A)*u.*(abs(u).^2)','u','A'); %nonlinear function
%========================= GENERIC SET UP =========================
x=(lam/N)*(1:N)'; [X,Y]=ndgrid(x,x);      %{[X,Y,Z]=ndgrid(x,x,x)}
u=randn(N,N)*0.1; v=fftn(u);              %random IC.{randn(N,N,N)}
k=[0:N/2-1 0 -N/2+1:-1]'/(lam/(2*pi)); %wave numbers
[xi,eta]=ndgrid(k,k); %2D wave numbers. {[xi,eta,zeta]=ndgrid(k,k,k)}
L=-D*(eta.^2+xi.^2);  %2D Laplacian.    {-D*(eta.^2+xi.^2+zeta.^2)}
Fr=logical(zeros(N,1)); %High frequencies for de-aliasing
Fr(N/2+1-round(N/6):N/2+round(N/6))=1;
[alxi,aleta]=ndgrid(Fr,Fr); %{[alxi,aleta,alzeta]=ndgrid(Fr,Fr,Fr)}
ind=alxi | aleta; %de-aliasing index.{alxi | aleta | alzeta}
%=============== PRECOMPUTING ETDRK4 COEFFS =====================
E=exp(h*L); E2=exp(h*L/2);
M=16; % no. of points for complex mean
r=exp(1i*pi*((1:M)-0.5)/M); % roots of unity
L=L(:); LR=h*L(:,ones(M,1))+r(ones(N^2,1),:); %{r(ones(N^3,1),:)}
Q=h*real(mean( (exp(LR/2)-1)./LR,2));
f1=h*real(mean( (-4-LR+exp(LR).*(4-3*LR+LR.^2))./LR.^3 ,2));
f2=h*real(mean( (4+2*LR+exp(LR).*(-4+2*LR))./LR.^3 ,2));
f3=h*real(mean( (-4-3*LR-LR.^2+exp(LR).*(4-LR))./LR.^3 ,2));
f1=reshape(f1,N,N); f2=reshape(f2,N,N); f3=reshape(f3,N,N);
L=reshape(L,N,N); Q=reshape(Q,N,N); clear LR %{reshape(*,N,N,N)}
%==================== TIME STEPPING LOOP =======================
tmax=100; nmax=round(tmax/h);
for n = 1:nmax
    t=n*h; %**Nonlinear terms are evaluated in physical space**
    Nv=fftn( g(ifftn(v),A) ); %Nonlinear evaluation. g(u,*)
    a=E2.*v + Q.*Nv; %Coefficient ?a? in ETDRK formula
    Na=fftn( g(ifftn(a),A) ); %Nonlinear evaluation. g(a,*)
    b=E2.*v + Q.*Na; %Coefficient ?b? in ETDRK formula
    Nb=fftn( g(ifftn(b),A) ); %Nonlinear evaluation. g(b,*)
    c=E2.*a + Q.*(2*Nb-Nv); %Coefficient ?c? in ETDRK formula
    Nc=fftn( g(ifftn(c),A) ); %Nonlinear evaluation. g(c,*)
    v=E.*v + Nv.*f1 + (Na+Nb).*f2 + Nc.*f3;%update
    v(ind) = 0; % High frequency removal --- de-aliasing
end
u=real(ifftn(v));
%figure(1)
%imagesc(u)