% VP 1d

clear all
a = []; % coefficient
prob = 1; % problem


T = 40;

Nx = 64;
Nv = 128;

% for truncation
opts.max_rank = min(Nx,Nv);
opts.rel_eps = 1e-6;

problem = [];
if prob==1
    problem = 'weak1d';
    a = 0.01;
elseif prob ==2
    problem = 'strong1d';
    a = 0.5;
elseif prob == 3
    problem = 'twostream1d';
elseif prob == 4
    problem = 'twostream1dii';
    a = 0.01;
end

isfilter = 0;
alpha = 35; %
etac = 0.2;
p = 2;





kx = 0.5; % wavenumber in x
vmax = 6;
if prob ==3
    kx = 0.2;
    vmax = 9;
end
Lx = 2*pi/kx;
Lv = vmax * 2; % [-6,6]

kv = 2*pi/Lv;

hx = Lx/Nx;
hv = Lv/Nv;

dt = min(hx,hv)/20;



Nt = T/dt;

% differentiation matrix
% column = [0 kx*5*(-1).^(1:Nx-1).*cot(kx*(1:Nx-1)*hx/2)];
% Dx = toeplitz(column,column([1 Nx:-1:2]));
% 
% column = [0 kv*.5*(-1).^(1:Nv-1).*cot(kv*(1:Nv-1)*hv/2)];
% Dv = toeplitz(column,column([1 Nv:-1:2]));

% for computing E via fft
Nk = [0:Nx/2-1 0 -Nx/2+1:-1]';
Nk(1) = 1;
Nk(Nx/2+1) = 1;
Nk = Nk*1i*kx;


Nkx = [0:Nx/2-1 0 -Nx/2+1:-1]'; % for first derivative
Nnx = abs(Nkx)*2/Nx;


etax = zeros(Nx,1);

ii  = (1:Nx)';

i1 = find(abs(Nnx)<= etac);

i2 = setdiff(ii,i1);

etax(i1) = 1.;

etax(i2) = exp(-alpha*((Nnx(i2)-etac)/(1-etac)).^p);


Nkv = [0:Nv/2-1 0 -Nv/2+1:-1]'; % for first derivative
Nnv = abs(Nkv)*2/Nv;


etav = zeros(Nv,1);

ii  = (1:Nv)';

i1 = find(abs(Nnv)<= etac);

i2 = setdiff(ii,i1);

etav(i1) = 1.;

etav(i2) = exp(-alpha*((Nnv(i2)-etac)/(1-etac)).^p);

if isfilter
    DNx = kx*1i*[0:Nx/2-1 0 -Nx/2+1:-1]'.*etax; % for first derivative
    DNv = kv*1i*[0:Nv/2-1 0 -Nv/2+1:-1]'.*etav; % for first derivative
else
    DNx = kx*1i*[0:Nx/2-1 0 -Nx/2+1:-1]'; % for first derivative
    DNv = kv*1i*[0:Nv/2-1 0 -Nv/2+1:-1]'; % for first derivative
end

DNx2 = DNx.^2;
DNv2 = DNv.^2;



x = hx*(1:Nx)';
v = -vmax + hv*(1:Nv)';


% for computing rho via contraction

htrho =  htensor({ones(Nv,1)});

htj =  htensor({v});

hte =  htensor({v.^2});

% for advection in x

Av = spdiags(v,0,Nv,Nv); % v in htd form


%landau damping
%landau damping
if prob == 1 || prob==2
    f = htensor({1/sqrt(2*pi)*(1+ a*cos(kx*x)), exp(-v.^2/2)}); % create a htd tensor from CP format
elseif prob ==3
    a1 = 1.; b1 = 2.4; eps1 = 1e-3; f = htensor({1/(2*sqrt(2*pi))*(1+ eps1*cos(kx*x)), exp(-((v-b1)/(a1*sqrt(2))).^2) + exp(-((v+b1)/(a1*sqrt(2))).^2)});
elseif prob ==4
    f = htensor({2/(7*sqrt(2*pi))*(1+ a*((cos(2*kx*x) + cos(3*kx*x))/1.2 + cos(kx*x))), (1+5*v.^2).*exp(-v.^2/2)}); % create a htd tensor from CP format
end

if ~f.is_orthog()
    f = orthog(f);
end

tmass = 1;
if prob == 4
    tmass = 12/7;
end
    


figure;
he = animatedline;

e_his = [];
e_rank = [];
mass_his = [];
ener_his = [];



tic;
% forward euler for first step

ft = f;



rho = ttt(f, htrho, 2, 1);

Kinetic = ttt(f, hte, 2, 1);

ener_int = sum(full(Kinetic))*hv*hx;



rhof = hv*full(rho) - tmass;



%Jf = hv*full(J); % may need to subtract the initial momentum;

%     plot(x, Jf);
%
%     return;

r_hat = fft(rhof);

h_hat = r_hat./Nk;
h_hat(1) = 0;
h_hat(Nx/2+1) = 0;

E = real(ifft(h_hat));

% plot(x,E)
% return

ener_int = ener_int + dot(E,E)*hx;




fx = ttm( ttm(f, @(y)fft(y), 1), @(y) real(ifft(DNx.*y)), 1);
fv = ttm( ttm(f, @(y)fft(y), 2), @(y) real(ifft(DNv.*y)), 2);


%     fx2 = ttm( ttm(f, @(y)fft(y), 1), @(y) real(ifft(DNx2.*y)), 1);
%     fv2 = ttm( ttm(f, @(y)fft(y), 2), @(y) real(ifft(DNv2.*y)), 2);
%
%     fxv = ttm( ttm(fx, @(y)fft(y), 2), @(y) real(ifft(DNv.*y)), 2);


% f1

f1x1 = -dt * ttm(fx, @(y) v.*y, 2); % v fx
f1v1 = -dt * ttm(fv, @(y) E.*y, 1); % E fv + J  fv
%   f1x2 = 0.5 *dt^2 * ttm(fx, @(y) E.*y, 1); % E fx
%    f1v1 = ttm(fv, @(y) (-dt*E + 0.5*dt^2*Jf).*y, 1); % E fv + J  fv
%     f1v2 = 0.5 *dt^2 * ttm(fv, {@(y) (rhof).*y, @(y) v.*y}, [1 2]); % v Ex fv

% f2

%     f2x1 =  0.5 *dt^2 * ttm(fx2, @(y) v.^2.*y, 2);  % v^2 fxx
%     f2xv1 = dt^2 * ttm(fxv, {@(y) E.*y, @(y) v.*y}, [1 2]); % E v fxv
%     f2v1 =  0.5 *dt^2 * ttm(fv2, @(y) E.^2.*y, 1);  % E^2 fvv


fsum = {f, f1x1, f1v1};



f = htensor.truncate_sum(fsum, opts);


for i=1:Nt
    
    
    
    rho = ttt(f, htrho, 2, 1);
    
    
    rhof = hv*full(rho) - tmass;
    
    mass = sum(rhof)*hx;
    
    
    Kinetic = ttt(f, hte, 2, 1);

    k_ener = sum(full(Kinetic))*hv*hx; % kinetic energy
    
    
    r_hat = fft(rhof);
    
    h_hat = r_hat./Nk;
    h_hat(1) = 0;
    h_hat(Nx/2+1) = 0;
    
    E = real(ifft(h_hat));
    
    
    
    
    fx = ttm( ttm(f, @(y)fft(y), 1), @(y) real(ifft(DNx.*y)), 1);
    fv = ttm( ttm(f, @(y)fft(y), 2), @(y) real(ifft(DNv.*y)), 2);
    
    

    
    
    % f1
    
    f1x1 = -(2*dt) * ttm(fx, @(y) v.*y, 2); % v fx
    f1v1 = -(2*dt) * ttm(fv, @(y) E.*y, 1); % E fv + J  fv

    

    fsum = {ft, f1x1, f1v1};
   
    ft = f;
    
    f = htensor.truncate_sum(fsum, opts);
    
    
    
    e_ener = dot(E,E)*hx; % electric energy
    
    ener = (k_ener + e_ener - ener_int)/ener_int;
    
    
    e_his = [e_his; e_ener];
    
    mass_his = [mass_his; mass];
    
    ener_his = [ener_his; ener];
    
    rk = rank(f);
    
    e_rank = [e_rank; rk(2)];
    
%     addpoints(he, i*dt, log(dot(E,E)*hx));
%     drawnow
    
    
    disp([i*dt, rk]);
    

    

    
    
end

toc;

%disp(toc);

[xx, vv] = meshgrid(x,v);

figure;

colormap(jet);



if prob==2
    cont = 0:0.01:0.5;
    zlim = [0, 0.5];
elseif prob==4
    cont = 0:0.01:0.5;
    zlim = [0, 0.5];
end
[C, H] = contourf(xx,vv, max(0,full(f)'),cont);

set(H,'LineColor','none');

colorbar;
caxis(zlim);

xlabel('x');
ylabel('v');
title(['T=',num2str(T)]);
set(gcf,'renderer','zbuffer');

filename = strcat(problem,'_fourier_t',num2str(T),'_',num2str(Nx),'_',num2str(Nv),'_contour_',num2str(opts.max_rank),'.eps');
saveas(gcf, filename, 'epsc2');
%export_fig(filename, '-eps');
%print(gcf,'-depsc','-painters',strcat(problem,'_fourier_t',num2str(T),'_',num2str(Nx),'_',num2str(Nv),'_contour.eps'));
%epsclean(strcat(problem,'_fourier_t',num2str(T),'_',num2str(Nx),'_',num2str(Nv),'_contour.eps'));
figure;

semilogy([1:Nt]*dt, e_his);




figure;

plot((1:Nt)'*dt, e_rank);
name1 = strcat(problem,'_','rank','_',num2str(Nx),'_',num2str(Nv),'_',num2str(opts.rel_eps),'_',num2str(opts.max_rank),'.mat');
e_rank = [(1:Nt)'*dt, e_rank];
save(name1, 'e_rank');

name2 = strcat(problem,'_','elec','_',num2str(Nx),'_',num2str(Nv),'_',num2str(opts.rel_eps),'_',num2str(opts.max_rank),'.mat');
e_elec = [(1:Nt)'*dt, e_his];
save(name2, 'e_elec');

name3 = strcat(problem,'_','ener','_',num2str(Nx),'_',num2str(Nv),'_',num2str(opts.rel_eps),'_',num2str(opts.max_rank),'.mat');
e_elec = [(1:Nt)'*dt, ener_his];
save(name3, 'e_elec');

name4 = strcat(problem,'_','mass','_',num2str(Nx),'_',num2str(Nv),'_',num2str(opts.rel_eps),'_',num2str(opts.max_rank),'.mat');
e_elec = [(1:Nt)'*dt, mass_his];
save(name4, 'e_elec');

disp(name1);

