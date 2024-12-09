function homoARZ_Elisa()
%% Homogeneous ARZ Model
% \rho_t+(\rho u)_x=0
%  z_t+(z u)_x=0
%where z=\rho.*(u+\rho.^gamma)=\rho w (in the overleaf)

%% with Lax Friedrichs scheme & Periodic boundary conditions
%clc;clear all;close all

global dx dt

%initial parameters
L=1;                    %road length
T=1;
dx= 0.01;
dt= 1 *dx;               % to respect CFL
Nx=round(L/dx);
Nt=round(T/dt);

gamma=1;                 % parameter for w

%% initialization

rho=zeros(Nx,Nt);
u=zeros(Nx,Nt);

%% initial data

init_bool = -1;
%0: Riemann
%1: Bell

%rho(:,1)=sin(x(:));
%u(:,1)=cos(x(:));

if init_bool == -1
% Small perturbation
rho(1:round(Nx/3),1)=0.55; rho(round(Nx/3)+1:round(2*Nx/3),1)=0.95; rho(round(2*Nx/3)+1:Nx,1)=0.2;
%rho(1:round(Nx/2),1)=0.4;  rho(round(Nx/2)+1:round(Nx/2)+4,1)=0.8; rho(round(Nx/2)+5:end,1)=0.4;
u(1:round(Nx/3),1)=0.4;    u(round(Nx/3)+1:round((2*Nx/3),1))=0.1;   u(round(2*Nx/3)+1:Nx,1)=0.8;
%u(1:round(Nx/3),1)=0.4 ; u(round(Nx/3)+1:round(2*Nx/3),1)=0.1 ; u(round(2*Nx/3)+1:Nx,1)=0.7;
%u(1:round(Nx/2),1)=0.5;    u(round(Nx/2)+1:round(Nx/2)+4,1)=0.1;   u(round(Nx/2)+5:end,1)=0.5;

elseif init_bool == 0
    %Riemann Problem
    rho(1:round(Nx/2),1)=0.7;  rho(round(Nx/2)+1:end,1)=0.2;
    u(1:round(Nx/2),1)=0.2;    u(round(Nx/2)+1:end,1)=0.7;
elseif init_bool == 1 
    rho_jam = 1;
    u_free = 1;

    rmin = 0.55;  rmax = 0.95;
    sigma = 0.1;
    %   user-specified parameters
    x = linspace(-L/2+dx/2, L/2-dx/2, Nx);
    rho(:,1) = rmin + (rmax-rmin) * exp(-x.^2/sigma^2);
    %u(:,1) = 1.1* u_free *(1-rho(:,1)./rho_jam); %ARZ
    u(:,1) = u_free *(1-rho(:,1)./rho_jam); %GSOM
    %u(1:round(Nx/2),1)=0.2;    u(round(Nx/2)+1:end,1)=0.9; %GSOM
    
elseif init_bool == 2
    nw = 2;
    rho(:,1) = (1 + sin(linspace(0, nw*2*pi, Nx))) / 2;
    u(:,1) = (1 + cos(linspace(0, nw*2*pi, Nx))) / 2;
end

%w=u+rho.^gamma; just to remember the relation between the variables
z(:,1) = rho(:,1).*(u(:,1)+rho(:,1).^gamma);

%% Time Evolution

for n = 1:Nt
    for i = 1:Nx
%         l = i-1;
%         r = i+1;
        if i>1	
            l = i - 1;  
        else l = Nx-1; 
        end
        if i < Nx   
            r = i + 1; 
        else r = 1+1; 
        end
        
    rho(i,n+1)=rho(i,n)-dt/dx*(num_flux(rho(r,n),rho(i,n), rho(i,n)*u(i,n), rho(r,n)*u(r,n))...
                              -num_flux(rho(i,n),rho(l,n), rho(l,n)*u(l,n), rho(i,n)*u(i,n)));
    z(i,n+1)=z(i,n)-dt/dx*(num_flux(z(r,n),z(i,n), z(i,n)*u(i,n), z(r,n)*u(r,n))...
                              -num_flux(z(i,n),z(l,n), z(l,n)*u(l,n), z(i,n)*u(i,n)));
    end
    
    %update of the velocity
    u(:,n+1)=z(:,n+1)./rho(:,n+1)-rho(:,n+1).^gamma;
end

% for n=1:Nt
%     %left cell (boundary)
%     rho(1,n+1)=rho(1,n)-dt/dx*(num_flux(rho(2,n),rho(1,n), rho(1,n)*u(1,n), rho(2,n)*u(2,n))...
%                               -num_flux(rho(1,n),rho(Nx-1,n), rho(Nx-1,n)*u(Nx-1,n), rho(1,n)*u(1,n)));
%     z(1,n+1)=z(1,n)-dt/dx*(num_flux(z(2,n),z(1,n), z(1,n)*u(1,n), z(2,n)*u(2,n))...
%                               -num_flux(z(1,n),z(Nx-1,n), z(Nx-1,n)*u(Nx-1,n), z(1,n)*u(1,n))); 
% 
%     %inner cells
%     for i=2:Nx-1
%     rho(i,n+1)=rho(i,n)-dt/dx*(num_flux(rho(i+1,n),rho(i,n), rho(i,n)*u(i,n), rho(i+1,n)*u(i+1,n))...
%                               -num_flux(rho(i,n),rho(i-1,n), rho(i-1,n)*u(i-1,n), rho(i,n)*u(i,n)));
%     z(i,n+1)=z(i,n)-dt/dx*(num_flux(z(i+1,n),z(i,n), z(i,n)*u(i,n), z(i+1,n)*u(i+1,n))...
%                               -num_flux(z(i,n),z(i-1,n), z(i-1,n)*u(i-1,n), z(i,n)*u(i,n)));
%     end
%     
%     %right cell (boundary)
%     rho(Nx,n+1)=rho(Nx,n)-dt/dx*(num_flux(rho(2,n),rho(Nx,n), rho(Nx,n)*u(Nx,n), rho(2,n)*u(2,n))...
%                               -num_flux(rho(Nx,n),rho(Nx-1,n), rho(Nx-1,n)*u(Nx-1,n), rho(Nx,n)*u(Nx,n)));
%     z(Nx,n+1)=z(Nx,n)-dt/dx*(num_flux(z(2,n),z(Nx,n), z(Nx,n)*u(Nx,n), z(2,n)*u(2,n))...
%                               -num_flux(z(Nx,n),z(Nx-1,n), z(Nx-1,n)*u(Nx-1,n), z(Nx,n)*u(Nx,n)));
%                           
% 
%     %update of the velocity
%     u(:,n+1)=z(:,n+1)./rho(:,n+1)-rho(:,n+1).^gamma;
%     
% end

w=u+rho.^gamma;


save rho z u w

%% Plot

figure;
for n=1:10:Nt-1
  subplot(3,1,1)   % density
  plot(rho(:,n),'bo-');
  hold on;  
  axis([1 Nx -0.1 1]);
  ylabel('rho (density)');
  drawnow
  title(['time = ', num2str((n-1)*dt, '%-1.2f')]);
hold off
subplot(3,1,2)  %velocity
  plot(u(:,n),'ro-');
  hold on;  
  axis([1 Nx -0.1 2]);
  ylabel('velocity');
  drawnow
hold off
subplot(3,1,3)   % w
  plot(w(:,n),'ko-');
  hold on;  
  axis([1 Nx -0.1 2]);
  ylabel('Lagrangian marker');
  drawnow
hold off
pause(0.2)
end

% figure;
% plotfig_3d(tra);
x = linspace(dx/2, L-dx/2, Nx);
t = linspace(0, T-dt, Nt);
[T X] = meshgrid(t, x);
    
figure;
    rho = mesh(T, X, rho(:,1:Nt));
    colormap('jet');
    rho.FaceColor='interp';
    xlabel('time');
    ylabel('space');
    zlabel('Traffic density');

end

function [flux] = num_flux(flux_r,flux_l,rho_l,rho_r) %num_flux(rho_r,rho_l,flux_l,flux_r) 
% Numerical flux for Lax Friedrichs schemeW
global dx dt
flux=0.5*(rho_l+rho_r)-0.5*dx/dt*(flux_r-flux_l);
end
