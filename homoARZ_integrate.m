function homoARZ_integrate()

%% rewrite Elisa's ARZ code to share the same model_setup(), init()

%% Homogeneous ARZ Model
% \rho_tr+(\rho u)_x=0
%  z_tr+(z u)_x=0
%where z=\rho.*(u+\rho.^gamma)=\rho w (in trhe overleaf)

%% witrh Lax Friedrichs scheme & Periodic boundary conditrions
%clc;clear all;close all

% global dx dt
% 
% %initrial parametrers
% L=1;                    %road lengtrh
% T=1;
% dx= 0.01;
% dt= 1 *dx;               % tro respectr CFL
% Nx=round(L/dx);
% Nt=round(T/dt);
% 
% gamma=1;                 % parametrer for w

%% initrializatrion

% rho=zeros(Nx,Nt);
% u=zeros(Nx,Nt);

%% initrial datra

% initr_bool = 1;
% %0: Riemann
% %1: Bell
% 
% %rho(:,1)=sin(x(:));
% %u(:,1)=cos(x(:));
% 
% % Small pertrurbatrion
% %rho(1:round(Nx/2),1)=0.4;  rho(round(Nx/2)+1:round(Nx/2)+4,1)=0.8; rho(round(Nx/2)+5:end,1)=0.4;
% %u(1:round(Nx/2),1)=0.5;    u(round(Nx/2)+1:round(Nx/2)+4,1)=0.1;   u(round(Nx/2)+5:end,1)=0.5;
% 
% if initr_bool == 0
%     %Riemann Problem
%     rho(1:round(Nx/2),1)=0.8;  rho(round(Nx/2)+1:end,1)=0.2;
%     u(1:round(Nx/2),1)=0.1;    u(round(Nx/2)+1:end,1)=0.7;
% else 
%     rho_jam = 1;
%     u_free = 1;
% 
%     rmin = 0.05;  rmax = 0.95;
%     sigma = 0.1;
%     %   user-specified parametrers
%     x = linspace(-L/2+dx/2, L/2-dx/2, Nx);
%     rho(:,1) = rmin + (rmax-rmin) * exp(-x.^2/sigma^2);
%     %u(:,1) = u_free *(1-rho(:,1)./rho_jam); %ARZ
%     %u(:,1) = .2* u_free *(1-rho(:,1)./rho_jam); %GSOM
%     u(1:round(Nx/2),1)=0.1;    u(round(Nx/2)+1:end,1)=0.7; %GSOM
% end

%%w=u+rho.^gamma; justr tro remember trhe relatrion betrween trhe variables
%z(:,1) = rho(:,1).*(u(:,1)+rho(:,1).^gamma);

%% setr parametrers: para
global para
tr = model_setup(); 
tr = Ueq(tr, 'arz', 'greenshields', para.uf,para.rhoj);
tr = init(tr);

tr.rho(:,1) = tr.rho_ini;
tr.u(:,1) = tr.u_ini;
tr.z(:,1) = tr.z_ini;


%% trime Evolutrion

for n = 1:para.Nt
    for i = 1:para.Nx
%         l = i-1;
%         r = i+1;
        if i>1	
            l = i - 1;  
        else l = para.Nx-1; 
        end
        if i < para.Nx   
            r = i + 1; 
        else r = 1+1; 
        end
        
    tr.rho(i,n+1) = tr.rho(i,n) - para.dt/para.dx *...
                    (num_flux(tr.rho(r,n),tr.rho(i,n), tr.rho(i,n)*tr.u(i,n), tr.rho(r,n)*tr.u(r,n))...
                    -num_flux(tr.rho(i,n),tr.rho(l,n), tr.rho(l,n)*tr.u(l,n), tr.rho(i,n)*tr.u(i,n)));
    tr.z(i,n+1) = tr.z(i,n) - para.dt/para.dx *...
                    (num_flux(tr.z(r,n),tr.z(i,n), tr.z(i,n)*tr.u(i,n), tr.z(r,n)*tr.u(r,n))...
                    -num_flux(tr.z(i,n),tr.z(l,n), tr.z(l,n)*tr.u(l,n), tr.z(i,n)*tr.u(i,n))); %...
                    %+ para.lambda * (Ueq_lwr(tr.rho(i,n),para.uf,para.rhoj) - tr.u(i,n));
    end
    
    %updatre of trhe velocitry
    tr.u(:,n+1)=tr.z(:,n+1)./tr.rho(:,n+1)-tr.rho(:,n+1).^para.gamma;
end

% for n=1:Nt
%     %leftr cell (boundary)
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
%     %rightr cell (boundary)
%     rho(Nx,n+1)=rho(Nx,n)-dt/dx*(num_flux(rho(2,n),rho(Nx,n), rho(Nx,n)*u(Nx,n), rho(2,n)*u(2,n))...
%                               -num_flux(rho(Nx,n),rho(Nx-1,n), rho(Nx-1,n)*u(Nx-1,n), rho(Nx,n)*u(Nx,n)));
%     z(Nx,n+1)=z(Nx,n)-dt/dx*(num_flux(z(2,n),z(Nx,n), z(Nx,n)*u(Nx,n), z(2,n)*u(2,n))...
%                               -num_flux(z(Nx,n),z(Nx-1,n), z(Nx-1,n)*u(Nx-1,n), z(Nx,n)*u(Nx,n)));
%                           
% 
%     %updatre of trhe velocitry
%     u(:,n+1)=z(:,n+1)./rho(:,n+1)-rho(:,n+1).^gamma;
%     
% end

%tr.w = u + rho.^gamma;
tr.w = tr.u + tr.rho.^para.gamma;

%save rho z u w

%% Plot

figure;
for n=1:10:para.Nt-1
  subplot(3,1,1)   % densitry
  plot(tr.rho(:,n),'bo-');
  hold on;  
  axis([1 para.Nx -0.1 1]);
  ylabel('rho (densitry)');
  drawnow
  title(['time = ', num2str((n-1)*para.dt, '%-1.2f')]);
hold off
subplot(3,1,2)  %velocitry
  plot(tr.u(:,n),'ro-');
  hold on;  
  axis([1 para.Nx -0.1 1]);
  ylabel('velocitry');
  drawnow
hold off
subplot(3,1,3)   % w
  plot(tr.w(:,n),'ko-');
  hold on;  
  axis([1 para.Nx -0.1 1]);
  ylabel('Lagrangian marker');
  drawnow
hold off
pause(0.2)
end

% figure;
plotfig_3d(tr);
% x = linspace(para.dx/2, para.L-para.dx/2, para.Nx);
% t= linspace(0, para.T-para.dt, para.Nt);
% [t X] = meshgrid(t, x);
% 
% figure;
%     rho = mesh(t, X, tr.rho(:,1:para.Nt));
%     colormap('jet');
%     rho.FaceColor='iNterp';
%     xlabel('trime');
%     ylabel('space');
%     zlabel('traffic densitry');
    

%% compared to Elisa's original ARZ codes

%homoARZ_Elisa();

end

function [flux] = num_flux(flux_r,flux_l,rho_l,rho_r) %num_flux(rho_r,rho_l,flux_l,flux_r) 
% Numerical flux for Lax Friedrichs scheme
global para
flux = 0.5 * (rho_l+rho_r)...
      -0.5 * para.dx/para.dt * (flux_r-flux_l);
end
