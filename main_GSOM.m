function main_GSOM()
clc;clear all;close all;

% figure 2: u_init and w_init should be from 'init.m'.
% figure 6: u_init has the same issue
% figure 11: same issue for u_init
% figure 12: w should within [0,1].




%% set parameters: para
global para

% define rho, u, z, w
tra = model_setup_w(); 
%   set Ueq (LWR or ARZ) & FD type
tra = Ueq(tra, 'arz', 'greenshields', para.uf,para.rhoj);


%  cost function
%tra.vlf = @(rho, w, Vx) arz(tra, rho, w, Vx);
tra.vlf = @(rho, w, Vx, Vw) gsom(tra, rho, w, Vx, Vw);
%tra.vlf_arz=@

%   initial and terminal condition
tra = init(tra);

%% solve FP
tic;  
maxiter = 200;
[tra err u_hist V_hist] = solve(tra, maxiter, 1e-2);
toc;  
save tra;
%save equilibrium
rhoe_=tra.rho;
ue_=tra.u;
Ve_=tra.V;
uwe_=tra.uw;
%writematrix(tra.rho, 'zm_bellshape_arz_GSOM_rho.csv')
%writematrix(tra.u, 'zm_bellshape_arz_GSOM_u.csv')
rhoe_w = repmat(rhoe_, [1, 20, 1]);
rhoe_w = reshape(rhoe_w, [20,20,21]);


%
%% surface
x = linspace(0, para.L, para.Nx);
t = linspace(0, para.T, para.Nt);
[T X] = meshgrid(t, x);

uw_0 = squeeze(tra.uw(:,1,:));
save uw_0;
uw_50 = squeeze(tra.uw(:,end/2, :));
save uw_50;
uw_100 = squeeze(tra.uw(:,end, :));
save uw_100;

figure;
u = mesh(T, X, uw_0(:,1:para.Nt));
colormap('jet');
u.FaceColor='interp';
% xlabel('time');
% ylabel('space');
zlabel('u');
pause(0.2);

figure;
u = mesh(T, X, uw_50(:,1:para.Nt));
colormap('jet');
u.FaceColor='interp';
% xlabel('time');
% ylabel('space');
zlabel('u');
pause(0.2);

figure;
u = mesh(T, X, uw_100(:,1:para.Nt));
colormap('jet');
u.FaceColor='interp';
% xlabel('time');
% ylabel('space');
zlabel('u');
pause(0.2);

%% surface
x = linspace(0, para.L, para.Nx);
t = linspace(0, para.T, para.Nt);
[T X] = meshgrid(t, x);

V_0 = squeeze(tra.V(:,1,:));
V_50 = squeeze(tra.V(:,end/2, :));
V_100 = squeeze(tra.V(:,end, :));

figure;
u = mesh(T, X, V_0(:,1:para.Nt));
colormap('jet');
u.FaceColor='interp';
% xlabel('time');
% ylabel('space');
zlabel('V');
pause(0.2);

figure;
u = mesh(T, X, V_50(:,1:para.Nt));
colormap('jet');
u.FaceColor='interp';
% xlabel('time');
% ylabel('space');
zlabel('V');
pause(0.2);

figure;
u = mesh(T, X, V_100(:,1:para.Nt));
colormap('jet');
u.FaceColor='interp';
% xlabel('time');
% ylabel('space');
zlabel('V');
pause(0.2);




%% Plot
figure;nfig=3;
for n=1:2:para.Nt-1
  subplot(nfig,1,1)   % density  
  plot(tra.rho(:,n),'bo-');
  hold on;  
  axis([1 para.Nx -0.1 1]);
  ylabel('\rho'); %(Density)');
  drawnow
  title(['time = ', num2str((n-1)*para.dt, '%-1.2f')]); 
hold off
subplot(nfig,1,2)  %velocity
  plot(tra.u(:,n),'ro-');
  hold on;  
  axis([1 para.Nx -0.1 1]);
  ylabel('u'); % (Velocity)');
  drawnow
hold off
subplot(nfig,1,3)   
  plot(tra.w(:,n),'ko-');
  hold on;  
  axis([1 para.Nx -0.1 1]);
  ylabel('w'); % (Lagrangian marker)');
  drawnow
hold off
% subplot(nfig,1,4)   % V
%   plot(tra.V(:,n),'g.-');
%   hold on;  
%   axis([1 para.Nx -0.1 1]);
%   ylabel('V'); % (Value)');
%   drawnow
%   hold off
% subplot(5,1,5)   % z
%   plot(tra.z(:,n),'c*');
%   hold on;  
%   axis([1 para.Nx -0.1 2]);
%   ylabel('z (Traffic flux)');
%   drawnow
% hold off
%pause(0.2)
end
pause(0.2);

nfig=3;
nvec = [1 ceil(para.Nt/2) para.Nt];
for n=1:size(nvec,2)
  figure;
  subplot(nfig,1,1)   % density  
  plot(tra.rho(:,nvec(n)),'bo-');
  %hold on;  
  axis([1 para.Nx -0.1 1]);
  ylabel('\rho'); %(Density)');
  title(['time = ', num2str((nvec(n)-1)*para.dt, '%-1.2f')]); 
subplot(nfig,1,2)  %velocity
  plot(tra.u(:,nvec(n)),'ro-');
  %hold on;  
  axis([1 para.Nx -0.1 1]);
  ylabel('u'); % (Velocity)');
subplot(nfig,1,3)   
  plot(tra.w(:,nvec(n)),'ko-');
  %hold on;  
  axis([1 para.Nx -0.1 1]);
  ylabel('w'); % (Lagrangian marker)');
  pause(0.2);
end
pause(0.2);
%figure;plotfig_3d(tra);
% x = linspace(0, para.L, para.Nx);
%     t = linspace(0, para.T, para.Nt);
%     [T X] = meshgrid(t, x);
%     
%     figure
%     rho = mesh(T, X, tra.rho(:,1:para.Nt));
%     colormap('jet');
%     rho.FaceColor='interp';
%     xlabel('time');
%     ylabel('space');
%     zlabel('Density');

% plotfig(tra);
plotfig_3d(tra);
pause(0.2);
% x = linspace(para.dx/2, para.L-para.dx/2, para.Nx);
% t= linspace(0, para.T-para.dt, para.Nt);
% [t X] = meshgrid(t, x);
% figure;
%     rho = mesh(t, X, tra.rho(:,1:para.Nt));
%     colormap('jet');
%     rho.FaceColor='iNterp';
%     xlabel('trime');
%     ylabel('space');
%     zlabel('traffic densitry');
% 
% figure;
%     u = mesh(t, X, tra.u(:,1:para.Nt));
%     colormap('jet');
%     u.FaceColor='iNterp';
%     xlabel('trime');
%     ylabel('space');
%     zlabel('traffic velocity');
% 
% figure;
%     w = mesh(t, X, tra.w(:,1:para.Nt));
%     colormap('jet');
%     w.FaceColor='iNterp';
%     xlabel('trime');
%     ylabel('space');
%     zlabel('traffic marker');


%save tra 
figure;%nfig=1;
  %subplot(nfig,1,1)   % err
  plot(err,'b*-');
  hold on;  
  %axis([1 maxiter]);
  xlabel('iteration');
  ylabel('Gap'); 

  % subplot(nfig,1,2)   % err
  % plot(err(:,1:maxiter),'b*-');
  % hold on;  
  % axis([1 maxiter]);
  % xlabel('iteration');
  % ylabel('error'); 

%% ARZ Elisa's code to compare
pause(0.2);
homoARZ_integrate_w();

%homoARZ_Elisa();
%}
end