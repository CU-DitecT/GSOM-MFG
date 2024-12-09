%  Set up a traffic flow model:  ring road, single class
function traffic = model_setup_w() % function traffic = model_setup(L, T, uf, rhoj, Nx, Nt)
global para
%% define parameters
 %   time of simulation
    para.T = 1;
    %   length of the road
    para.L = 1; %* para.T;
%     para.L = 1;
    %   free flow speed
    para.uf = 1;
    %   jam density
    para.rhoj = 1;  
    %   parameter for w
    para.gamma = 1; 
    
% length of the Lagrangian marker
	para.W = para.uf;
    

    %   coefficient for Nx Nt   
    para.multiplier = 100;
    
    ind_para = 2; 
if ind_para == 0
    %   spatial and time mesh size
%     para.Nx = 2 * para.multiplier;
%     para.Nt = 2 * para.multiplier; 
    para.Nx = 1 * para.multiplier;
    para.Nt = 1 * para.multiplier;

    para.dx = para.L / para.Nx;
    para.dt = para.T / para.Nt;
    para.c  = para.dt / para.dx;   %CFL number
elseif ind_para == 1      
    %c = tra.dt / tra.dx;      
    para.c  = 1;   %CFL number c<=1: to respect CFL 
    para.dx = 0.005;
    para.dt = para.c * para.dx;               
    para.Nx = round(para.L/para.dx);
    para.Nt = round(para.T/para.dt);  
    
elseif ind_para == 2 
    para.c  = 1;   %CFL number c<=1: to respect CFL 
    para.Nx = 20;
    para.Nt = 20* (para.T / para.L);    
    para.dx = para.L/para.Nx; 
    para.dt = para.T/para.Nt;               
    
end
if para.c * para.uf > 1  
	disp('CFL condition not satisfied');
	disp('Solver stopped');
	return
end

para.dw = para.dx;
para.Nw = para.W/para.dw;

%   parameter for inhomo ARZ: relaxation time
    para.tau = .1; 
    para.lambda = 0.1 * para.dt / para.tau;
    
%% define variables
    
    %   discretized density and velocity
    traffic.rho_ini = zeros(para.Nx, 1);
    
    traffic.rho = zeros(para.Nx, para.Nt+1);
    traffic.uw = zeros(para.Nx, para.Nw, para.Nt+1);
    traffic.u = zeros(para.Nx, para.Nt+1);    
    traffic.V = zeros(para.Nx, para.Nw, para.Nt+1);
    %equilibrium results
    traffic.uwe = zeros(para.Nx, para.Nw, para.Nt+1);
    traffic.rhoe = zeros(para.Nx, para.Nt+1);
    traffic.we = zeros(para.Nx, para.Nt+1);
    traffic.Ve = zeros(para.Nx, para.Nw, para.Nt+1);
    %arz factor
    traffic.V_arz = zeros(para.Nx, para.Nw, para.Nt+1);
    traffic.u_arz = zeros(para.Nx, para.Nw, para.Nt+1);
    traffic.c_arz=0.5;
    traffic.c_=0;
    traffic.u_ = zeros(para.Nx, para.Nw, para.Nt+1);
    
    % ARZ
    traffic.uw_ini = zeros(para.Nx, para.Nw, 1);
    traffic.u_ini = zeros(para.Nx, 1);
    traffic.z_ini = zeros(para.Nx, 1);
    
    traffic.z = zeros(para.Nx, para.Nt+1);
    traffic.w = zeros(para.Nx, para.Nt+1);
    
    
    % HJB
    traffic.V_ter = zeros(para.Nx, para.Nw, 1);
    %traffic.Vt = 0;
end