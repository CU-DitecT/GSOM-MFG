%  Set up a traffic flow model:  ring road, single class
function traffic = model_setup() % function traffic = model_setup(L, T, uf, rhoj, Nx, Nt)
global para
%% define parameters
 %   time of simulation
    para.T = 1;
    %   length of the road
    para.L = 1 * para.T;
%     para.L = 1;
    %   free flow speed
    para.uf = 1;
    %   jam density
    para.rhoj = 1;  
    %   parameter for w
    para.gamma = 1; 
    
    %   coefficient for Nx Nt   
    para.multiplier = 100;
    
    ind_para = 0; 
if ind_para == 0
    %   spatial and time mesh size
%     para.Nx = 2 * para.multiplier;
%     para.Nt = 2 * para.multiplier; 
    para.Nx = 1 * para.multiplier;
    para.Nt = 3 * para.multiplier;

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
    para.Nx = 16;
    para.Nt = 64;    
    para.dx = para.L/para.Nx; 
    para.dt = para.T/para.Nt;               
    
end
if para.c * para.uf > 1  
	disp('CFL condition not satisfied');
	disp('Solver stopped');
	return
end

%   parameter for inhomo ARZ: relaxation time
    para.tau = .1; 
    para.lambda = para.dt / para.tau;
    
%% define variables
    
    %   discretized density and velocity
    traffic.rho_ini = zeros(para.Nx, 1);
    
    traffic.rho = zeros(para.Nx, para.Nt+1);
    traffic.u = zeros(para.Nx, para.Nt+1);    
    traffic.V = zeros(para.Nx, para.Nt+1);
    
    % ARZ
    traffic.u_ini = zeros(para.Nx, 1);
    traffic.z_ini = zeros(para.Nx, 1);
    
    traffic.z = zeros(para.Nx, para.Nt+1);
    traffic.w = zeros(para.Nx, para.Nt+1);
    
    % HJB
    traffic.V_ter = zeros(para.Nx, 1);
    %traffic.Vt = 0;
end