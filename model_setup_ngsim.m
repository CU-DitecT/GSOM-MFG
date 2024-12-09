%  Set up a traffic flow model:  ring road, single class
function traffic = model_setup_ngsim() % function traffic = model_setup(L, T, uf, rhoj, Nx, Nt)
global para
%% define parameters
 %   time of simulation
    load rho_ng
    load u_ng
    %toy example
    rho_ng=rho_ng(1:5,1:5);
    u_ng=u_ng(1:5,1:5);
    [row,col]=size(rho_ng);


    para.T = col*1.5;
    %   length of the road
    para.L = row*30; %* para.T;
%     para.L = 1;
    %   free flow speed
    para.uf = 18;
    %   jam density
    para.rhoj = 0.45;  
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
    para.c  = 1/20;   %CFL number c<=1: to respect CFL 
    para.Nx = row;
    para.Nt = col;    
    para.dx = 30; 
    para.dt = 1.5; 
    
end
if para.c * para.uf > 1  
	disp('CFL condition not satisfied');
	disp('Solver stopped');
	return
end

%para.Nw = para.W/para.dw;
para.dw = 2;
para.Nw = para.W/para.dw;

%   parameter for inhomo ARZ: relaxation time
    para.tau = .1; 
    para.lambda = 0 * para.dt / para.tau;
    
%% define variables
    w1=zeros(para.Nx, para.Nt+1);
    for t=1:col
        %w-rho*u_max/rhoj=u
        for i=1:row
            %rho_ng(i,t)=max(0, min(rho_ng(i,t), para.rhoj));
            %u_ng(i,t)=max(0, min(u_ng(i,t), para.uf));
            w1(i,t)=u_ng(i,t)+rho_ng(i,t)*para.uf/para.rhoj;
            %w1(i,t) = max(0, min(w1(i,t), para.uf));
        end
    end
    
    u_ng(:,para.Nt+1)=zeros(para.Nx,1);
    rho_ng(:,para.Nt+1)=zeros(para.Nx,1);
       
    uwe=zeros(para.Nx, para.Nw, para.Nt+1)+0.5;
    for w_id=1:para.Nw
        uwe(:,w_id,:)=u_ng; %should be modified
    end

    %   discretized density and velocity
    traffic.rho_ini = zeros(para.Nx, 1);
    
    traffic.rho = zeros(para.Nx, para.Nt+1);
    traffic.uw = zeros(para.Nx, para.Nw, para.Nt+1);
    traffic.u = zeros(para.Nx, para.Nt+1);    
    traffic.V = zeros(para.Nx, para.Nw, para.Nt+1);
    
    %equilibrium results
    traffic.uwe = uwe;
    traffic.rhoe = rho_ng;
    traffic.we = w1;
    traffic.ue=u_ng;
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