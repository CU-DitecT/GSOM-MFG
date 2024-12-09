function tra = init(tra)
global para
% define intial/terminal conditions
%     rhob = 0.4; %default: 0.4
%     tra = set_rho_ini(tra, 'unifpert', rhob); 
    
    %tra = set_rho_ini(tra, 'rectwave', .8,.2);  %tra = set_rho_ini(tra, 'Riemann'); %, .8,.2,.1,.7);  
    %tra = set_rho_ini(tra, 'bellshape');
    %tra = set_rho_ini(tra, 'sinewave');    

    %tra = set_rho_ini(tra,'bellshape',0.05,0.95,0.15);
    % tra = set_rho_ini(tra,'bellshape',0.05,0.95,0.25);
    tra = set_rho_ini(tra,'bellshape',0.05,0.95,0.35);
    %tra = set_rho_ini(tra,'bellshape',0.05,0.95,0.45);
    
    tra.V_ter = 0 * ones(para.Nx, para.Nw, 1); %tra.V_ter = zeros(para.Nx, 1);
    for i = 1:para.Nw
        all_w = linspace(0, 1, para.Nw);
        for j = 1:para.Nx
            tra.V_ter(j, i, 1) = all_w(i) * j;%(i+j)/(para.Nx+para.Nw); %all_w(i)*j;
        end
    end

    % for i = 1:para.Nw
    %     all_w = linspace(0, 1, para.Nw);
    %     all_x = linspace(0, 1, para.Nx);
    %     for j = 1:para.Nx
    %         tra.V_ter(j, i, 1) = all_w(i)*all_x(j);
    %     end
    % end
    
    %tra.V_ter(1:round(para.Nx/2),1)=0; tra.V_ter(round(para.Nx/2)+1:end,1)=1;
    %tra.V_ter = linspace(0, para.L, para.Nx);
end

function tra = set_rho_ini(tra, type, varargin)
global para

    if strcmpi(type, 'bellshape')
        %   default parameters
        rmin = 0.05;  rmax = 0.95;
        sigma = 0.1;
        %   user-specified parameters
        if nargin > 2
            rmin = varargin{1};
            rmax = varargin{2};
        end
        if nargin > 4
            sigma = varargin{3};
        end
        
        x = linspace(-para.L/2 + para.dx/2, para.L/2-para.dx/2, para.Nx);
        tra.rho_ini = rmin + (rmax-rmin) * exp(-x.^2/sigma^2);
        
        pert = 0.05 * 2;
        nw = 1;
        %% u is Rieman problem
        a=zeros(1,para.Nx);
        a(:,1:round(para.Nx/2))=0.1; a(:,round(para.Nx/2)+1:end)=0.8;
        %a(:,1:round(para.Nx/3))=0.2; a(:,round(para.Nx/3)+1:round((2*para.Nx/3)))=0.5; a(:,round(2*para.Nx/3)+1:para.Nx)=0.9;
        tra.u_ini = a; %tra.u_ini(1:round(para.Nx/2),1)=0.1; tra.u_ini(round(para.Nx/2)+1:end,1)=0.7;
        
        %tra.u_ini = Ueq_lwr(tra.rho_ini, para.uf, para.rhoj); %...
                  %+ pert * sin(linspace(0, nw*2*pi, para.Nx));
        %tra.u_ini = 1.1* tra.fU_lwr(tra.rho_ini);
        
        tra.w_ini = tra.u_ini + tra.rho_ini.^para.gamma; %...
                  %- pert * sin(linspace(0, nw*2*pi, para.Nx));
    end
    
    if strcmpi(type, 'sinewave')
        %   default parameters
        nw = 2;
        %   user-specified parameters
        if nargin > 2
            nw = varargin{1};
        end
        tra.rho_ini = (1 + sin(linspace(0, nw*2*pi, para.Nx))) / 2;
        tra.u_ini = (1 + cos(linspace(0, nw*2*pi, para.Nx))) / 2;
        tra.w_ini = tra.u_ini + tra.rho_ini.^para.gamma;
    end
    %rho(:,1)=sin(x(:));
    %u(:,1)=cos(x(:));

    
    if strcmpi(type, 'unifpert')
        %   default parameters
        rhob = 0.5;
        pert = 0.05;
        nw = 1;
        %   user-specified parameters
        if nargin == 3
            rhob = varargin{1};
            pert = rhob / 10;
        end
        if nargin > 3
            rhob = varargin{1};
            pert = varargin{2};
        end
        if nargin > 4
            nw = varargin{3};
        end
        tra.rho_ini = rhob + pert * sin(linspace(0, nw*2*pi, para.Nx));
        tra.u_ini = Ueq_lwr(tra.rho_ini, para.uf, para.rhoj);
        
        tra.w_ini = tra.u_ini + tra.rho_ini.^para.gamma;
    end

%     if strcmpi(type, 'Riemann')
%         % Riemann Problem
%         tra.rho_ini(1:round(para.Nx/2))=0.8 -.0; tra.rho_ini(round(para.Nx/2)+1:end)=0.2 +.0;
%         tra.u_ini(1:round(para.Nx/2),1)=0.1; tra.u_ini(round(para.Nx/2)+1:end,1)=0.7;
%         
%         tra.w_ini = tra.u_ini + tra.rho_ini.^para.gamma;  %.65 * ones(para.Nx, 1);    
%         %tra.u_ini = Ueq_arz(tra.rho_ini, tra.w_ini, para.uf, para.rhoj);
%     end    
    if strcmpi(type, 'rectwave') % Riemann == rectwave
        %   default parameters
        rmin = 0.5;
        rmax = 0.75;
        %   user-specified parameters
        if nargin > 2
            rmin = varargin{1};
            rmax = varargin{2};
        end
        if mod(para.Nx, 2) ~= 0
            disp('Error : Nx should be even!')
            return
        end
        tra.rho_ini = [rmin * ones(para.Nx/2, 1);...
                       rmax * ones(para.Nx/2, 1)]; 
        tra.u_ini = [.3 * ones(para.Nx/2, 1);...
                     .8 * ones(para.Nx/2, 1)]; %[1 * rmin * ones(para.Nx/2, 1);1 * rmax * ones(para.Nx/2, 1)];           
        %tra.u_ini = Ueq_lwr(tra.rho_ini, para.uf, para.rhoj);
        
        tra.w_ini = tra.u_ini + tra.rho_ini.^para.gamma;  
        %tra.w_ini = .65 * ones(para.Nx, 1);  
    end

    
    tra.rho_ini = max(0, min(tra.rho_ini, para.rhoj));
    tra.u_ini = max(0, min(tra.u_ini, para.uf));

    %tra.w_ini = tra.rho_ini;
    %tra.w_ini = tra.w_ini + .1 * ones(1,para.Nx);
    tra.w_ini = max(0, min(tra.w_ini, para.uf));
    
    tra.z_ini = tra.rho_ini .* tra.w_ini;
    %tra.z_ini = tra.rho_ini .* (tra.u_ini + tra.h(tra.rho_ini));
end