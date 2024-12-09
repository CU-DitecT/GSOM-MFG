%   Solve the backward HJB equation
%   by an upwind scheme
function tra = solve_backward_w(tra)
global para

    tra.V(:,:, para.Nt+1) = tra.V_ter;
    %tra.V_arz(:,:, para.Nt+1) = tra.V_ter;
    
    for n = para.Nt : -1 : 1
        for j = 1 : para.Nx
	        for k = 1 : para.Nw

            if j < para.Nx-1 
                r = j + 2; 
            elseif j == para.Nx-1
                r = 1;
            else
                r = 2; 
            end

	        if k < para.Nw-1 
                rw = k + 2; 
            elseif k == para.Nw-1 
                rw = 1;
            else
                rw = 2;
            end	

            Vx = (tra.V(r,k, n+1) - tra.V(j,k, n+1)) / para.dx / 2;
	        Vw = (tra.V(j,rw, n+1) - tra.V(j,k, n+1)) / para.dw / 2;

            %V_arz_x = (tra.V_arz(r,k, n+1) - tra.V_arz(j,k, n+1)) / para.dx / 2;
	        %V_arz_w = (tra.V_arz(j,rw, n+1) - tra.V_arz(j,k, n+1)) / para.dw / 2;
            
            %upwind
            % if j < para.Nx r = j + 1; else r = 1; end
	        % if k < para.Nw rw = k + 1; else rw = 1; end	
            % Vx = (tra.V(r,k, n+1) - tra.V(j,k, n+1)) / para.dx;
	        % Vw = (tra.V(j,rw, n+1) - tra.V(j,k, n+1)) / para.dw;
            
            
            % %central difference
            % Vx = (tra.V(r,k, n+1) - tra.V(j,k, n+1)) / para.dx;
	        % Vw = (tra.V(j,rw, n+1) - tra.V(j,k, n+1)) / para.dw;
            %tra.u(j,n)=lwr_f(tra.rho(j,n));
            %tra.V(j,n)=0;%tra.dt*0.5*(1-tra.rho(j,n)-tra.u(j,n))^2+(1-tra.u(j,n))*tra.V(j,n+1)+tra.u(j,n)*tra.V(r,n+1);
            
            
            [tra.uw(j,k, n), Vt] = tra.vlf(tra.rho(j,n), tra.w(j,n), Vx, Vw);
            %tra.u(j,n) = tra.vlf(tra.rho(j,n), tra.w(j,n), Vx);
            

            tra.V(j,k,n) = tra.V(j,k,n+1) - Vt * para.dt;
            %disp(sprintf('k=%d, vx=%f, vw=%f', k, Vx, Vw));
            % stablize convergence
            if tra.V(j,k,n)<1e-7  
                tra.V(j,k,n)=0;
            end
            
        end
    end
end