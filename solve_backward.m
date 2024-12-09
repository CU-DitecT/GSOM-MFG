%   Solve the backward HJB equation
%   by an upwind scheme
function tra = solve_backward(tra)
global para

    tra.V(:, para.Nt+1) = tra.V_ter;
    
    for n = para.Nt : -1 : 1
        for j = 1 : para.Nx
            if j < para.Nx r = j + 1; else r = 1; end
            
            Vx =(tra.V(r,n+1) - tra.V(j,n+1)) / para.dx;
            %tra.u(j,n)=lwr_f(tra.rho(j,n));
            %tra.V(j,n)=0;%tra.dt*0.5*(1-tra.rho(j,n)-tra.u(j,n))^2+(1-tra.u(j,n))*tra.V(j,n+1)+tra.u(j,n)*tra.V(r,n+1);
            
            
            [tra.u(j,n), Vt] = tra.vlf(tra.rho(j,n), tra.w(j,n), Vx);
            %tra.u(j,n) = tra.vlf(tra.rho(j,n), tra.w(j,n), Vx);
            
            
            tra.V(j,n) = tra.V(j,n+1) - Vt * para.dt;
            
            % stablize convergence
            if tra.V(j,n)<1e-3  
                tra.V(j,n)=0;
            end
            
        end
    end
end