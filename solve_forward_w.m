%   Solve the forward Fokker-Planck equation
%   by a Lax-Friedrichs scheme
function tra = solve_forward_w(tra)
global para

    %   time matching
    tra.rho(:, 1) = tra.rho_ini;
    tra.z(:, 1) = tra.z_ini;
    tra.u = mean(tra.uw, 2);
    tra.u(:, 1) = tra.u_ini;
    tra.w(:,1) = tra.w_ini;
    
    %tra.u = mean(tra.uw, 2);  %tra.uw(para.Nx, para.Nw, para.Nt+1)


for n = 1 : para.Nt
        for j = 1 : para.Nx
            if j > 1 l = j - 1; else l = para.Nx; end
            if j < para.Nx r = j + 1; else r = 1; end
            
            tra.rho(j,n+1) = 0.5 * (tra.rho(l,n) + tra.rho(r,n)) -...
                0.5 * para.c * (tra.rho(r,n) * tra.u(r,n) -...
                tra.rho(l,n) * tra.u(l,n));
            
            tra.rho(j,n+1) = max(0, min(tra.rho(j,n+1), para.rhoj)); 
            
            tra.z(j,n+1) = 0.5 * (tra.z(l,n) + tra.z(r,n))...
                   - 0.5 * para.c * (tra.z(r,n) * tra.u(r,n) - tra.z(l,n) * tra.u(l,n))...
                   + para.lambda * (tra.fU_lwr(tra.rho(j,n)) - tra.u(j,n));


            tra.w(j,n+1) = tra.z(j,n+1)/ tra.rho(j,n+1);
            tra.w(j,n+1) = max(0, min(tra.w(j,n+1), para.uf));
                
                
%            tra.u(j,n+1) = z_temp / tra.rho(j,n+1) - tra.rho(j,n+1)^para.gamma;
%             if tra.rho(j,n+1) > 1e-2
% %             tra.u(j,n+1) = (para.tau * tra.fU_lwr(tra.rho(j,n+1)) + 
% %                              z_temp / tra.rho(j,n+1) - ...
% %                              tra.h(tra.rho(j,n+1))) / (1 + para.tau);                       
%                 tra.u(j,n+1) = z_temp / tra.rho(j,n+1) - tra.rho(j,n+1)^para.gamma;
%                 %%(para.lambda * tra.fU_lwr(tra.rho(j,n+1)) + ...
%             else
% %                 tra.u(j,n+1) = tra.fU_lwr(tra.rho_ini);
%                 tra.u(j,n+1) = Ueq_lwr(tra.rho(j,n+1),para.uf,para.rhoj);
%             end
%            tra.u(j,n+1) = max(0, min(tra.u(j,n+1), para.uf));
                         
%            tra.z(j,n+1) = tra.rho(j,n+1)*(tra.u(j,n+1) + tra.rho(j,n+1)^para.gamma);               
            %tra.z(j,n+1) = tra.rho(j,n+1) * (tra.u(j,n+1) + tra.h(tra.rho(j,n+1)));                      
%        end    
        
%        tra.w = tra.u + tra.rho.^para.gamma;    
        end
end