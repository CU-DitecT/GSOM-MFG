tr = model_setup_w(); 

global para
tr = Ueq(tr, 'arz', 'greenshields', para.uf,para.rhoj);
tr = init(tr);

tr.rho(:,1) = tr.rho_ini;
tr.u(:,1) = tr.u_ini;
tr.z(:,1) = tr.z_ini;
%tr.z(:,1) = tr.z_ini;


%% trime Evolutrion
rho = []; % Initialize rho as an empty array
u = []; % Initialize u as an empty array

for lb = 0.05:0.1:0.45
    for ub = 0.95:-0.1:0.55
        for exp = 0.15:0.1:0.55
            tr = init_bellshape_forloop(tra, lb, ub, exp);
            tr.rho(:,1) = tr.rho_ini;
            tr.u(:,1) = tr.u_ini;
            tr.z(:,1) = tr.z_ini;
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
                                    -num_flux(tr.z(i,n),tr.z(l,n), tr.z(l,n)*tr.u(l,n), tr.z(i,n)*tr.u(i,n)))...
                                    + para.lambda * (Ueq_lwr(tr.rho(i,n),para.uf,para.rhoj) - tr.u(i,n));
                end
                
                %updatre of trhe velocitry
                tr.u(:,n+1)=tr.z(:,n+1)./tr.rho(:,n+1)-tr.rho(:,n+1).^para.gamma;
            end
            rho = [rho; tr.rho]; % Concatenate tra.rho to rho
            u = [u; tr.u]; % Concatenate tra.u to u
        end
    end
end
writematrix(rho, "for_fd/bellshape/zm_bellshape_homoarz_rho.csv")
writematrix(u, "for_fd/bellshape/zm_bellshape_homoarz_u.csv")


function [flux] = num_flux(flux_r,flux_l,rho_l,rho_r) %num_flux(rho_r,rho_l,flux_l,flux_r) 
% Numerical flux for Lax Friedrichs scheme
global para
flux = 0.5 * (rho_l+rho_r)...
      -0.5 * para.dx/para.dt * (flux_r-flux_l);
end
