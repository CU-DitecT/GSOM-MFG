%   Solve the MFG system by fixed-point iterations
function [tra err u_hist V_hist] = solve(tra, maxiter, tol)
global para

%   iterate
    err = zeros(maxiter, 1);
    u_hist=zeros(para.Nx, para.Nt+1, maxiter);
    V_hist=zeros(para.Nx, para.Nw, para.Nt+1, maxiter);
    uw_hist=zeros(para.Nx, para.Nw, para.Nt+1, maxiter);
    disp('Iteration    Error')
    
    for iter = 1 : maxiter
        tra_ = tra;
        
        tra = solve_forward_w(tra);
        tra = solve_backward_w(tra);
        
        
%       %FP
        %
        u_hist(1:para.Nx, 1:para.Nt+1, iter) = tra.u;
        tra.u = squeeze(sum(u_hist, 3))/iter;
        
        V_hist(1:para.Nx, 1:para.Nw, 1:para.Nt+1, iter) = tra.V;
        tra.V = squeeze(sum(V_hist, 4))/iter;
        %}
        uw_hist(1:para.Nx, 1:para.Nw, 1:para.Nt+1, iter) = tra.uw;
        tra.uw = squeeze(sum(uw_hist, 4))/iter;
        
        
        err(iter) = norm(tra.rho - tra_.rho, inf); 
        %+...
        %norm(tra.u - tra_.u, inf);  %norm returns NaN if the input contains NaN values.      
        disp([num2str(iter), '    ', num2str(err(iter))])
        
        if err(iter) <= tol
            return;
        end
    end
    if err(iter) > tol
        disp(['Fixed point iteration does not converge by ', ...
            num2str(maxiter), ' steps']);
        disp(['Error: ', num2str(err(maxiter))]);
        disp('Solver failed');
    end
end