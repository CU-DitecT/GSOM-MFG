%   initial and terminal condition
%% set parameters: para
global para

% Define rho, u, z, w
tra = model_setup_w();

% Set Ueq (LWR or ARZ) & FD type
tra = Ueq(tra, 'arz', 'greenshields', para.uf, para.rhoj);

% Cost function
% tra.vlf = @(rho, w, Vx) arz(tra, rho, w, Vx);
tra.vlf = @(rho, w, Vx, Vw) gsom(tra, rho, w, Vx, Vw);
% tra.vlf_arz = @

% Initial and terminal condition
tra = init(tra);

maxiter = 200;
exp_num = 0;
rho = zeros(para.Nx, para.Nt+1, exp_num); % Initialize rho as an empty array
u = zeros(para.Nx, para.Nt+1, exp_num);
uw = zeros(para.Nx, para.Nw, para.Nt+1, exp_num);
for lb = 0.05:0.1:0.45
    for ub = 0.95:-0.1:0.55
        for exp = 0.15:0.1:0.55
            exp_num = exp_num + 1;
            tra = init_bellshape_forloop(tra, lb, ub, exp);
            [tra, err, u_hist, V_hist] = solve(tra, maxiter, 1e-2);
            rho(1:end, 1:end,exp_num) =tra.rho; % Concatenate tra.rho to rho
            u(1:end, 1:end,exp_num) = tra.u; % Concatenate tra.u to u
            uw(1:end, 1:end, 1:end, exp_num) = tra.uw;
        end
    end
end
a=1;
% Write the results to CSV files
%writematrix(rho, 'for_fd/bellshape/zm_bellshape_arz_rho.csv');
%writematrix(u, 'for_fd/bellshape/zm_bellshape_arz_u.csv');

%writematrix(tra.rho, 'for_fd/onebellshape/zm_bellshape_gsom_rho.csv');
%writematrix(tra.u, 'for_fd/onebellshape/zm_bellshape_gsom_u.csv');
%writematrix(tra.uw, 'for_fd/onebellshape/zm_bellshape_gsom_uw.csv');
%save("for_fd/bellshape_w/rho_u_uw.mat", "rho", "u", "uw");

%% solve FP