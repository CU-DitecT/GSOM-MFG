function plot_fd()

MARKER_SIZE = 1;
PLOT_FD = true;
PLOT_w = false;
NUM_FD = 7;

% Jitter amount
jitterAmount = 0.01;



%rho = csvread('for_fd/bellshape/zm_bellshape_arz_GSOM_rho.csv');
%rho = rho(:,1:300);

rho = load("rho_ng.mat");
rho = rho.rho_ng;
u = load("u_ng.mat");
u = u.u_ng;

% Add jittering
rng(1,"twister");
rho = rho + jitterAmount * randn(size(rho));



%u = csvread("for_fd/bellshape/zm_bellshape_arz_GSOM_u.csv");
%u = u(:,1:300);
q = rho.*u;
max_q = max(q, [], 'all');

%plot the scatter
figure;
scatter(rho(:, 2:end-1), q(:, 2:end-1), 0.3, 'k', 'filled','MarkerFaceAlpha',.0,'MarkerEdgeAlpha',.0); %。8 
hold on;
h1 = scatter(rho(1, 1), q(1, 1), 0.3, 'k', 'filled','MarkerFaceAlpha',.0,'MarkerEdgeAlpha',.0); % .2
xlabel('density','FontSize', 12);                
ylabel('flow','FontSize', 12);                   
%xlim([0,1]);
ylim([0, 1.2*max_q]);
set(gca, 'FontSize', 12); % Set tick label font size
%grid on;    
axis off;

%plot the fd
if PLOT_FD 
    rho_in = linspace(0, 1, 100);
    % fit the Q_eq
    for beta=[0.5]
    %for beta = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8] 
        params0 = [0.5, 20]; % rhomax and umax
        fittedParams = lsqnonlin(@(params) myObjective(params, rho, q, beta) ...
            , params0);
        
        % arbitrarily set the parameters
        %fittedParams = [1, 1];
        %fittedParams(1) =  0.4524;
    
        if PLOT_w
        % plot Q with different w
        u_max = fittedParams(2); % maximum value of w
        lb = 0.85;
        ub = 1.15;
        % for w = lb*u_max: (ub-lb)*u_max / NUM_FD : ub*u_max
        %     q_pred = Q_eq(fittedParams, rho_in);
        % 
        %     % below is eq. 3.8 of Shimao's dissertation
        %     % https://www.dropbox.com/s/fzybjd7i3ct1epg/tudiss_shimao.pdf
        %     %q_pred = q_pred + rho_in.*(w - u_max);
        %     hold on;
        %     plot(rho_in, q_pred, 'b-', 'LineWidth', 1);
        % end
        end
        % plot the equilibrium one
        q_pred = Q_eq(fittedParams, rho_in);
        hold on;
        hold on;
        if beta == 0.5
            h2=plot(rho_in, q_pred, 'r-', 'LineWidth', 3);
        else
            h3=plot(rho_in, q_pred, 'b-', 'LineWidth', 1);
        end
    end

end
ylim([0, 4]);
hold off;
%legend([h1 h2 h3],{'sensor data','equilibrium curve','family of flow rate curves'})
%ylim([0, 0.5]);
xticks([0, 0.2,  0.4, 0.6]);
yticks([0, 1, 2, 3, 4]);
set(gca, 'FontSize', 16);
f = gcf;
set(f, 'Units', 'inches', 'Position', [1, 1, 5.8, 4.3]);
exportgraphics(f,'new_FD.png','Resolution',300)
%saveas(gcf, 'scatter_withFD_w_gsom.png');
end



function q = Q_eq(params, rho)

% green shield
RHO_MAX = params(1);
U_MAX = params(2);
q = U_MAX .* rho .* (1-rho/RHO_MAX);

end



function modified_residuals = myObjective(params, rho, q, beta)
q_pred = Q_eq(params, rho);
residuals = q_pred - q;  
modified_residuals = residuals;

% Find the indices of positive elements
positive_indices = residuals > 0;

% Find the indices of negative elements
negative_indices = residuals < 0;

% Multiply the positive elements by 0.1
modified_residuals(positive_indices) = residuals(positive_indices) * beta;

% Multiply the negative elements by 0.9
modified_residuals(negative_indices) = residuals(negative_indices) * (1-beta);
end