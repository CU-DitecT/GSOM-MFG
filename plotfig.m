%   plot the evolution of density and velocity
%   as a movie
function mov = plotfig(tra)
global para

    %   parameters
%     x = linspace(para.dx/2, para.L-para.dx/2, para.Nx);
    x = linspace(0, para.L, para.Nx);
    
    hfig = figure;
    %   full screen
    set(gcf,'outerposition',get(0,'screensize'));
    %   speed of movie
    fps = 50;
    
    for n = 1 : fps : para.Nt
        plot(x, tra.rho(:,n) / para.rhoj, 'LineWidth', 2, 'Color', 'b');
        xlim([0, para.L]);
        ylim([0, 1]);
        legend('Traffic density');
        title(['time = ', num2str((n-1)*para.dt, '%-1.2f')]);  
        mov(n) = getframe(hfig);
    end    
%     for n = 1 : fps : para.Nt
%         plot(x, tra.rho(:,n) / para.rhoj, 'LineWidth', 2, 'Color', 'b');
%         hold on;
%         plot(x, tra.u(:,n) / para.uf, 'LineWidth', 2, 'Color', 'r');
%         hold off;
%         plot(x, tra.w(:,n) / para.uf, 'LineWidth', 2, 'Color', 'g');
%         hold off;
%         xlim([0, para.L]);
%         ylim([0, 1]);
%         legend('Traffic density', 'Traffic velocity', 'Lagrangian marker');
%         title(['time = ', num2str((n-1)*para.dt, '%-1.2f')]);  
%         mov(n) = getframe(hfig);
%     end
end