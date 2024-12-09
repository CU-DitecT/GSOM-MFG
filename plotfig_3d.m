%   Show the density evolution in 3D
%   space-time-density plot
function plotfig_3d(tra)
global para

%     x = linspace(para.dx/2, para.L-para.dx/2, para.Nx);
%     t = linspace(0, para.T-para.dt, para.Nt);
    x = linspace(0, para.L, para.Nx);
    t = linspace(0, para.T, para.Nt);
    [T X] = meshgrid(t, x);
    
    %subplot(4,1,1);
    figure;
    rho = mesh(T, X, tra.rho(:,1:para.Nt));
    colormap('jet');
    rho.FaceColor='interp';
    % xlabel('time');
    % ylabel('space');
    zlabel('\rho');
    %zlabel('Density');
    
    %subplot(4,1,2);
    figure;
    u = mesh(T, X, tra.u(:,1:para.Nt));
    colormap('jet');
    u.FaceColor='interp';
    % xlabel('time');
    % ylabel('space');
    zlabel('u');
    %zlabel('Velocity');
    
    %subplot(4,1,3);
    figure;
    w = mesh(T, X, tra.w(:,1:para.Nt));
    colormap('jet');
    w.FaceColor='interp';
    % xlabel('time');
    % ylabel('space');
    % zlabel('Lagrangian marker');
    zlabel('w');
    
    % subplot(4,1,4);
    % V = mesh(T, X, tra.V(:,1:para.Nt));
    % colormap('jet');
    % V.FaceColor='interp';
    % % xlabel('time');
    % % ylabel('space');
    % % zlabel('Value');
    % xlabel('t');
    % ylabel('x');
    % zlabel('V');
end