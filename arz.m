function [u Vt] = arz(tra, rho, w, Vx)
global para

    u = tra.fU(rho, w) -  para.uf^2 * Vx;
    u= max(0, min(u, para.uf));
    
    % w.o. w
    Vt = -u*Vx - (1-rho/para.rhoj-u/para.uf)^2/2; % + u/uf - u*rho/uf/rhoj; %LWR
    %Vt = -u*Vx -1/2* (u/para.uf)^2+u/para.uf - u*rho/(para.uf*para.rhoj);     %non-separable
    
    %Vt = -u*Vx - 1/2 * (w/para.uf - rho/para.rhoj - u./para.uf).^2; % + u/para.uf - u*rho/para.uf/para.rhoj;
    %Vt = -u*Vx - 1/2 * (w/para.uf - rho/para.rhoj - u./para.uf).^2 + w*rho/para.uf/para.rhoj - u*rho/para.uf/para.rhoj;
    
    %Vt = @(rho, w, Vx) -u*Vx - 1/2 * (w/para.uf - rho/para.rhoj - u./para.uf).^2; % + u/para.uf - u*rho/para.uf/para.rhoj;
    % other cost functionals
    %Vt = -u*Vx - 1/2 * (w/para.uf - rho/para.rhoj)^2 - 1/2 * (rho/para.rhoj)^2;
    
    %Vt = -u*Vx - 1/2 * (w/para.uf)^2 + (u *rho)/(para.uf * para.rhoj) - u./para.uf;
end


% function [u Vt] = arz(rho, w, Vx)
% global para
% 
%     u = Ueq_arz(rho, w, para.uf, para.rhoj) -  para.uf^2 * Vx;
%     %u= max(0, min(u, para.uf));
%     
%     Vt = -u*Vx - 1/2 * (w/para.uf - rho/para.rhoj - u/para.uf)^2; % + u/para.uf - u*rho/para.uf/para.rhoj;
%     
%     % other cost functionals
%     %Vt = -u*Vx - 1/2 * (w/para.uf - rho/para.rhoj)^2 - 1/2 * (rho/para.rhoj)^2;
% end