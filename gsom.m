function [uw Vt] = gsom(tra, rho, w, Vx, Vw)
global para
    c=0.5;
    uw = tra.fU(rho, w) -  para.uf^2 * (Vx - para.lambda * Vw)/(2*c);
    uw = max(0, min(uw, para.uf));
    
    % w.o. w
    %Vt = -uw*Vx - (1-rho/para.rhoj-u/para.uf)^2/2; % + u/uf - u*rho/uf/rhoj; %LWR
    %Vt = -u*Vx -1/2* (u/para.uf)^2+u/para.uf - u*rho/(para.uf*para.rhoj);     %non-separable
    
    %Vt = -uw * Vx - 1/2 * (w/para.uf - rho/para.rhoj - uw./para.uf).^2; % + u/para.uf - u*rho/para.uf/para.rhoj;
    %Vt = -u*Vx - 1/2 * (w/para.uf - rho/para.rhoj - u./para.uf).^2 + w*rho/para.uf/para.rhoj - u*rho/para.uf/para.rhoj;
    
    %Vt = @(rho, w, Vx) -u*Vx - 1/2 * (w/para.uf - rho/para.rhoj - u./para.uf).^2; % + u/para.uf - u*rho/para.uf/para.rhoj;
    % other cost functionals
    %Vt = -u*Vx - 1/2 * (w/para.uf - rho/para.rhoj)^2 - 1/2 * (rho/para.rhoj)^2;
    
    % with both Vx & Vw (nonseparable)
    %Vt = 1/(2*para.uf) * (tra.fU(rho, w) - para.uf^2 * (Vx - para.lambda*Vw))^2 - para.lambda * tra.fU(rho, w) * Vw...
        %- 1/2*(1-w/para.uf)^2 + 1/2*(1-rho/para.rhoj)^2;
    %with both Vx & Vw (arz)
    Vt=-uw*(Vx-para.lambda*Vw)- para.lambda * tra.fU(rho, w) * Vw...
        -c/(para.uf^2)*(tra.fU(rho, w)-uw)^2; 
    
end