function [u Vt] = lwr(rho, Vx)
global para

    u = Ueq(rho, para.uf, para.rhoj) -  para.uf^2 * Vx;
    
    Vt = -u * Vx - 1/2 * (1 - rho/para.rhoj - u/para.uf)^2; % + u/para.uf - u*rho/para.uf/para.rhoj;
end