function tra = Ueq(tra, type, typeFD, uf,rhoj)

if strcmpi(type, 'lwr')    
    if strcmpi(typeFD, 'greenshields')
        tra.fU = @(rho) uf * (1 - rho / rhoj);
        tra.fQ = @(rho) rho .* tra.fU(rho);
    end    
    if strcmpi(typeFD, 'smnewell-daganzo')
        b = 1.0 / 3.0;
        lambda = 0.1;
        c = 0.078 * rhoj * uf;
        g = @(y) sqrt(1 + ((y-b)/lambda).^2);
        tra.fQ = @(rho) c * (g(0) + (g(1) - g(0)) * rho/rhoj - ...
                   g(rho/rhoj));
        tra.fU = @(rho) tra.fQ(rho) ./ rho;
    end
elseif strcmpi(type, 'arz')    
    if strcmpi(typeFD, 'greenshields')
        tra.fU = @(rho,w) max(0, min(uf * (w / uf - rho / rhoj), uf));
%     u = uf * (w ./ uf - rho ./ rhoj);
%     u = max(0, min(u, uf));  

%     % value function
%     u = Ueq_lwr(rho, uf,rhoj) -  para.uf^2 * Vx;
%     u = max(0, min(u, para.uf));    
    
    %tra.Vt = @(rho, w, Vx) -u*Vx - 1/2 * (w/uf - rho/rhoj - u./uf).^2; % + u/para.uf - u*rho/para.uf/para.rhoj;
    tra.Vt = @(rho, w, Vx) -u*Vx - 1/2 * (w/uf - rho/rhoj - u./uf).^2 + w/uf - u*rho/uf/rhoj;
    % other cost functionals
    %tra.Vt = @(rho, w, Vx) -u*Vx - 1/2 * (w/uf - rho/rhoj)^2 - 1/2 * (rho/rhoj)^2 + w/uf - u*rho/uf/rhoj;
    end
    
    tra.fU_lwr = @(rho) uf * (1 - rho / rhoj);
    
    %   hesitation function
    tra.h = @(rho) 0.3 * uf * sqrt((rho/rhoj) ./ (1 - rho/rhoj));
end


end

% function [u] = Ueq_lwr(rho, uf,rhoj)
%     u = uf * (1 - rho ./ rhoj);
%     u = max(0, min(u, uf));
% end

% function [u] = Ueq_arz(rho,w, uf,rhoj)
%     u = uf * (w ./ uf - rho ./ rhoj);
%     u = max(0, min(u, uf)); 
% end