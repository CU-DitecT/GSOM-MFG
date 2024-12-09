function [u] = Ueq_lwr(rho, uf,rhoj)
    u = uf * (1 - rho ./ rhoj);
    u = max(0, min(u, uf));
end
