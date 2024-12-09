%mean field RL
global para;
tra1 = model_setup_w();

%   set Ueq (LWR or ARZ) & FD type
tra1 = Ueq(tra1, 'arz', 'greenshields', para.uf,para.rhoj);
tra1 = init(tra1);
para.lambda = 0.07 * para.dt / para.tau;
tic;
% total = 5*5*5;
total = 2 * 2 * 2;
assemble.Rhoe = zeros(para.Nx, para.Nt+1,total);
assemble.We = zeros(para.Nx, para.Nt+1,total);
assemble.Uwe = zeros(para.Nx, para.Nw, para.Nt+1,total);
assemble.Ue = zeros(para.Nx, para.Nt+1,total);
assemble.Ve = zeros(para.Nx, para.Nw, para.Nt+1,total);
count=0;
for lb = 0.05:0.3:0.45
    for ub = 0.95:-0.3:0.55
        for exp = 0.15:0.3:0.55
            count = count + 1;
            maxiter=200;
            tra1 = init_bellshape_forloop(tra1, lb, ub, exp);
            uw_hist=zeros(para.Nx, para.Nw, para.Nt+1, maxiter);
            tra1.V(:,:, para.Nt+1) = tra1.V_ter;
            err = zeros(maxiter, 1);
            for iter=1:maxiter
                tra1_=tra1;
                for t=1:para.Nt+1
                    for i =1:para.Nx
                        %forward
                        if t==1
                            tra1.rho(i,t) = tra1.rho_ini(i);
                            tra1.z(i,t) = tra1.z_ini(i);
                            tra1.u(i,t) = tra1.u_ini(i);
                            tra1.w(i,t) = tra1.w_ini(i);
                        else
                            tra1.u(:,t) = mean(tra1.uw(:,:,t), 2);
                            if i > 1 l = i - 1; else l = para.Nx; end
                            if i < para.Nx r = i + 1; else r = 1; end

                            tra1.rho(i,t) = 0.5 * (tra1.rho(l,t-1) + tra1.rho(r,t-1)) -...
                                0.5 * para.c * (tra1.rho(r,t-1) * tra1.u(r,t-1) -...
                                tra1.rho(l,t-1) * tra1.u(l,t-1));

                            tra1.rho(i,t) = max(0, min(tra1.rho(i,t), para.rhoj));

                            tra1.z(i,t) = 0.5 * (tra1.z(l,t-1) + tra1.z(r,t-1))...
                                - 0.5 * para.c * (tra1.z(r,t-1) * tra1.u(r,t-1) - tra1.z(l,t-1) * tra1.u(l,t-1))...
                                + para.lambda * (tra1.fU_lwr(tra1.rho(i,t-1)) - tra1.u(i,t-1));


                            tra1.w(i,t) = tra1.z(i,t)/ tra1.rho(i,t);
                            tra1.w(i,t) = max(0, min(tra1.w(i,t), para.uf));

                        end
                        if t<para.Nt+1
                            for k = 1 : para.Nw
                                if i < para.Nx-1
                                    r = i + 2;
                                elseif i == para.Nx-1
                                    r = 1;
                                else
                                    r = 2;
                                end
                                if k < para.Nw-1
                                    rw = k + 2;
                                elseif k == para.Nw-1
                                    rw = 1;
                                else
                                    rw = 2;
                                end
                                Vx = (tra1.V(r,k, t+1) - tra1.V(i,k, t+1)) / para.dx / 2;
                                Vw = (tra1.V(i,rw, t+1) - tra1.V(i,k, t+1)) / para.dw / 2;
                                %arz_cost
                                %{
                %c=0.5; % if only use one feature, then any c (>0) can be the coefficient. (if separate into different terms)
                c1=0.5; c2=1; c3=0.5;
                tra1.uw(i,k,t) = c2*tra1.fU(tra1.rho(i,t), tra1.w(i,t))/(2*c3) -  para.uf^2 * (Vx - para.lambda * Vw)/(2*c3);
                tra1.uw(i,k,t) = max(0, min(tra1.uw(i,k,t), para.uf)); %fp buffer
                uw_hist(i,k,t, iter) = tra1.uw(i,k,t);
                %
                Vt=-tra1.uw(i,k,t)*(Vx-para.lambda*Vw)- para.lambda * tra1.fU(tra1.rho(i,t), tra1.w(i,t)) * Vw...
                    -c1/(para.uf^2)*tra1.fU(tra1.rho(i,t), tra1.w(i,t))^2 ...
                    +c2/(para.uf^2)*tra1.fU(tra1.rho(i,t), tra1.w(i,t))*tra1.uw(i,k,t) ...
                    -c3/(para.uf^2)*tra1.uw(i,k,t)^2;
                                %}
                                %{
                Vt=-tra1.uw(i,k,t)*(Vx-para.lambda*Vw)- para.lambda * tra1.fU(tra1.rho(i,t), tra1.w(i,t)) * Vw...
                    -c1/(para.uf^2)*tra1.fU(tra1.rho(i,t), tra1.w(i,t)-tra1.uw(i,k,t))^2;
                                %}
                                %}
                                %arz_nonsep
                                %
                                c1=0.5; c2=0.5; c3=0.5; %1. equilibrium speed 2. w-dependency 3. safety
                                tra1.uw(i,k,t) = tra1.fU(tra1.rho(i,t), tra1.w(i,t)) -  para.uf^2 * (Vx - para.lambda * Vw)/(2*c1);
                                tra1.uw(i,k,t) = max(0, min(tra1.uw(i,k,t), para.uf));
                                uw_hist(i,k,t, iter) = tra1.uw(i,k,t); %fp buffer
                                Vt=-tra1.uw(i,k,t)*(Vx-para.lambda*Vw)- para.lambda * tra1.fU(tra1.rho(i,t), tra1.w(i,t)) * Vw...
                                    -c1/(para.uf^2)*(tra1.fU(tra1.rho(i,t), tra1.w(i,t))-tra1.uw(i,k,t))^2 ...
                                    -c2*(1-tra1.w(i,t)/para.uf)^2 ...
                                    +c3*(1-tra1.rho(i,t)/para.rhoj)^2; %cost function
                                %}

                                tra1.V(i,k,t) = tra1.V(i,k,t+1) - Vt * para.dt;
                                %if tra1.V(i,k,t)<1e-7 tra1.V(i,k,t)=0; end
                            end
                        end
                    end
                end
                tra1.uw = squeeze(sum(uw_hist, 4))/iter;
                err(iter) = norm(tra1.rho - tra1_.rho, inf);
                %disp([num2str(iter), '    ', num2str(err(iter))])
            end
            disp([num2str(count), '    ', num2str(err(iter))])
            toc;
            tra1.rhoe=tra1.rho;
            tra1.we=tra1.w;
            tra1.uwe=tra1.uw;
            tra1.ue=tra1.u;
            tra1.Ve=tra1.V;
            assemble.Rhoe(:,:,count) = tra1.rhoe;
            assemble.We(:,:,count) = tra1.we;
            assemble.Uwe(:,:,:,count) = tra1.uwe;
            assemble.Ue(:,:,count) = tra1.ue;
            assemble.Ve(:,:,:,count) = tra1.Ve;
        end
    end
end
% save("arz_nonsep_fd.mat","assemble");