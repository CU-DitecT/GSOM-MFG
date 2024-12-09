function [mu_U_,mu_U_u_,mu_u_2_]=mu_inter_arz_ngsim(tra,u_)
global para
mu_U_=zeros(para.Nx, para.Nt+1);
mu_U_(:,para.Nt+1) = tra.V_ter;

mu_U_u_=zeros(para.Nx, para.Nt+1);
mu_U_u_(:,para.Nt+1) = tra.V_ter;

mu_u_2_=zeros(para.Nx, para.Nt+1);
mu_u_2_(:,para.Nt+1) = tra.V_ter;


for t=para.Nt:-1:1
    for i = 1 : para.Nx-1
            %c1
            mu_U_x = (mu_U_(i+1,t+1) - mu_U_(i,t+1)) / para.dx;
            mu_U_t=-u_(i,t)*(mu_U_x) ...
               -1/(para.uf^2)*(para.uf*(tra.we(i,t)/para.uf-tra.rhoe(i,t)/para.rhoj)-u_(i,t))^2; %cost 1
            mu_U_(i,t) = mu_U_(i,t+1) - mu_U_t * para.dt;
            if mu_U_(i,t)<1e-7 mu_U_(i,t)=0; end
            %c2
            mu_U_u_x = (mu_U_u_(i+1,t+1) - mu_U_u_(i,t+1)) / para.dx ;
            mu_U_u_t=-u_(i,t)*(mu_U_u_x) ...
               -(1-tra.we(i,t)/para.uf)^2;
            mu_U_u_(i,t) = mu_U_u_(i,t+1) - mu_U_u_t * para.dt;
            if mu_U_u_(i,t)<1e-7 mu_U_u_(i,t)=0; end
            %c3
            mu_u_2_x = (mu_u_2_(i+1,t+1) - mu_u_2_(i,t+1)) / para.dx;
            mu_u_2_t=-u_(i,t)*(mu_u_2_x) ...
               +(1-tra.rhoe(i,t)/para.rhoj)^2;
            mu_u_2_(i,t) = mu_u_2_(i,t+1) - mu_u_2_t * para.dt;
            if mu_u_2_(i,t)<1e-7 mu_u_2_(i,t)=0; end
            %{
            V_test_x = (V_test_(r,k, t+1) - V_test_(i,k, t+1)) / para.dx / 2;
	        V_test_w = (V_test_(i,rw, t+1) - V_test_(i,k, t+1)) / para.dw / 2;

            V_test_t=-uw_(i,k,t)*(V_test_x-para.lambda*V_test_w)- para.lambda * tra.fU(tra.rhoe(i,t), tra.we(i,t)) * V_test_w...
                     -0.5*1/(para.uf^2)*(tra.fU(tra.rhoe(i,t), tra.we(i,t))-uw_(i,k,t))^2 ...
                     -0.5*(1-tra.we(i,t)/para.uf)^2 ...
                     +0.5*(1-tra.rhoe(i,t)/para.rhoj)^2;
              %cost 1
            V_test_(i,k,t) = V_test_(i,k,t+1) - V_test_t * para.dt;
            if V_test_(i,k,t)<1e-7 V_test_(i,k,t)=0; end
            %}

    end
end