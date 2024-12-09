global para
[row,col]=size(rho_ng);
c=30/1.5;
rho_max=0.45;
u_max=18;

rho_ini=rho_ng(:,1);
u_ini=u_ng(:,1);
w_ini =u_ini + rho_ini.^para.gamma;
z_ini = rho_ini.*w_ini;


w1=zeros(row,col);
for t=1:col
    %w-rho*u_max/rhoj=u
    for i=1:row
        w1(i,t)=u_ng(i,t)+rho_ng(i,t)*u_max/rho_max;
        %w1(i,t) = max(0, min(w1(i,t), u_max));
    end
end

w=zeros(row,col);
z=zeros(row,col);
rho_test_ng=zeros(row,col);
for t=1:col
    if t==1
    z(:,t)=z_ini;
    rho_test_ng(:,t)=rho_ini;
    w(:,t)=w_ini;
    else
        for i =1:row
            %l=i-1;r=i+1;
            %tra1.u(:,t) = mean(tra1.uw(:,:,t), 2);            
            if i > 1 l = i - 1; else l = 1; end
            if i < row r = i + 1; else r = row; end
               
            rho_test_ng(i,t) = 0.5 * (rho_test_ng(l,t-1) + rho_test_ng(r,t-1)) -...
                0.5 * c * (rho_test_ng(r,t-1) * u_ng(r,t-1) -...
                rho_test_ng(l,t-1) * u_ng(l,t-1));
            
            rho_test_ng(i,t) = max(0, min(rho_test_ng(i,t), rho_max));
            
            z(i,t) = 0.5 * (z(l,t-1) + z(r,t-1))...
                - 0.5 * c * (z(r,t-1) * u_ng(r,t-1) - z(l,t-1) * u_ng(l,t-1))...
                + para.lambda * (tra1.fU_lwr(rho_ng(i,t-1)) - u_ng(i,t-1));
    
    
            w(i,t) = z(i,t)/rho_ng(i,t);
            w(i,t) = max(0, min(w(i,t), u_max));
            %w(i,t) = max(0, w(i,t));
        end
    end
end
error_rho=norm(rho_ng(2:r-1,:) - rho_test_ng(2:r-1,:), inf);
error_w=norm(w - w1, inf);
w_=reshape(w1,[row*col,1]);
max_w1=max(w_);
min_w1=min(w_);

uwe=zeros(para.Nx, para.Nw, para.Nt+1)+0.5;
for w_id=1:para.Nw
    uwe(:,w_id,:)=u_ng;
end

