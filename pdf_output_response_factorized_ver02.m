function pd_y = pdf_output_response_factorized_ver02(pd_in_opt,file_s_name,Finv,MCS_bnd)
load(file_s_name)

rho=0.284;W=5.0;T=5.0;Fx=500.0;
N=10^6;

nMCS=10^6;
for k=1:3
    rand_nMCS=(MCS_bnd(2,k)-MCS_bnd(1,k))*rand(nMCS,1)+MCS_bnd(1,k);
    z_MCS(:,k) = feval(Finv{k,:},rand_nMCS);
end

L=z_MCS(:,1);
E=z_MCS(:,2);
Fy=z_MCS(:,3);

xx_set_exp=[W*ones(N,1) T*ones(N,1) L E Fx*ones(N,1) Fy rho*ones(N,1)];
[y(:,1) y(:,2) y(:,3)]=canti_SDF(xx_set_exp); % hist(y_exp(:,3))

sigma_y=std(y,1);
bw_y=(4/(3*N))^(1/5)*sigma_y;
for kk=1:3
pd_y{kk,:} = fitdist(y(:,kk),'Kernel','kernel','normal','support','unbounded','Width',bw_y(kk));
end







