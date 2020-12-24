function V = nonlcon_correct_ver03_new(x_dc,pd_in_opt,y_nom,y0,file_s_name,MCS_bnd)

load(file_s_name)
rho=0.284;W=5.0;T=5.0;Fx=500.0;
n=1.0; % Index of validation metric

nMCS=10^6;
Finv=kde_correct_ver02(x_dc,pd_in_opt);
for k=1:3
    rand_nMCS=(MCS_bnd(2,k)-MCS_bnd(1,k))*rand(nMCS,1)+MCS_bnd(1,k);
    z_MCS(:,k) = feval(Finv{k,:},rand_nMCS);
end

L=z_MCS(:,1);
E=z_MCS(:,2);
Fy=z_MCS(:,3);

xx_set=[W*ones(nMCS,1) T*ones(nMCS,1) L E Fx*ones(nMCS,1) Fy rho*ones(nMCS,1)];
[y(:,1) y(:,2) y(:,3)]=canti_SDF(xx_set);
r=(y-y_nom)./y;

sigma_y=std(y,1);
bw_y=(4/(3*nMCS))^(1/5)*sigma_y;
sigma_y0=std(y0,1);
bw_y0=(4/(3*nMCS))^(1/5)*sigma_y0;
for kk=1:3
pd_y{kk,:} = fitdist(y(:,kk),'Kernel','kernel','normal','support','unbounded','Width',bw_y(kk));
pd_y0{kk,:} = fitdist(y0(:,kk),'Kernel','kernel','normal','support','unbounded','Width',bw_y0(kk));
end

ci_int(:,1)=[4*10^3 0 2]';
ci_int(:,2)=[12*10^3 0.8 12]';
for kk=1:3
q_int(:,kk)=linspace(ci_int(kk,1),ci_int(kk,2),1000);
pdf_y(:,kk)=pdf(pd_y{kk,:},q_int(:,kk));
pdf_y0(:,kk)=pdf(pd_y0{kk,:},q_int(:,kk));
pdf_min(:,kk)=min(pdf_y(:,kk),pdf_y0(:,kk));
Af(kk,:)=trapz(q_int(:,kk),pdf_min(:,kk),1);
end
V=Af.^n;






