function pd_y = pdf_output_response_true(file_s_name)
load(file_s_name)

rho=0.284;W=5.0;T=5.0;Fx=500.0;
N=10^6;

% 시험데이터 생성
E=random('Weibull',2.91*10^7,5,[N 1])+7.36*10^6; % wblrnd(2.91*10^7,7.36*10^6) 
L=random('Normal',100.0,10.0,[N 1]);
Fy=random('Lognormal',6.933,0.082,[N 1]);

xx_set_exp=[W*ones(N,1) T*ones(N,1) L E Fx*ones(N,1) Fy rho*ones(N,1)];
[y(:,1) y(:,2) y(:,3)]=canti_SDF(xx_set_exp); % hist(y_exp(:,3))

sigma_y=std(y,1);
bw_y=(4/(3*N))^(1/5)*sigma_y;
for kk=1:3
pd_y{kk,:} = fitdist(y(:,kk),'Kernel','kernel','normal','support','unbounded','Width',bw_y(kk));
end







