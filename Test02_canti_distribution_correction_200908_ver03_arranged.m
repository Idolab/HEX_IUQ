%% Section 3.2 & 3.3
clc
close all
clear

%% Distribution correction
load([pwd,'\Results_section_4_2\results_trial1\','UQ_results_case4_exp_set18.mat']) % Section 3.2 결과 불러오기

% % seed 생성 - MCS
% mkdir([pwd,'\s_save_200720'])
% for trial=1
%     file_name=[pwd,'\s_save_200720','\',sprintf('s_save_trial%d_s_save_200720.mat',trial)];
%     s=rng('shuffle');% 함수 밖에서
%     save(file_name,'s')
% end

%% Nominal Output Response
file_s_name=[pwd,'\s_save_200720','\',sprintf('s_save_trial%d_s_save_200720.mat',1)];
rho=0.284;W=5.0;T=5.0;Fx=500.0;
E_nom=2.94*10^7;L_nom=100.0;Fy_nom=1029;
[y_nom(:,1) y_nom(:,2) y_nom(:,3)]=canti_SDF([W T L_nom E_nom Fx Fy_nom rho]);

%% Initial Weight Normal Function (Inverse CDF 계산)
x0_dc=[L_nom L_nom*0.4 E_nom E_nom*0.6 Fy_nom Fy_nom*0.4];

lbx=[50.0 1.0*10^7 700.0];ubx=[150.0 5.0*10^7 1400.0];
[Finv0_fac Finv0_fac_bnd] = kde_correct_ver02_factorized(pd_in_opt); % Factorized Distribution
[Finv0 Finv0_bnd] = kde_correct_ver02(x0_dc,pd_in_opt); % Factorized Distribution x Initial Weight Normal Function
close all
for k=1:3
    figure
    hold on
    plot(linspace(0.00,1.00,100)',feval(Finv0{k,:},linspace(0.00,1.00,100)'),'r--','LineWidth',2.0);
    box on
    grid on
    ylim([lbx(k),ubx(k)])
    xlim([0 1])
end

%% 시험데이터 생성
rho=0.284;W=5.0;T=5.0;Fx=500.0;N=10^6;
E=random('Weibull',2.91*10^7,5,[N 1])+7.36*10^6; % wblrnd(2.91*10^7,7.36*10^6) 
L=random('Normal',100.0,10.0,[N 1]);
Fy=random('Lognormal',6.933,0.082,[N 1]);
xx_set_exp=[W*ones(N,1) T*ones(N,1) L E Fx*ones(N,1) Fy rho*ones(N,1)];
[y_true(:,1) y_true(:,2) y_true(:,3)]=canti_SDF(xx_set_exp); % hist(y_exp(:,3))
r_true=(y_true-y_nom)./y_true;

%% UQ 결과 (Factorized, Factorized + Weight)
nMCS=10^6;
for k=1:3
    MCS_bnd(1,k)=max(Finv0_bnd(1,k),Finv0_fac_bnd(1,k));
    MCS_bnd(2,k)=min(Finv0_bnd(2,k),Finv0_fac_bnd(2,k));
    rand_nMCS=(MCS_bnd(2,k)-MCS_bnd(1,k))*rand(nMCS,1)+MCS_bnd(1,k);
    z0_MCS(:,k) = feval(Finv0{k,:},rand_nMCS);

    z_low(k,:) = feval(Finv0{k,:},0.025);
    z_up(k,:) = feval(Finv0{k,:},0.975);
    
    z0_MCS_fac(:,k) = feval(Finv0_fac{k,:},rand_nMCS);
    
    z_low_fac(k,:) = feval(Finv0_fac{k,:},0.025);
    z_up_fac(k,:) = feval(Finv0_fac{k,:},0.975);
end
S_a_nom_fac=z_up_fac-z_low_fac;
S_a_nom=z_up-z_low;

% UQ 결과 - Output Sample Set (Factorized Distribution)
L0=z0_MCS_fac(:,1);E0=z0_MCS_fac(:,2);Fy0=z0_MCS_fac(:,3);
xx_set0=[W*ones(nMCS,1) T*ones(nMCS,1) L0 E0 Fx*ones(nMCS,1) Fy0 rho*ones(nMCS,1)];
[y0(:,1) y0(:,2) y0(:,3)]=canti_SDF(xx_set0);
r0=(y0-y_nom)./y0;

%% Optimization (Factorized => Factorized + Optimal Weight)
lb_dc=[L_nom*0.5 L_nom*0.01 E_nom*0.5 E_nom*0.01 Fy_nom*0.5 Fy_nom*0.01];
ub_dc=[L_nom*1.5 L_nom*4.00 E_nom*1.5 E_nom*4.00 Fy_nom*1.5 Fy_nom*4.00];
obj2=@(x) obj_correct_ver02(x,pd_in_opt,S_a_nom);
nonlcon2=@(x) nonlcon_correct_ver03(x,pd_in_opt,y_nom,y0,file_s_name,MCS_bnd);

opts = optimoptions('fmincon','Display','iter','Algorithm','sqp','SpecifyObjectiveGradient',false,'StepTolerance',1e-15,'OptimalityTolerance',1e-15,'ConstraintTolerance',1e-15,'MaxFunctionEvaluations',1000);
[xopt_correct,fval_correct,exitflag,output] = fmincon(obj2,x0_dc,[],[],[],[],lb_dc,ub_dc,nonlcon2,opts);

%% Check Results
[~,S_a_opt] = obj_correct_ver02(xopt_correct,pd_in_opt,S_a_nom);
V_opt = nonlcon_correct_ver03_new(xopt_correct,pd_in_opt,y_nom,y_true,file_s_name,MCS_bnd);
V_ini = nonlcon_correct_ver03_new(x0_dc,pd_in_opt,y_nom,y_true,file_s_name,MCS_bnd);
V_0 = nonlcon_correct_ver04_factorized(pd_in_opt,y_nom,y_true,file_s_name,Finv0_fac,MCS_bnd);

ratio_S=S_a_opt./S_a_nom_fac;
varia_V=V_opt-V_0;

save('xopt_correct_sqp_201116_trial8.mat','xopt_correct','fval_correct','ratio_S','V_opt','V_0') % x0_dc=[L_nom L_nom*0.4 E_nom E_nom*0.6 Fy_nom Fy_nom*0.4];

Finv0_opt = kde_correct_ver02(xopt_correct,pd_in_opt);
% UQ 결과
nMCS=10^6;
for k=1:3
    rand_nMCS=(MCS_bnd(2,k)-MCS_bnd(1,k))*rand(nMCS,1)+MCS_bnd(1,k);
    z0_MCS_opt(:,k) = feval(Finv0_opt{k,:},rand_nMCS);

    z_low_opt(k,:) = feval(Finv0_opt{k,:},0.025);
    z_up_opt(k,:) = feval(Finv0_opt{k,:},0.975);
end
S_a_nom_opt=z_up_opt-z_low_opt;

%% 최적화 결과 - Graph Plot
close all
graph_name2_x={'\it L (in)','\it E (psi)','\it F_{y}(lbf)'};
sub_n=[2 1 3];
Fig2=figure;
set(Fig2,'pos',[200 541 2260 543]);
for k=1:3
    ax=subplot(1,3,sub_n(k));
    cdf_output=linspace(0.00,1.00,10000)';
    cdf_input_fac(:,k)=feval(Finv0_fac{k,:},linspace(0.00,1.00,10000)');
    s1=plot(cdf_input_fac(:,k),cdf_output,'-','LineWidth',2.0,'Color',[0.4660 0.6740 0.1880]);
    hold on
    s2=scatter([z_low_fac(k,1) z_up_fac(k,1)],[0.025 0.975],50,'filled','LineWidth',1.5,'MarkerEdgeColor',[0.4660 0.6740 0.1880],'MarkerFaceColor',[1 1 1]);
    
    hold on
    cdf_input_opt(:,k)=feval(Finv0_opt{k,:},linspace(0.00,1.00,10000)');
    s3=plot(cdf_input_opt(:,k),cdf_output,'-','LineWidth',2.0,'Color',[0.8500 0.3250 0.0980]);
    hold on
    s4=scatter([z_low_opt(k,1) z_up_opt(k,1)],[0.025 0.975],50,'filled','LineWidth',1.5,'MarkerEdgeColor',[0.8500 0.3250 0.0980],'MarkerFaceColor',[1 1 1]);
    
    ax = gca;
    set(gca,'fontsize', 14);
    
    legend([s1 s3 s2 s4],{'CDF (Factorized)','CDF (Corrected)','S_{0}','S_{\alpha}'},'fontsize',14,'location','best')
    box on
    grid on
    xlim([lbx(k),ubx(k)])
    if k==2
        xlim([1*10^7,5*10^7])
    end
    ylim([0 1])
    xlabel(graph_name2_x{k},'fontsize',15)
    ylabel('CDF','fontsize',15)
    set(gca,'fontname','times')  % Set it to times
end
file_save_name2=[pwd,'\Results_section_3_3','\','CDF_of_corrected_input_variables'];
saveas(Fig2,[file_save_name2,'.png'])
close all

pdf_corrected = pdf_correct(xopt_correct,pd_in_opt);
graph_name2_x={'\it L (in)','\it E (psi)','\it F_{y}(lbf)'};
sub_n=[2 1 3];
Fig2=figure;
set(Fig2,'pos',[200 541 2260 543]);
ubx=[150.0 6.0*10^7 1400.0];
for k=1:3
    ax=subplot(1,3,sub_n(k));
    s1=plot(q_int(:,k),pdf_true(:,k),'k--','LineWidth',2.0);
    hold on
    pdf_factorized(:,k)=pdf(pd_in_opt{k,:},q_int(:,k));
    s2=plot(q_int(:,k),pdf_factorized(:,k),'-','LineWidth',2.0,'Color',[0.4660 0.6740 0.1880]);
    hold on
    s3=plot(q_int(:,k),pdf_corrected(:,k),'-','LineWidth',2.0,'Color',[0.8500 0.3250 0.0980]);
    
    ax = gca;
    set(gca,'fontsize', 14);
    
    legend([s1 s2 s3],{'Target','Factorized','Corrected'},'fontsize',14,'location','northeast')
    box on
    grid on
    xlim([lbx(k),ubx(k)])
    xlabel(graph_name2_x{k},'fontsize',15)
    ylabel('CDF','fontsize',15)
    set(gca,'fontname','times')  % Set it to times
end
file_save_name2=[pwd,'\Results_section_3_3','\','PDF_of_corrected_input_variables'];
saveas(Fig2,[file_save_name2,'.png'])
close all

ci_resp(:,1)=[4*10^3 0 2]';
ci_resp(:,2)=[12*10^3 0.8 12]';
pd_y_corrected = pdf_output_response_correct_ver02(xopt_correct,pd_in_opt,file_s_name,MCS_bnd);
pd_y_true = pdf_output_response_true(file_s_name);
pd_y_factorized = pdf_output_response_factorized_ver02(pd_in_opt,file_s_name,Finv0_fac,MCS_bnd);
graph_name2_x={'\it \tau (psi)','\it \delta (in)','\it \omega (rad/s)'};
sub_n=[1 2 3];
Fig2=figure;
set(Fig2,'pos',[200 541 2260 543]);
for k=1:3
    y_resp=linspace(ci_resp(k,1),ci_resp(k,2),100);
    Y_resp(:,k)=y_resp';
    pdf_resp_true(:,k)=pdf(pd_y_true{k,:},y_resp);
    pdf_resp_factorized(:,k)=pdf(pd_y_factorized{k,:},y_resp);
    pdf_resp_corrected(:,k)=pdf(pd_y_corrected{k,:},y_resp);

    ax=subplot(1,3,sub_n(k));
    s1=plot(y_resp,pdf_resp_true(:,k),'k--','LineWidth',2.0);
    hold on
    s2=plot(y_resp,pdf_resp_factorized(:,k),'-','LineWidth',2.0,'Color',[0.4660 0.6740 0.1880]);
    hold on
    s3=plot(y_resp,pdf_resp_corrected(:,k),'-','LineWidth',2.0,'Color',[0.8500 0.3250 0.0980]);
    
    ax = gca;
    set(gca,'fontsize', 14);
    if k==2
    legend([s1 s2 s3],{'Target','Propagated before distribution correction','Propagated after distribution correction'},'fontsize',14,'location','best')
    end
    box on
    grid on
    xlim([ci_resp(k,1),ci_resp(k,2)])
    xlabel(graph_name2_x{k},'fontsize',15)
    ylabel('PDF','fontsize',15)
    set(gca,'fontname','times')  % Set it to times
end
file_save_name2=[pwd,'\Results_section_3_3','\','PDF_of_corrected_output_responses'];
saveas(Fig2,[file_save_name2,'.png'])
close all
