%% Section 3.2 & 3.3
clc
close all
clear

for case_num=1:4
    for trial=1:10
        folder_path=pwd;
        cd(folder_path)
        folder_name=[pwd,'\Results_section_4_2'];
        mkdir(folder_name)
        
        folder_name_trial=[folder_name,'\',sprintf('results_trial%d',trial)];
        folder_name_figure=[folder_name_trial,'\graph'];
        mkdir(folder_name_trial)
        mkdir(folder_name_figure)
        for case_exp=1:18
            err_cond=[0 0;0.01 0;0 0.01;0.01 0.01];
            n=1.0; % Index of validation metric
            graph_name2_x={'L (in)','E (psi)','F_{y}(lbf)'};
            
            rho=0.284;
            W=5.0;
            T=5.0;
            Fx=500.0;
            p=[W T Fx rho];
            
            N_set=[25 50 100 200 300 400 500 1000 2000 3000 4000 5000 10000 20000 30000 40000 50000 100000]';
            N=N_set(case_exp);
            
            % Experiment data
            E=random('Weibull',2.91*10^7,5,[N 1])+7.36*10^6; % wblrnd(2.91*10^7,7.36*10^6)
            L=random('Normal',100.0,10.0,[N 1]);
            Fy=random('Lognormal',6.933,0.082,[N 1]);
            dist_type={'Normal','Weibull','Lognormal'};
            custompdf = @(x,a,b,c) (x>c).*(b/a).*(((x-c)/a).^(b-1)).*exp(-((x-c)/a).^b); % Weibull three-parameter
            p_stat=[100.0,10.0;2.91*10^7,5;6.933,0.082];
            
            xx_set_exp=[W*ones(N,1) T*ones(N,1) L E Fx*ones(N,1) Fy rho*ones(N,1)];
            [y_exp(:,1) y_exp(:,2) y_exp(:,3)]=canti_SDF(xx_set_exp); % hist(y_exp(:,3))
            
            CoV=err_cond(case_num,2);
            err_m=random('Uniform',-err_cond(case_num,1),err_cond(case_num,1),[N 1]).*y_exp;
            err_sim=random('Normal',zeros(N,size(y_exp,2)),y_exp*CoV);
            y_exp=(y_exp+err_sim)+err_m;
            
%% Step1: Error-lumped inverse UQ  
            % Optimization
            tic
            for ei=1:N
                y_num=[ei];
                obj=@(x) func_e_UF_ver02(x,p,y_exp,y_num);
                
                lbx=[50.0 1.0*10^7 700.0];
                ubx=[150.0 6.0*10^7 1400.0];
                x0=[100.0 2.94*10^7 1029.0];
                
                lbx_n=zeros(1,size(lbx,2));
                ubx_n=ones(1,size(ubx,2));
                x0_n=(x0-lbx)./(ubx-lbx);
                
                obj_n = @(x) obj_norm(x,obj,lbx,ubx);
                
                % Optimization1
                opts = optimoptions('fmincon','Display','off','Algorithm','sqp','SpecifyObjectiveGradient',false,'StepTolerance',1e-15,'OptimalityTolerance',1e-15,'ConstraintTolerance',1e-15,'MaxFunctionEvaluations',1000);
                [xopt_n,fval_sqp,exitflag,output] = fmincon(obj_n,x0_n,[],[],[],[],lbx_n,ubx_n,[],opts);
                xopt_sqp(ei,:)=xopt_n.*(ubx-lbx)+lbx;
                
                if mod(ei,N/5)==0
                    disp(sprintf('Experiment data %d Test: %d-th complete',N,ei))
                    toc
                end
            end
            
            %%%%%%%%% Validation metric %%%%%%%%%
            % PDF of estimated input distribution
            mean_xopt=mean(xopt_sqp,1);
            sigma_xopt=std(xopt_sqp,1);
            bw_opt=(4/(3*N))^(1/5)*sigma_xopt;
            for kk=1:3
                pd_in_opt{kk,:} = fitdist(xopt_sqp(:,kk),'Kernel','kernel','normal','support','unbounded','Width',bw_opt(kk));
             end
            
            ci_int(:,1)=lbx';
            ci_int(:,2)=ubx';
            for kk=1:3
                q_int(:,kk)=linspace(ci_int(kk,1),ci_int(kk,2),100);
                pdf_est(:,kk)=pdf(pd_in_opt{kk,:},q_int(:,kk));
                if kk==2
                    pdf_true(:,kk)=custompdf(q_int(:,kk),p_stat(kk,1),p_stat(kk,2),7.36*10^6);
                else
                    pdf_true(:,kk)=pdf(dist_type{kk},q_int(:,kk),p_stat(kk,1),p_stat(kk,2));
                end
                pdf_int(:,kk)=min(pdf_est(:,kk),pdf_true(:,kk));
                Af(kk,:)=trapz(q_int(:,kk),pdf_int(:,kk),1);
            end
            V(case_exp,:)=Af.^n;
            
            file_save_name=[folder_name_trial,'\',sprintf('UQ_results_case%d_exp_set%d',case_num,case_exp)];
            save(file_save_name)
            
            file_save_name1=[folder_name_figure,'\',sprintf('UQ_results_case%d_exp_set%d',case_num,case_exp)];
            Fig=figure;
            set(Fig,'pos',[500 341 860 443]);
            a1=semilogx(N_set([1:case_exp]),V(:,1),'k:s','LineWidth',2.0,'MarkerSize',7);
            hold on
            a2=semilogx(N_set([1:case_exp]),V(:,2),'b-.d','LineWidth',2.0,'MarkerSize',7);
            hold on
            a3=semilogx(N_set([1:case_exp]),V(:,3),'r--o','LineWidth',2.0,'MarkerSize',7);
            ylim([0.75 1.00])
            xlim([10 10^5])
            legend([a1 a2 a3],{'E','L','Fy'},'fontsize',11,'location','best')
            xlabel('No. of samples','fontsize',13)
            ylabel('V','fontsize',13)
            title(['Validation metric (',sprintf('Case%d',case_num),')'],'fontsize',15)
            grid on
            saveas(Fig,[file_save_name1,'.png'])
            close all
            
            for kk=1:3
                file_save_name2=[folder_name_figure,'\',sprintf('UQ_results_input_PDF_case%d_exp_set%d_input%d',case_num,case_exp,kk)];
                Fig2=figure;
                p1=plot(q_int(:,kk),pdf_true(:,kk),'k-.','LineWidth',2.0);
                hold on
                p2=plot(q_int(:,kk),pdf_est(:,kk),'b-','LineWidth',2.0);
                xlim(ci_int(kk,:))
                grid on
                legend([p1 p2],{'Target','Estimated'},'fontsize',11,'location','best')
                xlabel(graph_name2_x{kk},'fontsize',13)
                ylabel('PDF','fontsize',13)
                saveas(Fig2,[file_save_name2,'.png'])
                close all
            end
            
            clear y_exp
            
        end
        clearvars -except case_num
    end
end


%% Graph - 10 Trials
for case_num=1:4   % case_num=1;
    V_total=[];
    for trial=1:10
        folder_path=pwd;
        folder_name=[pwd,'\Results_section_4_2'];
        folder_name_trial=[folder_name,'\',sprintf('results_trial%d',trial)];
        folder_name_figure=[folder_name_trial,'\graph'];
        for case_exp=18    % case_exp=18;
            file_save_name=[folder_name_trial,'\',sprintf('UQ_results_case%d_exp_set%d',case_num,case_exp)];
            load(file_save_name);
            V_total(:,:,trial)=V;
            if trial==1
                pdf_est_total{case_num,:}=pdf_est;
            end
        end
        %         clearvars -except case_num V_total V_case V_min
    end
    V_case{case_num,:}=V_total;
    V_min{case_num,:}=min(min(V_total,[],3),[],2);
    V_max{case_num,:}=max(max(V_total,[],3),[],2);
    V_mean{case_num,:}=mean(V_total,3);
end

Fig=figure;
set(Fig,'pos',[100 141 1860 1043]);
sub_n=[3 4 1 2];
for case_num=1:4
    ax=subplot(2,2,sub_n(case_num));
    a1=semilogx(N_set,V_mean{case_num,:}(:,1),'-s','LineWidth',2.0,'MarkerSize',7,'Color',[0.8500 0.3250 0.0980],'MarkerFaceColor',[0.8500 0.3250 0.0980]);
    hold on
    a2=semilogx(N_set,V_mean{case_num,:}(:,2),'-d','LineWidth',2.0,'MarkerSize',7,'Color',[0 0.4470 0.7410],'MarkerFaceColor',[0 0.4470 0.7410]);
    hold on
    a3=semilogx(N_set,V_mean{case_num,:}(:,3),'-^','LineWidth',2.0,'MarkerSize',7,'Color',[0.4660 0.6740 0.1880],'MarkerFaceColor',[0.4660 0.6740 0.1880]);
    hold on
    m1=semilogx(N_set,V_min{case_num,:}(:,1),':','LineWidth',2.0,'MarkerSize',7,'Color',[0.3010 0.7450 0.9330]);
    hold on
    m2=semilogx(N_set,V_max{case_num,:}(:,1),':','LineWidth',2.0,'MarkerSize',7,'Color',[0.3010 0.7450 0.9330]);
    
    ax = gca;
    set(gca,'fontsize', 14);
    
    ylim([0.75 1.00])
    % ylim([0.50 1.00])
    xlim([10 10^5])
    legend([a2 a1 a3 m1 m2],{'\it E','\it L','\it Fy','\it Min.','\it Max.'},'fontsize',14,'location','best')
    xlabel('No. of samples','fontsize',15)
    ylabel('V','fontsize',15)
    title(['Validation Metric (',sprintf('Case%d',case_num),')'],'fontsize',16)
    set(gca,'fontname','times')  % Set it to times
    grid on
end
file_save_name1=[folder_name,'\','UQ_results_total_case'];
saveas(Fig,[file_save_name1,'.png'])
close all

graph_name2_x={'\it L (in)','\it E (psi)','\it F_{y}(lbf)'};
graph_color={[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.4660 0.6740 0.1880],[0.4940 0.1840 0.5560]};
sub_n=[2 1 3];
Fig2=figure;
set(Fig2,'pos',[200 541 2260 543]);
for kk=1:3
    ax=subplot(1,3,sub_n(kk));
    p1=plot(q_int(:,kk),pdf_true(:,kk),'k--','LineWidth',2.0);
    for case_num=1:4
        hold on
        plot(q_int(:,kk),pdf_est_total{case_num}(:,kk),'-','LineWidth',1.5,'Color',graph_color{case_num});
    end
    
    ax = gca;
    set(gca,'fontsize', 14);
    
    xlim(ci_int(kk,:))
    grid on
    legend({'Target','Case 1','Case 2','Case 3','Case 4'},'fontsize',14,'location','best')
    xlabel(graph_name2_x{kk},'fontsize',15)
    ylabel('PDF','fontsize',15)
    set(gca,'fontname','times')  % Set it to times
end
file_save_name2=[folder_name,'\','UQ_results_input_PDF_trial1_total_case'];
saveas(Fig2,[file_save_name2,'.png'])
close all





