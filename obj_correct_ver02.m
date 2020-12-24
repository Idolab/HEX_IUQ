function [S_a_total S_a_opt] = obj_correct_ver02(x_dc,pd_in_opt,S_a_nom)

Finv = kde_correct_ver02(x_dc,pd_in_opt);
for k=1:3
    z_low(k,:) = feval(Finv{k,:},0.025);
    z_up(k,:) = feval(Finv{k,:},0.975);
end
S_a_total=mean(abs(z_up-z_low)./S_a_nom);
S_a_opt=z_up-z_low;
save('Finv_temp.mat','Finv','S_a_total')




