function [Finv Finv_bnd] = kde_correct_ver02(x_dc,pd_in_opt)

% Inverse CDF °è»ê
lbx=[50.0 1.0*10^7 700.0];
ubx=[150.0 5.0*10^7 1400.0];
for k=1:3
    x(:,k) = linspace(lbx(k),ubx(k),50);
    pdf_unnorm=@(x) normpdf(x,x_dc(2*k-1),x_dc(2*k)).*pdf(pd_in_opt{k,:},x);
    F_norm(k,:)=integral(pdf_unnorm,lbx(k),ubx(k));
    pdf_norm=@(x) normpdf(x,x_dc(2*k-1),x_dc(2*k)).*pdf(pd_in_opt{k,:},x)/F_norm(k,:);
    for j=1:length(x)
        F(j,k) = integral(pdf_norm,-Inf,x(j,k));
    end
    
     Finv{k,:} = @(u) func_interp(u,F(:,k),x(:,k));
     Finv_bnd(:,k)=F([1 end],k);
end
