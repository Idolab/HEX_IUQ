function pdf_corrected = pdf_correct(x_dc,pd_in_opt)

% Inverse CDF °è»ê
lbx=[50.0 1.0*10^7 700.0];
% ubx=[150.0 5.0*10^7 1400.0];
ubx=[150.0 6.0*10^7 1400.0];
for k=1:3
    x(:,k) = linspace(lbx(k),ubx(k),100);
    pdf_unnorm=@(x) normpdf(x,x_dc(2*k-1),x_dc(2*k)).*pdf(pd_in_opt{k,:},x);
    F_norm(k,:)=integral(pdf_unnorm,lbx(k),ubx(k));
    pdf_norm=@(x) normpdf(x,x_dc(2*k-1),x_dc(2*k)).*pdf(pd_in_opt{k,:},x)/F_norm(k,:);

    pdf_corrected(:,k)=pdf_norm(x(:,k));
end
