function [e_UF] = func_e_UF_ver02(xx,p,y_exp,y_num)

y_len=length(y_num);
x(1:y_len,[3 4 6])=repmat(xx,[y_len 1]);
x(1:y_len,[1 2 5 7])=repmat(p,[y_len 1]);

[y_sim(:,1) y_sim(:,2) y_sim(:,3)]=canti_SDF(x);
err_s=y_exp(y_num,:)-(y_sim);
e_UF=sqrt(sum(((err_s)./y_exp(y_num,:)).^2,2));
