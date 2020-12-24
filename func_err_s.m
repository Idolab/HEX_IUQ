function [err_s] = func_err_s(xx,p,N,y_exp)

x(:,[3 4 6])=repmat(xx,[N 1]);
x(:,[1 2 5 7])=repmat(p,[N 1]);

[y_sim(:,1) y_sim(:,2) y_sim(:,3)]=canti_SDF(x);
err_s=y_exp-y_sim;


