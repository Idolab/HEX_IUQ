function [z] = icdf_correct_1dim(u,x,pd)

% dim=size(x,2);
mu=x(1);
sig=x(2);

inp=pd.InputData.data;
bw=pd.BandWidth;

mu_new=(repmat(mu,[size(inp,1),1])*bw^2+inp*sig^2)./repmat(sig^2+bw^2,[size(inp,1),1]);
sig_new=sqrt(repmat((bw^2*sig^2)/(bw^2+sig^2),[size(inp,1),1]));

z=mean(icdf('Normal',u,mu_new,sig_new),1)';
% z=[mu_new sig_new];







