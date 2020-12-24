%% Likelihood Function Metric
function [ff]=obj_norm(x_n,obj,lbx,ubx)

x=(ubx-lbx).*x_n+lbx;
ff=obj(x);