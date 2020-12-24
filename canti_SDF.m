function [S D F] = canti_SDF(xx)

%% Input variables
% High sensitivity
W = xx(:,1); 
T = xx(:,2);
L = xx(:,3); % Unknown variable (to be estimated)

% Low sensitivity
E = xx(:,4); % Unknown variable (to be estimated)
Fx = xx(:,5);
Fy = xx(:,6); % Unknown variable (to be estimated)
rho = xx(:,7);

% Parameter
ck=1.875; % k=1¿œ ∂ß

%% System responses
S=6*L./(W.*T).*(Fx./T+Fy./W);
D=3*L.^3./(E.*W.*T).*sqrt((Fy./T.^2).^2+(Fx./W.^2).^2);
F=ck^2*T./(L.^2).*sqrt(E./(12*rho));


end
