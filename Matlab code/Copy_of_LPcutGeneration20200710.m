function [alpha,constant] = LPcutGeneration(ei0,gi0,x_hat,A,b,r)
[ncons,nvar] = size(A);
nvar_ei = size(ei0,1)*size(ei0,2);
nvar_gi = size(gi0,1)*size(gi0,2);
nvar_x = nvar - nvar_ei - nvar_gi;
alpha_ei = sdpvar(nvar_ei,1);
alpha_gi = sdpvar(nvar_gi,1);
alpha_x = sdpvar(nvar_x,1);
u = sdpvar(ncons,1); 
v = sdpvar(ncons,1);
sdpvar u0 v0 
sdpvar beta01 % beta is occupied as a system function name
 
alpha = [alpha_ei;alpha_gi;alpha_x];
opts = sdpsettings('verbose',2,'solver','gurobi','usex0',0,'debug',1,'gurobi.InfUnbdInfo',1,...
    'savesolverinput',1,'savesolveroutput',1);
obj = alpha'*[ei0;gi0;x_hat] - beta01;
cons = [
    alpha' == u'*A - u0*generateUnitVector(r,nvar)';
    alpha' == v'*A + v0*generateUnitVector(r,nvar)';
    beta01 <= u'*b;
    beta01 <= v'*b + v0;
%     sum(u) + sum(v) + u0 + v0 == 1;
    u >= 0;
    v >= 0;
%     [u0,v0] >= 0;
    ];
diag = optimize(cons,obj,opts);
[alpha_ei,alpha_gi] = deal(value(alpha_ei), value(alpha_gi));
alpha_x = value(alpha_x);
alpha = [alpha_ei;alpha_gi;alpha_x];
beta01 = value(beta01);
constant = alpha_ei'*ei0 + alpha_gi'*gi0 - beta01;
end
function unitVector = generateUnitVector(j,size)
unitVector = zeros(size,1);
unitVector(j) = 1;
end