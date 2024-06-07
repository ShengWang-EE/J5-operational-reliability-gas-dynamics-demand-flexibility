function [alpha,constant] = LPcutGeneration(x_hat,A,b,r)
[ncons,nvar] = size(A);

alpha = sdpvar(nvar,1);
u = sdpvar(ncons,1); 
v = sdpvar(ncons,1);
sdpvar u0 v0 
sdpvar beta01 
% beta is occupied as a system function name
 
opts = sdpsettings('verbose',2,'solver','gurobi','usex0',0,'debug',1,'gurobi.InfUnbdInfo',1,...
    'savesolverinput',1,'savesolveroutput',1);
obj = alpha'*x_hat - beta01;
cons = [
    alpha' == u'*A - u0*generateUnitVector(r,nvar)';
    alpha' == v'*A + v0*generateUnitVector(r,nvar)';
    beta01 <= u'*b;
    beta01 <= v'*b + v0;
    sum(u) + sum(v) + u0 + v0 == 1;
    u >= 0;
    v >= 0;
    [u0,v0] >= 0;
    ];
diag = optimize(cons,obj,opts);

alpha = value(alpha);
beta01 = value(beta01);
constant = - beta01;
end
function unitVector = generateUnitVector(j,size)
unitVector = zeros(size,1);
unitVector(j) = 1;
end