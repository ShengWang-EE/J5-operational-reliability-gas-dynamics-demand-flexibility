function [nextTimeState] = transientGasPowerFlowForOneTimeStep(initialCondition, boundaryCondition, mpc, gtd,nx,dx,dt,para)
% reform initial condition to x
nGl = size(mpc.Gline,1);
nVar = sum(nx+1);
x0_P = zeros(nVar,1);x0_Q = zeros(nVar,1);
cumnx = cumsum(nx+1);

for i = 1:nGl
    if i==1
        x0_P(1:nx(i)+1) = initialCondition.P{i};
        x0_Q(1:nx(i)+1) = initialCondition.Q{i};
    else
        x0_P(cumnx(i-1)+1:cumnx(i)) = initialCondition.P{i};
        x0_Q(cumnx(i-1)+1:cumnx(i)) = initialCondition.Q{i};
    end
end
% ฑ๊็Ûปฏ
x0_P = x0_P / 1e5;
%
x0 = [x0_P;x0_Q];


% options = optimoptions('fsolve','OptimalityTolerance',1e-6,'FunctionTolerance',1e-6,'StepTolerance',1e-10);
% options = optimset('TolFun',1e-1,'TolX',1e-1);
options = optimoptions('fsolve','Display','iter');
fun = @(x)fcn_transientGasFlowModelForOneTimeStep2(x,initialCondition,boundaryCondition,mpc,gtd,nx,dx,dt,para);
% [x,~,~,~,~] = fsolve(fun,x0);
[x,~,~,~,~] = fsolve(fun,x0,options);
% [x,~,~,~,~] = fsolve(fun,x0);
[nextTimeState] = formatXtoPQcell(x,nGl,nx);
for i = 1:size(nextTimeState.P)
    nextTimeState.P{i} = nextTimeState.P{i} * 1e5;
end
end
