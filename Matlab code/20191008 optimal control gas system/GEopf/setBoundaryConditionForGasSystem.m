function boundaryCondition = setBoundaryConditionForGasSystem(result,mpc,nx)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
nGl = size(mpc.Gline,1); nGb = size(mpc.Gbus,1);
boundaryCondition.P = zeros(nGb,1); boundaryCondition.Q = zeros(nGl,2);

for i = 1:nGl
    boundaryCondition.Q(i,1) = result.Gline(i,6); % Mm3/day
    boundaryCondition.Q(i,2) = boundaryCondition.Q(i,1);
end             
for i = 1:nGb
    boundaryCondition.P(i) = result.Gbus(i,7) * 100000; % Pa
end
boundaryCondition.gasSource = result.Gsou;
end

