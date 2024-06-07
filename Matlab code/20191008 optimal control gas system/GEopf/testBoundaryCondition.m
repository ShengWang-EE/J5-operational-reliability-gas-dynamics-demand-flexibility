function [boundaryCondition] = testBoundaryCondition(para)
%TESTBOUNDARYCONDITION Summary of this function goes here
%   Detailed explanation goes here
[~,~,~,A,~,~,B,~,~,F,D,rhon,dx,nx] = unpackPara(para);

p(1,1) = 60 * 100000;%Pa
Q(1,1:nx) = 6*1000000 / 86400; %m3/s at t=1 and x=1(t=0 and x=0);
dp2dx = -4 * B^2 * rhon^2 * Q(1,1) * abs(Q(1,1)) /( F^2 * A^2 * D );
p(1,2:nx) = sqrt(p(1,1)^2 + dp2dx * (1:nx-1) * dx);

boundaryCondition.P = p(1,:);
boundaryCondition.Q = Q(1,:)*86400/1000000; % at initial state, Q is a constant
boundaryCondition.P_T = p(1,end);
% boundaryCondition.P_T = p(1,1);
boundaryCondition.Q_T = 0;

boundaryCondition.flag = 1;
end

