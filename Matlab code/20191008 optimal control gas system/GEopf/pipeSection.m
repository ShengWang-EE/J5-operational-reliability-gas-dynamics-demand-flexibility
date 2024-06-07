function [P_T,Q_T] = pipeSection(boundaryCondition,dt,para)
%UNTITLED a single pipe between two nodes
%   Detailed explanation goes here

%% parameters
[~,~,~,~,~,~,~,~,~,~,~,~,~,nx] = unpackPara(para);
%% test input
% [boundaryCondition] = testBoundaryCondition();
%% solve the euqation
% x0 = [boundaryCondition.P(1:nx-1),boundaryCondition.Q(2:nx)];
x0 = [boundaryCondition.P(2:nx),boundaryCondition.Q(2:nx)];
fun = @(x)pipeSectionForOneTimeStep(x,boundaryCondition,dt,para);
[x,~,~,~,~] = fsolve(fun,x0);

P_T = [x(1:nx-1),boundaryCondition.P_T];
% P_T = [boundaryCondition.P_T, x(1:nx-1)];
Q_T = [boundaryCondition.Q_T,x(nx:2*(nx-1))];
P = [boundaryCondition.P;P_T];
Q = [boundaryCondition.Q;Q_T];


end

