function [ newmpc ] = mpcUpdateBinary(componentState,mpc,nGen)
%UNTITLED modified from original mpcFailureGeneratingWithoutBranch.m in
%MPCE paper
%   Detailed explanation goes here
newmpc = mpc;
%% update gencost
newmpc.gencost(mpc.originalGenNumber+1:end,6) = 10000; % assumes the same CDF for all bus
%%
genState = componentState(1:nGen); GsouState = componentState(nGen+1:end);
failureGen = genState==0; failureGsou = GsouState==0;
newmpc.gen(failureGen,[2 3 4 5 9 10]) = 0;%让gen失效
% if GsouState(12) ~= 0 % 最大气源故障不考虑
    newmpc.Gsou(failureGsou,[2 3 4]) = 0;
% end

end

