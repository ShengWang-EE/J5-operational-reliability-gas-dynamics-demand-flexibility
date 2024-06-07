function [ newmpc ] = mpcUpdate(Info_components,interuptionTime,mpc,rts)
%UNTITLED modified from original mpcFailureGeneratingWithoutBranch.m in
%MPCE paper
%   Detailed explanation goes here
newmpc = mpc;
numberOfGen = size(rts.gen,1);numberOfGsou = size(rts.Gsou,1);
%% update failure tfu (binary-state)
failureGen = Info_components(:,4:(3+numberOfGen))==0;
newmpc.gen(failureGen,[2 3 4 5 9 10]) = 0;%»√gen ß–ß
%% update failure gfu (multi-state)
gfuIndex = rts.gfu(:,1);
gfuUnitNumber = rts.gfu(:,4);
gfuActualUnitNumber = Info_components(:,3+rts.gfu(:,1));
for i = 1:size(gfuIndex,1)
    newmpc.gen(gfuIndex(i),[2 3 4 5 9 10]) = mpc.gen(gfuIndex(i),[2 3 4 5 9 10]) * gfuActualUnitNumber(i) / gfuUnitNumber(i);
end
%% update failure Gas source (multi-state)
gsUnitNumber = rts.Gsou(:,4);
gsActualUnitNumber = Info_components(:,(3+numberOfGen+1):3+numberOfGen+numberOfGsou);
for i = 1:size(mpc.Gsou,1)
    newmpc.Gsou(i,[2 3 4]) = mpc.Gsou(i,[2 3 4]) * gsActualUnitNumber(i) / gsUnitNumber(i);
end
%% update gencost
% customerType = 1; % industry
% CDF = calculateElectricityCDF(mpc,interuptionTime,customerType);
newmpc.gencost(mpc.originalGenNumber+1:end,6) = 10000; % assumes the same CDF for all bus
%% fix slack bus
% newmpc.gen(15,9) = 900; 
end

