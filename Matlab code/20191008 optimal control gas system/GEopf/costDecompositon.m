function [ genCost, gasCost, LCeCost, LCgCost] = costDecompositon( result, mpc )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% 提取变量
PGs = result.Gsou(:,5);
PgAndLCe = result.gen(:,2);
indexGen = 1:mpc.originalGenNumber;
indexLCe = mpc.originalGenNumber+1 : size(mpc.gen);
LCg = result.Gbus(mpc.Gbus(:,3)>0,10);

% calculate gas CDF
gasCDF = calculateGasCDF(mpc);

%%
gasCost = PGs'*mpc.Gcost;
LCgCost = sum(gasCDF .* LCg);

genAndLCeCost = totcost(mpc.gencost, PgAndLCe);

LCeCost = sum(genAndLCeCost(indexLCe));
genCost = sum(genAndLCeCost(indexGen));


end

