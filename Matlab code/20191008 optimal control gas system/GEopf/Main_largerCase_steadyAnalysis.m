%% this is a larger case study, with steady state analysis only
clear
clc
tic
% load case
[mpc, gtd] = case24GEv3();
mpc0 = mpc; %back up
rts = Case24ReliabillityDatav2();
% set MCS parameters
simulationTimes = 500;
missionTime = 12;
numberOfPeriods = 144;
nGl = size(mpc.Gline,1);
%set transient parameters
[result0, success] = GErunopf(mpc);
[para] = initializeParameters2();
%
[systemReliabilityForOneSimulation, nodalElectricityReliabilityForOneSimulation, ...
    nodalGasReliabilityForOneSimulation] = deal(cell(simulationTimes,1));
for i = 1:simulationTimes
    %% MCS
    [Info_components,failureComponents] = MCSformingScenarioV5(rts,missionTime);
    numberOfStates = size(Info_components,1);
    duration = Info_components(:,3) - Info_components(:,2);
    steadyResult = cell(numberOfStates,1);
    for j = 1:numberOfStates
        %% run IOPF 
        % the units used in IOPF are bar and Mm3/day
        newmpc = mpcUpdate(Info_components(j,:),duration(j),mpc0,rts);% update CDF by the way
        % test-----------------------------
        newmpc = testmpc_case24GEv3();
        % test---------------------------------
        opfOptions = mpoption('out.all',0);
        [steadyResult{j},flag] = GErunopf(newmpc,opfOptions);
        if flag ~= 1 % didn't converge£¬optimize without constraints
            [steadyResult{j}] = opfWithoutConstraint(steadyResult{j});
        end
        [ genCost, gasCost, LCeCost, LCgCost] = costDecompositon( steadyResult{j}, newmpc );
    end
    % get the reliability indices in one simualtion
    [systemReliabilityForOneSimulation{i}, nodalElectricityReliabilityForOneSimulation{i},...
        nodalGasReliabilityForOneSimulation{i}] = getReliabilityIndice_steady(steadyResult,Info_components(:,[2 3]),numberOfStates,numberOfPeriods,missionTime);
end
sumOfSystemReliability = zeros(size(systemReliabilityForOneSimulation{1}));
sumOfNodalElectricityReliability = zeros(size(nodalElectricityReliabilityForOneSimulation{1}));
sumOfNodalGasReliability = zeros(size(nodalGasReliabilityForOneSimulation{1}));
% std format: [electricityLOLP,electricityEENS,gasLOLP,gasEENS,totalCost];
stdSystemReliability = zeros(numberOfPeriods,11,simulationTimes);
allSystemReliability = zeros(numberOfPeriods,11,simulationTimes);
for i = 1:simulationTimes
    sumOfSystemReliability = systemReliabilityForOneSimulation{i} + sumOfSystemReliability;
    sumOfNodalElectricityReliability = nodalElectricityReliabilityForOneSimulation{i} + sumOfNodalElectricityReliability;
    sumOfNodalGasReliability = nodalGasReliabilityForOneSimulation{i} + sumOfNodalGasReliability;
    % std
    allSystemReliability(:,:,i) = systemReliabilityForOneSimulation{i}(:,:); % fit two demenssion into three demenssion    
    stdSystemReliability(:,:,i) = std(allSystemReliability(:,:,1:i),0,3)./mean(allSystemReliability(:,:,1:i),3)/sqrt(i);  
end
% reliability indices
systemReliability = sumOfSystemReliability / simulationTimes;
nodalElectricityReliability = sumOfNodalElectricityReliability / simulationTimes;
nodalGasReliability = sumOfNodalGasReliability / simulationTimes;
save    
elapseTime = toc  