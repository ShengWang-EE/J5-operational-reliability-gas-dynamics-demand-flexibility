function [systemReliability, nodalElectricityReliability,nodalGasReliability] = getReliabilityIndice_steady(steadyResult,fromTimeToTime,numberOfStates,numberOfPeriods,missionTime)
%% parameters
gasLCthreshold = 0.001;
electricityLCthreshold = 0.2;
nd = size(find(steadyResult{1}.bus(:,3)>0),1);
nGd = size(find(steadyResult{1}.Gbus(:,3)>0),1);
gasLoadIndex = find(steadyResult{1}.Gbus(:,3)>0);
% system reliability indices:
% [1 from time, 2 to time, 3 elelctricityLOLP, 4 electricityEENS, 5 gasLOLP,
% 6 gasEENS, 7 electricity gen cost, 8 gas purchasing cost, 9 LCe cost, 10 LCg cost,
% 11 total cost]; 11 rows
% nodal reliability indice
% nodel reliability:
% [1 from time, 2 to time, 3 electricity/gas LOLP, 4 electricity/gas EENS]
%
%% get reliability info
Info_systemReliability = zeros(numberOfStates,11);
Info_nodalElectricityReliability = zeros(numberOfStates,11,nd);
Info_nodalGasReliability = zeros(numberOfStates,11,nGd);

Info_systemReliability(:,1:2) = fromTimeToTime;
Info_nodalElectricityReliability(:,1:2,:) = repmat(fromTimeToTime,1,1,nd);
Info_nodalGasReliability(:,1:2,:) = repmat(fromTimeToTime,1,1,nGd);
for j = 1:numberOfStates
    nodalGasLC = steadyResult{j}.Gbus(gasLoadIndex,10);
    systemGasLC = sum(nodalGasLC);
    nodalElectricityLC = steadyResult{j}.gen(steadyResult{j}.originalGenNumber+1:end,2);
    systemELectricityLC = sum(nodalElectricityLC);
    %
    Info_systemReliability(j,3) = (systemELectricityLC > electricityLCthreshold);
    Info_systemReliability(j,5) = (systemGasLC > gasLCthreshold);
    Info_systemReliability(j,4) = systemELectricityLC;
    Info_systemReliability(j,6) = systemGasLC;
    Info_systemReliability(j,11) = steadyResult{j}.f;
    % cost decomposition (not doing this now)
    % nodal reliability
    Info_nodalElectricityReliability(j,3,:) = (nodalElectricityLC > electricityLCthreshold);
    Info_nodalGasReliability(j,3,:) = (nodalGasLC > gasLCthreshold);
    Info_nodalElectricityReliability(j,4,:) = nodalElectricityLC;
    Info_nodalGasReliability(j,4,:) = nodalGasLC;
end
%% 
% system reliabiliyt
systemElectricityLOLP = adjustTimeScaleofReliability(Info_systemReliability(:,[1 2 3]),numberOfPeriods,missionTime);
systemELectricityEENS = adjustTimeScaleofReliability(Info_systemReliability(:,[1 2 4]),numberOfPeriods,missionTime);
systemGasLOLP = adjustTimeScaleofReliability(Info_systemReliability(:,[1 2 5]),numberOfPeriods,missionTime);
systemGasEENS = adjustTimeScaleofReliability(Info_systemReliability(:,[1 2 6]),numberOfPeriods,missionTime);
systemGenCost = adjustTimeScaleofReliability(Info_systemReliability(:,[1 2 7]),numberOfPeriods,missionTime);
systemGasCost = adjustTimeScaleofReliability(Info_systemReliability(:,[1 2 8]),numberOfPeriods,missionTime);
systemElectricityLCCost = adjustTimeScaleofReliability(Info_systemReliability(:,[1 2 9]),numberOfPeriods,missionTime);
systemGasLCCost = adjustTimeScaleofReliability(Info_systemReliability(:,[1 2 10]),numberOfPeriods,missionTime);
systemTotalCost = adjustTimeScaleofReliability(Info_systemReliability(:,[1 2 11]),numberOfPeriods,missionTime);
systemReliability = [systemElectricityLOLP, systemELectricityEENS(:,3), systemGasLOLP(:,3), systemGasEENS(:,3),...
    systemGenCost(:,3), systemGasCost(:,3), systemElectricityLCCost(:,3), systemGasLCCost(:,3), systemTotalCost(:,3)];

% nodal reliability
nodalElectricityReliability = zeros(numberOfPeriods,4,nd);
nodalGasReliability = zeros(numberOfPeriods,4,nGd);
for i = 1:nd
    nodalElectricityLOLP = adjustTimeScaleofReliability(Info_nodalElectricityReliability(:,[1 2 3],i),numberOfPeriods,missionTime);
    nodalElectricityEENS = adjustTimeScaleofReliability(Info_nodalElectricityReliability(:,[1 2 4],i),numberOfPeriods,missionTime);
   
    nodalElectricityReliability(:,1:3,i) = nodalElectricityLOLP;
    nodalElectricityReliability(:,4,i) = nodalElectricityEENS(:,3);

end
for i = 1:nGd
    nodalGasLOLP = adjustTimeScaleofReliability(Info_nodalGasReliability(:,[1 2 5],i),numberOfPeriods,missionTime);
    nodalGasEENS = adjustTimeScaleofReliability(Info_nodalGasReliability(:,[1 2 6],i),numberOfPeriods,missionTime);
    
    nodalGasReliability(:,1:3,i) = nodalGasLOLP;
    nodalGasReliability(:,4,i) = nodalGasEENS(:,3);
end
end
%% 
