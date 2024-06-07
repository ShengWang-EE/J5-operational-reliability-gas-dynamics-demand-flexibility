function gasReliability  = calculateGasReliabilityFromTransient(gasLC,nt,dt,numberOfPeriods,missionTime,energyType)

NK = size(gasLC,1);
nGd = size(gasLC,2);
switch energyType
    case 'electricity'
        gasLCthreshold = 1; % electricityLC threshold actually
    case 'gas'
        gasLCthreshold = 0.02;
end
% electricityLCthreshold = 0.2;

% system reliability column: [from time, to time, gas LOLP, gas EENS]
Info_gasReliability.system = zeros(NK,4);
gasReliability.system = zeros(numberOfPeriods,4);
% nodal reliability column: [each gas load bus LOLP, each gas load bus EENS]
Info_gasReliability.nodal = zeros(NK,2*nGd);
gasReliability.nodal = zeros(numberOfPeriods,2*nGd);

% 
cumnt = cumsum([0;nt]);

for k = 1:NK
    j = max(find(cumnt<k));
    if k>1
        Info_gasReliability.system(k,1) = Info_gasReliability.system(k-1,1) + dt(j);
    end
    Info_gasReliability.system(k,2) = Info_gasReliability.system(k,1) + dt(j);
end
Info_gasReliability.system(:,1:2) = Info_gasReliability.system(:,1:2) / 3600;
Info_gasReliability.system(end,2) = missionTime;
% 
for k = 1:NK
    nodalGasLC = gasLC(k,:);
    systemGasLC = sum(nodalGasLC);
    Info_gasReliability.system(k,3) = (systemGasLC > gasLCthreshold);
    Info_gasReliability.system(k,4) = systemGasLC;
    Info_gasReliability.nodal(k,1:nGd) = (nodalGasLC > gasLCthreshold);
    Info_gasReliability.nodal(k,nGd+1:2*nGd) = nodalGasLC;
end
%% 
systemGasLOLP = adjustTimeScaleofReliability(Info_gasReliability.system(:,[1 2 3]),numberOfPeriods,missionTime);
systemGasEENS = adjustTimeScaleofReliability(Info_gasReliability.system(:,[1 2 4]),numberOfPeriods,missionTime);
gasReliability.system(:,1:3) = systemGasLOLP;
gasReliability.system(:,4) = systemGasEENS(:,3);
% nodal reliability, does not contain the from time and to time information
fromTimeToTime = Info_gasReliability.system(:,[1 2]);
for i = 1:nGd
    nodalGasLOLP = adjustTimeScaleofReliability([fromTimeToTime, Info_gasReliability.nodal(:,i)],numberOfPeriods,missionTime);
    nodalGasEENS = adjustTimeScaleofReliability([fromTimeToTime, Info_gasReliability.nodal(:,nGd+i)],numberOfPeriods,missionTime);
    
    gasReliability.nodal(:,i) = nodalGasLOLP(:,3);
    gasReliability.nodal(:,nGd+i) = nodalGasEENS(:,3);
end

end

