%%
% this is an illustrative example for the transient of one pipeline
% gas source------gas fired unit-----------electricity load

%% 
clear
clc
tic
% parameters
[para] = initializeParameters();
[C_q,R,T,A,Z_T,Z,B,g,alpha,F,D,rhon,dx,nx] = unpackPara(para);
missionTime = 24;
numberOfPeriods = 144;
simulationTimes = 100;
% initialize the reliability indices
reliabilityForOnesimulation = cell(simulationTimes,1);
allReliabilityIndices = zeros(numberOfPeriods,6,simulationTimes);
allReliabilityIndices(:,1,:) = [0:1:numberOfPeriods-1]' * missionTime/numberOfPeriods .* ones(numberOfPeriods,1,simulationTimes); 
allReliabilityIndices(:,2,:) = allReliabilityIndices(:,1,:) + missionTime/numberOfPeriods;
std_reliabilityindices = zeros(numberOfPeriods,4,simulationTimes);
% for observation
ob.countsForTransientSimulationTimes = 0;
%%
parfor i = 1:simulationTimes
    % MCS
    Info_gasSource = MCS_gasSource(missionTime);
    % test
    % Info_gasSource = [1,0,3,6;
    %                     2,3,12,4;
    %                     1,12,24,6];
    % predetermine the time step for each scenario
    numberOfStates = size(Info_gasSource,1);
    nt = zeros(1, numberOfStates); dt = nt;
    for j = 1:numberOfStates
        duration = Info_gasSource(j,3) - Info_gasSource(j,2);
        nt(j) = floor(duration*3600 / 600) + 1;%keep the time steps around 1h
        dt(j) = duration / nt(j) * 3600;
    end
    cumsumnt = cumsum(nt);
    % the steady initial state
    [boundaryCondition] = testBoundaryCondition(para);
    % start the transient state
    P = zeros(sum(nt),nx); Q = P;
    P(1,:) = boundaryCondition.P;
    Q(1,:) = boundaryCondition.Q;% just for scenario 1 (normal state). save the trouble

%
    counter = nt(1); % time step counter
    P(2:nt(1),:) = repmat(P(1,:),nt(1)-1,1);
    Q(2:nt(1),:) = repmat(Q(1,:),nt(1)-1,1); % there is no need to calculate the transient
    if numberOfStates == 1  % there is only one normal state
        disp('no transient')
    else       
        % ob
%         ob.countsForTransientSimulationTimes = ob.countsForTransientSimulationTimes + 1;
        %
        for j = 2:numberOfStates
                % convert MCS infomation into boundary conditions
                boundaryCondition.P_T = boundaryCondition.P_T;% the same
                boundaryCondition.Q_T = Info_gasSource(j,4);

            for k = 1:nt(j)
                % update initial condition according to the last time P,Q
                boundaryCondition.P = P(counter,:);
                boundaryCondition.Q = Q(counter,:);

                % calculate the newP,newQ
                [P_T,Q_T] = pipeSection(boundaryCondition,dt(j),para);
                P(counter+1,:) = P_T; Q(counter+1,:) = Q_T;
                counter = counter + 1;
            end
        end
    end
    
    % the P,Q for one simulation over whole period are obtained. calculate
    % the gas curtailment, electricity curtailment
    % remember to re-scele time
    % InfoReliability = [from time, to time, P at bus1, Q at bus2, gas LC, electricity LC]
    numberOfTimestep = sum(nt,2);
    InfoReliability = zeros(numberOfTimestep, 6);
    gasLoad = 5; electricityLoad = 200; GFUefficiency = 200; 
    
    pointer = 1; timePointer = 0;
    for j = 1:numberOfStates
        InfoReliability(pointer:pointer+nt(j),1) = timePointer + (0:1:nt(j)) * dt(j); % from time
        InfoReliability(pointer:pointer+nt(j)-1,2) = InfoReliability(pointer+1:pointer+nt(j),1); % to time
        pointer = pointer + nt(j);
        timePointer = timePointer + nt(j)*dt(j);        
    end
    InfoReliability(end,:) = [];
    InfoReliability(:,1:2) = InfoReliability(:,1:2) / 3600; % convert into hour
    InfoReliability(:,3) = P(:,1);
    InfoReliability(:,4) = Q(:,end);
    for k = 1:numberOfTimestep
        if Q(k,end) >= gasLoad % only curtail GFU
            gasLCforGFU = electricityLoad / GFUefficiency - (Q(k,end) - gasLoad);
            gasLoadLC = 0;
        else
            gasLCforGFU = electricityLoad / GFUefficiency;
            gasLoadLC = gasLoad - Q(k,end);
        end
        InfoReliability(k,5) = gasLoadLC;
        InfoReliability(k,6) = gasLCforGFU * GFUefficiency;
    end
    % obtain the reliability indices for one simulation
    % (need to adjust time scale)
    [reliabilityForOnesimulation{i}] = obtainReliabilityIndicesForOneSimulation(InfoReliability,numberOfPeriods,missionTime);


end
% results
for i = 1:simulationTimes
        % counts the reliability indices for bigger loop
    allReliabilityIndices(:,3:6,i) = reliabilityForOnesimulation{i}(:,3:6); % fit two demenssion into three demenssion    
    std_reliabilityindices(:,:,i) = std(allReliabilityIndices(:,3:6,1:i),0,3)./mean(allReliabilityIndices(:,3:6,1:i),3)/sqrt(i);  
end
reliabilityIndice = mean(allReliabilityIndices,3);





% 
ElapsedTime = toc
% 
% if numberOfStates ~= 1 
%     plot(Q(:,end))
% end
%% compare group
% t = 1:1:168;
% lamda = 1/960; mu = 1/40;
% y = mu/(lamda+mu) + lamda/(lamda+mu) * exp(-(lamda+mu)*t);
% plot(t,y)
%% function
function    [operationalReliability] = obtainReliabilityIndicesForOneSimulation(InfoReliability,numberOfPeriods,missionTime)
    % operationalReliability = [from time, to time, gasLOLP, gas EENS, electricity LOLP, electricity EENS]
    % LOLP 和EENS不同，LOLP是用checkpoint来，只检查某一时间点上是否失负荷，而EENS需要统计一个区间上的积累量
    % LOLP从1开始好了，而不是从0~1，因为t=0，LOLP肯定为0
    operationalReliability = zeros(numberOfPeriods,6);
    operationalReliability(:,1) = [0:1:numberOfPeriods-1]' * missionTime/numberOfPeriods; 
    operationalReliability(:,2) = operationalReliability(:,1) + missionTime/numberOfPeriods;
    for i = 1:numberOfPeriods
        % about LOLP
        indexj = max(find(InfoReliability(:,1)<=operationalReliability(i,2)));
        if InfoReliability(indexj,5)>0.005 % gas LC
            operationalReliability(i,3) = 1;
        end
        if InfoReliability(indexj,6)>1 %electricity LC
            operationalReliability(i,5) = 1;
        end
        % about EENS
        ft = operationalReliability(i,1);% from time
        tt = operationalReliability(i,2);% to time
        indexj_ft = max(find(InfoReliability(:,1)<=ft));
        indexj_tt = max(find(InfoReliability(:,1)<=tt));
        % if the period is completely contained in time step (not likely to happen)
        if indexj_ft == indexj_tt
            gasENS = InfoReliability(indexj_ft,5) * (tt - ft); % LC accumulated by hour
            electricityENS = InfoReliability(indexj_ft,6) * (tt - ft); 
        elseif indexj_tt - indexj_ft == 1
            gasENS = (InfoReliability(indexj_ft,2) - ft) * InfoReliability(indexj_ft,5)...
                + (tt - InfoReliability(indexj_tt,1)) * InfoReliability(indexj_ft,5);
            electricityENS = (InfoReliability(indexj_ft,2) - ft) * InfoReliability(indexj_ft,6)...
                + (tt - InfoReliability(indexj_tt,1)) * InfoReliability(indexj_ft,6);
        elseif indexj_tt - indexj_ft > 1
            gasENS = (InfoReliability(indexj_ft,2) - ft) * InfoReliability(indexj_ft,5)...
                + (tt - InfoReliability(indexj_tt,1)) * InfoReliability(indexj_ft,5);
            electricityENS = (InfoReliability(indexj_ft,2) - ft) * InfoReliability(indexj_ft,6)...
                + (tt - InfoReliability(indexj_tt,1)) * InfoReliability(indexj_ft,6);
            for k = (indexj_ft+1):(indexj_tt-1)
                gasENS = gasENS + (InfoReliability(k,2) - InfoReliability(k,1)) * InfoReliability(k,5);
                electricityENS = electricityENS + (InfoReliability(k,2) - InfoReliability(k,1)) * InfoReliability(k,6);
            end                
        else
            warnning('!!!!!!');
        end
        
        operationalReliability(i,4) = gasENS/((missionTime/numberOfPeriods));
        operationalReliability(i,6) = electricityENS/((missionTime/numberOfPeriods));  
    end
end
