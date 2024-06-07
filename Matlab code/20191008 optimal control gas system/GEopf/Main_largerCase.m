%% 
% this is a larger case study
clear
clc
tic
% load case
[mpc, gtd] = case24GEv3();
mpc0 = mpc; %back up
rts = Case24ReliabillityDatav2();
% set MCS parameters
simulationTimes = 1;
missionTime = 12;
numberOfPeriods = 144;
nGl = size(mpc.Gline,1);
nGd = size(find(mpc.Gbus(:,3) > 0),1);
nEd = size(mpc.gen,1)-mpc.originalGenNumber;
[electricityReliabilityForOneSimulation, gasReliabilityForOneSimulation] = deal(cell(simulationTimes,1));
fromTimeToTime = zeros(numberOfPeriods,2);
fromTimeToTime(:,1) = [0:1:numberOfPeriods-1]' * missionTime/numberOfPeriods; 
fromTimeToTime(:,2) = fromTimeToTime(:,1) + missionTime/numberOfPeriods;
%set transient parameters
[result0, success] = GErunopf(mpc);
[para] = initializeParameters2();

for i = 1:simulationTimes
    %% MCS
    [Info_components,failureComponents] = MCSformingScenarioV5(rts,missionTime);
    numberOfStates = size(Info_components,1);

    % first scenario there is no failure, and it is used to determine the initial state only
    interruptionTime = Info_components(1,3) - Info_components(1,2);
    newmpc0 = mpcUpdate(Info_components(1,:),interruptionTime,mpc0,rts);% update CDF by the way
    [steadyResult0] = GErunopf(newmpc0);
    % decide pipe section and time step for all scenarios
    [dx, nx] = decidePipelineCell(gtd);% maximun value of 3 miles per cell according to the reference
    initialCondition = setInitialConditionForGasSystem(steadyResult0,gtd,nx,dx);
    duration = Info_components(:,3) - Info_components(:,2);
    [dt, nt] = decideTimeStep(duration);% 5min approximately
    cumnt = cumsum(nt);
    %% initialization
    timeStepCounter = nt(1)+1;
    transientResult = struct();
    transientResult.P = cell(nGl,1); transientResult.Q = cell(nGl,1);
    for m = 1:nGl
    transientResult.P{m} = [repmat(initialCondition.P{m},nt(1),1);
                            zeros(sum(nt)-nt(1),nx(m)+1);]; 
    transientResult.Q{m} = [repmat(initialCondition.Q{m},nt(1),1);
                            zeros(sum(nt)-nt(1),nx(m)+1);]; 
    end
    steadyResult = cell(numberOfStates,1);
    newmpc = cell(numberOfStates,1);
    steadyResult{1} = steadyResult0;
    newmpc{1} = newmpc0;
    
    if numberOfStates == 1
        disp('no transient')
        electricityReliabilityForOneSimulation{i}.system = zeros(numberOfPeriods,4);
        electricityReliabilityForOneSimulation{i}.nodal = zeros(numberOfPeriods,2*nEd);
        gasReliabilityForOneSimulation{i}.system = zeros(numberOfPeriods,4);
        gasReliabilityForOneSimulation{i}.nodal = zeros(numberOfPeriods,2*nGd);
        electricityReliabilityForOneSimulation{i}.system(:,[1 2]) = fromTimeToTime;
        gasReliabilityForOneSimulation{i}.system(:,[1 2]) = fromTimeToTime;
    else
        for j = 2:numberOfStates
            %% run IOPF to determine the LC in steady state
            % the units used in IOPF are bar and Mm3/day
            newmpc{j} = mpcUpdate(Info_components(j,:),duration(j),mpc0,rts);% update CDF by the way
            % test-----------------------------
            newmpc{j} = testmpc_case24GEv3();
            % test---------------------------------
            [steadyResult{j},flag] = GErunopf(newmpc{j});
            if flag == 0 % didn't converge£¬optimize without constraints
                [steadyResult{j}] = opfWithoutConstraint(steadyResult{j});
            end
    %         [LCe,LCg,Pg_gfu,totalCost] = GEopfResultDecomposition(result1);
            %% gas transient
            % the units in transient gas flow are Pa and Mm3/day
            boundaryCondition = setBoundaryConditionForGasSystem(steadyResult{j},newmpc{j});

            for k = 1:nt(j)
                [nextTimeState] = transientGasPowerFlowForOneTimeStep(initialCondition,boundaryCondition, newmpc{j},gtd, nx,dx,dt(j),para);
                
                % renew initial condition
                initialCondition = nextTimeState;
                % output
                for m = 1:nGl
                    transientResult.P{m}(timeStepCounter,:) = nextTimeState.P{m}; 
                    transientResult.Q{m}(timeStepCounter,:) = nextTimeState.Q{m};
                end
                % determine whether should skip the remainning time steps
                YESskip = 0;
                if timeStepCounter > 2
                    [YESskip] = skipRemainingTimeStep(transientResult,nGl,timeStepCounter);
                end
                if YESskip
                    for m = 1:nGl
                        transientResult.P{m}(timeStepCounter+1:cumnt(j),:) = repmat(nextTimeState.P{m},cumnt(j)-timeStepCounter,1); 
                        transientResult.Q{m}(timeStepCounter+1:cumnt(j),:) = repmat(nextTimeState.Q{m},cumnt(j)-timeStepCounter,1);
                    end
                    timeStepCounter = cumnt(j)+1;
                    break                
                else
                    timeStepCounter = timeStepCounter + 1;
                end

            end 
        end
        % diminish the over shoot
%         transientResult = diminishOvershoot(transientResult,nGl);
        %% calculate the gas curtailment £¨all time in one simulation)
        [gasLC,GFUcapacity] = calculateGasLC(transientResult,steadyResult,newmpc{j},nt);

        gasReliabilityForOneSimulation{i}  = calculateGasReliabilityFromTransient(gasLC,nt,dt,numberOfPeriods,missionTime,'gas');
        % run electricity power flow to determine the electricity load
        % curtailement
        electricityLC = zeros(cumnt(end),nEd); 
        for k = 1:cumnt(end)     
            [newmpcForElectricityPowerFlow] = updatempcFromTransient(GFUcapacity(k,:),newmpc,k,nt);
            %
            steadyElectricityResult = runopf(newmpcForElectricityPowerFlow);
            electricityLC(k,:) = steadyElectricityResult.gen(newmpc{j}.originalGenNumber+1:end,2)';
        end
        % diminish the overshoot of electricity directly to speed up the convergence
        for k = 1:size(electricityLC,2)
            electricityLC(:,i) = diminishOvershootOfaCurve(electricityLC(:,k));
        end
        % using the same function as in the gas reliability, but with
        % differernt threshold
        electricityReliabilityForOneSimulation{i}  = ...
            calculateGasReliabilityFromTransient(electricityLC,nt,dt,numberOfPeriods,missionTime,'electricity');
    end    
end
% deal with the reliability over whole simulations
[allElectricityReliability.system, allGasReliability.system, allElectricityReliability.nodal, allGasReliability.nodal] = ...
    deal(zeros(numberOfPeriods,4,simulationTimes),zeros(numberOfPeriods,4,simulationTimes), ...
    zeros(numberOfPeriods,2*nEd,simulationTimes),zeros(numberOfPeriods,2*nGd,simulationTimes));
stdSystemReliability = zeros(numberOfPeriods,4,simulationTimes);
for i = 1:simulationTimes
    allElectricityReliability.system(:,:,i) = electricityReliabilityForOneSimulation{i}.system;
    allElectricityReliability.nodal(:,:,i) = electricityReliabilityForOneSimulation{i}.nodal;
    allGasReliability.system(:,:,i) = gasReliabilityForOneSimulation{i}.system;
    allGasReliability.nodal(:,:,i) = gasReliabilityForOneSimulation{i}.nodal;
    
    if mod(i,1000) == 0
        stdSystemReliability(:,:,i) = std(allElectricityReliability.system(:,:,1:i),0,3)./mean(allElectricityReliability.system(:,:,1:i),3)/sqrt(i); 
    end
end
% 
electricityReliability.system = mean(allElectricityReliability.system,3);
electricityReliability.nodal = mean(allElectricityReliability.nodal,3);
gasReliability.system = mean(allGasReliability.system,3);
gasReliability.nodal = mean(allGasReliability.nodal,3);

elapseTime = toc  
save