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

nGl = size(mpc.Gline,1);
nGd = size(find(mpc.Gbus(:,3) > 0),1);
nEd = size(mpc.gen,1)-mpc.originalGenNumber;
[electricityReliabilityForOneSimulation, gasReliabilityForOneSimulation] = deal(cell(simulationTimes,1));

%set transient parameters
mpc.bus(:,3) = mpc.bus(:,3)*1.1;
[result0, success] = GErunopf(mpc);
[para] = initializeParameters2();

for i = 1:simulationTimes
    %% MCS
    [Info_components,failureComponents] = MCSformingScenarioV5(rts,missionTime);
    % -------------test-----------------
    Info_components = testInfo_components();
    numberOfStates = size(Info_components,1);

    % first scenario there is no failure, and it is used to determine the initial state only
    interruptionTime = Info_components(1,3) - Info_components(1,2);
    newmpc0 = mpcUpdate(Info_components(1,:),interruptionTime,mpc0,rts);% update CDF by the way
    % decide pipe section and time step for all scenarios
    [dx, nx] = decidePipelineCell(gtd);% maximun value of 3 miles per cell according to the reference
    duration = Info_components(:,3) - Info_components(:,2);
    [dt, nt] = decideTimeStep(duration);% 5min approximately
    cumnt = cumsum(nt);
    
    % calculate the steady result for each state
    newmpc = cell(numberOfStates,1); steadyResult = cell(numberOfStates,1);
    flag = zeros(numberOfStates,1);
    for j = 1:numberOfStates
        newmpc{j} = mpcUpdate(Info_components(j,:),duration(j),mpc0,rts);% update CDF by the way
        [steadyResult{j},flag(j)] = GErunopf(newmpc{j});
        if flag(j) == 0 % didn't converge£¬optimize without constraints
            [steadyResult{j}] = opfWithoutConstraint(steadyResult{j});
        end
    end


    %% initialization
    timeStepCounter = nt(1)+1;
    transientResult = struct();
    transientResult.P = cell(nGl,1); transientResult.Q = cell(nGl,1);
    initialCondition = setInitialConditionForGasSystem(steadyResult{1},gtd,nx,dx);
    initialSolutionPoint = setInitialConditionForGasSystem(steadyResult{2},gtd,nx,dx);
    for m = 1:nGl
    transientResult.P{m} = [repmat(initialCondition.P{m},nt(1),1);
                            zeros(sum(nt)-nt(1),nx(m)+1);]; 
    transientResult.Q{m} = [repmat(initialCondition.Q{m},nt(1),1);
                            zeros(sum(nt)-nt(1),nx(m)+1);]; 
    end
    save copy1.mat
    %%
    load copy1.mat
    for j = 2:numberOfStates
        % obtain the probabilities of state 3
        [possibleResults, probability] = predictNextState(Info_components(j,4:end),rts,mpc);
        optimalControlResult = lookaheadContingencyManagement_yalmip(initialCondition,initialSolutionPoint,newmpc{j},gtd,possibleResults, probability,dx,nx,dt(j),nt(j));
    end
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
%     end    

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