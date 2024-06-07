function [results, probability] = predictNextState(currentComponentInfo,rts,mpc)
%% calculate possible state and probability
possibleState = []; counter = 0;
numberOfGenAndGfu = size(rts.gen,1);  numberOfGsou = size(rts.Gsou,1);
numberOfComponents = numberOfGenAndGfu +  numberOfGsou;
for i = 1:numberOfComponents
    % transition rate of this component
    if i <= numberOfGenAndGfu && ~any(rts.gfu(:,1)==i) % is tfu
        lamda = 1/rts.gen(i,3); mu = 1/rts.gen(i,4);
    elseif i <= numberOfGenAndGfu && any(rts.gfu(:,1)==i)
        gfuIndex = find(rts.gfu(:,1)==i);
        lamda = 1/rts.gfu(gfuIndex,2); mu = 1/rts.gfu(gfuIndex,3);
    else 
        lamda = 1/rts.Gsou(i-numberOfGenAndGfu,2); mu = 1/rts.Gsou(i-numberOfGenAndGfu,3);
    end
    
    state = currentComponentInfo(i);
    if i <= numberOfGenAndGfu
        if ismember(i,mpc.gfuIndex)
            % is gfu
            Nstate = rts.gfu(find(mpc.gfuIndex==i),4);
                     
        else
            Nstate = 1;
        end
    else
        Nstate = rts.Gsou(i-numberOfGenAndGfu,4);
    end
    
    if state == 0
        addPossibleState = currentComponentInfo;
        addPossibleState(i) = addPossibleState(i) + 1;       
        counter = counter + 1;
        transitionRate(counter) = mu;
    elseif state == Nstate
        addPossibleState = currentComponentInfo;
        addPossibleState(i) = addPossibleState(i) - 1;
        counter = counter + 1;
        transitionRate(counter) = lamda;
    else
        addPossibleState1 = currentComponentInfo; addPossibleState2 = currentComponentInfo;
        addPossibleState1(i) = addPossibleState(i) + 1;
        addPossibleState2(i) = addPossibleState(i) - 1;
        counter = counter + 1;
        transitionRate(counter) = mu;
        counter = counter + 1;
        transitionRate(counter) = lamda;
        addPossibleState = [addPossibleState1; addPossibleState2];
    end
    possibleState = [possibleState; addPossibleState];
end
probability = transitionRate' / sum(transitionRate);
%% calculate the pressure in each scenario
Nscenario = size(possibleState,1);
mpc0 = mpc;
results = cell(Nscenario,1);
for i = 1:Nscenario
    newmpc = mpcUpdate([i,0,3600,possibleState(i,:)],3600,mpc0,rts);
    results{i} = GErunopf(newmpc);
end
end


