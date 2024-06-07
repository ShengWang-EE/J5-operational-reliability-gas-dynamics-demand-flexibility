function [f,jacobi] = fcn_transientGasFlowModelForOneTimeStep2(x,initialCondition,boundaryCondition,mpc,gtd,nx,dx,dt,para)
% version 2 use bus to calculate bus gas conservation and pressure euqation
%% grab dimenssions
nGl = size(mpc.Gline,1);
nGb = size(mpc.Gbus,1);
% nVar = size(x,1);
% cumnx = [0;cumsum(nx+1)];
%% sort x
deltat = formatXtoPQcell(x,nGl,nx);
for i = 1:size(deltat.P)
    deltat.P{i} = deltat.P{i} * 1e5;
end
%% f consist of fcontinuity and motion, boundary conditions   
f = [];
% jacobi = zeros(nVar,nVar);
% jacobi_continuity_motion = cell(nGl,1);

f_continuity_motion = cell(nGl,1);

for i = 1:nGl
% fcontinuity and motion    
    % variables parameters
    allP = initialCondition.P{i}; allQ = initialCondition.Q{i};
    allP_T = deltat.P{i}; allQ_T = deltat.Q{i};
    diameter = gtd.Gline(i,3);
    F = gtd.Gline(i,5);
    [f_continuity_motion{i}] = pipeSectionForOneTimeStep2(allP,allQ,allP_T,allQ_T,diameter,F,nx(i),dx(i),dt,para);
 
    f = [f; f_continuity_motion{i}];  
end
% jacobi = jacobi';
%% for cosntrainsts concernning every bus
f_equalP = zeros(2*nGl-nGb,1);
nGb1 = size(find(mpc.Gbus(:,2)==1),1);nGb2 = size(find(mpc.Gbus(:,2)==2),1);nGb3 = size(find(mpc.Gbus(:,2)==3),1);
f_specifyP = zeros(nGb2,1);
f_conservationQ = zeros(nGb1+nGb3,1);
counter1 = 1; counter2 = 1; counter3 = 1;
for i = 1:nGb
    % equal P
    connectLineAsFromBus = find(mpc.Gline(:,1)==i);
    connectLineAsToBus = find(mpc.Gline(:,2)==i);
    connectLine = [connectLineAsFromBus; connectLineAsToBus];

    if size(connectLine,1) > 1
        for j = 2:size(connectLineAsFromBus,1) %take the first P in from line as a base
            f_equalP(counter1) = deltat.P{connectLineAsFromBus(1)}(1) - deltat.P{connectLineAsFromBus(j)}(1);
            counter1 = counter1 + 1;
        end
        for j = 1:size(connectLineAsToBus,1)
            f_equalP(counter1) = deltat.P{connectLineAsFromBus(1)}(1) - deltat.P{connectLineAsToBus(j)}(end);
            counter1 = counter1 + 1;
        end  
    end
    % specify P
    switch mpc.Gbus(i,2)
        case 2
            % if is a junction of multiple bus, then only have to specify
            % one bus, for the P is already equal
            firstConnectLine = min([find(mpc.Gline(:,1)==i); find(mpc.Gline(:,2)==i)]);
            if i==mpc.Gline(firstConnectLine,1) % is the from bus
                f_specifyP(counter2) = deltat.P{firstConnectLine}(1) - boundaryCondition.P(i);
                counter2 = counter2 + 1;
            else % is the to bus
                f_specifyP(counter2) = deltat.P{firstConnectLine}(end) - boundaryCondition.P(i);
                counter2 = counter2 + 1;
            end
        case {1,3}

        for j = 1:size(connectLineAsFromBus,1)
            f_conservationQ(counter3) = f_conservationQ(counter3) + deltat.Q{connectLineAsFromBus(j)}(1);  
        end
        for j = 1:size(connectLineAsToBus,1)
            f_conservationQ(counter3) = f_conservationQ(counter3) - deltat.Q{connectLineAsToBus(j)}(end);
        end
        % add gas source to the sum, and there is only one gas source
        % at each bus
         gasSourceIndex = find(boundaryCondition.gasSource(:,1)==i);
         gasSource = boundaryCondition.gasSource(gasSourceIndex,5);
         if ~isempty(gasSource)             
             f_conservationQ(counter3) = f_conservationQ(counter3) - gasSource;       
         end
         counter3 = counter3 + 1;
    end
end
f_specifyP = f_specifyP / 1e5;
f = [f; f_equalP; f_specifyP; f_conservationQ]; 

end
