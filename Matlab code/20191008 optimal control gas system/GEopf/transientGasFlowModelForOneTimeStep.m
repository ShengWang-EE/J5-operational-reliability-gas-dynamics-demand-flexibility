function [f,jacobi] = transientGasFlowModelForOneTimeStep(x,initialCondition,boundaryCondition,mpc,gtd,nx,dx,nt,dt)
%% grab dimenssions
nGl = size(mpc.Gline,1);
nVar = size(x,1);
cumnx = [0;cumsum(nx+1)];
%% sort x
deltat = formatXtoPQcell(x,nGl,nx);
    
%% f consist of fcontinuity and motion, boundary conditions   
f = [];
jacobi = zeros(nVar,nVar);
jacobi_continuity_motion = cell(nGl,1);
j_fbBoundary = cell(nGl,1);
j_tbBoundary = cell(nGl,1);
f_continuity_motion = cell(nGl,1);
f_fbBoundary_ob = 999*ones(nGl,1);
f_tbBoundary_ob = 999*ones(nGl,1);
for i = 1:nGl
% fcontinuity and motion    
    % variables parameters
    allP = initialCondition.P{i}; allQ = initialCondition.Q{i};
    allP_T = deltat.P{i}; allQ_T = deltat.Q{i};
    diameter = gtd.Gline(i,3);
    F = gtd.Gline(i,5);
    [f_continuity_motion{i}, jacobi_continuity_motion{i}] = pipeSectionForOneTimeStep2(allP,allQ,allP_T,allQ_T,diameter,F,nx(i),dx(i),nt,dt);

    
% specify boundary conditions as equations for pipelines 
% only have to specify one value in each bus?
    j_fbBoundary{i} = zeros(nVar,1); j_tbBoundary{i} = zeros(nVar,1); % x order:[P(1~n+1); Q(1~n+1)]
    fb = mpc.Gline(i,1); tb = mpc.Gline(i,2);
    switch mpc.Gbus(fb,2)
        case 1 % fist type: gas system-in bus, constant Q
            f_fbBoundary = allQ_T(1) - boundaryCondition.Q(i,1); 
            j_fbBoundary{i}(nVar/2 + cumnx(i) + 1) = 1;
        case 2 % gas load, constant P
            f_fbBoundary = allP_T(1) - boundaryCondition.P(i,1);
            j_fbBoundary{i}(cumnx(i) + 1) = 1;
        case 3 % connection bus
            fbConnectBusAsFromBusLine = find(mpc.Gline(:,1)==fb);
            fbConnectBusAsFromBus = mpc.Gline(fbConnectBusAsFromBusLine,2);
            fbConnectBusAsToBusLine = find(mpc.Gline(:,2)==fb);
            fbConnectBusAsToBus = mpc.Gline(fbConnectBusAsToBusLine,1);         
            %
%             f_fbBoundary = allQ_T(1); % out-flow is assumed positive in fb
            f_fbBoundary = 0;
            for j = 1:size(fbConnectBusAsFromBusLine,1)
                Q_fbConnectBusAsFromBus = deltat.Q{fbConnectBusAsFromBusLine(j)};
                f_fbBoundary = f_fbBoundary + Q_fbConnectBusAsFromBus(1);
                j_fbBoundary{i}(nVar/2 + cumnx(fbConnectBusAsFromBusLine(j)) + 1) = ...
                    j_fbBoundary{i}(nVar/2 + cumnx(fbConnectBusAsFromBusLine(j)) + 1) + 1;            
            end
            for j = 1:size(fbConnectBusAsToBusLine,1)
                Q_fbConnectBusAsToBus = deltat.Q{fbConnectBusAsToBusLine(j)};
                f_fbBoundary = f_fbBoundary - Q_fbConnectBusAsToBus(end);
                j_fbBoundary{i}(nVar/2 + cumnx(fbConnectBusAsToBusLine(j)+1)) = ...
                    j_fbBoundary{i}(nVar/2 + cumnx(fbConnectBusAsToBusLine(j)+1)) - 1;
            end
            % add gas source to the sum, and there is only one gas source
            % at each bus
             gasSourceIndex = find(boundaryCondition.gasSource(:,1)==fb);
             gasSource = boundaryCondition.gasSource(gasSourceIndex,5);
             if ~isempty(gasSource)             
                 f_fbBoundary = f_fbBoundary - gasSource;       
             end

             
    end
    
    switch mpc.Gbus(tb,2)
        case 1 % fist type: gas source bus, constant Q
            f_tbBoundary = allQ_T(nx(i)+1) - boundaryCondition.Q(i,2); 
            j_tbBoundary{i}(nVar/2 + cumnx(i+1)) = 1;
        case 2 % gas load, constant P
            f_tbBoundary = allP_T(nx(i)+1) - boundaryCondition.P(i,2);
            j_tbBoundary{i}(nVar/2 + cumnx(i+1)) = 1;
        case 3 % connection bus
            tbConnectBusAsFromBusLine = find(mpc.Gline(:,1)==tb);
            tbConnectBusAsFromBus = mpc.Gline(tbConnectBusAsFromBusLine,2);
            tbConnectBusAsToBusLine = find(mpc.Gline(:,2)==tb);
            tbConnectBusAsToBus = mpc.Gline(tbConnectBusAsToBusLine,1);

%             f_tbBoundary = allQ_T(end); % in-flow is assumed positive
            f_tbBoundary = 0;
            for j = 1:size(tbConnectBusAsFromBusLine,1)
                Q_tbConnectBusAsFromBus = deltat.Q{tbConnectBusAsFromBusLine(j)};
                f_tbBoundary = f_tbBoundary - Q_tbConnectBusAsFromBus(1);
                j_tbBoundary{i}(nVar/2 + cumnx(tbConnectBusAsFromBusLine(j)) + 1) = ...
                    j_tbBoundary{i}(nVar/2 + cumnx(tbConnectBusAsFromBusLine(j)) + 1) - 1;  
            end
            for j = 1:size(tbConnectBusAsToBusLine,1)
                Q_tbConnectBusAsToBus = deltat.Q{tbConnectBusAsToBusLine(j)};
                f_tbBoundary = f_tbBoundary + Q_tbConnectBusAsToBus(end);
                j_tbBoundary{i}(nVar/2 + cumnx(tbConnectBusAsToBusLine(j)+1)) = ...
                    j_tbBoundary{i}(nVar/2 + cumnx(tbConnectBusAsToBusLine(j)+1)) + 1;       
            end
            % add gas source to the sum, and there is only one gas source
            % at each bus
             gasSourceIndex = find(boundaryCondition.gasSource(:,1)==tb);
             gasSource = boundaryCondition.gasSource(gasSourceIndex,5);
             if ~isempty(gasSource)
                f_tbBoundary = f_tbBoundary + gasSource;
             end
    end
 
    f = [f; f_continuity_motion{i}; f_fbBoundary; f_tbBoundary];
    
    
    jacobi_continuity_motion_extend = zeros(nVar,2*nx(i));
    jacobi_fbBoundary_extend = zeros(nVar,1); 
    jacobi_tbBoundary_extend = zeros(nVar,1);
    
    jacobi_continuity_motion_extend(cumnx(i)+1:cumnx(i+1),:) = jacobi_continuity_motion{i}(1:(nx(i)+1),:); % the half of dP
    jacobi_continuity_motion_extend((nVar/2 + cumnx(i)+1):(nVar/2+cumnx(i+1)),:) = jacobi_continuity_motion{i}((nx(i)+1+1):(2*(nx(i)+1)),:); % the half of dQ
    jacobi_fbBoundary_extend(cumnx(i)+1:cumnx(i+1)) = j_fbBoundary{i}(1:(nx(i)+1)); % the half of dP
    jacobi_tbBoundary_extend(cumnx(i)+1:cumnx(i+1)) = j_tbBoundary{i}(1:(nx(i)+1));
    jacobi(:,2*cumnx(i)+1:2*cumnx(i+1)) = [jacobi_continuity_motion_extend, jacobi_fbBoundary_extend,jacobi_tbBoundary_extend];
    
    
    
    % ob
    f_fbBoundary_ob(i) = f_fbBoundary;
    f_tbBoundary_ob(i) = f_tbBoundary;
    ob = max(f);
end
jacobi = jacobi';
% some of the constraints in bus type 3 are linear correlated, therefore add constraints of P    

end

