function [MP_IEGSresult, solverTime,MPformulation] = MP_IEGSdispatch_solver(addCutsCoefficient,exitflagSP,SP_EHresult,...
    dayahead_IEGSresult_basicLoad,dayaheadEHschedule,EHreserve,mpc,gtd,nx,dx,dt,NH,nEH,NK,iK,nodalPrice,iter,MPformulation)
% 根据每小时日前运行的结果，作为暂态分析每小时及其间每个时间片段的初值 
% 注意：原来dayahead中所设定的电力、天然气的节点负荷中包含了ei，gi
nGb  = size(mpc.Gbus,1); % number of gas bus
nGl  = size(mpc.Gline,1); % number of gas line
counter = 0;
for k = iK.studyPeriod(1):iK.studyPeriod(end)
    % 判断k所在的小时
    h = k2h_h2k(k,4,'k2h');
    counter = counter +1 ;
    
    PGs0(:,counter) = dayahead_IEGSresult_basicLoad{h}.Gsou(:,5);% 这个不是变量，而是定值

end
%% add cut 
[Prs,Gf,Va,Pg,ei,gi,EHcost_hat] = deal(MPformulation.x{:});
addCutEquations = [];addCutRHSs = [];
addCutCons = [];
if ~isempty(SP_EHresult) && ~isempty(addCutsCoefficient)


    for i = 1:nEH
%         dual_ei = addCutsCoefficient{i}.dual_ei; dual_gi = addCutsCoefficient{i}.dual_gi;
%         ei_hat = addCutsCoefficient{i}.ei_hat;   gi_hat = addCutsCoefficient{i}.gi_hat;
%         EHcost_hat_hat = addCutsCoefficient{i}.EHcost_hat;   
%         addCutEquations = [addCutEquations; 
%             EHcost_hat_hat - dual_ei'*(mpc.baseMVA*ei(:,i)-ei_hat) - dual_gi'*(200*gi(:,i)-gi_hat);]; % 注意gi前面的系数需要折算
%         if exitflagSP(iter-1,i) == 0 % SP infeasible   
%             addCutRHSs = [addCutRHSs;0;];
%         elseif exitflagSP(iter-1,i) == 1 % SP feasible
%             addCutRHSs = [addCutRHSs;EHcost_hat(i);];
%         end

        A.ei = addCutsCoefficient{i}.ei.A; A.gi = addCutsCoefficient{i}.gi.A; 
        rhs = addCutsCoefficient{i}.rhs;            
        addCutEquations = [addCutEquations; 
             A.ei * mpc.baseMVA * ei(:,i) + A.gi * 200 * gi(:,i);]; % 注意gi前面的系数需要折算
        if exitflagSP(iter-1,i) == 0 % SP infeasible   
            addCutRHSs = [addCutRHSs;rhs;];
        elseif exitflagSP(iter-1,i) == 1 % SP feasible
            addCutRHSs = [addCutRHSs;rhs+EHcost_hat(i);];
        end
        
    end

addCutCons = [addCutEquations <= addCutRHSs];
end



%%
objfcn = MPformulation.obj;
constraints = [MPformulation.cons;addCutCons;];
options = MPformulation.opts;
% options.gurobi.Method = 1;
diagnostics  = optimize(constraints,objfcn,options);
if diagnostics.problem == 4
    flag = 1;
end
% diagnostics  = optimize(constraints,[],options);
solverTime = diagnostics.solvertime;
%% results

IEGScost = objfcn_IEGSdispatch(Pg,mpc,NK.studyPeriod,nodalPrice);
[~,sumtfuCost,sumgfuCost] = objfcn_IEGSdispatch(Pg,mpc,NK.studyPeriod,nodalPrice);
constraintViolation = nan;
[Pg1,Va1,PrsRaw,GfRaw,ei1,gi1,EHcost_hat1,IEGScost1,sumtfuCost1,sumgfuCost1,constraintViolation1,f_MP] = deal(value(Pg),value(Va),value(Prs),value(Gf),...
    value(ei),value(gi),value(EHcost_hat),value(IEGScost),value(sumtfuCost),value(sumgfuCost),value(constraintViolation),value(objfcn)); 

PrsPipe = mat2cell(PrsRaw,NK.studyPeriod,nx+1); Gf1 = mat2cell(GfRaw,NK.studyPeriod,nx+1);
PrsNodal = zeros(nGb,1);
for i = 1:nGl
    PrsNodal(mpc.Gline(i,1)) = PrsPipe{i}(1);
    PrsNodal(mpc.Gline(i,2)) = PrsPipe{i}(end);
end


[MP.Pg,MP.PrsPipe,MP.PrsNodal,MP.Gf,MP.PGs,MP.ei,MP.gi,MP.EHcost_hat,MP.IEGScost,MP.sumtfuCost, ...
    MP.sumgfuCost,MP.constraintViolation,MP.f_MP] = ...
    deal(Pg1,PrsPipe,PrsNodal,Gf1,PGs0,ei1,gi1,EHcost_hat1,IEGScost1,sumtfuCost1, ...
    sumgfuCost1,constraintViolation1,f_MP);
MP_IEGSresult = MP;
MPformulation.cons = constraints;
end

