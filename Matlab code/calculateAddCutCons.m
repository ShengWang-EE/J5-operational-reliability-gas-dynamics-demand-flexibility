function addCutCons = calculateAddCutCons(SP_EHresult,addCutsCoefficient,nEH,iter,mpc)
addCutEquations = [];addCutRHSs = [];

    for i = 1:nEH
        dual_ei = addCutsCoefficient{iter-1,i}.dual_ei; dual_gi = addCutsCoefficient{iter-1,i}.dual_gi;
        ei_hat = addCutsCoefficient{iter-1,i}.ei_hat;   gi_hat = addCutsCoefficient{iter-1,i}.gi_hat;
        EHcost_hat_hat = addCutsCoefficient{iter-1,i}.EHcost_hat;   
        addCutEquations = [addCutEquations; 
            EHcost_hat_hat - dual_ei'*(mpc.baseMVA*ei(:,i)-ei_hat) - dual_gi'*(200*gi(:,i)-gi_hat);]; % 注意gi前面的系数需要折算
        if exitflagSP(iter,i) == 0 % SP infeasible   
            addCutRHSs = [addCutRHSs;0;];
        elseif exitflagSP(iter,i) == 1 % SP feasible
            addCutRHSs = [addCutRHSs;EHcost_hat(i);];
        end
    end

addCutCons = [addCutEquations <= addCutRHSs];
end