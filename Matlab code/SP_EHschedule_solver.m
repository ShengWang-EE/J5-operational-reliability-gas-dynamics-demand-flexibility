function [sol,addCutsCoefficient,solverTime,exitflag] = SP_EHschedule_solver(SPformulation,ei_hat,gi_hat,EHpara,NK,iK)
%% unpack
optimalityObj = SPformulation.optimalityObj;
SPconstraints = SPformulation.SPconstraints;
feasibilityCheckObj = SPformulation.feasibilityCheckObj;
opts = SPformulation.opts;
[ei,gi,eee,ee3,e13,e1h,gg1,gg2,h1h,h14,h24,h2h,c3c,c4c,so,si,lc,...
    slackE1,slackE2,slackG1,slackG2] = deal(SPformulation.x{:});
%%
optimalityCons = [
     SPconstraints; % all the EH schedule constraints
    (ei == ei_hat):'eiBalance';
    (gi == gi_hat):'giBalance';
];    
feasibilityCheckCons = [
    SPconstraints; % all the EH schedule constraints
    (ei + slackE1 == ei_hat + slackE2):'eiBalance';
    (gi + slackG1 == gi_hat + slackG2):'giBalance';
    [slackE1 slackE2 slackG1 slackG2] >= 0;
    ];  
% check the feasibility of SP
% solve the optimality problem (SP)
diagnostics = optimize(optimalityCons,optimalityObj,opts);

if diagnostics.problem == 0 % feasible
    exitflag = 1;
    dual_ei = dual(optimalityCons('eiBalance'));
    dual_gi = dual(optimalityCons('giBalance'));
    EHcost_hat = value(optimalityObj);
else % infeasible, get other kind of dual
    exitflag = 0;
    diagnostics = optimize(feasibilityCheckCons,feasibilityCheckObj,opts);  
    % get dual and generate feasibility cut
    dual_ei = dual(feasibilityCheckCons('eiBalance'));
    dual_gi = dual(feasibilityCheckCons('giBalance'));
    EHcost_hat = value(feasibilityCheckObj);  
end


addCutsCoefficient.dual_ei = dual_ei; addCutsCoefficient.dual_gi = dual_gi;
addCutsCoefficient.ei_hat = ei_hat; addCutsCoefficient.gi_hat = gi_hat;
addCutsCoefficient.EHcost_hat = EHcost_hat;    





solverTime = diagnostics.solvertime;

h3h = deal(zeros(NK.studyPeriod,1));
[sol.eee, sol.ee3, sol.e13, sol.e1h, sol.gg1, sol.gg2, ...
    sol.h3h, sol.h1h, sol.h14, sol.h24, sol.h2h, sol.c3c, sol.c4c,...
    sol.so, sol.si, sol.lc] = deal ...
    (value(eee),value(ee3),value(e13),value(e1h),value(gg1),...
    value(gg2),value(h3h),value(h1h),value(h14),value(h24),value(h2h),value(c3c),value(c4c),...
    value(so), value(si), value(lc));
sol.EHcost = addCutsCoefficient.EHcost_hat;

if value(so)>0
    flag=1;
end

end
