function [sol,addCutsCoefficient,solverTime,exitflag] = SP_EHschedule_lp_solver(SPformulation,ei_hat,gi_hat,EHpara,NK,iK)
%% unpack
obj = SPformulation.obj;
generalEquations = SPformulation.generalEquations;
opts = SPformulation.opts;
[ei,gi,eee,ee3,e13,e1h,gg1,gg2,h1h,h14,h24,h2h,c3c,c4c,so,si,lc] = deal(SPformulation.x{:});
SPequations = replace(generalEquations,[ei,gi],[ei_hat,gi_hat]);
SPconstraints = [SPequations <= 0];
%%
diagnostics = optimize(SPconstraints,obj,opts);
solverTime = diagnostics.solvertime;

h3h = deal(zeros(NK.studyPeriod,1));
[sol.eee, sol.ee3, sol.e13, sol.e1h, sol.gg1, sol.gg2, ...
    sol.h3h, sol.h1h, sol.h14, sol.h24, sol.h2h, sol.c3c, sol.c4c,...
    sol.so, sol.si, sol.lc] = deal ...
    (value(eee),value(ee3),value(e13),value(e1h),value(gg1),...
    value(gg2),value(h3h),value(h1h),value(h14),value(h24),value(h2h),value(c3c),value(c4c),...
    value(so), value(si), value(lc));
sol.fval = value(obj);

if value(so)>0
    flag=1;
end


if diagnostics.problem ~= 0 %unbounded or infeasible

    exitflag = 0;
    u = diagnostics.solveroutput.result.farkasdual;
    sol.cost_hat = 0;
elseif diagnostics.problem == 0 %feasible
    exitflag = 1;
    u = dual(SPconstraints);
    sol.cost_hat = sol.fval;
else
    error('error');
end
addCutsEquation = [u'*generalEquations];%不管添加什么cut，这一部分的式子不变
% replace x(subproblem variables) with 0, only left with ei,gi
addCutsEquation = replace(addCutsEquation,...
    [eee, ee3, e13, e1h, gg1, gg2, h1h, h14, h24, h2h, c3c, c4c], ...
    [zeros(NK.studyPeriod,1),zeros(NK.studyPeriod,1),zeros(NK.studyPeriod,1),zeros(NK.studyPeriod,1),...
     zeros(NK.studyPeriod,1),zeros(NK.studyPeriod,1),zeros(NK.studyPeriod,1),zeros(NK.studyPeriod,1),...
     zeros(NK.studyPeriod,1),zeros(NK.studyPeriod,1),zeros(NK.studyPeriod,1),zeros(NK.studyPeriod,1)]);
addCutsEquation = replace(addCutsEquation,so,zeros(NK.DR,3));
% cannot replace multi-dimenssional sdpvars at once
addCutsEquation = replace(addCutsEquation,reshape(si,[],3),reshape(zeros(NK.DR,NK.nonDR,3),[],3));
addCutsEquation = replace(addCutsEquation,lc,zeros(NK.DR,3));
% 如果直接export，只能输出系数不为0的项的系数，这里做一点小处理: 
% 找一个几乎不可能出现的系数，手动加上，得到系数矩阵后再去掉
% coefficient of ei
addCutsEquation_wrt_ei = addCutsEquation + pi*sum(ei);
addCutsEquation_wrt_ei = replace(addCutsEquation_wrt_ei,gi,zeros(size(gi)));
addCut = [addCutsEquation_wrt_ei <=0];
[model,recoverymodel] = export(addCut,[],opts);
A.ei = model.A - pi; rhs = model.rhs; % rhs 是一样的
% coefficient of gi
addCutsEquation_wrt_gi = addCutsEquation + pi*sum(gi);
addCutsEquation_wrt_gi = replace(addCutsEquation_wrt_gi,ei,zeros(size(ei)));
addCut = [addCutsEquation_wrt_gi <=0];
[model,recoverymodel] = export(addCut,[],opts);
A.gi = model.A - pi;
%
[addCutsCoefficient.ei.A, addCutsCoefficient.rhs, ...
    addCutsCoefficient.gi.A] = ...
    deal(A.ei, rhs, A.gi);

end


function totalCost = objfcn_SP(lc,si,CDF,NK,iK)
curtailmentCost = 0; shiftCost = 0;
for l = 1:3 % 3 energies
    for ift = 1:NK.DR % from time
        curtailmentCost = curtailmentCost + lc(ift,l) * CDF;
        ft = iK.DR(ift);
        for itt = 1:NK.nonDR % to time
            tt = iK.nonDR(itt);
            shiftCost = shiftCost + si(ift,itt,l) * CDF/ 24 / 4 * abs(ft-tt);
        end
    end
end
totalCost = curtailmentCost + shiftCost; % sum all
end