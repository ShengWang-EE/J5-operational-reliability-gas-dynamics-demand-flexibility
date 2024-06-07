function [sol,addCutsCoefficient,solverTime,exitflag] = SP_EHschedule(ei_hat,gi_hat,EHpara,NK,iK)
yalmip('clear');
% data
capacityCell = num2cell(EHpara.capacity);
efficiencyCell = num2cell(EHpara.efficiency);
scheduleCell = num2cell(EHpara.schedule);
[HA,EA,HB,EB,HC,EC,HD,ED,HO2MAX,HO2MIN,HO3MAX,HO3MIN,CO3MAX,CO3MIN,CO4MAX,CO4MIN] = deal(capacityCell{:});
[eta1e,eta1h,eta2,COP3hSummer,COP3hWinter,COP3c,COP4] = deal(efficiencyCell{:});
[alpha,beta,CDF] = deal(scheduleCell{:});
load0 = EHpara.load;
COP3h = COP3hSummer;
gamma = 0;%summer,EHP cooling
NL = 3;

%------------------test---------------
beta = 1;
%-------------------------------
iK.nonDR = iK.studyPeriod; iK.nonDR(ismember(iK.studyPeriod,iK.DR)) = [];
%% vars
% [ei; gi; eee, ee3, e13, e1h, gg1, gg2, h3h, h1h, h14, h24, h2h, c3c, c4c] = deal(sdpvar(NK.studyPeriod,1));
[ei] = sdpvar(NK.studyPeriod,1);[gi] = sdpvar(NK.studyPeriod,1);[eee] = sdpvar(NK.studyPeriod,1);
[ee3] = sdpvar(NK.studyPeriod,1);[e13] = sdpvar(NK.studyPeriod,1);[e1h] = sdpvar(NK.studyPeriod,1);
[gg1] = sdpvar(NK.studyPeriod,1);[gg2] = sdpvar(NK.studyPeriod,1);[h3h] = sdpvar(NK.studyPeriod,1);
[h1h] = sdpvar(NK.studyPeriod,1);[h14] = sdpvar(NK.studyPeriod,1);[h24] = sdpvar(NK.studyPeriod,1);
[h2h] = sdpvar(NK.studyPeriod,1);[c3c] = sdpvar(NK.studyPeriod,1);[c4c] = sdpvar(NK.studyPeriod,1);
so = sdpvar(NK.DR,NL); % shift out for DR period for three loads
si = sdpvar(NK.DR,NK.nonDR,NL); % shift in quantity for each shiftout period and shiftin period
lc = sdpvar(NK.DR,NL); % load curtailment for each DR period
% start value
[ei0, gi0, eee0, ee30, e130, e1h0, gg10, gg20, h3h0, h1h0, h140, h240, h2h0, c3c0, c4c0] = deal(zeros(NK.studyPeriod,1));
so0 = zeros(NK.DR,3); % shift out for DR period for three loads
si0 = zeros(NK.DR,NK.nonDR,3); % shift in quantity for each shiftout period and shiftin period
lc0 = zeros(NK.DR,3); % load curtailment for each DR period
% some constants
h3h = h3h0;

assign([ei, gi, eee, ee3, e13, e1h, gg1, gg2, h1h, h14, h24, h2h, c3c, c4c],...
    [ei_hat, gi_hat, eee0, ee30, e130, e1h0, gg10, gg20, h1h0, h140, h240, h2h0, c3c0, c4c0]);
assign(so,so0); assign(si,si0); assign(lc,lc0);
%%
% calculate new load curve
newLoad(:,iK.DR) = load0(:,iK.DR) - so' - lc';
newLoad(:,1:(iK.DR(1)-1)) = load0(:,1:(iK.DR(1)-1));
newLoad = [newLoad load0(:,(iK.DR(end)+1):NK.all)];
newLoad(:,iK.nonDR) = newLoad(:,iK.nonDR) + reshape(sum(si,1),size(si,2),size(si,3))';

[el,hl,cl] = deal(newLoad(1,iK.studyPeriod)',newLoad(2,iK.studyPeriod)',newLoad(3,iK.studyPeriod)');
% schedule capacity
scheduleEquations = [
    -so;
    so - alpha * load0(:,iK.DR)';% limits for shiftout capability
    -lc;
    lc - beta * load0(:,iK.DR)'; % limits for curtailment
    so - reshape(sum(si,2),size(si,1),size(si,3)) ; % shiftout should equal shiftin (here relax the equalty constraints)
    -reshape(si,[],3);
    ];
scheduleEquations = reshape(scheduleEquations,[],1);
scheduleConstraints = [scheduleEquations <= 0]:'scheduleConstraints';
% energy conversion constraints (automatically formulated for each period)
% energyConversionConstraints.inner = [
%     eee+e1h >= el;
%     (ee3+e13)*COP3h*gamma >= h3h;
%     (ee3+e13)*COP3h*(1-gamma) >= c3c;
%     gg1*eta1e >= e13+e1h;
%     gg1*eta1h >= h1h + h14;
%     gg2*eta2 >= h24+h2h;
%     (h14+h24)*COP4 >= c4c;
%     h3h+h1h+h2h >= hl;
%     c3c+c4c >= cl;
% ]:'energyConversionConstraints.inner';
% energyConversionConstraints.outter = [
%     ei_hat >= eee+ee3;
%     gi_hat >= gg1+gg2;
%     ]:'energyConversionConstraints.outter';
energyConversionEquations = [
     eee+ee3 - ei;
     gg1+gg2 - gi;
    -(eee+e1h) + el;
%     h3h - ((ee3+e13)*COP3h*gamma);
    c3c - ((ee3+e13)*COP3h*(1-gamma));
    e13+e1h - gg1*eta1e;
    h1h + h14 - gg1*eta1h;
    h24+h2h - gg2*eta2;
    c4c - (h14+h24)*COP4;
    hl - (h3h+h1h+h2h);
    cl - (c3c+c4c);
    -eee;
    -ee3;
    -e13;
    -e1h;
    -gg1;
    -gg2;
    -h1h;
    -h14;
    -h24;
    -h2h;
    -c3c;
    -c4c;
    ];
energyConversionConstraints = [energyConversionEquations<=0]:'energyConversionConstraints';

% capacity constraints
% capacityConstraints = [
%     h1h+h14 >= 0;
%     e13+e1h-EA-(EA-EB)/(HA-HB)*(h1h+h14) <= 0;
%     e13+e1h-EB-(EB-EC)*(HB-HC)*(h1h+h14-HB) >=0;
%     e13+e1h-ED-(EC-ED)*(HC-HD)*(h1h+h14) >= 0;
%     HO2MIN <= h24+h2h <= HO2MAX;
%     CO4MIN <= c4c <= CO4MAX;
%     gamma*HO3MIN <= h3h <= gamma*HO3MAX;
%     (1-gamma)*CO3MIN <= c3c <= (1-gamma)*CO3MAX;
%     ]:'capacityConstraints';
capacityEquations = [
    -(h1h+h14);
    e13+e1h-EA-(EA-EB)/(HA-HB)*(h1h+h14);
    -( e13+e1h-EB-(EB-EC)*(HB-HC)*(h1h+h14-HB) );
    -( e13+e1h-ED-(EC-ED)*(HC-HD)*(h1h+h14) );
    -(h24+h2h) + HO2MIN;
    h24+h2h - HO2MAX;
    -c4c + CO4MIN;
    c4c - CO4MAX;
%     -h3h + gamma*HO3MIN;  
%     h3h - gamma*HO3MAX;
    -c3c + (1-gamma)*CO3MIN; 
    c3c - (1-gamma)*CO3MAX;
    ];
capacityConstraints = [capacityEquations <= 0]:'capacityConstraints';
% DR constraint 这个约束考虑在了IEGS主网中，还是考虑在这里？
% DRconstraint = [
%     ei <= ei0 - reserve;
%     ];
% in this general form, MP vars (ei,gi) is considered vars. in this SP,
% replace them as constants
generalEquations = [scheduleEquations;  energyConversionEquations;  capacityEquations;]; 
energyConversionEquations1 = replace(energyConversionEquations,[ei,gi],[ei_hat,gi_hat]);
SPequations = [scheduleEquations;  energyConversionEquations1;  capacityEquations;];
% generalEquations = [scheduleEquations;  energyConversionEquations;  capacityEquations;]; 
% SPequations = generalEquations;
% SPequations = replace(generalEquations([1:2]),[ei,gi],[ei_hat,gi_hat]);
SPconstraints = [SPequations <= 0];
%%
status = 0;
% solve
obj = objfcn_SP(lc,si,CDF,NK,iK);
% obj = [];
opts = sdpsettings('verbose',2,'solver','gurobi','usex0',0,'debug',1,'gurobi.InfUnbdInfo',1,...
    'savesolverinput',1,'savesolveroutput',1);
diagnostics = optimize(SPconstraints,obj,opts);
solverTime = diagnostics.solvertime;


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

% 先进行运算，再把x代为0
% % replace x(subproblem variables) with 0, only left with ei,gi
% cutConstraints = replace(generalConstraints,...
%     [eee, ee3, e13, e1h, gg1, gg2, h3h, h1h, h14, h24, h2h, c3c, c4c],...
%     [zeros(NK.study,1),zeros(NK.study,1),zeros(NK.study,1),zeros(NK.study,1),zeros(NK.study,1),zeros(NK.study,1),zeros(NK.study,1),zeros(NK.study,1),zeros(NK.study,1),zeros(NK.study,1),zeros(NK.study,1),zeros(NK.study,1),zeros(NK.study,1)]);
% cutConstraints = replace(cutConstraints,so,zeros(NK.DR,3));
% cutConstraints = replace(cutConstraints,si,zeros(NK.DR,NK.nonDR,3));
% cutConstraints = replace(cutConstraints,lc,zeros(NK.DR,3));

% u = dual(SPconstraints); % 虽然ei，gi为定值，但这里SPconstraint和generalconstraint的维度应该是一致的，因为没有只有ei，gi的约束

% obtain the coeeficient of the linear constraints (cut)
% 用recoverymodel里的，应该会保留变量的字符信息
% [model,recoverymodel] = export(generalConstraints);
% A = model.A; b = model.rhs;
% if diagnostics.problem ~= 0 
%     %unbounded or infeasible, （确定目标函数以及自变量约束，原问题不会unbounded，所以肯定是infeasible）
%     % relax the problem
%     nSlack = size(SPequations,1);
%     slack_SP = sdpvar(nSlack,1);
%     SPcons_relaxed = [
%         (SPequations - slack_SP <= 0):'SPcons_relaxed';
%         (55555555>=slack_SP >= 0):'slackCons';
% %         -10e6 <= lc <= 10e6;
%         -10e10 <= si <= 10e10;
%         ];
%     diag_relaxed = optimize(SPcons_relaxed,obj,opts); % this must be feasible? (or unbounded)
%     if diag_relaxed.problem ~= 0
%         error('still not feasible?');
%     else % consider it solved, obtain the dual vars
%         u = dual(SPcons_relaxed('SPcons_relaxed'));
%     end      
% [model,recoverymodel] = export(generalConstraints);
% A = model.A; b = model.rhs;    
if diagnostics.problem ~= 0 %unbounded or infeasible
    % 20200601 no need for this trouble
%     % obtain the dual problem


%     [model,recoverymodel] = export(SPconstraints,obj,opts);
%     dualVar = sdpvar(1,1032);
%     obj_dual = dualVar * model.rhs;
%     SPcons_dual = [
%         dualVar * model.A <= model.obj';
%         dualVar >= 0    ];

%     [SPcons_dual,obj_dual,primalConeVars,primalFreeVars,err,complexInfo] = dualize(SPconstraints,obj,0,0);
%     opts_dual = sdpsettings('verbose',2,'solver','gurobi','usex0',0,'debug',1,'gurobi.InfUnbdInfo',1,...
%     'savesolverinput',1,'savesolveroutput',1);
%     diagnostics_dual = optimize(SPcons_dual,-obj_dual,opts_dual);
%     % unbounded ray by gurobi
%     u = diagnostics_dual.solveroutput.result.unbdray;   
    % add feasibility cut 先保留变量，等会一起替换
%     addConstraint = u*A*recover(recoverymodel.used_variables) <= u*b;
    exitflag = 0;
    u = diagnostics.solveroutput.result.farkasdual;
    sol.cost_hat = 0;
elseif diagnostics.problem == 0 %feasible
    exitflag = 1;
    u = dual(SPconstraints);
    % add optimality cut
%     addConstraint = u*A*recover(recoverymodel.used_variables) <= u*b + t - f; % f is the optimal value of the objective function(feasible)
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