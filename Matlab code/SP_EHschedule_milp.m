function [sol,addCutsCoefficient,solverTime,exitflag] = SP_EHschedule_milp(ei_hat,gi_hat,EHpara,NK,iK)
bigM = 100;
yalmip('clear');
% data
capacityCell = num2cell(EHpara.capacity);
efficiencyCell = num2cell(EHpara.efficiency);
scheduleCell = num2cell(EHpara.schedule);
[HA,EA,HB,EB,HC,EC,HD,ED,HO2MAX,HO2MIN,HO3MAX,HO3MIN,CO3MAX,CO3MIN,CO4MAX,CO4MIN] = deal(capacityCell{:});
[eta1e,eta1h,eta2,COP3hSummer,COP3hWinter,COP3c,COP4] = deal(efficiencyCell{:});
[alpha,beta,CDF] = deal(scheduleCell{:});
% test
% beta = 1;
%
load0 = EHpara.load;
COP3h = COP3hSummer;
gamma = 0;%summer,EHP cooling
NL = 3;

iK.nonDR = iK.studyPeriod; iK.nonDR(ismember(iK.studyPeriod,iK.DR)) = [];
%% vars
% first strategy:
[ei] = sdpvar(NK.studyPeriod,1);[gi] = sdpvar(NK.studyPeriod,1);[eee] = sdpvar(NK.studyPeriod,1);
[ee3] = sdpvar(NK.studyPeriod,1);[e13] = sdpvar(NK.studyPeriod,1);[e1h] = sdpvar(NK.studyPeriod,1);
[gg1] = sdpvar(NK.studyPeriod,1);[gg2] = sdpvar(NK.studyPeriod,1);[h3h] = sdpvar(NK.studyPeriod,1);
[h1h] = sdpvar(NK.studyPeriod,1);[h14] = sdpvar(NK.studyPeriod,1);[h24] = sdpvar(NK.studyPeriod,1);
[h2h] = sdpvar(NK.studyPeriod,1);[c3c] = sdpvar(NK.studyPeriod,1);[c4c] = sdpvar(NK.studyPeriod,1);
% second strategy
so = sdpvar(NK.DR,3); % shift out for DR period for three loads
si = sdpvar(NK.DR,NK.nonDR,NL); % shift in quantity for each shiftout period and shiftin period
sw = intvar(NK.DR,NK.nonDR,NL); % weather the si has been deployed
kai = sdpvar(NK.DR,NK.nonDR,NL); % axulary variable
% third strategy
lc = sdpvar(NK.DR,3); % load curtailment for each DR period
% start value
[ei0, gi0, eee0, ee30, e130, e1h0, gg10, gg20, h3h0, h1h0, h140, h240, h2h0, c3c0, c4c0] = deal(zeros(NK.studyPeriod,1));
so0 = zeros(NK.DR,3); % shift out for DR period for three loads
si0 = zeros(NK.DR,NK.nonDR,3); % shift in quantity for each shiftout period and shiftin period
sw0 = zeros(NK.DR,NK.nonDR,3);
kai0 = zeros(NK.DR,NK.nonDR,3);
lc0 = zeros(NK.DR,3); % load curtailment for each DR period
lc1 = load0(:,iK.DR)';
% some constants
h3h = h3h0;

assign([ei, gi, eee, ee3, e13, e1h, gg1, gg2, h1h, h14, h24, h2h, c3c, c4c],...
    [ei_hat, gi_hat, eee0, ee30, e130, e1h0, gg10, gg20, h1h0, h140, h240, h2h0, c3c0, c4c0]);
assign(so,so0); assign(si,si0); assign(sw,sw0);assign(kai,kai0);assign(lc,lc0);
%%
% calculate new load curve 一个数组变量要先定义变量项再定义常数项
% newLoad = load0; newLoad(:,iK.DR) = newLoad(:,iK.DR) - so' - lc';
% newLoad(:,iK.nonDR) = newLoad(:,iK.nonDR) + reshape(sum(kai,1),size(kai,2),size(kai,3))';

newLoad(:,iK.DR) = load0(:,iK.DR) - so' - lc';
newLoad(:,1:(iK.DR(1)-1)) = load0(:,1:(iK.DR(1)-1));
newLoad = [newLoad load0(:,(iK.DR(end)+1):NK.all)];
newLoad(:,iK.nonDR) = newLoad(:,iK.nonDR) + reshape(sum(kai,1),size(kai,2),size(kai,3))';

[el,hl,cl] = deal(newLoad(1,iK.studyPeriod)',newLoad(2,iK.studyPeriod)',newLoad(3,iK.studyPeriod)');
% schedule capacity
scheduleEquations = [
    -so;
    so - alpha * load0(:,iK.DR)';% limits for shiftout capability
    -lc;
    lc - beta * load0(:,iK.DR)'; % limits for curtailment
    so - reshape(sum(kai,2),size(kai,1),size(kai,3)) ; % shiftout should equal shiftin (here relax the equalty constraints)
    -reshape(si,[],3);
    -reshape(sw,[],3);
    reshape(sw,[],3)-1;
    ];
scheduleEquations = reshape(scheduleEquations,[],1);
% reformulation constraints(eliminate the bilinears)
reformulationEquations = [
    -kai;
    kai-sw*bigM;
    kai-si;
    si - bigM*(1-sw) - kai;
    ];

reformulationEquations = reshape(reformulationEquations,[],1);
% reformulationEquations = [];
% energy conversion constraints (automatically formulated for each period)
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

% in this general form, MP vars (ei,gi) is considered vars. in this SP,
% replace them as constants
generalEquations = [scheduleEquations; reformulationEquations; energyConversionEquations;  capacityEquations;]; 
generalConstraints = [generalEquations <= 0];
energyConversionEquations1 = replace(energyConversionEquations,[ei,gi],[ei_hat,gi_hat]);
SPequations = [scheduleEquations;  reformulationEquations; energyConversionEquations1;  capacityEquations;];
SPconstraints = [SPequations <= 0];
%% check the feasibility of the relaxed SP
% solve
obj = objfcn_SP(lc,si,CDF,NK,iK);
optsRelaxed = sdpsettings('verbose',2,'solver','gurobi','usex0',0,'debug',1,'gurobi.InfUnbdInfo',1,...
    'savesolverinput',1,'savesolveroutput',1,'relax',2);
diagnostics = optimize(SPconstraints,obj,optsRelaxed);
% get the solution of the relaxed SP for further use in LP cut generation
[SPmodel,SPrecoverymodel] = export(SPconstraints,obj,optsRelaxed);
x_hat_relaxed = value(recover(SPrecoverymodel.used_variables));
sw_relaxed = value(sw);
% [generalmodel,generalrecoverymodel] = export(generalConstraints,[],optsRelaxed);
useLP = 0;
if diagnostics.problem ~= 0 %unbounded or infeasible
    exitflag = 0;
    lambda = diagnostics.solveroutput.result.farkasdual;
elseif diagnostics.problem == 0 %feasible
    % solve energy substitution SP
    energySubSPCons = [SPconstraints;
                        so == so0;
                        si == si0;
                        sw == sw0;
                        lc == lc0;];
    diagnostics = optimize(energySubSPCons,0,optsRelaxed);
    if diagnostics.problem == 0 % energy substitution without shifting or curtailing is feasible
        %add optimality cut
        exitflag = 1;
        lambda = dual(SPconstraints);
    elseif diagnostics.problem ~= 0 % means we need shifting or curtailing
        % lift-and-project
        addLPcutEquation = [];
        for r1 = 1:size(sw_relaxed,1)
            for r2 = 1:size(sw_relaxed,2)
                for r3 = 1:size(sw_relaxed,3)
                    if value(sw_relaxed(r1,r2,r3)) ~= fix(value(sw_relaxed(r1,r2,r3)))
                        useLP = 1;
                        % index of the integer varible, the constraint is irrelevant
                        r = getvariables(sw(r1,r2,r3));
                        % L&P cut generation
                        [alpha,beta0] = LPcutGeneration(x_hat_relaxed,-SPmodel.A,-SPmodel.rhs,r);
                        % cut
                        addLPcutEquation = [
                                    addLPcutEquation;
                                    -alpha'*recover(SPrecoverymodel.used_variables) + beta0;
                                    ];
                    end
                end
            end
        end
        % solve the relaxed SP with LP cut attached
        addLPcutCons = [addLPcutEquation <= 0];
        SPconstraintsWithLP = [SPconstraints;addLPcutCons];
        diagnostics = optimize(SPconstraintsWithLP,obj,optsRelaxed);
        if diagnostics.problem == 0 % feasible
            exitflag = 1;
            lambda = dual(SPconstraintsWithLP);
        elseif diagnostics.problem ~= 0 % infeasible
            exitflag = 0;
            lambda = diagnostics.solveroutput.result.farkasdual;
        end       
    end
end
% assign values according to the exitflag
if exitflag==1 %
    % get solution for SP
    [sol.eee, sol.ee3, sol.e13, sol.e1h, sol.gg1, sol.gg2, ...
    sol.h3h, sol.h1h, sol.h14, sol.h24, sol.h2h, sol.c3c, sol.c4c,...
    sol.so, sol.si,sol.kai,sol.sw, sol.lc] = deal ...
    (value(eee),value(ee3),value(e13),value(e1h),value(gg1),...
    value(gg2),value(h3h),value(h1h),value(h14),value(h24),value(h2h),value(c3c),value(c4c),...
    value(so), value(si), value(kai), value(sw), value(lc));
    sol.fval = value(obj);
    sol.cost_hat = sol.fval;
elseif exitflag == 0
    sol.cost_hat = 0;
end
%% generate benders cut
if useLP == 1 % lift and project is used, so the constraint numbers and lambda are different
    addCutsEquation = lambda'*SPgeneralEquationsWithLP;
elseif useLP == 0 % lift and project are not used
    addCutsEquation = [lambda'*generalEquations];
end
% replace x(subproblem variables) with 0, only left with ei,gi
addCutsEquation = replace(addCutsEquation,...
    [eee, ee3, e13, e1h, gg1, gg2, h1h, h14, h24, h2h, c3c, c4c], ...
    [zeros(NK.studyPeriod,1),zeros(NK.studyPeriod,1),zeros(NK.studyPeriod,1),zeros(NK.studyPeriod,1),...
     zeros(NK.studyPeriod,1),zeros(NK.studyPeriod,1),zeros(NK.studyPeriod,1),zeros(NK.studyPeriod,1),...
     zeros(NK.studyPeriod,1),zeros(NK.studyPeriod,1),zeros(NK.studyPeriod,1),zeros(NK.studyPeriod,1)]);
addCutsEquation = replace(addCutsEquation,so,zeros(NK.DR,3));
% cannot replace multi-dimenssional sdpvars at once
addCutsEquation = replace(addCutsEquation,reshape(si,[],3),reshape(zeros(NK.DR,NK.nonDR,3),[],3));
addCutsEquation = replace(addCutsEquation,reshape(kai,[],3),reshape(zeros(NK.DR,NK.nonDR,3),[],3));
addCutsEquation = replace(addCutsEquation,lc,zeros(NK.DR,3));
% 如果直接export，只能输出系数不为0的项的系数，这里做一点小处理: 
% 找一个几乎不可能出现的系数，手动加上，得到系数矩阵后再去掉
% coefficient of ei
addCutsEquation_wrt_ei = addCutsEquation + pi*sum(ei);
addCutsEquation_wrt_ei = replace(addCutsEquation_wrt_ei,gi,zeros(size(gi)));
addCut = [addCutsEquation_wrt_ei <=0];
[model,recoverymodel] = export(addCut,[],optsRelaxed);
A.ei = model.A - pi; rhs = model.rhs; % rhs 是一样的
% coefficient of gi
addCutsEquation_wrt_gi = addCutsEquation + pi*sum(gi);
addCutsEquation_wrt_gi = replace(addCutsEquation_wrt_gi,ei,zeros(size(ei)));
addCut = [addCutsEquation_wrt_gi <=0];
[model,recoverymodel] = export(addCut,[],optsRelaxed);
A.gi = model.A - pi;
%
[addCutsCoefficient.ei.A, addCutsCoefficient.rhs, ...
    addCutsCoefficient.gi.A] = ...
    deal(A.ei, rhs, A.gi);
solverTime = diagnostics.solvertime;
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