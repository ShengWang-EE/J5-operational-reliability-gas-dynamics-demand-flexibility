function [sol,addCutsCoefficient,solverTime,exitflag] = SP_EHschedule(ei_hat,gi_hat,EHpara,NK,iK)
% yalmip('clear');
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
    ];
energyConversionConstraints = [energyConversionEquations==0]:'energyConversionConstraints';

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
capacityConstraints = [capacityEquations <= 0]:'capacityConstraints';

SPconstraints = [scheduleConstraints;energyConversionConstraints;capacityConstraints];
%%
status = 0;
% check the feasibility of SP
% solve the optimality problem (SP)
optimalityObj = objfcn_SP(lc,si,CDF,NK,iK);
optimalityCons = [
     SPconstraints; % all the EH schedule constraints
    (ei == ei_hat):'eiBalance';
    (gi == gi_hat):'giBalance';
];    
opts = sdpsettings('verbose',2,'solver','gurobi','usex0',0,'debug',1,'gurobi.InfUnbdInfo',1,...
    'savesolverinput',1,'savesolveroutput',1);
diagnostics = optimize(optimalityCons,optimalityObj,opts);
if diagnostics.problem == 0 % feasible
    exitflag = 1;
    dual_ei = dual(optimalityCons('eiBalance'));
    dual_gi = dual(optimalityCons('giBalance'));
    EHcost_hat = value(optimalityObj);
else % infeasible, get other kind of dual
    exitflag = 0;
    % create slack vars
    slackE1 = sdpvar(NK.studyPeriod,1); slackE2 = sdpvar(NK.studyPeriod,1); 
    slackG1 = sdpvar(NK.studyPeriod,1); slackG2 = sdpvar(NK.studyPeriod,1);
  
    feasibilityCheckObj = sum(slackE1+slackE2+slackG1+slackG2);
    feasibilityCheckCons = [
        SPconstraints; % all the EH schedule constraints
        (ei + slackE1 == ei_hat + slackE2):'eiBalance';
        (gi + slackG1 == gi_hat + slackG2):'giBalance';
        [slackE1 slackE2 slackG1 slackG2] >= 0;
        ];    

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


function totalCost = objfcn_SP(lc,si,CDF0,NK,iK)
curtailmentCost = 0; shiftCost = 0;
for l = 1:3 % 3 energies
    if l == 1 % electricity
        CDF = CDF0;
    else
        CDF = 0.3 * CDF0;% ¥Ÿ πshift cooling load
    end
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