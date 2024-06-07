function [SPformulation] = SP_EHschedule_lp_optimizer(ei_hat,gi_hat,EHpara,NK,iK)
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
    [ei0, gi0, eee0, ee30, e130, e1h0, gg10, gg20, h1h0, h140, h240, h2h0, c3c0, c4c0]);
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

%%
status = 0;
% solve
obj = objfcn_SP(lc,si,CDF,NK,iK);
% obj = [];
opts = sdpsettings('verbose',2,'solver','gurobi','usex0',0,'debug',1,'gurobi.InfUnbdInfo',1,...
    'savesolverinput',1,'savesolveroutput',1);

SPformulation.obj = obj;
SPformulation.generalEquations = generalEquations;
SPformulation.opts = opts;
SPformulation.x = {ei,gi,eee,ee3,e13,e1h,gg1,gg2,h1h,h14,h24,h2h,c3c,c4c,so,si,lc};

end

