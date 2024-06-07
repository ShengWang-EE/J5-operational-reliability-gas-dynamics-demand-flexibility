function [MPformulation] = MP_IEGSdispatch_optimizer(addCutsCoefficient,exitflagSP,SP_EHresult,...
    dayahead_IEGSresult_basicLoad,dayaheadEHschedule,EHreserve,mpc,gtd,nx,dx,dt,NH,nEH,NK,iK,nodalPrice,Niter)
% 根据每小时日前运行的结果，作为暂态分析每小时及其间每个时间片段的初值 
% 注意：原来dayahead中所设定的电力、天然气的节点负荷中包含了ei，gi
for h = 1:NH
    terminalCondition{h} = setInitialConditionForGasSystem(dayahead_IEGSresult_basicLoad{h},gtd,nx,dx);
end
initialCondition = terminalCondition{k2h_h2k((iK.studyPeriod(1)-1),4,'k2h')}; % studyPeriod 开始时间点的前一个时间点日前运行的状态
% test: simplified:
% terminalCondition = terminalCondition{1};
%---------------------------
mpc.basePrs = 1e6;
exitflag = nan;
%% define named indices into data matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
[PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST, COST] = idx_cost;

% create (read-only) copies of individual fields for convenience\
mpopt = [];
[baseMVA, bus, gen, branch, gencost, Au, lbu, ubu, mpopt, ...
    N, fparm, H, Cw, z0, zl, zu, userfcn] = opf_args(mpc, mpopt);
[para] = initializeParameters2();
%% data dimensions
nb   = size(mpc.bus, 1);    %% number of buses
nl   = size(mpc.branch, 1); %% number of branches
ngen   = size(mpc.gen, 1);    %% number of dispatchable injections

%add GFU,LCe are included in the mpc
nGb  = size(mpc.Gbus,1); % number of gas bus
nGl  = size(mpc.Gline,1); % number of gas line
nGs  = size(mpc.Gsou,1); % number of gas source
nLCg = size(find(mpc.Gbus(:,3)~=0),1);
%%  define the decision variable
Prs = sdpvar(NK.studyPeriod,sum(nx+1));                                                    % gas pressure along the pipeline
Gf = sdpvar(NK.studyPeriod,sum(nx+1));                                                    % gas flow volume
Va = sdpvar(NK.studyPeriod,nb,'full'); 
Pg = sdpvar(NK.studyPeriod,ngen);
%
ei = sdpvar(NK.studyPeriod,nEH);
gi = sdpvar(NK.studyPeriod,nEH);
EHcost_hat = sdpvar(nEH,1);% 所有时间段总的cost
%% set up initial variables and bounds
nodalGasPrice = zeros(NK.studyPeriod,nEH);
counter = 0;
for k = iK.studyPeriod(1):iK.studyPeriod(end)
    % 判断k所在的小时
    h = k2h_h2k(k,4,'k2h');
    counter = counter +1 ;
    Pg0(:,counter) = dayahead_IEGSresult_basicLoad{h}.gen(:, PG) / baseMVA;
    Qg0(:,counter)   = dayahead_IEGSresult_basicLoad{h}.gen(:, QG) / baseMVA;
    Va0(:,counter)   = dayahead_IEGSresult_basicLoad{h}.bus(:, VA) * (pi/180);
    Vm0(:,counter)   = dayahead_IEGSresult_basicLoad{h}.bus(:, VM);
    V0(:,counter) = Vm0(:,counter) .* exp(1j*Va0(:,counter));
%     Vr0(:,k) = real(V0);
%     Vi0(:,k) = imag(V0);

    Prs0(:,counter) = cell2mat(terminalCondition{h}.P')'/mpc.basePrs;% bar
    Gf0(:,counter) =  cell2mat(terminalCondition{h}.Q')';%Mm3/day
    PGs0(:,counter) = dayahead_IEGSresult_basicLoad{h}.Gsou(:,5);% 这个不是变量，而是定值
    for i = 1:nEH
        ei0(i,counter) = dayaheadEHschedule{i}(h,1)/baseMVA;
        gi0(i,counter) = dayaheadEHschedule{i}(h,2)/200;
    end
    nodalGasPrice(counter,:) = nodalPrice.gas(h,mpc.EHlocation(:,1)');
end
%这两项初值等会根据上次SP结果赋予
EHcost_hat0 = zeros(nEH,1);
% if ~isempty(SP_EHresult)
%     for i = 1:nEH      
%         EHcost_hat0(i) = SP_EHresult{end,i}.EHcost_hat;
%     end
% end

Pgmin = gen(:, PMIN) / baseMVA;
Pgmax = gen(:, PMAX) / baseMVA;
Qgmin = gen(:, QMIN) / baseMVA;
Qgmax = gen(:, QMAX) / baseMVA;

Vclim = 1.1 * bus(:, VMAX);
Vrmin = -Vclim; Vrmax = Vclim; Vimin = -Vclim; Vimax = Vclim; 
refs = find(bus(:, BUS_TYPE) == REF);
Vau = Inf(nb, 1);       %% voltage angle limits
Val = -Vau;
Vau(refs) = Va0(refs);   %% voltage angle reference constraints
Val(refs) = Va0(refs);
Vml = bus(:, VMIN);
Vmu = bus(:, VMAX);

% Prsmin = min(mpc.Gbus(:,5)) * ones(sum(nx+1),1)/10;
% Prsmax = max(mpc.Gbus(:,6)) * ones(sum(nx+1),1)/10;
% Prsmin = Prs0*0.96;
% Prsmax = Prs0*1.04;
Prsmin = Prs0*0.95;
Prsmax = Prs0*1.05;

% Gfmin = -999*ones(sum(nx+1),1);
% Gfmax = 999*ones(sum(nx+1),1);
Gfmax = [];
for i = 1:nGl
    addGfmax = mpc.Gline(i,5) * ones(nx(i)+1,1);
    Gfmax = [Gfmax; addGfmax];
end
Gfmin = -Gfmax;

% PGsmin = mpc.Gsou(:,3);
% PGsmax = mpc.Gsou(:,4);
% LCgmin = zeros(nLCg,1);
% LCgmax = mpc.Gbus(mpc.Gbus(:,3)~=0,3).*0.2;  

% Pg0(Pg0<=Pgmin) = Pgmin(Pg0<=Pgmin);

eimin = zeros(NK.studyPeriod,nEH); gimin = zeros(NK.studyPeriod,nEH);
% 设置nonDR的ei不限制？还是限制ei只能变得更低，这样就不能shiftin了？
eimax = ei0'- EHreserve(iK.studyPeriod(1):iK.studyPeriod(end),:)/baseMVA; 
% eimax = inf * ones(NK.studyPeriod,nEH);
% eimax((iK.DR(1)-iK.studyPeriod(1)+1):(iK.DR(end)-iK.studyPeriod(1)+1),:) = ...
%     ei0(:,(iK.DR(1)-iK.studyPeriod(1)+1):(iK.DR(end)-iK.studyPeriod(1)+1))'- EHreserve(iK.DR(1):iK.DR(end),:)/baseMVA; 
gimax = inf * ones(NK.studyPeriod,nEH);

assign(Pg,Pg0'); assign(Va,Va0');
assign(Prs,Prs0'); assign(Gf,Gf0');
assign(ei,ei0'); assign(gi,gi0'); 
assign(EHcost_hat,EHcost_hat0);

%% define constraints
% upper and lower limits
PrsBoxCons = [Prsmin' <= Prs <= Prsmax'];
% PrsBoxCons = [repmat(Prsmin',NK.studyPeriod,1) <= Prs <= repmat(Prsmax',NK.studyPeriod,1)];
GfBoxCons = [repmat(Gfmin',NK.studyPeriod,1) <= Gf <= repmat(Gfmax',NK.studyPeriod,1)];
VaBoxCons = [repmat(Val',NK.studyPeriod,1) <= Va <= repmat(Vau',NK.studyPeriod,1)];
PgBoxCons = [repmat(Pgmin',NK.studyPeriod,1) <= Pg <= repmat(Pgmax',NK.studyPeriod,1)];




%% nonlinear constraints  
% preparations
il = find(branch(:, RATE_A) ~= 0 & branch(:, RATE_A) < 1e10);
[Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch);
cumnx = cumsum(nx+1);
pp = ones(nGl,2);
pp(2:nGl,1) = cumnx(1:(nGl-1))+1;
pp(1:nGl,2) = cumnx(1:nGl);
%
vv.i1.Pg = 1; vv.iN.Pg = ngen; vv.N.Pg = ngen;
vv.i1.Qg = vv.iN.Pg + 1; vv.iN.Qg = vv.iN.Pg + ngen; vv.N.Qg = ngen;
vv.i1.Va = vv.iN.Qg + 1; vv.iN.Va = vv.iN.Qg + nb; vv.N.Va = nb;
vv.i1.Vm = vv.iN.Va + 1; vv.iN.Vm = vv.iN.Va + nb; vv.N.Vm = nb;
vv.i1.Prs = vv.iN.Vm + 1; vv.iN.Prs = vv.iN.Vm + sum(nx+1); vv.N.Prs = sum(nx+1);
vv.i1.Gf = vv.iN.Prs + 1; vv.iN.Gf = vv.iN.Prs + sum(nx+1); vv.N.Gf = sum(nx+1);
vv.i1.PGs = vv.iN.Gf + 1; vv.iN.PGs = vv.iN.Gf + nGs; vv.N.PGs = nGs;
vv.i1.LCg = vv.iN.PGs + 1; vv.iN.LCg = vv.iN.PGs + nLCg; vv.N.LCg = nLCg;
vv.N.all = vv.N.Pg + vv.N.Qg + vv.N.Pg + vv.N.Va + vv.N.Vm + vv.N.Prs + ...
    vv.N.Gf + vv.N.PGs + vv.N.LCg;

% for i = 1:nGl % 第一次用的，接下来第二次第三次就用上一次求解结果的平均
    avrgQP = Gf0 ./ Prs0;
%     0.5 * (initialCondition.Q{i}./initialCondition.P{i} + terminalCondition.Q{i}./terminalCondition.P{i});
%     avrgQP{i} = repmat(avrgQP{i},NK.studyPeriod+1,1);
% end

% 
% ei_inbus = zeros(NK.studyPeriod,nb); gi_inbus = zeros(NK.studyPeriod,nGb);
gi_connect = zeros(nEH,nGb);ei_connect = zeros(nEH,nb);
for i = 1:nEH
    gi_connect(i,mpc.EHlocation(i,1)) = 1;
    ei_connect(i,mpc.EHlocation(i,2)) = 1;
end
gi_inbus = gi * gi_connect ;
ei_inbus = ei * ei_connect ;


%%
% objfcn = sum(slack);
% objfcn = [];
dynamicGasFlowContinuityCons = [optimalcontrol_dynamicGasflow_fcn_continuity(Prs,Gf,avrgQP,...
    dayahead_IEGSresult_basicLoad,gtd,initialCondition,para,vv,dx,nx,dt,NK.studyPeriod) ==0]:'dynamicGasFlowContinuityCons';
dynamicGasFlowMotionCons = [optimalcontrol_dynamicGasflow_fcn_motion(Prs,Gf,avrgQP,...
    dayahead_IEGSresult_basicLoad,gtd,initialCondition,para,vv,dx,nx,dt,NK.studyPeriod) <= 10]:'dynamicGasFlowMotionCons';


electricPowerBalanceConsDC = [optimalcontrol_power_balance_fcn2_dc(Va,Pg,ei_inbus, ...
    dayahead_IEGSresult_basicLoad, NK.studyPeriod,iK) == 0]:'electricPowerBalanceConsDC';
electricBranchFlowConsDC = [optimalcontrol_branch_flow_fcn_dc(Va, dayahead_IEGSresult_basicLoad, ...
    il,NK.studyPeriod,iK) <= 0]:'electricBranchFlowConsDC';
gasBalanceCons = [optimalcontrol_gas_balance_fcn_yalmip(Pg,Gf,PGs0',gi_inbus,pp,vv,...
    dayahead_IEGSresult_basicLoad,nx,NK.studyPeriod,iK)==0]:'gasBalanceCons';
gasPressureCons = [optimalcontrol_gas_pressure_fcn_yalmip(Prs,pp,dayahead_IEGSresult_basicLoad,NK.studyPeriod,iK)==0]:'gasPressureCons';

%%
constraints = [ 
                PrsBoxCons;
                GfBoxCons;
                VaBoxCons;
                PgBoxCons;
                electricPowerBalanceConsDC;
                electricBranchFlowConsDC;
                gasBalanceCons;
                gasPressureCons;
                dynamicGasFlowContinuityCons;  
                dynamicGasFlowMotionCons;
                %%
                eimin <= ei <= eimax;
                gimin <= gi <= gimax;
%                 addCutCons;
                EHcost_hat >= 0;  
%                 slack1 >= 0;
%                 slack2 >= 0;
%                 -slack2 <= Prs-Prs0' <= slack2;
                ];

IEGScost = objfcn_IEGSdispatch(Pg,mpc,NK.studyPeriod,nodalPrice);
giCost = sum(sum(nodalGasPrice .* gi));


objfcn = IEGScost + sum(EHcost_hat) + giCost;

options = sdpsettings('verbose',2,'solver','gurobi','usex0',0,'debug',1);


MPformulation.obj = objfcn;
MPformulation.cons = constraints;
MPformulation.opts = options;
MPformulation.x = {Prs,Gf,Va,Pg,ei,gi,EHcost_hat};

end

