function [scenarioConstraints,scenarioCost] = generateFormulationsForEachScenario(Prs,Gf,PGs,LCg,Va,Pg,ei,gi,PrsVioPercentge,... % IEGS vars
    eeeAll,ee3All,e13All,e1hAll,gg1All,gg2All,h3hAll,h1hAll,h14All,h24All,h2hAll,c3cAll,c4cAll,lcAll,... %EH vars
    PGsmax,Pgmax,initialP,initialQ,electricityLoad,gasLoad,EHloadAll,initialCondition,... % parameter vars
    time0,mpc,nodalEHpara,gtd,dx,nx,dt,NL,NK,nEH)
%% initialization

[para] = initializeParameters2();
avrgQP = initialP./initialQ;
% define named indices into data matrices
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
% data dimensions
nb   = size(mpc.bus, 1);    %% number of buses
nl   = size(mpc.branch, 1); %% number of branches
ngen   = size(mpc.gen, 1);    %% number of dispatchable injections

%add GFU,LCe are included in the mpc
nGb  = size(mpc.Gbus,1); % number of gas bus
nGl  = size(mpc.Gline,1); % number of gas line
nGs  = size(mpc.Gsou,1); % number of gas source
nLCg = size(find(mpc.Gbus(:,3)~=0),1);
%% IEGS upper and lower limits
Pgmin = mpc.gen(:, PMIN) / baseMVA *0; %Pgmin is set to zero

refs = find(mpc.bus(:, BUS_TYPE) == REF);
Vau = Inf(nb, 1);       %% voltage angle limits
Val = -Vau;
Vau(refs) = 1;   %% voltage angle reference constraints
Val(refs) = 1;

Prsmin = initialP*(1-PrsVioPercentge); Prsmax = initialP*(1+PrsVioPercentge); %bar/10
% Prs约束是有意义的，因为第二调度时段的prs约束导致第一时段不能用太多prs
Gfmax = [];
for i = 1:nGl
    addGfmax = mpc.Gline(i,5) * ones(nx(i)+1,1);
    Gfmax = [Gfmax; addGfmax];
end
Gfmin = -Gfmax;
PGsmin = mpc.Gsou(:,3)*0; %test: PGs下限设为0

LCgmin = zeros(nLCg,1);
LCgmax = mpc.Gbus(mpc.Gbus(:,3)~=0,3).*1;  
eimin = zeros(NK,nEH); gimin = zeros(NK,nEH);
eimax = inf * ones(NK,nEH); % no upper limits 
gimax = inf * ones(NK,nEH);
% ------------------
PrsBoxCons = [Prsmin <= Prs <= Prsmax];
GfBoxCons = [repmat(Gfmin',NK,1) <= Gf <= repmat(Gfmax',NK,1)];
PGsBoxCons = [repmat(PGsmin',NK,1) <= PGs <= repmat(PGsmax,NK,1)];
LCgBoxCons = [repmat(LCgmin',NK,1) <= LCg <= repmat(LCgmax',NK,1)];
VaBoxCons = [repmat(Val',NK,1) <= Va <= repmat(Vau',NK,1)];
PgBoxCons = [repmat(Pgmin',NK,1) <= Pg <= repmat(Pgmax,NK,1)];
eiBoxCons = [eimin <= ei <= eimax];
giBoxCons = [gimin <= gi <= gimax];
%% IEGS other limtis
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
% 
% ei_inbus = zeros(NK.studyPeriod,nb); gi_inbus = zeros(NK.studyPeriod,nGb);
gi_connect = zeros(nEH,nGb);ei_connect = zeros(nEH,nb);
for i = 1:nEH
    gi_connect(i,mpc.EHlocation(i,1)) = 1;
    ei_connect(i,mpc.EHlocation(i,2)) = 1;
end
gi_inbus = gi * gi_connect ;
ei_inbus = ei * ei_connect ;

% ----- other constraints -----
dynamicGasFlowContinuityCons = [optimalcontrol_dynamicGasflow_fcn_continuity3(Prs,Gf,avrgQP,...
    mpc,gtd,initialCondition,para,vv,dx,nx,dt,NK) ==0]:'dynamicGasFlowContinuityCons';
dynamicGasFlowMotionCons = [optimalcontrol_dynamicGasflow_fcn_motion3(Prs,Gf,avrgQP,...
    mpc,gtd,initialCondition,time0,para,vv,dx,nx,dt,NK) == 0]:'dynamicGasFlowMotionCons';
electricPowerBalanceConsDC = [optimalcontrol_power_balance_fcn_dc3(Va,Pg,ei_inbus, electricityLoad,...
    mpc, NK) == 0]:'electricPowerBalanceConsDC';
electricBranchFlowConsDC = [optimalcontrol_branch_flow_fcn_dc3(Va, mpc, ...
    il,NK) <= 0]:'electricBranchFlowConsDC';
gasBalanceCons = [optimalcontrol_gas_balance_fcn_yalmip3(Pg,Gf,PGs,LCg,gi_inbus,gasLoad,pp,vv,...
    mpc,nx,NK)==0]:'gasBalanceCons';
gasPressureCons = [optimalcontrol_gas_pressure_fcn_yalmip3(Prs,pp,mpc,NK)==0]:'gasPressureCons';
%% EH constraints
EHconstraints = []; EHobjective = 0;
for i = 1:nEH
    eiNodal = ei(:,i)*baseMVA; % MW
    giNodal = gi(:,i) * 200; % MW
    eee = eeeAll(:,i);ee3 = ee3All(:,i);e13 = e13All(:,i);e1h = e1hAll(:,i);gg1 = gg1All(:,i);
    gg2 = gg2All(:,i);h3h = h3hAll(:,i);h1h = h1hAll(:,i); h14 = h14All(:,i); h24 = h24All(:,i);
    h2h = h2hAll(:,i); c3c = c3cAll(:,i); c4c = c4cAll(:,i); lc = reshape(lcAll(:,i,:),[NK,1,NL]);
    EHload = reshape(EHloadAll(:,i,:),[NK,1,NL]); % MW
    [EHconstraintsEach{i},EHobjectiveEach{i}] = EHscheduleFormulation(eiNodal,giNodal,eee,ee3,e13,e1h,gg1,gg2,h3h,h1h,h14,h24,h2h,c3c,c4c,lc...
        ,nodalEHpara{i},EHload,NK);
    EHconstraints = [EHconstraints,EHconstraintsEach{i}];
    EHobjective = EHobjective + EHobjectiveEach{i};
end

%% formulation
IEGSconstraints = [ 
                PrsBoxCons;
                GfBoxCons;
                PGsBoxCons;
                LCgBoxCons;
                VaBoxCons;
                PgBoxCons;
                electricPowerBalanceConsDC;
                electricBranchFlowConsDC;
                gasBalanceCons;
                gasPressureCons;
                dynamicGasFlowContinuityCons;  
                dynamicGasFlowMotionCons;
                %%
                eiBoxCons;
                giBoxCons;
                ];
linepackCost = ((sum(gasLoad)-sum(LCg)) - sum(PGs)) * max(mpc.Gcost)*2; %linepack补丁，防止过度消耗linepack
scenarioConstraints = [IEGSconstraints;EHconstraints];
scenarioCost = objfcn_IEGSdispatch(Pg,PGs,LCg,mpc,NK) + EHobjective +linepackCost;

    
end