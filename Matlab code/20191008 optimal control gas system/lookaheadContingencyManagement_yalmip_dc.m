% function optimalControlResult = lookaheadContingencyManagement(initialCondition,initialSolutionPoint,mpc,gtd,possibleResults, probability,dx,nx,dt,NK)
% save f
clc
clear
load f
[mpc, gtd] = case24GEv3();
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
Prs = sdpvar(NK,sum(nx+1));                                                    % gas pressure along the pipeline
Gf = sdpvar(NK,sum(nx+1));                                                    % gas flow volume
PGs = sdpvar(NK,nGs);                                                      % gas source production
LCg = sdpvar(NK,nLCg);                                                      % gas load curtailment

Va = sdpvar(NK,nb,'full');  
% Vi = sdpvar(NK,nb,'full'); 
% V = sdpvar(NK,nb,'full','complex');
% V = sdpvar(NK,nb)+sqrt(-1)*sdpvar(NK,nb);
Pg = sdpvar(NK,ngen);                                                      % active power of a generator
% Qg = sdpvar(NK,ngen);  

%% set up initial variables and bounds
Pg0   = gen(:, PG) / baseMVA;
Qg0   = gen(:, QG) / baseMVA;
Va0   = bus(:, VA) * (pi/180)+0.01;
Vm0   = bus(:, VM);
V0 = Vm0 .* exp(1j*Va0);
Vr0 = real(V0);
Vi0 = imag(V0);

Prs0 = cell2mat(initialSolutionPoint.P')'/1e5;% bar
Gf0 =  cell2mat(initialSolutionPoint.Q')';%Mm3/day
PGs0 = mpc.Gsou(:,2);
LCg0 = zeros(nLCg,1); 

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



Prsmin = min(mpc.Gbus(:,5)) * ones(sum(nx+1),1);
Prsmax = max(mpc.Gbus(:,6)) * ones(sum(nx+1),1);
% Gfmin = -999*ones(sum(nx+1),1);
% Gfmax = 999*ones(sum(nx+1),1);
Gfmax = [];
for i = 1:nGl
    addGfmax = mpc.Gline(i,5) * ones(nx(i)+1,1);
    Gfmax = [Gfmax; addGfmax];
end
Gfmin = -Gfmax;

PGsmin = mpc.Gsou(:,3);
PGsmax = mpc.Gsou(:,4);
LCgmin = zeros(nLCg,1);
LCgmax = mpc.Gbus(mpc.Gbus(:,3)~=0,3).*0.2;  

Pg0(Pg0<=Pgmin) = Pgmin(Pg0<=Pgmin);
assign(Pg,repmat(Pg0',NK,1));
assign(Qg,repmat(Qg0',NK,1));
assign(Vr,repmat(Vr0',NK,1));
assign(Vi,repmat(Vi0',NK,1));
% assign(V,repmat(V0',NK,1));

assign(Prs,repmat(Prs0',NK,1));
assign(Gf,repmat(Gf0',NK,1));
assign(PGs,repmat(PGs0',NK,1));
assign(LCg,repmat(LCg0',NK,1)); 

 
%% define constraints
% upper and lower limits
PrsBoxCons = [repmat(Prsmin',NK,1) <= Prs <= repmat(Prsmax',NK,1)];
GfBoxCons = [repmat(Gfmin',NK,1) <= Gf <= repmat(Gfmax',NK,1)];
PGsBoxCons = [repmat(PGsmin',NK,1) <= PGs <= repmat(PGsmax',NK,1)];
LCgBoxCons = [repmat(LCgmin',NK,1) <= LCg <= repmat(LCgmax',NK,1)];

% VrBoxCons = [repmat(Vrmin',NK,1) <= Vr <= repmat(Vrmax',NK,1)];
% ViBoxCons = [repmat(Vimin',NK,1) <= Vi <= repmat(Vimax',NK,1)];
VUpperBoxCons = [(Vr.^2+Vi.^2) <= repmat(Vmu',NK,1).^2];
VLowerBoxCons = [(Vr.^2+Vi.^2) >= repmat(Vml',NK,1).^2 ];
% VRealUpperBoxCons = [real(V) <= repmat(Vmu',NK,1)];
% VRealLowerBoxCons = [repmat(Vml',NK,1) <= real(V)];
% VImagUpperBoxCons = [-repmat(Vmu',NK,1) <= imag(V) <= repmat(Vmu',NK,1)];
% VImagLowerBoxCons = [repmat(Vml',NK,1) <= imag(V)];
% VRealBoxCons = [repmat(Vml',NK,1).^2<= real(V).^2 <= repmat(Vmu',NK,1).^2];
% VImagBoxCons = [repmat(Vml',NK,1).^2<= imag(V).^2 <= repmat(Vmu',NK,1).^2];
% VBoxCons = [Vcons(V,Vmu,Vml,NK) <= 0];

PgBoxCons = [repmat(Pgmin',NK,1) <= Pg <= repmat(Pgmax',NK,1)];
QgBoxCons = [repmat(Qgmin',NK,1) <= Qg <= repmat(Qgmax',NK,1)];


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
%
% fcn_PQmis = @(V, Pg, Qg)optimalcontrol_power_balance_fcn2(V, Pg, Qg, mpc, Ybus, mpopt,NK);
% fcn_flow = @(V)optimalcontrol_branch_flow_fcn(V, mpc, Yf(il, :), Yt(il, :), il, mpopt,NK);
% 
% fcn_Gmis = @(Pg,Prs,PGs,LCg)optimalcontrol_gas_balance_fcn(Pg,Gf,PGs,LCg,pp,vv,mpc,nx,NK);
% fcn_PrsEqual = @(Prs)optimalcontrol_gas_pressure_fcn(Prs,pp,mpc,NK);
% fcn_Gasflow = @(Prs,Gf)optimalcontrol_dynamicGasflow_fcn(Prs,Gf,mpc,gtd,initialCondition,para,dx,nx,dt,NK);

electricPowerBalanceCons = [optimalcontrol_power_balance_fcn2(Vr,Vi, Pg, Qg, mpc, Ybus, mpopt,NK) == 0];
electricBranchFlowCons = [optimalcontrol_branch_flow_fcn(Vr,Vi, mpc, Yf(il, :), Yt(il, :), il, mpopt,NK) <= 0];
gasBalanceCons = [optimalcontrol_gas_balance_fcn_yalmip(Pg,Gf,PGs,LCg,pp,vv,mpc,nx,NK) == 0];
gasPressureCons = [optimalcontrol_gas_pressure_fcn_yalmip(Prs,pp,mpc,NK) == 0];
dynamicGasFlowCons = [optimalcontrol_dynamicGasflow_fcn_yalmip(Prs,Gf,mpc,gtd,initialCondition,para,vv,dx,nx,dt,NK) == 0];



constraints = [ 
                PrsBoxCons;
                GfBoxCons;
                PGsBoxCons;
                LCgBoxCons;
                VUpperBoxCons;
                VLowerBoxCons;
                PgBoxCons;
                QgBoxCons;
                electricPowerBalanceCons;
%                 electricBranchFlowCons;
%                 gasBalanceCons;
%                 gasPressureCons;
%                 dynamicGasFlowCons;               
                ];

    
objfcn = objfcn_LaCMS_yalmip(Pg,Prs,PGs,LCg,mpc,vv,possibleResults, probability,nx,NK,baseMVA);
options = sdpsettings('solver','IPOPT','verbose',2,'usex0',1,'debug',1);
% options.ipopt.max_iter = 1e8;
% options.ipopt.max_cpu_time = 1e4;
% options.ipopt.display = 'iter';
% options = sdpsettings('debug','1');
diagnostics  = optimize(constraints,[],options);
% diagnostics  = optimize(constraints,[],options);
%% results
[Pg1,Qg1,Vr1,Vi1,PrsRaw,GfRaw,PGs1,LCg1,fval] = deal(value(Pg),value(Qg),value(Vr),value(Vi),value(Prs),value(Gf),value(PGs),value(LCg),value(objfcn)); 
PrsPipe = mat2cell(PrsRaw,NK,nx+1); Gf1 = mat2cell(GfRaw,NK,nx+1);
PrsNodal = zeros(nGb,1);
for i = 1:nGl
    PrsNodal(mpc.Gline(i,1)) = PrsPipe{i}(1);
    PrsNodal(mpc.Gline(i,2)) = PrsPipe{i}(end);
end
[J_dynamic,J_endpoint,genLCeCost,gasCurtailmentCost,gasPurchasingCost] = LaCMScost(Pg1,LCg1,PGs1,PrsNodal,mpc,possibleResults,probability,NK);
optimalControlResults = struct('Pg',Pg1,'Qg',Qg1,'Vr',Vr1,'Vi',Vi1,'PrsPipe',PrsPipe,'PrsNodal',PrsNodal,...
    'Gf',Gf1,'PGs',PGs1,'LCg',LCg1,'f',fval,'J_dynamic',J_dynamic,'J_endpoint',J_endpoint,'genLCeCost',genLCeCost,...
    'gasCurtailmentCost',gasCurtailmentCost,'gasPurchasingCost',gasPurchasingCost);

