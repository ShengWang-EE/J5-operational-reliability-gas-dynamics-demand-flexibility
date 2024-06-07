function om = GEopc_setup(mpc,gtd,mpopt,steadyResults,NK,initialSolutionPoint,initialCondition,dx,nx,dt,possibleResults,probability)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
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
[baseMVA, bus, gen, branch, gencost, Au, lbu, ubu, mpopt, ...
    N, fparm, H, Cw, z0, zl, zu, userfcn] = opf_args(mpc);
[para] = initializeParameters2();
%% data dimensions
nb   = size(mpc.bus, 1);    %% number of buses
nl   = size(mpc.branch, 1); %% number of branches
ng   = size(mpc.gen, 1);    %% number of dispatchable injections

%add GFU,LCe are included in the mpc
nGb  = size(mpc.Gbus,1); % number of gas bus
nGl  = size(mpc.Gline,1); % number of gas line
nGs  = size(mpc.Gsou,1); % number of gas source
nLCg = size(find(mpc.Gbus(:,3)~=0),1);
%% set up initial variables and bounds
Pg0   = repmat(steadyResults.gen(:, PG) / baseMVA,NK,1);
Qg0   = repmat(steadyResults.gen(:, QG) / baseMVA,NK,1);
Va0   = repmat(steadyResults.bus(:, VA) * (pi/180) / baseMVA,NK,1);
Vm0   = repmat(steadyResults.bus(:, VM),NK,1);
V0 = Vm0 .* exp(1j*Va0);
Vr0 = real(V0);
Vi0 = imag(V0);
Prs0 = repmat(cell2mat(initialSolutionPoint.P')'/1e5,NK,1);% bar
Gf0 =  repmat(cell2mat(initialSolutionPoint.Q')',NK,1);%Mm3/day
PGs0 = repmat(steadyResults.Gsou(:,5),NK,1);
LCg0 = repmat(steadyResults.Gbus(steadyResults.Gbus(:,3)>0,10),NK,1); 

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

% 
Pgmin = repmat(Pgmin,NK,1);     Pgmax = repmat(Pgmax,NK,1);
Qgmin = repmat(Qgmin,NK,1);     Qgmax = repmat(Qgmax,NK,1);
Val = repmat(Val,NK,1);     Vau = repmat(Vau,NK,1);
Vml = repmat(Vml,NK,1);     Vmu = repmat(Vmu,NK,1);
Prsmin = repmat(Prsmin,NK,1);     Prsmax = repmat(Prsmax,NK,1);
Gfmin = repmat(Gfmin,NK,1);     Gfmax = repmat(Gfmax,NK,1);
PGsmin = repmat(PGsmin,NK,1);     PGsmax = repmat(PGsmax,NK,1);
LCgmin = repmat(LCgmin,NK,1);     LCgmax = repmat(LCgmax,NK,1); 

% preparations
il = find(branch(:, RATE_A) ~= 0 & branch(:, RATE_A) < 1e10);
[Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch);
cumnx = cumsum(nx+1);
pp = ones(nGl,2);
pp(2:nGl,1) = cumnx(1:(nGl-1))+1;
pp(1:nGl,2) = cumnx(1:nGl);
%% define variables and constraints
fcn_Emis = @(x)optimalcontrol_power_balance_fcn_mips(x, mpc, Ybus, mpopt,NK);
% hess_Emis = @(x, lam)opf_power_balance_hess(x, lam, mpc, Ybus, mpopt);
% electricity branch limit    
fcn_Eflow = @(x)optimalcontrol_branch_flow_fcn_mips(x, mpc, Yf(il, :), Yt(il, :), il, mpopt,NK);
% hess_Eflow = @(x, lam)opf_branch_flow_hess(x, lam, mpc, Yf(il, :), Yt(il, :), il, mpopt);
% gas nodal balance
fcn_Gmis = @(x)optimalcontrol_gas_balance_fcn_mips(x,pp,mpc,nx,NK);
% pressure equation
fcn_PrsEqual = @(x)optimalcontrol_gas_pressure_fcn_mips(x,pp,mpc,NK,nx);
% 
fcn_gasDynamic = @(x)optimalcontrol_dynamicGasflow_fcn_mips(x,mpc,gtd,initialCondition,para,dx,nx,dt,NK);
% objective function
Jcost = @(x)objfcn_LaCMS_mips(x,mpc,possibleResults, probability,nx,NK,baseMVA);

%% construct opf model
om = opf_model(mpc);
% user data
[Apqh, ubpqh, Apql, ubpql, Apqdata] = makeApq(baseMVA, gen);
  iang  = [];
om.userdata.Apqdata = Apqdata;
om.userdata.iang = iang;
% optimization variables
om.add_var('Va', NK*nb, Va0, Val, Vau);
om.add_var('Vm', NK*nb, Vm0, Vml, Vmu);
om.add_var('Pg', NK*ng, Pg0, Pgmin, Pgmax);
om.add_var('Qg', NK*ng, Qg0, Qgmin, Qgmax);
om.add_var('Prs',NK*sum(nx+1),Prs0,Prsmin,Prsmax);
om.add_var('Gf',NK*sum(nx+1),Prs0,Gfmin,Gfmax);
om.add_var('PGs',NK*nGs,PGs0,PGsmin,PGsmax);
om.add_var('LCg',NK*nLCg,LCg0,LCgmin,LCgmax);
%add constraints
JcostVar = {'Pg','Prs','PGs','LCg'};
nodal_balance_vars = {'Va', 'Vm', 'Pg', 'Qg'};
flow_lim_vars = {'Va', 'Vm'};
gas_nodal_balance_vars = {'Pg','Gf','PGs','LCg'};
PrsEqualVars = {'Prs'};
gasDynamicVars = {'Prs','Gf'};
% electricity nodal balance      
Emis_cons = {'Pmis', 'Qmis'};

om.add_nln_constraint(Emis_cons, NK*[nb;nb], 1, fcn_Emis, [], nodal_balance_vars); % active and reactive nodal power balance
om.add_nln_constraint({'Sf', 'St'}, NK*[nl;nl], 0, fcn_Eflow, [], flow_lim_vars);% branch limits
% om.add_nln_constraint('Gmis',NK*nGb,1,fcn_Gmis,[],gas_nodal_balance_vars);%节点气守恒
% om.add_nln_constraint('PrsEqual',NK*(2*nGl-nGb),1,fcn_PrsEqual,[],PrsEqualVars);%(0表示不等式约束）
% om.add_nln_constraint('gasDynamic',2*NK*sum(nx),1,fcn_gasDynamic,[],gasDynamicVars);

om.add_nln_cost('Jcost', 1, Jcost,JcostVar);

% [vv, ~, ~, ~] = om.get_idx();


end

