function optimalControlResult = lookaheadContingencyManagement_fmincon(initialCondition,mpc,gtd,possibleResults, probability,dx,nx,dt,NK)
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
%%
prob = optimproblem('ObjectiveSense','minimize');
%%  define the decision variable
Prs = optimvar('Prs',NK,sum(nx+1));                                                    % gas pressure along the pipeline
Gf = optimvar('Gf',NK,sum(nx+1));                                                    % gas flow volume
PGs = optimvar('PGs',NK,nGs);                                                      % gas source production
LCg = optimvar('LCg',NK,nLCg);                                                      % gas load curtailment

Va = optimvar('Va',NK,nb);  
Vm = optimvar('Vm',NK,nb); 
Pg = optimvar('Pg',NK,ngen);                                                      % active power of a generator
Qg = optimvar('Qg',NK,ngen);  
%% set up initial variables and bounds
Pg0   = gen(:, PG) / baseMVA;
Qg0   = gen(:, QG) / baseMVA;
Va0   = bus(:, VA) * (pi/180);
Vm0   = bus(:, VM);

Prs0 = cell2mat(initialCondition.P')';
Gf0 =  cell2mat(initialCondition.Q')';
PGs0 = mpc.Gsou(:,2);
LCg0 = zeros(nLCg,1); 

% assign(Pg,repmat(Pg0',NK,1));
% assign(Qg,repmat(Qg0',NK,1));
% assign(Va,repmat(Va0',NK,1));
% assign(Vm,repmat(Vm0',NK,1));
% 
% assign(Prs,repmat(Prs0',NK,1));
% assign(Gf,repmat(Gf0',NK,1));
% assign(PGs,repmat(PGs0',NK,1));
% assign(LCg,repmat(LCg0',NK,1));


Pgmin = gen(:, PMIN) / baseMVA;
Pgmax = gen(:, PMAX) / baseMVA;
Qgmin = gen(:, QMIN) / baseMVA;
Qgmax = gen(:, QMAX) / baseMVA;

refs = find(bus(:, BUS_TYPE) == REF);
Vau = Inf(nb, 1);       %% voltage angle limits
Val = -Vau;
Vau(refs) = Va0(refs);   %% voltage angle reference constraints
Val(refs) = Va0(refs);
Vml = bus(:, VMIN);
Vmu = bus(:, VMAX);

% [Prsmin,Prsmax] = deal(zeros(sum(nx+1),1));

Prsmin = min(mpc.Gbus(:,5)) * ones(sum(nx+1),1);
Prsmax = max(mpc.Gbus(:,6)) * ones(sum(nx+1),1);
Gfmin = -999*ones(sum(nx+1),1);
Gfmax = 999*ones(sum(nx+1),1);
PGsmin = mpc.Gsou(:,3);
PGsmax = mpc.Gsou(:,4);
LCgmin = zeros(nLCg,1);
LCgmax = mpc.Gbus(mpc.Gbus(:,3)~=0,3).*0.1;  
 
%% define the objective function

%% define constraints
% upper and lower limits
Prs.LowerBound = repmat(Prsmin',NK,1); Prs.UpperBound = repmat(Prsmax',NK,1);
Gf.LowerBound = repmat(Gfmin',NK,1); Gf.UpperBound = repmat(Gfmax',NK,1);
PGs.LowerBound = repmat(PGsmin',NK,1); PGs.UpperBound = repmat(PGsmax',NK,1);
LCg.LowerBound = repmat(LCgmin',NK,1); LCg.UpperBound = repmat(LCgmax',NK,1);

Va.LowerBound = repmat(Val',NK,1); Va.UpperBound = repmat(Vau',NK,1);
Vm.LowerBound = repmat(Vml',NK,1); Vm.UpperBound = repmat(Vmu',NK,1);
Pg.LowerBound = repmat(Pgmin',NK,1); Pg.UpperBound = repmat(Pgmax',NK,1);
Qg.LowerBound = repmat(Qgmin',NK,1); Qg.UpperBound = repmat(Qgmax',NK,1);


%% nonlinear constraints  
% preparations
il = find(branch(:, RATE_A) ~= 0 & branch(:, RATE_A) < 1e10);
[Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch);
cumnx = cumsum(nx+1);
pp = ones(nGl,2);
pp(2:nGl,1) = cumnx(1:(nGl-1))+1;
pp(1:nGl,2) = cumnx(1:nGl);
%
fcn_PQmis = @(Va, Vm, Pg, Qg)optimalcontrol_power_balance_fcn(Va, Vm, Pg, Qg, mpc, Ybus, mpopt,NK);
fcn_flow = @(Va, Vm)optimalcontrol_branch_flow_fcn(Va, Vm, mpc, Yf(il, :), Yt(il, :), il, mpopt,NK);

fcn_Gmis = @(Pg,Prs,PGs,LCg)optimalcontrol_gas_balance_fcn(Pg,Prs,PGs,LCg,pp,mpc,NK);
fcn_PrsEqual = @(Prs)optimalcontrol_gas_pressure_fcn(Prs,pp,mpc,NK);
fcn_Gasflow = @(Prs,Gf)optimalcontrol_dynamicGasflow_fcn(Prs,Gf,mpc,gtd,para,dx,nx,dt,NK);

prob.Constraints.PQmis = optimalcontrol_power_balance_fcn(Va, Vm, Pg, Qg, mpc, Ybus, mpopt,NK) == 0;
% prob.Constraints.electricityFlow = fcn_flow;
% prob.Constraints.Gmis = fcn_Gmis;
% prob.Constraints.PrsEqual = fcn_PrsEqual;
% prob.Constraints.Gasflow = fcn_Gasflow;


objfcn = @(Pg,LCg,Prs,PGs)objfcn_LaCMS(Pg,LCg,Prs,PGs,possibleResults, probability);
prob.Objective = objfcn;

sol = solve(prob);

    


end

function J = objfcn_LaCMS(Pg,LCg,PrsRaw,PGs,possibleResults, probability)
sumGenLCeCost = 0; 
sumGasCurtailmentCost = 0; sumGasPurchasingCost = 0;
% obtain the nodal pressure
nGb = size(mpc.Gbus);
Prs = zeros(nGb,1);
p = mat2cell(PrsRaw(NK,:),nx);
for i = 1:nGl
    Prs(mpc.Gline(i,1)) = p{i}(1);
    Prs(mpc.Gline(i,2)) = p{i}(end);
end
    
for k = 1:NK
    sumGenLCeCost = sumGenLCeCost + sum(totcost(gencost, Pg(k,:)));
    sumGasCurtailmentCost = sumGasCurtailmentCost + sum(LCg(k,:)) * CDFg;
    sumGasPurchasingCost = sumGasPurchasingCost + PGs(k,:) * Gcost;
end
J_dynamic = sumGenLCeCost + sumGasCurtailmentCost + sumGasPurchasingCost;
J_endpoint = 0;
penaltyFactor = 0.1;
nScenario = size(probability,1);
for j = 1:nScenario
    for i = 1:nGb
        J_endpoint = J_endpoint + probability(j) * ((Prs(i)-possibleResults.Gbus(i,7))/possibleResults.Gbus(i,7))^2;
    end
end
J_endpoint = J_endpoint * penaltyFactor;

J = J_dynamic + J_endpoint;
end

