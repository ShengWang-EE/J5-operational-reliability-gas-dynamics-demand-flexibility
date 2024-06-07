function optimalControlResults = lookaheadContingencyManagement_fminconv2(initialCondition,initialSolutionPoint,mpc,gtd,possibleResults, probability,dx,nx,dt,NK)
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

%% set up initial variables and bounds
Pg0   = gen(:, PG) / baseMVA;
Qg0   = gen(:, QG) / baseMVA;
Va0   = bus(:, VA) * (pi/180);
Vm0   = bus(:, VM);

Prs0 = cell2mat(initialSolutionPoint.P')'/1e5;% bar
Gf0 =  cell2mat(initialSolutionPoint.Q')';%Mm3/day
PGs0 = mpc.Gsou(:,2);
LCg0 = zeros(nLCg,1); 

x0 = repmat([Pg0;Qg0;Va0;Vm0;Prs0;Gf0;PGs0;LCg0],NK,1);

Pgmin = gen(:, PMIN) / baseMVA;
Pgmax = gen(:, PMAX) / baseMVA;
Qgmin = gen(:, QMIN) / baseMVA;
Qgmax = gen(:, QMAX) / baseMVA;

refs = find(bus(:, BUS_TYPE) == REF);
Vamax = Inf(nb, 1);       %% voltage angle limits
Vamin = -Vamax;
% Vamax(refs) = Va0(refs);   %% voltage angle reference constraints
% Vamin(refs) = Va0(refs);
Vmmin = bus(:, VMIN);
Vmmax = bus(:, VMAX);

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
LCgmax = mpc.Gbus(mpc.Gbus(:,3)~=0,3).*0.1;  
 
LB = repmat([Pgmin;Qgmin;Vamin;Vmmin;Prsmin;Gfmin;PGsmin;LCgmin],NK,1);
UB = repmat([Pgmax;Qgmax;Vamax;Vmmax;Prsmax;Gfmax;PGsmax;LCgmax],NK,1);
%% define the objective function
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


% preparations
il = find(branch(:, RATE_A) ~= 0 & branch(:, RATE_A) < 1e10);
[Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch);
cumnx = cumsum(nx+1);
pp = ones(nGl,2);
pp(2:nGl,1) = cumnx(1:(nGl-1))+1;
pp(1:nGl,2) = cumnx(1:nGl);
%% nonlinear constraints  
% fcn_nonlinearCons_LaCMS(x0,pp,vv,mpc,gtd,initialCondition,Ybus,Yf,Yt,il,mpopt,para,dx,nx,dt,NK);
% objfcn_LaCMS(x0,mpc,vv,possibleResults, probability,nx,NK);
nonlinearCons = @(x)fcn_nonlinearCons_LaCMS(x,pp,vv,mpc,gtd,initialCondition,Ybus,Yf,Yt,il,mpopt,para,dx,nx,dt,NK);

objfcn = @(x)objfcn_LaCMS(x,mpc,vv,possibleResults, probability,nx,NK,baseMVA);

options = optimoptions('fmincon','Display','iter');
% options = optimoptions('fmincon','UseParallel',true);
options.MaxFunctionEvaluations = 3.000000e+08;
options.MaxIterations = 1.000000e+08;
options.ConstraintTolerance = 1e-5;
options.OptimalityTolerance = 1e-5;
options.StepTolerance = 1e-8;

[x1,fval,exitflag,output] = fmincon(objfcn,x0,[],[],[],[],LB,UB,nonlinearCons,options);

%% results
[Pg1,Qg1,Va1,Vm1,PrsRaw,GfRaw,PGs1,LCg1] = unpackVariables(x1,vv,NK); 
PrsPipe = mat2cell(PrsRaw,NK,nx+1); Gf1 = mat2cell(GfRaw,NK,nx+1);
PrsNodal = zeros(nGb,1);
for i = 1:nGl
    PrsNodal(mpc.Gline(i,1)) = PrsPipe{i}(1);
    PrsNodal(mpc.Gline(i,2)) = PrsPipe{i}(end);
end
[J_dynamic,J_endpoint,genLCeCost,gasCurtailmentCost,gasPurchasingCost] = LaCMScost(Pg1,LCg1,PGs1,PrsNodal,mpc,possibleResults,probability,NK);
optimalControlResults = struct('Pg',Pg1,'Qg',Qg1,'Va',Va1,'Vm',Vm1,'PrsPipe',PrsPipe,'PrsNodal',PrsNodal,...
    'Gf',Gf1,'PGs',PGs1,'LCg',LCg1,'f',fval','J_dynamic',J_dynamic,'J_endpoint',J_endpoint,'genLCeCost',genLCeCost,...
    'gasCurtailmentCost',gasCurtailmentCost,'gasPurchasingCost',gasPurchasingCost);


end



