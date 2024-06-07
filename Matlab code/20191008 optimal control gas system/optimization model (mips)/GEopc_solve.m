function [results, success, raw] = GEopc_solve(om)
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
%% unpack data
mpc = om.get_mpc();
[baseMVA, bus, gen, branch, gencost] = ...
    deal(mpc.baseMVA, mpc.bus, mpc.gen, mpc.branch, mpc.gencost);
mpopt = [];
[vv, ll, nne, nni] = om.get_idx();
%% problem dimensions
nb = size(bus, 1);          %% number of buses
nl = size(branch, 1);       %% number of branches
ny = om.getN('var', 'y');   %% number of piece-wise linear costs
%add
nGb  = size(mpc.Gbus,1);
nGl  = size(mpc.Gline,1);
nGs  = size(mpc.Gsou,1);
nLCg = size(find(mpc.Gbus(:,3)~=0),1);
%add
%% linear constraints
[A, l, u] = om.params_lin_constraint();

%% bounds on optimization vars
[x0, xmin, xmax] = om.params_var();

%% build admittance matrices
[Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch);

%% try to select an interior initial point

%% find branches with flow limits
il = find(branch(:, RATE_A) ~= 0 & branch(:, RATE_A) < 1e10);
nl2 = length(il);           %% number of constrained lines
%%-----  run opf  -----
opt = [];
f_fcn = @(x)GEopc_costfcn(x, om);
gh_fcn = @(x)GEopf_consfcn(x, om, Ybus, Yf(il,:), Yt(il,:), mpopt, il);
hess_fcn = @(x, lambda, cost_mult)GEopf_hessfcn(x, lambda, cost_mult, om, Ybus, Yf(il,:), Yt(il,:), mpopt, il);
[x, f, info, Output, Lambda] = ...
  mips(f_fcn, x0, A, l, u, xmin, xmax, gh_fcn, hess_fcn, opt);
success = (info > 0);

end