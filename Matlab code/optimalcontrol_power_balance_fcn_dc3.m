function [g] = optimalcontrol_power_balance_fcn_dc3(Va0,Pg0,ei_inbus0,electricityLoad0,mpc0,NK)
for k = 1:NK
    mpc = mpc0;% 这里就涉及到取负荷值，不涉及故障约束，因此可以取basicLoad，但时间要对应
    Va = Va0(k,:)';
    Pg = Pg0(k,:)';
    electricityLoad = electricityLoad0(k,:)';
    ei_inbus = ei_inbus0(k,:)';
    [g(k,:)] = opf_power_balance_fcn(Va, Pg,ei_inbus,electricityLoad, mpc)';
%     [~,dg{k}] = opf_power_balance_fcn(Va, Vm, Pg, Qg, mpc, Ybus, mpopt);
end
end
%OPF_POWER_BALANCE_FCN  Evaluates AC power balance constraints and their gradients.
%   [G, DG] = OPF_POWER_BALANCE_FCN(X, OM, YBUS, MPOPT)
%
%   Computes the active or reactive power balance equality constraints for
%   AC optimal power flow. Computes constraint vectors and their gradients.
%
%   Inputs:
%     X : optimization vector
%     MPC : MATPOWER case struct
%     YBUS : bus admittance matrix
%     MPOPT : MATPOWER options struct
%
%   Outputs:
%     G  : vector of equality constraint values (active/reactive power balances)
%     DG : (optional) equality constraint gradients
%
%   Examples:
%       g = opf_power_balance_fcn(x, mpc, Ybus, mpopt);
%       [g, dg] = opf_power_balance_fcn(x, mpc, Ybus, mpopt);
%
%   See also OPF_POWER_BALANCE_HESS

%   MATPOWER
%   Copyright (c) 1996-2017, Power Systems Engineering Research Center (PSERC)
%   by Carlos E. Murillo-Sanchez, PSERC Cornell & Universidad Nacional de Colombia
%   and Ray Zimmerman, PSERC Cornell
%   and Baljinnyam Sereeter, Delft University of Technology
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.
%% 
function [g] = opf_power_balance_fcn(Va, Pg,ei_inbus, electricityLoad,mpc)
%%----- initialize -----
%% define named indices into data matrices
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
nb = length(Va);             %% number of buses
ng = length(Pg);            %% number of dispatchable injections

%% unpack data
[baseMVA, bus, gen, branch] = deal(mpc.baseMVA, mpc.bus, mpc.gen, mpc.branch);
%%

% gen(:, PG) = Pg * baseMVA;  %% active generation in MW
  [B, Bf, Pbusinj, Pfinj] = makeBdc(baseMVA, bus, branch);
  neg_Cg = sparse(gen(:, GEN_BUS), 1:ng, -1, nb, ng);   %% Pbus w.r.t. Pg
  Amis = [B neg_Cg];
  ei_inbus = ei_inbus*baseMVA;
  electricityLoad = electricityLoad*baseMVA;%MW
  %-----modified-------
  bmis = -(electricityLoad+ei_inbus + bus(:, GS)) / baseMVA - Pbusinj;
%     bmis = -(bus(:, PD)+ei_inbus + bus(:, GS)) / baseMVA - Pbusinj;
  g = Amis*[Va;Pg]-bmis;%MW/baseMVA
end
