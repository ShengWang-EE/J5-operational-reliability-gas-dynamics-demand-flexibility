function [h] = optimalcontrol_branch_flow_fcn_dc3(Va0, mpc0, il,NK)
for k = 1:NK

    mpc = mpc0;
    Va = Va0(k,:)';
    [h(k,:)] = opf_branch_flow_fcn(Va, mpc, il);
end
end
%OPF_BRANCH_FLOW_FCN  Evaluates AC branch flow constraints and Jacobian.
%   [H, DH] = OPF_BRANCH_FLOW_FCN(X, OM, YF, YT, IL, MPOPT)
%
%   Active power balance equality constraints for AC optimal power flow.
%   Computes constraint vectors and their gradients.
%
%   Inputs:
%     X : optimization vector
%     MPC : MATPOWER case struct
%     YF : admittance matrix for "from" end of constrained branches
%     YT : admittance matrix for "to" end of constrained branches
%     IL : vector of branch indices corresponding to branches with
%          flow limits (all others are assumed to be unconstrained).
%          YF and YT contain only the rows corresponding to IL.
%     MPOPT : MATPOWER options struct
%
%   Outputs:
%     H  : vector of inequality constraint values (flow limits)
%          where the flow can be apparent power, real power, or
%          current, depending on the value of opf.flow_lim in MPOPT
%          (only for constrained lines), normally expressed as
%          (limit^2 - flow^2), except when opf.flow_lim == 'P',
%          in which case it is simply (limit - flow).
%     DH : (optional) inequality constraint gradients, column j is
%          gradient of H(j)
%
%   Examples:
%       h = opf_branch_flow_fcn(x, mpc, Yf, Yt, il, mpopt);
%       [h, dh] = opf_branch_flow_fcn(x, mpc, Yf, Yt, il, mpopt);
%
%   See also OPF_BRANCH_FLOW_HESS.

%   MATPOWER
%   Copyright (c) 1996-2018, Power Systems Engineering Research Center (PSERC)
%   by Carlos E. Murillo-Sanchez, PSERC Cornell & Universidad Nacional de Colombia
%   and Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.
function [h] = opf_branch_flow_fcn(Va, mpc, il)
%%----- initialize -----
%% define named indices into data matrices
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;

[B, Bf, Pbusinj, Pfinj] = makeBdc(mpc.baseMVA, mpc.bus, mpc.branch);
    upf = mpc.branch(il, RATE_A) / mpc.baseMVA - Pfinj(il);
    upt = mpc.branch(il, RATE_A) / mpc.baseMVA + Pfinj(il);
%% unpack data
h1 = Bf(il,:)*Va - upf;
h2 = -upt - Bf(il,:)*Va;
h = [h1;h2];
end
