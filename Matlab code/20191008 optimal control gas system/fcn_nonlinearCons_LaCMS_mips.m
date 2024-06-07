function [c,ceq,dc,dceq] = fcn_nonlinearCons_LaCMS_mips(x,pp,vv,mpc,gtd,initialCondition,Ybus,Yf,Yt,il,mpopt,para,dx,nx,dt,NK)
Nx = size(x,1) / NK;
Pg = zeros(NK,vv.N.Pg);Qg = zeros(NK,vv.N.Qg);Va = zeros(NK,vv.N.Va);Vm = zeros(NK,vv.N.Vm);
Prs = zeros(NK,vv.N.Prs);Gf = zeros(NK,vv.N.Gf);PGs = zeros(NK,vv.N.PGs);LCg = zeros(NK,vv.N.LCg);
for k = 1:NK
    xk = x((k-1)*Nx+1:k*Nx);
    
    Pg(k,:) = xk(vv.i1.Pg:vv.iN.Pg)';Qg(k,:) = xk(vv.i1.Qg:vv.iN.Qg)';
    Va(k,:) = xk(vv.i1.Va:vv.iN.Va)';Vm(k,:) = xk(vv.i1.Vm:vv.iN.Vm)';
    Prs(k,:) = xk(vv.i1.Prs:vv.iN.Prs)';Gf(k,:) = xk(vv.i1.Gf:vv.iN.Gf)';
    PGs(k,:) = xk(vv.i1.PGs:vv.iN.PGs)';LCg(k,:) = xk(vv.i1.LCg:vv.iN.LCg)';
end
f_PQmis = optimalcontrol_power_balance_fcn(Va, Vm, Pg, Qg, mpc, Ybus, mpopt,NK);
f_flow = optimalcontrol_branch_flow_fcn(Va, Vm, mpc, Yf(il, :), Yt(il, :), il, mpopt,NK);

f_Gmis = optimalcontrol_gas_balance_fcn(Pg,Prs,PGs,LCg,pp,mpc,NK);
f_PrsEqual = optimalcontrol_gas_pressure_fcn(Prs,pp,mpc,NK);
f_Gasflow = optimalcontrol_dynamicGasflow_fcn(Prs,Gf,mpc,gtd,initialCondition,para,dx,nx,dt,NK) / 1e5;

c0=f_flow;
ceq0 = [f_PQmis,f_Gmis,f_PrsEqual,f_Gasflow];

c = reshape(c0,[],1);
ceq = reshape(ceq0,[],1);

%%
dc = zeros(size(c));
dceq = zeros(size(ceq));
end