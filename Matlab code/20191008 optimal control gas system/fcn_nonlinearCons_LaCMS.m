function [c,ceq,dc,dceq] = fcn_nonlinearCons_LaCMS(x,pp,vv,mpc,gtd,initialCondition,Ybus,Yf,Yt,il,mpopt,para,dx,nx,dt,NK)
[Pg,Qg,Va,Vm,Prs,Gf,PGs,LCg] = unpackVariables(x,vv,NK);
f_PQmis = optimalcontrol_power_balance_fcn_mips(Va, Vm, Pg, Qg, mpc, Ybus, mpopt,NK);
f_flow = optimalcontrol_branch_flow_fcn_mips(Va, Vm, mpc, Yf(il, :), Yt(il, :), il, mpopt,NK);

f_Gmis = optimalcontrol_gas_balance_fcn_mips(Pg,Gf,PGs,LCg,pp,vv,mpc,nx,NK);
f_PrsEqual = optimalcontrol_gas_pressure_fcn_mips(Prs,pp,mpc,NK);
f_Gasflow = optimalcontrol_dynamicGasflow_fcn_mips(Prs,Gf,mpc,gtd,initialCondition,para,dx,nx,dt,NK) / 1e5;

c = reshape(f_flow,1,[])';
ceq = [reshape(f_PQmis,1,[])';reshape(f_Gmis,1,[])';reshape(f_PrsEqual,1,[])';reshape(f_Gasflow,1,[])'];



end