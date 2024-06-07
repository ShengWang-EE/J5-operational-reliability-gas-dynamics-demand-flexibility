% function optimalControlResult = lookaheadContingencyManagement(initialCondition,initialSolutionPoint,mpc,gtd,possibleResults, probability,dx,nx,dt,NK)
% save f
clc
clear
load f
load mpopt
% [mpc, gtd] = case24GEv3();
steadyResults = GErunopf(mpc);
%%
om = GEopc_setup(mpc,gtd,mpopt,steadyResults,NK,initialSolutionPoint,initialCondition,dx,nx,dt,possibleResults,probability);


[results, success, raw] = GEopc_solve(om);

    
options = sdpsettings('solver','IPOPT','verbose',2,'usex0',1,'debug',1);
% options.ipopt.max_iter = 1e8;
% options.ipopt.max_cpu_time = 1e4;
options.ipopt.tol = 1e-3;
% options.fmincon.TolCon = 1e-4;
% options.fmincon.TolFun = 1e-4;
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
a=1;
% end
