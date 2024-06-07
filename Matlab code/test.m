clear
clc
sdpvar x y z 
intvar zz
cons = [
    [x y z zz] <= 2;
    x+y <= z+zz;
%     norm([x,z]) <= zz
    ];
obj = x+2*y;
opt = sdpsettings('verbose',2,'solver','gurobi','usex0',0,'debug',1,'gurobi.InfUnbdInfo',1,...
    'savesolverinput',1,'savesolveroutput',1);
diag = optimize(cons,obj,opt);
% [cons_d,obj_d,varCone,varFree] = dualize(cons,obj);
% diag_d = optimize(cons_d,obj_d,opt);
% x1 = value(x); y1 = value(y);
% u = dual(cons);
%%
clear
clc
x = sdpvar(2,1);
y = sdpvar(1,1);
p = [1 + y;sum(x)+y];
p1 = replace(p,[y ],[1 ]);
sdisplay(p1)
%%
clear
sdpvar a11 a12 a21 b1 b2 b3 b4 b5
obj = 8*a11+2*a12+5*a21;
cons = [
    a11 + b1 <= 1;
    a12 + b2 <= 1;
    a21 + b3 <= 1;
    a11 + a21 + b4 <= 1;
    a11 + a12 + b5 <= 1;
    [b1 b2 b3 b4 b5] <= 0;
    ];
opts = sdpsettings('verbose',2,'solver','gurobi','usex0',0,'debug',1,'gurobi.InfUnbdInfo',1,...
    'savesolverinput',1,'savesolveroutput',1);
diag = optimize(cons,obj,opts);
%%
    clc
    clear
    yalmip('clear')
    load MP20200615.mat
    LB(iter) = MP_IEGSresult{iter}.f_MP;
    for i = 1:nEH
        % 注意MP得到的ei，gi都是标幺化后的
        ei0 = MP_IEGSresult{iter}.ei(:,i) * mpc.baseMVA; gi0 = MP_IEGSresult{iter}.gi(:,i) * 200; 
        ei0(iK.DR-iK.studyPeriod(1)+1) = ei0(iK.DR-iK.studyPeriod(1)+1) - 30;
        gi0(iK.DR-iK.studyPeriod(1)+1) = gi0(iK.DR-iK.studyPeriod(1)+1) - 10;
        [SP_EHresult{i},addCutsCoefficient{iter,i},SP_solverTime(iter,i),exitflagSP(iter,i)] = SP_EHschedule_milp(ei0,gi0,nodalEHpara{i},NK,iK);
    end