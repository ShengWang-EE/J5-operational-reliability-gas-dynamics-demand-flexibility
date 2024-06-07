function [solution, diagnostics] = EHschedule(electricityPrice,gasPrice,EHpara,k)
% data
capacityCell = num2cell(EHpara.capacity);
efficiencyCell = num2cell(EHpara.efficiency);
[HA,EA,HB,EB,HC,EC,HD,ED,HO2MAX,HO2MIN,HO3MAX,HO3MIN,CO3MAX,CO3MIN,CO4MAX,CO4MIN] = deal(capacityCell{:});
[eta1e,eta1h,eta2,COP3hSummer,COP3hWinter,COP3c,COP4] = deal(efficiencyCell{:});

el = EHpara.load(1,k);
hl = EHpara.load(2,k);
cl = EHpara.load(3,k);
% electricityPrice = nodalPrice.electricity(time);
% gasPrice = nodalPrice.gas(time);
COP3h = COP3hSummer;
gamma = 0;%summer,EHP cooling
% vars
sdpvar ei gi eee ee3 e13 e1h gg1 gg2 h3h h1h h14 h24 h2h c3c c4c


% energy conversion constraints
energyConversionConstraints = [
    ei == eee+ee3;
    gi == gg1+gg2;
    eee+e1h == el;
    (ee3+e13)*COP3h*gamma == h3h;
    (ee3+e13)*COP3c*(1-gamma) == c3c;
    gg1*eta1e == e13+e1h;
    gg1*eta1h == h1h + h14;
    gg2*eta2 == h24+h2h;
    (h14+h24)*COP4 == c4c;
    h3h+h1h+h2h == hl;
    c3c+c4c == cl;
];
% capacity constraints
capacityConstraints = [
    [ei gi eee ee3 e13 e1h gg1 gg2 h3h h1h h14 h24 h2h c3c c4c] >= 0;
    h1h+h14 >= 0;
    e13+e1h-EA-(EA-EB)/(HA-HB)*(h1h+h14) <= 0;
    e13+e1h-EB-(EB-EC)*(HB-HC)*(h1h+h14-HB) >=0;
    e13+e1h-ED-(EC-ED)*(HC-HD)*(h1h+h14) >= 0;
    HO2MIN <= h24+h2h <= HO2MAX;
    CO4MIN <= c4c <= CO4MAX;
    gamma*HO3MIN <= h3h <= gamma*HO3MAX;
    (1-gamma)*CO3MIN <= c3c <= (1-gamma)*CO3MAX;
    ];

% solve
objective = electricityPrice * ei + gasPrice * gi;
constraints = [energyConversionConstraints;capacityConstraints]; 
opts = sdpsettings('solver','mosek');
diagnostics = optimize(constraints,objective,opts);
 
solution = [value(ei),value(gi),value(eee),value(ee3),value(e13),value(e1h),value(gg1),...
    value(gg2),value(h3h),value(h1h),value(h14),value(h24),value(h2h),value(c3c),value(c4c)];




end