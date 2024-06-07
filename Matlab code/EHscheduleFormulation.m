function [EHconstraints,EHobjective] = EHscheduleFormulation(ei,gi,eee,ee3,e13,e1h,gg1,gg2,h3h,h1h,h14,h24,h2h,c3c,c4c,lc,...
    nodalEHpara,NK,startk,endk)
% yalmip('clear');
% data

capacityCell = num2cell(nodalEHpara.capacity);
efficiencyCell = num2cell(nodalEHpara.efficiency);
scheduleCell = num2cell(nodalEHpara.schedule);
[HA,EA,HB,EB,HC,EC,HD,ED,HO2MAX,HO2MIN,HO3MAX,HO3MIN,CO3MAX,CO3MIN,CO4MAX,CO4MIN] = deal(capacityCell{:});
[eta1e,eta1h,eta2,COP3hSummer,COP3hWinter,COP3c,COP4] = deal(efficiencyCell{:});
[alpha,beta,CDF] = deal(scheduleCell{:});
EHload = nodalEHpara.load';

COP3h = COP3hSummer;
gamma = 0;%summer,EHP cooling
NL = 3;
[el,hl,cl] = deal(EHload(startk:endk,1),EHload(startk:endk,2),EHload(startk:endk,3));

%------------------test---------------
beta = 1;
%-------------------------------
%%   
% schedule capacity
scheduleEquations = [
    -lc;
    lc - beta * EHload(startk:endk,:); % limits for curtailment
    ];
scheduleConstraints = [scheduleEquations <= 0]:'scheduleConstraints';
% energy conversion constraints
energyConversionEquations = [
     eee+ee3 - ei;
     gg1+gg2 - gi;
    -(eee+e1h) + el;
%     h3h - ((ee3+e13)*COP3h*gamma);
    c3c - ((ee3+e13)*COP3h*(1-gamma));
    e13+e1h - gg1*eta1e;
    h1h + h14 - gg1*eta1h;
    h24+h2h - gg2*eta2;
    c4c - (h14+h24)*COP4;
    hl - (h3h+h1h+h2h);
    cl - (c3c+c4c);
    -eee;
    -ee3;
    -e13;
    -e1h;
    -gg1;
    -gg2;
    -h1h;
    -h14;
    -h24;
    -h2h;
    -c3c;
    -c4c;
    ];
energyConversionConstraints = [energyConversionEquations<=0;
    h3h == 0;]:'energyConversionConstraints';
capacityEquations = [
    -(h1h+h14);
    e13+e1h-EA-(EA-EB)/(HA-HB)*(h1h+h14);
    -( e13+e1h-EB-(EB-EC)/(HB-HC)*(h1h+h14-HB) );
    -( e13+e1h-ED-(EC-ED)/(HC-HD)*(h1h+h14) );
    -(h24+h2h) + HO2MIN;
    h24+h2h - HO2MAX;
    -c4c + CO4MIN;
    c4c - CO4MAX;
%     -h3h + gamma*HO3MIN;  
%     h3h - gamma*HO3MAX;
    -c3c + (1-gamma)*CO3MIN; 
    c3c - (1-gamma)*CO3MAX;
    ];
capacityConstraints = [capacityEquations <= 0]:'capacityConstraints';
%%
EHconstraints = [scheduleConstraints;energyConversionConstraints;capacityConstraints];
EHobjective = objfcn_SP(lc,CDF,NK);


end
