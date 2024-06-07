function [LCe,LCg1,PGs1,totalCost,sumGenAndLCeCost,sumGasPurchasingCost,sumGasCurtailmentCost,nextInitialCondition] = ...
    lookAheadContingencyManagement_Solver(LaCMS_optimizer,NK,NS,mpc,...
   mpc0,newmpc,initialCondition,Info_components,nx)
%% initialization
NL = 3;
mpc0.basePrs = 1e6;
[para] = initializeParameters2();
% data dimensions
nb   = size(mpc.bus, 1);    %% number of buses
nl   = size(mpc.branch, 1); %% number of branches
ngen   = size(mpc.gen, 1);    %% number of dispatchable injections

%add GFU,LCe are included in the mpc
nGb  = size(mpc.Gbus,1); % number of gas bus
nGl  = size(mpc.Gline,1); % number of gas line
nGs  = size(mpc.Gsou,1); % number of gas source
nLCg = size(find(mpc.Gbus(:,3)~=0),1);
% 
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

%% specify parameter variables
PGsmax = zeros(NK,nGs); %Mm3/day
Pgmax = zeros(NK,ngen); % 每个场景相同，2个调度时段不同，其中第NS个场景为元件状态不变化的

for s = 1:NS
    duration = Info_components(s,3)-Info_components(s,2);
    PGsmax(Info_components(s,2)+1:Info_components(s,3),:) = repmat(newmpc{s}.Gsou(:,4)',duration,1);
    Pgmax(Info_components(s,2)+1:Info_components(s,3),:) = repmat(newmpc{s}.gen(:,PMAX)',duration,1)/baseMVA;
end
initialP = cell2mat(initialCondition.P')/1e6;initialQ = cell2mat(initialCondition.Q');
paraVar = {PGsmax,Pgmax,initialP,initialQ};
[stateVar,errorCode] = LaCMS_optimizer(paraVar);
[Prs,Gf,PGs,LCg,Va,Pg,ei,gi,...
    eeeAll,ee3All,e13All,e1hAll,gg1All,gg2All,h3hAll,...
    h1hAll,h14All,h24All,h2hAll,c3cAll,c4cAll,lcAll]=deal(stateVar{:});
%% results processing
[PrsRaw1,GfRaw1,PGsRaw1,LCg1,Va1,Pg1,ei1,gi1] = ...
    deal(value(Prs),value(Gf),value(PGs),value(LCg),value(Va),value(Pg)*baseMVA,value(ei)*baseMVA,value(gi));
PrsPipe = mat2cell(PrsRaw1*10,NK,nx+1); GfPipe = mat2cell(GfRaw1,NK,nx+1); 
PrsNodal = zeros(NK,nGb);
for i = 1:nGl
    PrsNodal(:,mpc.Gline(i,1)) = PrsPipe{i}(:,1);%bar
    PrsNodal(:,mpc.Gline(i,2)) = PrsPipe{i}(:,end);
end
PGs1 = [sum(PGsRaw1(:,1:5),2),sum(PGsRaw1(:,6:9),2),sum(PGsRaw1(:,10:11),2),PGsRaw1(:,12:14)];
% EH
[eeeAll1,ee3All1,e13All1,e1hAll1,gg1All1,gg2All1,h3hAll1,h1hAll1,h14All1,h24All1,h2hAll1,c3cAll1,c4cAll1] ...
    = deal(value(eeeAll),value(ee3All),value(e13All),value(e1hAll),value(gg1All),value(gg2All),...
    value(h3hAll),value(h1hAll),value(h14All),value(h24All),value(h2hAll),value(c3cAll),value(c4cAll));
ho1 = h1hAll1+h14All1; eo1 = e13All1+e1hAll1;
ho2 = h24All1 + h2hAll1;
co3 = c3cAll1;
co4 = c4cAll1;
%% extract LC PGs cost
LCe = Pg1(:,34:end); 
[totalCost,sumGenAndLCeCost,sumGasPurchasingCost,sumGasCurtailmentCost] = objfcn_IEGSdispatch(Pg1/baseMVA,PGsRaw1,LCg1,mpc,NK);
%% next initial condition
nextInitialCondition.P = mat2cell(PrsRaw1(end,:)*1e6,1,nx+1)';
nextInitialCondition.Q = mat2cell(GfRaw1(end,:),1,nx+1)';
end

