function [g,dg] = optimalcontrol_gas_balance_fcn_mips(x,pp,mpc,nx,NK)
nGb = size(mpc.Gbus,1);
[Pg,Gf, PGs, LCg] = deal(x{:});
[nPg,nGf,nPGs,nLCg] = deal(size(Pg,1),size(Gf,1),size(PGs,1),size(LCg,1));
cumnvar = cumsum([nPg,nGf,nPGs,nLCg]);
nvar = sum([nPg,nGf,nPGs,nLCg]);
[Pg0,Gf0, PGs0, LCg0] = deal(reshape(Pg,[],NK),reshape(Gf,[],NK),reshape(PGs,[],NK),reshape(LCg,[],NK));%列是k
g = zeros(NK*nGb,1); dg = zeros(NK*nGb,nvar);[g0,dg0] = deal(cell(NK,1));
for k = 1:NK
    Pg = Pg0(:,k);Gf = Gf0(:,k);
    PGs = PGs0(:,k);LCg = LCg0(:,k);
    [g0{k},dg0{k}] = opf_gas_balance_fcn(Pg,Gf,PGs,LCg,pp,mpc,nx);
%     [g((k-1)*nGb+1:k*nGb)] = opf_gas_balance_fcn(Pg,Gf,PGs,LCg,pp,mpc,vv,nx);
%     [dg((k-1)*2*nGb+1:k*4*nGb,(k-1)*nGb+1:k*nGb)] = opf_gas_balance_fcn(Pg,Gf,PGs,LCg,pp,mpc,vv,nx);
    g((k-1)*nGb+1:k*nGb) = g0{k};    
end
for i = 1:nvar
    itype = min(find(cumnvar>=i));
    k = mod(i,NK);
    if k==0
        k=NK;
    end
    dg((k-1)*nGb+1:k*nGb,i) = dg0{k}(:,itype);
end

end
function [g, dg] = opf_gas_balance_fcn(Pg,Gf,PGs,LCg,pp,mpc,nx)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
[baseMVA, bus, gen,Gbus,Gline,Gsou,gfuIndex] = ...
    deal(mpc.baseMVA, mpc.bus, mpc.gen,mpc.Gbus,mpc.Gline,mpc.Gsou,mpc.gfuIndex);
[nPg,nGf,nPGs,nLCg] = deal(size(Pg,1),size(Gf,1),size(PGs,1),size(LCg,1));
Nx = nPg + nGf + nPGs + nLCg;
[iPg,iGf,iPGs,iLCg] = deal(1:nPg,nPg+1:nPg+nGf,nPg+nGf+1:nPg+nGf+nPGs,nPg+nGf+nPGs+1:nPg+nGf+nPGs+nLCg);
nGl = size(mpc.Gline,1); nGb = size(mpc.Gbus,1);
%% unpack variables
Gf1 = zeros(nGl,2);
Gf1(:,1) = Gf(pp(:,1));
Gf1(:,2) = Gf(pp(:,2));
iGf1(:,1) = iGf(pp(:,1));
iGf1(:,2) = iGf(pp(:,2));


Pg = 100 * Pg; %100MVA
%exclude the LCe (because of LCe, the order of original gen is mixed
LCeIndex = find(gen(:,22)==2);
genIndex = find(gen(:,22)~=2);
% Pg(LCeIndex,:)=[];
% mpc.gen(LCeIndex,:) = [];
Pgen = Pg(genIndex,:);
gen_woLCe = gen(genIndex,:); % exclude LCe gen
iPgen = iPg(genIndex);

%% 首先计算不涉及到线路注入潮流的几个分量
%计算节点本身净产气量，注意要用变量x里面的值

g = zeros(nGb,1);dg = zeros(Nx,nGb);
gtpEfficiency = 200;
for i=1:nGb
    g(i)=Gbus(i,3);
    if ismember(i,mpc.Gsou(:,1))
       g(i)=g(i)+PGs(find(i==Gsou(:,1)));%-负荷+气源
       dg(iPGs(find(i==Gsou(:,1))),i) = 1;
    end
end
% we need to consider the gas consumption of gfu
for i = 1:size(mpc.gfuIndex,1)
    electricityBus = gen_woLCe(gfuIndex(i),1);
    gasBus = mpc.GEcon(find(mpc.GEcon(:,2)==electricityBus),1);
    g(gasBus)=g(gasBus)-Pgen(mpc.gfuIndex(i)) / gtpEfficiency;%找到相应的Ggtp编号,再-GTP
    dg(iPgen(mpc.gfuIndex(i)),gasBus) = -1/gtpEfficiency;
end
%% 然后计算线路流量（流出为正）

for i = 1:nGl
    fb = mpc.Gline(i,1); tb = mpc.Gline(i,2);
    g(fb) = g(fb) - Gf1(i,1);
    g(tb) = g(tb) + Gf1(i,2);%流入为正
    dg(iGf1(i,1),fb) = -1;
    dg(iGf1(i,2),tb) = 1;
end   
%% 然后修改对应的LCg
LCg_bus = find(mpc.Gbus(:,3)~=0);
g(LCg_bus) = g(LCg_bus) + LCg;
dg(iLCg,LCg_bus) = 1;
dg = dg';   

end


