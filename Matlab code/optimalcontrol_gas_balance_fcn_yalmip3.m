function [g] = optimalcontrol_gas_balance_fcn_yalmip3(Pg0,q0,g0,gc0,gi_inbus0,gasLoad0,pp,vv,mpc0,nx,NK)
for k = 1:NK
    mpc = mpc0;
    Pg = Pg0(k,:)';
    Gf = q0(k,:)';
    PGs = g0(k,:)';
    LCg = gc0(k,:)';
    gi_inbus = gi_inbus0(k,:)';
    gasLoad = gasLoad0(k,:)';
    [g(k,:)] = opf_gas_balance_fcn(Pg,Gf,PGs,LCg,gi_inbus,gasLoad,pp,mpc,vv,nx)';
%     [~,dg{k}] = opf_gas_balance_fcn(Pg,Gf,PGs,LCg,pp,mpc,vv,nx);
end
end
function [g, dg] = opf_gas_balance_fcn(Pg,q,PGs,LCg,gi_inbus,gasLoad,pp,mpc,vv,nx)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
nGl = size(mpc.Gline,1);
%% unpack variables

% for i = 1:nGl
    gasflow(:,1) = q(pp(:,1));
    gasflow(:,2) = q(pp(:,2));
% end
 

Pgen = 100 * Pg; %100MVA
%exclude the LCe (because of LCe, the order of original gen is mixed
LCeIndex = find(mpc.gen(:,22)==2);
PgIndex = find(mpc.gen(:,22)~=2);
Pgen = Pgen(PgIndex,:);
mpc.gen = mpc.gen(PgIndex,:);
[g,dg] = calculate_g(Pgen,gasflow,PGs,LCg,gi_inbus,gasLoad,mpc,vv,nx);


end
%%
function [g,dg] = calculate_g(Pg,Gf,PGs,LCg,gi_inbus,gasLoad,mpc,vv,nx)
    % numbers
    Gnb=size(mpc.Gbus,1); Gnl = size(mpc.Gline,1);
    nGs = size(mpc.Gsou,1);
    Nx = vv.N.all; dg = zeros(Nx,Gnb);

    %% 首先计算不涉及到线路注入潮流的几个分量
    %计算节点本身净产气量，注意要用变量x里面的值
    PGsbus = mpc.Gsou(:,1); 
    Cgs = sparse(PGsbus, (1:nGs)', 1, Gnb, nGs); % connection matrix
    gtpEfficiency = 200;
    Pgd = Cgs*PGs - gasLoad - gi_inbus;
%     for i=1:Gnb
%         Pgd(i)=-mpc.Gbus(i,3);
%         if ismember(i,mpc.Gsou(:,1))
%            Pgd(i)=Pgd(i)+PGs(find(i==mpc.Gsou(:,1)));%-负荷+气源
% %            dg(iPGs(find(i==mpc.Gsou(:,1)))) = 1;
%         end
%     end
    % we need to consider the gas consumption of gfu
    for i = 1:size(mpc.gfuIndex,1)
        electricityBus = mpc.gen(mpc.gfuIndex(i),1);
        gasBus = mpc.GEcon(find(mpc.GEcon(:,2)==electricityBus),1);
        Pgd(gasBus)=Pgd(gasBus)-Pg(mpc.gfuIndex(i)) / gtpEfficiency;%找到相应的Ggtp编号,再-GTP
%         dg(iPg(mpc.gfuIndex(i)),gasBus) = -1/gtpEfficiency;
    end
    %% 然后计算线路流量
%     dGP = zeros(Gnb,1);
    for i = 1:Gnl
        fb = mpc.Gline(i,1); tb = mpc.Gline(i,2);
        Pgd(fb) = Pgd(fb) - Gf(i,1);
        Pgd(tb) = Pgd(tb) + Gf(i,2);%流入为正
%         dg(iGf{i}(1),fb) = 1;
%         dg(iGf{i}(end),tb) = -1;
    end
        
%     dPgd=Pgd-dGP;  
  %% 然后修改对应的LCg
    LCg_bus = find(mpc.Gbus(:,3)~=0);
    Pgd(LCg_bus) = Pgd(LCg_bus) + LCg;
%     dg(iLCg(LCg_bus)) = 1;    
    
    g = Pgd;

end


