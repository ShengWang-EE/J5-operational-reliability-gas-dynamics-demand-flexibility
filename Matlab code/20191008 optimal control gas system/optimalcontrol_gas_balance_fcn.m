function [g] = optimalcontrol_gas_balance_fcn(Pg0,q0,g0,gc0,pp,vv,mpc,nx,NK)
for k = 1:NK
    Pg = Pg0(k,:)';
    Gf = q0(k,:)';
    PGs = g0(k,:)';
    LCg = gc0(k,:)';
    [g(k,:)] = opf_gas_balance_fcn(Pg,Gf,PGs,LCg,pp,mpc,vv,nx)';
%     [~,dg{k}] = opf_gas_balance_fcn(Pg,Gf,PGs,LCg,pp,mpc,vv,nx);
end
end
function [g, dg] = opf_gas_balance_fcn(Pg,q,PGs,LCg,pp,mpc,vv,nx)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
nGl = size(mpc.Gline,1);
%% unpack variables
gasflow = zeros(nGl,2);
for i = 1:nGl
    gasflow(:,1) = q(pp(:,1));
    gasflow(:,2) = q(pp(:,2));
end


Pgen = 100 * Pg; %100MVA
%exclude the LCe (because of LCe, the order of original gen is mixed
LCeIndex = find(mpc.gen(:,22)==2);
Pgen(LCeIndex,:)=[];
mpc.gen(LCeIndex,:) = [];
[g,dg] = calculate_g(Pgen,gasflow,PGs,LCg,mpc,vv,nx);


end
%%
function [g,dg] = calculate_g(Pg,Gf,PGs,LCg,mpc,vv,nx)
    % numbers
    Gnb=size(mpc.Gbus,1); Gnl = size(mpc.Gline,1);
    Nx = vv.N.all; dg = zeros(Nx,Gnb);
%     iPGs = vv.i1.PGs:vv.iN.PGs; iPg = vv.i1.Pg:vv.iN.Pg;
%     iGf = mat2cell(vv.i1.Gf:vv.iN.Gf,1,nx+1);
    %% 首先计算不涉及到线路注入潮流的几个分量
    %计算节点本身净产气量，注意要用变量x里面的值
%     Pgd=zeros(Gnb,1); % nodal gas injection
    gtpEfficiency = 200;
    for i=1:Gnb
        Pgd(i)=-mpc.Gbus(i,3);
        if ismember(i,mpc.Gsou(:,1))
           Pgd(i)=Pgd(i)+PGs(find(i==mpc.Gsou(:,1)));%-负荷+气源
%            dg(iPGs(find(i==mpc.Gsou(:,1)))) = 1;
        end
    end
    % we need to consider the gas consumption of gfu
    for i = 1:size(mpc.gfuIndex,1)
        electricityBus = mpc.gen(mpc.gfuIndex(i),1);
        gasBus = mpc.GEcon(find(mpc.GEcon(:,2)==electricityBus),1);
        Pgd(gasBus)=Pgd(gasBus)-Pg(mpc.gfuIndex(i)) / gtpEfficiency;%找到相应的Ggtp编号,再-GTP
%         dg(iPg(mpc.gfuIndex(i)),gasBus) = -1/gtpEfficiency;
    end
    %% 然后计算线路流量（流出为正）
%     dGP = zeros(Gnb,1);
    for i = 1:Gnl
        fb = mpc.Gline(i,1); tb = mpc.Gline(i,2);
        dGP(fb) = dGP(fb) + Gf(i,1);
        dGP(tb) = dGP(tb) - Gf(i,2);%流出为正;%流出为正
%         dg(iGf{i}(1),fb) = 1;
%         dg(iGf{i}(end),tb) = -1;
    end
        
    dPgd=Pgd-dGP;  
    addg=dPgd;
    %% 然后修改对应的LCg
    LCg_bus = find(mpc.Gbus(:,3)~=0);
    addg(LCg_bus) = addg(LCg_bus) + LCg;
%     dg(iLCg(LCg_bus)) = 1;
    
    g = addg;
end


