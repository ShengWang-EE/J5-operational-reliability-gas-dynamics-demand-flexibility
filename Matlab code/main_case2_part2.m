clear
clc
%% 1 data input
% ----- 1.1 general paras -----
missionTime = 12; % 所仿真的运行时段时长
ND = 48; % 15min interval 联合调度时段时间间隔
dd = missionTime / ND;
NK = 48; % 15min time step 天然气暂态仿真时间步长
dt = missionTime / NK;
KK = 0:dt:missionTime;
simulationTimes = 10000; % TSMCS 次数

% ----- 1.2 IEGS para ------
[mpc0, gtd] = case24GEv5(); mpc = mpc0; % physical para
rts = Case24ReliabillityDatav3(); % reliability parameters
% preprocess: 要把发电机、天然气源的出力下限取消，不然负荷太小的时候不收敛（主要是电力侧）
mpc0.gen(:,10) = 0; 
%
nb   = size(mpc.bus, 1);    %% number of buses
nGb  = size(mpc.Gbus,1); % number of gas bus
nEH = size(mpc.EHlocation,1); % number of EHs
nGl = size(mpc.Gline,1);
nGen = sum(mpc.gen(:,22)==1)+sum(mpc.gen(:,22)==0);% all gen(TFU and GFU), excluded dispatchable loads
nGs = size(mpc.Gsou,1);
nComponent = nGen+nGs;
% --------------test-----------------------
% testmpc = mpc0; testmpc.bus(:,3) = mpc0.bus(:,3)*1.25; %经测试，最多负荷到1.27倍左右还能收敛
% testGEresult = GErunopf(testmpc);
%-----------------------
GEresult0 = GErunopf(mpc0);

% ----- 1.2 the load data for each energy hub -----
[loadCurve,EHpara] = EHdata();% original EH　data
% expand the time resolution of load profile
[load0.electricity, load0.heating, load0.cooling] = deal(interp(loadCurve(1,:),1/dt),interp(loadCurve(2,:),1/dt),interp(loadCurve(3,:),1/dt));
% 根据原始电力load最大值占原始EH各设备容量的比例，设置新的个节点的EH各设备容量。（按照mpc.EH排序）
[nodalEHpara] = scaleEH(mpc,[load0.electricity;load0.heating;load0.cooling],EHpara);% 96个点的EH负荷

% ----- 1.3 gas transient data -----
[dx, nx] = decidePipelineCell(gtd);% maximun value of 3 miles per cell according to the reference

% ----- 1.4 reliability data preprocess
for i = 1:nComponent
    if i<= nGen % is generator
        lamda(i) = 1 / rts.gen(i,1); mu(i) = 1 / rts.gen(i,2);
    else
        lamda(i) = 1 / rts.Gsou(i-nGen,1); mu(i) = 1/ rts.Gsou(i-nGen,2);
    end
    avaliability{i,1} = mu(i)/(lamda(i)+mu(i)) + lamda(i)/(lamda(i)+mu(i)) * exp(-(lamda(i)+mu(i))*KK);
    avaliability{i,2} = 1 - avaliability{i,1};
end

%%  这里仅用于获取节点能源价格，和负荷无关
ob.dayahead = zeros(NK,1);
[dayahead_mpc,dayahead_result] = deal(cell(NK,1));
[nodalPrice.electricity, nodalPrice.gas] = deal(zeros(NK,nb),zeros(NK,nGb));
for k = 1:NK
    % mpc for each period
    % 负荷波动不能太大，而且不现实。对一个用户而言可能波动很剧烈，但是大范围电网可能有很多固定负荷，比如工业负荷，波动不会那么大
    % 所以设定只有一半负荷随EH变化，另一半是固定的
    dayahead_mpc{k} = mpc0;
    dayahead_mpc{k}.bus(:,3:4) = mpc0.bus(:,3:4) * (0.6 + 0.4* load0.electricity(k)/max(load0.electricity)); 
    %把这个负荷变化趋势放缓一点
    dayahead_mpc{k}.Gbus(:,3) = mpc0.Gbus(:,3) *(0.6 + 0.4 * (load0.heating(k)+load0.cooling(k))/ max((load0.heating+load0.cooling)));
    % 天然气负荷就先这么这算，实时EH进行调节的时候在上面做差额计算，看这样行不行   
end
for k = 1:NK
    [dayahead_IEGSresult{k},flag] = GErunopf(dayahead_mpc{k});
    ob.dayahead(k) = flag;
    [nodalPrice.electricity(k,:),nodalPrice.gas(k,:)] = deal(dayahead_IEGSresult{k}.bus(:,14)',dayahead_IEGSresult{k}.Gbus(:,11)');
end
% 节点价格折算：原来的是$/MWh和$/(Mm3)，现在电不变，气换算成$/(MWh)
nodalPrice.gas = nodalPrice.gas/200;
%% calculate the day ahead optimal schedule of nodal EHs, according to the nodal energy price
for i = 1:nEH
    dayaheadEHschedule{i} = zeros(NK,15);%一个EH15个变量（包含ei，gi）
    for k = 1:NK
        Ebus = mpc.EHlocation(i,2); Gbus = mpc.EHlocation(i,1);
        [dayaheadEHschedule{i}(k,:),diagnostics] = EHschedule(nodalPrice.electricity(k,Ebus),nodalPrice.gas(k,Gbus),nodalEHpara{i},k);
        % attention: gi is in MW      
        ob.EH(i,k) = diagnostics.problem;
    end
end
% 修改dayaheadIEGS调度数据，去掉EH部分的负荷
dayahead_IEGSresult_basicLoad = dayahead_IEGSresult; % 这个也是在dayahead_mpc波动负荷的基础上算的
for k = 1:NK
    for i = 1:nEH
        EHgbus = mpc0.EHlocation(i,1); EHebus = mpc0.EHlocation(i,2); 
        dayahead_IEGSresult_basicLoad{k}.bus(EHebus,3:4) = dayahead_IEGSresult{k}.bus(EHebus,3:4) * ...
            (1-dayaheadEHschedule{i}(k,1)/dayahead_IEGSresult{k}.bus(EHebus,3));
        dayahead_IEGSresult_basicLoad{k}.Gbus(EHgbus,3) = dayahead_IEGSresult{k}.Gbus(EHgbus,3) * ...
            (1-dayaheadEHschedule{i}(k,2)/200/dayahead_IEGSresult{k}.Gbus(EHgbus,3));
        if (dayahead_IEGSresult_basicLoad{k}.bus(EHebus,3)<=0) || (dayahead_IEGSresult_basicLoad{k}.Gbus(EHgbus,3) <=0)
            ob.basicLoad = 1;%确保EH节点上的基本负荷不会变成负的
        end
    end
end
% 更新完EH的能源消耗后IEGS还调度吗？都可以
% EH的容量和所需要供给的负荷根据一定规则设计好了，那把剩余部分就认为是固定部分。
% 这样的前提下，节点的总负荷还是不变的，所以不用再计算一遍
save stop2.mat
%% stage2: MCS 
clear
clc
load stop2.mat

normalCondition = setInitialConditionForGasSystem(dayahead_IEGSresult{1},gtd,nx,dx);
%------------test-----------------
for k = 1:NK
    electricityLoad(k) = sum(dayahead_IEGSresult_basicLoad{k}.bus(:,3)) * 1 * 2850/2300;
    gasLoad(k) = sum(dayahead_IEGSresult_basicLoad{k}.Gbus(:,3)) * 46/45;
    dayahead_IEGSresult_basicLoad{k}.bus(:,3) = dayahead_IEGSresult_basicLoad{k}.bus(:,3) * 1 * 2850/2300;
    dayahead_IEGSresult_basicLoad{k}.Gbus(:,3) = dayahead_IEGSresult_basicLoad{k}.Gbus(:,3) * 46/45;
end
PrsVioPercentge = 0.1;terminalFactor = 0.00;startk = 1; endk = 24;
[LaCMS_optimizer1st] = lookAheadContingencyManagement_Optimizer(dt*3600,gtd,dx, NK/2,mpc,...
   nodalEHpara,dayahead_IEGSresult_basicLoad,mpc0,normalCondition,nEH,nx,startk,endk,PrsVioPercentge,terminalFactor);
PrsVioPercentge = 0.1;terminalFactor = 0.05;startk = 25; endk = 48;
[LaCMS_optimizer2nd] = lookAheadContingencyManagement_Optimizer(dt*3600,gtd,dx, NK/2,mpc,...
   nodalEHpara,dayahead_IEGSresult_basicLoad,mpc0,normalCondition,nEH,nx,startk,endk,PrsVioPercentge,terminalFactor);

% %%
clc
tic
simulationTimes = 3;
LCe = cell(simulationTimes,1); LCg = cell(simulationTimes,1);
for i = 1:simulationTimes
    % ----- 2.1 generate the system state sequence -----
    [Info_components] = MCSformingScenarioV6(rts,missionTime);
    Info_components(:,2:3) = round(Info_components(:,2:3)/dd); % 按照调度时段取整
    Info_components = deleteRepeatedLine(Info_components);
    % -------------- test ----------------------
    % 手动设置典型故障场景，先通过稳态OPF看看切负荷情况，好用就上暂态
    info1 = ones(1,nComponent);
    info2 = info1; info2(nGen+([1])) = 0;%气源1的1/5失效，触发紧急状态管理
    info3 = info2; info3(nGen+([2:4])) = 0;%气源1的3/5进一步失效，说明向前看的必要性，并变动参数，说明经济性和可靠性的权衡
    info4 = info1; info4([23,32]) = 0; % 创造电力系统故障，此时天然气作为边际机组给电力系统发电。变动气压范围，说明天然气系统对电力系统的支撑
    info5 = info4; info5(nGen+([1:4])) = 0; % 天然气系统故障，通过前一时段的不同的支撑说明，给太多支撑会给天然气系统可靠性带来问题
    Info_components =    [  1   0    8   info1;
                            2   8   16   info2;
                            3   16   24   info3;
                            4   24   32   info1;%恢复
                            5   32   40   info4;
                            6   40   48   info5;
                            ];
%     Info_components =    [  1   0    8   info1;
%                             2   8   16   info1;
%                             3   16   24   info1;
%                             4   24   32   info1;%恢复
%                             5   32   40   info4;
%                             6   40   48   info5;
%                             ];
% 比较不同气压限制对于电力系统故障负荷削减的影响
%     Info_components =    [  1   0    8   info1;
%                             2   8   16   info4;
%                             3   16   24   info1;
%                             4   24   32   info3;%恢复
%                             5   32   40   info3;
%                             6   40   48   info1;
%                             ];

% 拆成两部分（从第24个时段）
divideOrder = min(find(Info_components(:,2)>=NK/2));
Info_components1 = Info_components(1:divideOrder,:);
Info_components1(end,3) = NK/2;
Info_components2 = Info_components(divideOrder:end,:);
Info_components2(1,2) = NK/2;
Info_components2(:,2:3) = Info_components2(:,2:3) - NK/2;
Info_components1 = deleteRepeatedLine(Info_components1);
Info_components2 = deleteRepeatedLine(Info_components2);


    % ------------------------------------------
    % 1st 
    NS = size(Info_components1,1);
    newmpc = cell(NS,1);
    for s = 1:NS
        newmpc{s} = mpcUpdateBinary(Info_components1(s,4:end),mpc0,nGen);
    end

    [LCe1st{i},LCg1st{i},PGs1st{i},totalCost1st{i},GenAndLCeCost1st{i},GasPurchasingCost1st{i},GasCurtailmentCost1st{i},...
        nextInitialCondition] = lookAheadContingencyManagement_Solver(LaCMS_optimizer1st,NK/2,NS,mpc,...
       mpc0,newmpc,normalCondition,Info_components1,nx); 
   % 2nd
   NS = size(Info_components2,1);
    for s = 1:NS
        newmpc{s} = mpcUpdateBinary(Info_components2(s,4:end),mpc0,nGen);
    end
    [LCe2nd{i},LCg2nd{i},PGs2nd{i},totalCost2nd{i},GenAndLCeCost2nd{i},GasPurchasingCost2nd{i},GasCurtailmentCost2nd{i},...
        nextInitialCondition] = lookAheadContingencyManagement_Solver(LaCMS_optimizer2nd,NK/2,NS,mpc,...
       mpc0,newmpc,nextInitialCondition,Info_components2,nx); 
   
end
save stop3.mat
%%
% 注意我把天然气的故障率提升了10倍。按理来说拆成小气源故障率提升2倍就差不多了。
% 所以朴素来说，天然气相关的可靠性、切气成本等都除以5
nLCg = size(find(mpc.Gbus(:,3)~=0),1); nLCe = size(find(mpc.bus(:,3)~=0),1);
% 1st
[systemEDNSsum1st,systemEGNSsum1st] = deal(zeros(NK/2,1));
[sumPGs1st] = zeros(1,6);
[sumTotalCost1st,sumGenAndLCeCost1st,sumGasPurchasingCost1st,sumGasCurtailmentCost1st] = deal(zeros(1,1));
for i = 1:simulationTimes
    systemEDNSsum1st = systemEDNSsum1st + sum(LCe1st{i},2);systemEGNSsum1st = systemEGNSsum1st + sum(LCg1st{i},2);    
    sumPGs1st = sumPGs1st + sum(PGs1st{i});sumTotalCost1st = sumTotalCost1st + totalCost1st{i};
    sumGenAndLCeCost1st = sumGenAndLCeCost1st + GenAndLCeCost1st{i};
    sumGasPurchasingCost1st = sumGasPurchasingCost1st + GasPurchasingCost1st{i};
    sumGasCurtailmentCost1st = sumGasCurtailmentCost1st + GasCurtailmentCost1st{i};    
end
systemEDNS1st = systemEDNSsum1st/simulationTimes;systemEGNS1st = systemEGNSsum1st/simulationTimes;
expectedPGs1st = sumPGs1st/4/24/simulationTimes;expectedTotalCost1st = sumTotalCost1st/simulationTimes;
expectedGenAndLCeCost1st = sumGenAndLCeCost1st / simulationTimes; 
expectedGasPurchasingCost1st = sumGasPurchasingCost1st / simulationTimes;
expectedGasCurtailmentCost1st = sumGasCurtailmentCost1st / simulationTimes;
EENS1st = sum(systemEDNS1st)/4; EVNS1st = sum(systemEGNS1st)/24/4;
% 2nd
[systemEDNSsum2nd,systemEGNSsum2nd] = deal(zeros(NK/2,1));
[sumPGs2nd] = zeros(1,6);
[sumTotalCost2nd,sumGenAndLCeCost2nd,sumGasPurchasingCost2nd,sumGasCurtailmentCost2nd] = deal(zeros(1,1));
for i = 1:simulationTimes
    systemEDNSsum2nd = systemEDNSsum2nd + sum(LCe2nd{i},2);systemEGNSsum2nd = systemEGNSsum2nd + sum(LCg2nd{i},2);    
    sumPGs2nd = sumPGs2nd + sum(PGs2nd{i});sumTotalCost2nd = sumTotalCost2nd + totalCost2nd{i};
    sumGenAndLCeCost2nd = sumGenAndLCeCost2nd + GenAndLCeCost2nd{i};
    sumGasPurchasingCost2nd = sumGasPurchasingCost2nd + GasPurchasingCost2nd{i};
    sumGasCurtailmentCost2nd = sumGasCurtailmentCost2nd + GasCurtailmentCost2nd{i};    
end
systemEDNS2nd = systemEDNSsum2nd/simulationTimes;systemEGNS2nd = systemEGNSsum2nd/simulationTimes;
expectedPGs2nd = sumPGs2nd/4/24/simulationTimes;expectedTotalCost2nd = sumTotalCost2nd/simulationTimes;
expectedGenAndLCeCost2nd = sumGenAndLCeCost2nd / simulationTimes; 
expectedGasPurchasingCost2nd = sumGasPurchasingCost2nd / simulationTimes;
expectedGasCurtailmentCost2nd = sumGasCurtailmentCost2nd / simulationTimes;
EENS2nd = sum(systemEDNS2nd)/4; EVNS2nd = sum(systemEGNS2nd)/24/4;

% std
LCeAll = zeros(NK/2,simulationTimes);LCgAll = zeros(NK/2,simulationTimes);
for i = 1:simulationTimes
    LCeAll(:,i) = sum(LCe1st{i},2);LCgAll(:,i) = sum(LCg1st{i},2);
    if mod(i,1000) == 0
        stdSystemEDNS(:,i) = std(LCeAll(:,1:i),0,2)./mean(LCeAll(:,1:i),2)/sqrt(i); 
        stdSystemEGNS(:,i) = std(LCgAll(:,1:i),0,2)./mean(LCgAll(:,1:i),2)/sqrt(i); 
    end
end
toc    
%% convergence
for i = 1:simulationTimes
    if mod(i,1000) == 0
        EDNScvrg(:,i/1000) = stdSystemEDNS(40:48,i);
        EGNScvrg(:,i/1000) = stdSystemEGNS(40:48,i);
    end
end
EDNScvrg = EDNScvrg';EGNScvrg = EGNScvrg';
% nodal EENS and EGNS
nodalEENS = sum(nodalEDNS)'/4;%MWh
nodalEVNS = sum(nodalEGNS)'/24/4;%Mm3

%%
[totalCost1st,GenAndLCeCost1st,GasPurchasingCost1st,GasCurtailmentCost1st] = ...
    deal(totalCost1st,sumGenAndLCeCost1st,sumGasPurchasingCost1st,sumGasCurtailmentCost1st);
[totalCost2nd,GenAndLCeCost2nd,GasPurchasingCost2nd,GasCurtailmentCost2nd] = ...
    deal(totalCost2nd,sumGenAndLCeCost2nd,sumGasPurchasingCost2nd,sumGasCurtailmentCost2nd);