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
% [constraints,objective,variables] = 
% LaCMS_controller = lookAheadContingencyManagement_Optimizer(nComponent,dd,dt,nodalEHpara, gtd,...
%     nx,dx,nEH, mpc0);
% [model,recoverymodel,diagnostic,internalmodel] = export(constraints,objective);
% %%
clc
normalCondition = setInitialConditionForGasSystem(dayahead_IEGSresult{1},gtd,nx,dx);
%------------test-----------------
for k = 1:NK
%     dayahead_IEGSresult_basicLoad{k} = GEresult0;
end
[LaCMS_optimizer] = lookAheadContingencyManagement_Optimizer(dt*3600,gtd,dx, NK,mpc,...
   nodalEHpara,dayahead_IEGSresult_basicLoad,mpc0,normalCondition,nEH,nx);
%%
tic
simulationTimes = 1000;
LCe = cell(simulationTimes,1); LCg = cell(simulationTimes,1);
for i = 1:simulationTimes
    % ----- 2.1 generate the system state sequence -----
    [Info_components] = MCSformingScenarioV6(rts,missionTime);
    Info_components(:,2:3) = round(Info_components(:,2:3)/dd); % 按照调度时段取整
    nSystemStates = size(Info_components,1);
    deleteLine = [];
    for s = 1:nSystemStates
        if Info_components(s,2) == Info_components(s,3)
            deleteLine = [deleteLine s];
        end
    end
    Info_components(deleteLine,:) = []; % 删去四舍五入后持续时间过少的状态
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
%                             2   8   16   info1;
%                             3   16   24   info1;
%                             4   24   32   info1;%恢复
%                             5   32   40   info4;
%                             6   40   48   info1;
%                             ];
%     Info_components2 =    [  1   0    8   info3;
%                             2   8   16   info1;
%                             3   16   24   info1;
%                             4   24   32   info1;%恢复
%                             5   32   40   info1;
%                             6   40   48   info1;
%                             ];
%     for ii = 1:size(testInfoComponents,1)
%         testmpc = mpcUpdateBinary(testInfoComponents(ii,4:end),mpc0,nGen);
% %         testmpc.gencost(34:50,6)=2600;
%         testmpc.Gcost=mpc0.Gcost * 1000;
%         testmpc.branch(:,6:8) = 9999;
% %         testmpc.gen(mpc0.gfuIndex,9) = mpc0.gen(mpc0.gfuIndex,9)/3;
%         testmpc.gencost(12:14,6)=16;%降低这个非燃气机组的价格，使燃气机组成为边际机组
%         % 为啥怎么调燃气机组都会差不多满发？
%         testResultSteady{ii} = GErunopf(testmpc);
%     end
    % ------------------------------------------
    % initialization
    NS = size(Info_components,1);
    newmpc = cell(NS,1);
    for s = 1:NS
        newmpc{s} = mpcUpdateBinary(Info_components(s,4:end),mpc0,nGen);
    end

    [LCe{i},LCg{i}] = lookAheadContingencyManagement_Solver(LaCMS_optimizer,NK,NS,mpc,...
       mpc0,newmpc,normalCondition,Info_components,nx);   
end
save stop3.mat
%%
nLCg = size(find(mpc.Gbus(:,3)~=0),1); nLCe = size(find(mpc.bus(:,3)~=0),1);
[systemEDNSsum,systemEGNSsum,systemLOLPsum,systemLOGPsum] = deal(zeros(NK,1));
[nodalEDNSsum,nodalEGNSsum,nodalLOLPsum,nodalLOGPsum] = deal(zeros(NK,nLCe),zeros(NK,nLCg),zeros(NK,nLCe),zeros(NK,nLCg));
for i = 1:simulationTimes
    nodalEDNSsum = nodalEDNSsum + LCe{i};nodalLOLPsum = nodalLOLPsum + double(LCe{i}>1e-3);
    nodalEGNSsum = nodalEGNSsum + LCg{i};nodalLOGPsum = nodalLOGPsum + double(LCg{i}>1e-5);
    systemEDNSsum = systemEDNSsum + sum(LCe{i},2);systemLOLPsum = systemLOLPsum + double(sum(LCe{i},2)>1e-3);
    systemEGNSsum = systemEGNSsum + sum(LCg{i},2);systemLOGPsum = systemLOGPsum + double(sum(LCg{i},2)>1e-5);

end
nodalEDNS = nodalEDNSsum/simulationTimes;nodalEGNS = nodalEGNSsum/simulationTimes;
systemEDNS = systemEDNSsum/simulationTimes;systemEGNS = systemEGNSsum/simulationTimes;
nodalLOLP = nodalLOLPsum/simulationTimes;nodalLOGP = nodalLOGPsum/simulationTimes;
systemLOLP = systemLOLPsum/simulationTimes;systemLOGP = systemLOGPsum/simulationTimes;  
% std
LCeAll = zeros(NK,simulationTimes);LCgAll = zeros(NK,simulationTimes);
for i = 1:simulationTimes
    LCeAll(:,i) = sum(LCe{i},2);LCgAll(:,i) = sum(LCg{i},2);
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
EENS = sum(nodalEDNS)'/4;%MWh
EVNS = sum(nodalEGNS)'/24/4;%Mm3
%%
% 消除尖峰
[systemEGNSnew] = flatten_manual(EGNS);
