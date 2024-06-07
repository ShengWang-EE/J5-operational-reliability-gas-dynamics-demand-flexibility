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
simulationTimes = 100000; % TSMCS 次数

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
    dayahead_IEGSresult_basicLoad{k} = GEresult0;
end
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
    for s = 1:NS
        newmpc{s} = mpcUpdateBinary(Info_components(s,4:end),mpc0,nGen);
    end
    initialCondition = setInitialConditionForGasSystem(dayahead_IEGSresult{1},gtd,nx,dx);
    [solution,nextInitialCondition] = lookAheadContingencyManagement(dt*3600,gtd,dx, NK,NS,mpc,...
       nodalEHpara,dayahead_IEGSresult_basicLoad,dayaheadEHschedule,mpc0,newmpc,initialCondition,normalCondition,Info_components,...
       nEH,nx);
    NS = size(Info_components2,1);
    for s = 1:NS
        newmpc{s} = mpcUpdateBinary(Info_components2(s,4:end),mpc0,nGen);
    end
   [solution,nextInitialCondition] = lookAheadContingencyManagement(dt*3600,gtd,dx, NK,NS,mpc,...
       nodalEHpara,dayahead_IEGSresult_basicLoad,dayaheadEHschedule,mpc0,newmpc,nextInitialCondition,normalCondition,Info_components2,...
       nEH,nx);
   
    % 结果整理gu
    [a] = resultProcessing(stateVar,mpc,ND/NK,nComponent+1,nx);
    % 根据新的系统状态更新初始条件
    initialCondition = setInitialConditionForGasSystem(operationResults.IEGS{kst},gtd,nx,dx);

    % 整理出与结果评估相关的一些指标，便于并行计算。包括：
    for k = 1:NK
%         ec{i}(k,i) = operationResults.IEGS
    end
end
            
            
            
%%         
        
    % ----- 2.2 
iK.DR = (53+1):1:59; % 冷负荷大于700，觉得设置为DR时间比较合适
NK.DR = iK.DR(end)-iK.DR(1)+1;
notifyTime = 48; %中午12点通知
endTime = 68; %下午17点下班，最多调整到17点
iK.studyPeriod = notifyTime:endTime;
NK.studyPeriod = endTime-notifyTime+1;
NK.nonDR = NK.studyPeriod - NK.DR;
% 每个EH的DR提供运行备用的需求是多少？假设负荷在日前的基础上负荷*120%，然后只有EH上可以削负荷决定
EHreserve = zeros(NK.all,nEH);
counter = 0;
for k = iK.DR(1):iK.DR(end)
    counter = counter+1;
    h = k2h_h2k(k,4,'k2h');
    % electricity load only
    realtime_mpc{k} = dayahead_mpc{k};
%     realtime_mpc{k}.bus(:,3:4) = dayahead_mpc{k}.bus(:,3:4) * 1.29;
    realtime_mpc{k}.bus(:,3:4) = dayahead_mpc{k}.bus(:,3:4) * 1.25;
    % 去掉其他负荷的DR能力，只有EH负荷能DR
    load_igen = find(mpc0.gen(:,22)==2);
    EH_igen = [];
    for i = 1:size(load_igen)
        if ~ismember(realtime_mpc{k}.gen(load_igen(i),1),mpc0.EHlocation(:,2)) % not EH
            realtime_mpc{k}.gen(load_igen(i),9) = 0; % if not EH bus, don't curtail loads
        else % is EH
            % 修改DR的上限为目前EH的ei
            EHindex = find(mpc.EHlocation(:,2) == realtime_mpc{k}.gen(load_igen(i),1));
            realtime_mpc{k}.gen(load_igen(i),9) = dayaheadEHschedule{EHindex}(h,1);
        end
    end
    [reserve_result{counter},flag] = GErunopf(realtime_mpc{k});
    ob.reserve(counter) = flag;
    for i = 1:nEH
        EH_igen = [EH_igen; find(mpc.gen(load_igen,1) == mpc.EHlocation(i,2))+load_igen(1)-1];% 把EH对应的负荷在mpc.gen里面的编号储存在EH_igen里面
    end
        
    EHreserve(k,:) = reserve_result{counter}.gen(EH_igen,2)';
end

% knowing the reserve requirement at each time period, 
[dx, nx] = decidePipelineCell(gtd);
dt = 900; % 15min 
% iK = endTime - notifyTime;
iter = 1;
UB(iter) = 10e20; LB(iter) = 0;
addCutsCoefficient = [];SP_EHresult = [];exitflagSP = [];
% control variables: 
[MPformulation] = MP_IEGSdispatch_optimizer(addCutsCoefficient,exitflagSP,SP_EHresult,...
    dayahead_IEGSresult_basicLoad,dayaheadEHschedule,EHreserve,mpc,gtd,nx,dx,dt,NH,nEH,NK,iK,nodalPrice,iter);
for i = 1:nEH
    [SPformulation{i}] = SP_EHschedule_lp_optimizer([],[],nodalEHpara{i},NK,iK);
end
while abs(UB(iter)-LB(iter)) / abs(UB(iter)+LB(iter)) > 1e-2    
%     yalmip('clear');
    iter = iter + 1;
%     load periter.mat
    [MP_IEGSresult, solverTime,MPformulation] = MP_IEGSdispatch_solver(addCutsCoefficient,exitflagSP,SP_EHresult,...
    dayahead_IEGSresult_basicLoad,dayaheadEHschedule,EHreserve,mpc,gtd,nx,dx,dt,NH,nEH,NK,iK,nodalPrice,iter,MPformulation);


%     save MP.mat
%     %%%%%%%%%%%%%%%%
%     clear
%     yalmip('clear')
%     load MP20200615.mat
    LB(iter) = MP_IEGSresult.f_MP;
    for i = 1:nEH
        % 注意MP得到的ei，gi都是标幺化后的
        ei_hat = MP_IEGSresult.ei(:,i) * mpc.baseMVA; gi_hat = MP_IEGSresult.gi(:,i) * 200; 
        [SP_EHresult{i},addCutsCoefficient{i},SP_solverTime(iter,i),exitflagSP(iter,i)] = SP_EHschedule_lp_solver(SPformulation{i},ei_hat,gi_hat,nodalEHpara{i},NK,iK);
    end
    if min(exitflagSP(iter,:)) == 1 % if all EH schedule is feasible, then update upper bound
        UB(iter) = MP_IEGSresult.f_MP;
        for i = 1:nEH
            UB(iter) = UB(iter) + SP_EHresult{i}.EHcost;
        end
    else
        UB(iter) = UB(iter-1);
    end
%     UB(iter) = min(UB(iter),UB(iter-1));

% save periter.mat
end
totalSolverTime = sum(MP_solverTime) + sum(sum(SP_solverTime));

%% visualize the results
% 1 computation speed
% SP_solverTime = 10*SP_solverTime;
computationTime.MPsolver = sum(MP_solverTime);
computationTime.SPsolver = sum(SP_solverTime);
computationTime.SPwaiting = sum(max(SP_solverTime,[],2) - SP_solverTime);
computationTime.total = totalSolverTime;
UB1(1) = UB(1); LB1(1) = LB(1);
for i = 2:iter
    UB1(i) = min(UB(i),UB1(i-1));
    LB1(i) = max(LB(i),LB1(i-1));
end
% 2 self-schedule
load1 = nodalEHpara{1}.load(:,iK.DR)' - SP_EHresult{end,1}.so;
SP_EHresult{end,1}.h3h + SP_EHresult{end,1}.h1h + SP_EHresult{end,1}.h2h;
nodalEHpara{1,1}.load(3,iK.studyPeriod(1):iK.studyPeriod(end))' - SP_EHresult{end,1}.c3c - SP_EHresult{end,1}.c4c;
% operating condition

% 3 IEGS
% pressure:
prsInfo = cell(nGl,1);
for i = 1:nGl
    prsInfo{i} = [];
    for j = 1:size(MP_IEGSresult{end}.PrsPipe{i},2)
        prsInfo{i} = [prsInfo{i};max(MP_IEGSresult{end}.PrsPipe{i}(:,j)),min(MP_IEGSresult{end}.PrsPipe{i}(:,j)),mean(MP_IEGSresult{end}.PrsPipe{i}(:,j))];
    end
end
    
% cost

