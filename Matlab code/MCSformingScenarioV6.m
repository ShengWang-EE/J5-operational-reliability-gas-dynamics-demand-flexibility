function [allComponentsInfo, failureComponents] = MCSformingScenarioV6(rts,missionTime)
% 把大元件分解为小元件，全是两状态。
%% Info gen
numberOfGenAndGfu = size(rts.gen,1);  numberOfGsou = size(rts.Gsou,1);

numberOfComponents = numberOfGenAndGfu +  numberOfGsou;
allComponentsInfo = [];info1 = []; info2 = [];
%% conponents info order gen和gfu穿插，gas source
for i = 1:numberOfComponents
    % calculate lamda and mu according to different types of components
    %components order: gen(tfu,gfu), Gsou
    if i <= numberOfGenAndGfu
        lamda = 1/rts.gen(i,1); mu = 1/rts.gen(i,2);
        unitNumber = 1;%two states
        [Info_Component] = MCS_twoStates(missionTime,unitNumber,lamda,mu);
    else % is Gsou
        lamda = 1/rts.Gsou(i-numberOfGenAndGfu,1); mu = 1/rts.Gsou(i-numberOfGenAndGfu,2);
        unitNumber = 1; 
        [Info_Component] = MCS_twoStates(missionTime,unitNumber,lamda,mu);
    end
   

%% merge info1 and info2
    info2 = Info_Component;
    if i>1
        ok = 1;
        pointer1 = 1;  %% for the first unit, state pointer
        pointer2 = 1;  %% for the second unit, state pointer
        iter =1;
        StateNumber1 = size(info1,1);  %经历的状态个数   
        StateNumber2 = size(info2,1); 
        componentIndex = i+2;%info1的列数
        while ok 
          if iter ==1
             allComponentsInfo(iter,1)=1;  %% State number
             allComponentsInfo(iter,2)=0;  %% Starting time
          end   

          if info1(pointer1,3) <= info2(pointer2,3) %用于判定到底是那两个状态叠加，不是枚举，是根据unit实际模拟情况按时间顺序叠加

             if pointer1 == StateNumber1  %% Point of the first unit moves to the last state 

                  allComponentsInfo(iter,3)=info1(pointer1,3);  %% This is the state end time  
                  allComponentsInfo(iter,4:componentIndex)=info1(pointer1,4:componentIndex);
                  allComponentsInfo(iter,(componentIndex+1)) = info2(pointer2,4);  %% info1有很多行状态，info2只有一行
                  ok=0;  

             else 
                  allComponentsInfo(iter,3)=info1(pointer1,3);
                  allComponentsInfo(iter,4:componentIndex)=info1(pointer1,4:componentIndex);
                  allComponentsInfo(iter,(componentIndex+1)) = info2(pointer2,4);
                  pointer1 = pointer1+1;

                  iter = iter +1;
                  allComponentsInfo(iter,1)=iter;%状态等于下标
                  allComponentsInfo(iter,2)= info1(pointer1,2);
             end     


          else
             if pointer2 == StateNumber2  %% Point of the second unit moves to the last state

                 allComponentsInfo(iter,3)=info2(pointer2,3);   %% This is the state end time   
                 allComponentsInfo(iter,4:componentIndex)=info1(pointer1,4:componentIndex);
                 allComponentsInfo(iter,(componentIndex+1)) = info2(pointer2,4);
                 ok=0; 
              else
                 allComponentsInfo(iter,3)=info2(pointer2,3);
                 allComponentsInfo(iter,4:componentIndex)=info1(pointer1,4:componentIndex);
                 allComponentsInfo(iter,(componentIndex+1)) = info2(pointer2,4);
                 pointer2 = pointer2+1;


                 iter = iter +1;
                 allComponentsInfo(iter,1)=iter;
                 allComponentsInfo(iter,2)= info2(pointer2,2);
             end   
          end    
        end
        info1 = allComponentsInfo;
    end  
    if i == 1
        info1 = Info_Component;
    end
end
allComponentsInfo(end,3) = missionTime;
%% 
% failureComponents.genIndex = cell(size(allComponentsInfo,1),1);
% failureComponents.gasSourceIndex = cell(size(allComponentsInfo,1),1);
% allComponentsInfoForGen = allComponentsInfo(:,4:3+numberOfGenAndGfu);
% allComponentsInfoForGasSource = allComponentsInfo(:,4+numberOfGenAndGfu:end);
% for i = 1:size(allComponentsInfo,1)
%     failureComponents.genIndex{i} = find(allComponentsInfoForGen(i,:)< allComponentsInfoForGen(1,:));
%     failureComponents.gasSourceIndex{i} = find(allComponentsInfoForGasSource(i,:)<allComponentsInfoForGasSource(1,:));
% end
