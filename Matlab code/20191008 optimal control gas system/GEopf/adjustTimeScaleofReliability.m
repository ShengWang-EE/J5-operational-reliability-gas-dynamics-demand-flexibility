function [operationalReliability] = adjustTimeScaleofReliability(InfoReliability,numberOfPeriods,missionTime)
% LOLP 和EENS不同，LOLP是用checkpoint来，只检查某一时间点上是否失负荷，而EENS需要统计一个区间上的积累量
% LOLP从1开始好了，而不是从0~1，因为t=0，LOLP肯定为0
% the format of Info_reliability: [from time, to time, LOLP/LC]
% the format of operationalReliability: [new from time, new to time, newLOLP/newEENS]
%% v1
% operationalReliability = zeros(numberOfPeriods,3);
% operationalReliability(:,1) = [0:1:numberOfPeriods-1]' * missionTime/numberOfPeriods; 
% operationalReliability(:,2) = operationalReliability(:,1) + missionTime/numberOfPeriods;
% for i = 1:numberOfPeriods
%     switch LOLPorEENS
%         case 'LOLP'
%             % about LOLP
%             indexj = max(find(InfoReliability(:,1)<=operationalReliability(i,2)));
%             if InfoReliability(indexj,3)>threshold % gas LC
%                 operationalReliability(i,3) = 1;
%             end
%         case 'EENS'
%             % about EENS
%             ft = operationalReliability(i,1);% from time
%             tt = operationalReliability(i,2);% to time
%             indexj_ft = max(find(InfoReliability(:,1)<=ft));
%             indexj_tt = max(find(InfoReliability(:,1)<=tt));
%             % if the period is completely contained in time step (not likely to happen)
%             if indexj_ft == indexj_tt
%                 ENS = InfoReliability(indexj_ft,3) * (tt - ft); % LC accumulated by hour
%             elseif indexj_tt - indexj_ft == 1
%                 ENS = (InfoReliability(indexj_ft,2) - ft) * InfoReliability(indexj_ft,3)...
%                     + (tt - InfoReliability(indexj_tt,1)) * InfoReliability(indexj_ft,3);
%             elseif indexj_tt - indexj_ft > 1
%                 ENS = (InfoReliability(indexj_ft,2) - ft) * InfoReliability(indexj_ft,3)...
%                     + (tt - InfoReliability(indexj_tt,1)) * InfoReliability(indexj_ft,3);
%                 for k = (indexj_ft+1):(indexj_tt-1)
%                     ENS = ENS + (InfoReliability(k,2) - InfoReliability(k,1)) * InfoReliability(k,3);
%                 end                
%             else
%                 warnning('!!!!!!');
%             end
%     
%         operationalReliability(i,3) = ENS/((missionTime/numberOfPeriods));
%     end
% end
%% 20190709 v2
operationalReliability = zeros(numberOfPeriods,3);
operationalReliability(:,1) = [0:1:numberOfPeriods-1]' * missionTime/numberOfPeriods; 
operationalReliability(:,2) = operationalReliability(:,1) + missionTime/numberOfPeriods;
for i = 1:numberOfPeriods

            % about EENS
            ft = operationalReliability(i,1);% from time
            tt = operationalReliability(i,2);% to time
            indexj_ft = max(find(InfoReliability(:,1)<=ft));
            indexj_tt = max(find(InfoReliability(:,1)<=tt));
            % if the period is completely contained in time step (not likely to happen)
            if indexj_ft == indexj_tt
                ENS = InfoReliability(indexj_ft,3) * (tt - ft); % LC accumulated by hour
            elseif indexj_tt - indexj_ft == 1
                ENS = (InfoReliability(indexj_ft,2) - ft) * InfoReliability(indexj_ft,3)...
                    + (tt - InfoReliability(indexj_tt,1)) * InfoReliability(indexj_tt,3);
            elseif indexj_tt - indexj_ft > 1
                ENS = (InfoReliability(indexj_ft,2) - ft) * InfoReliability(indexj_ft,3)...
                    + (tt - InfoReliability(indexj_tt,1)) * InfoReliability(indexj_tt,3);
                for k = (indexj_ft+1):(indexj_tt-1)
                    ENS = ENS + (InfoReliability(k,2) - InfoReliability(k,1)) * InfoReliability(k,3);
                end                
            else
                warnning('!!!!!!');
            end
    
        operationalReliability(i,3) = ENS/((missionTime/numberOfPeriods));

end
end