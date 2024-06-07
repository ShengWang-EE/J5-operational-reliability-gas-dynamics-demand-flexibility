function [Info_twostateplant] = MCS_gasSource(mission_time)
% use original multiple two state plant 
% requirement:the same capacity and failure rate
%   Detailed explanation goes here

lamda_unit = 1/960; 
Mu_unit =    1/40;
twostateplant_rates = [3*lamda_unit, 0;
                       2*lamda_unit, Mu_unit;
                       lamda_unit, 2*Mu_unit;
                                  0, 3*Mu_unit;]; %ÏµÍ³²ãÃæ¿¼ÂÇ
twostateplant_state  =  [6,4,2,0];    %% Represent how many units are available for each state 


%Suppose it starts from the best state: state 1, 
Info_twostateplant(1,1)= 1;    %% State number
Info_twostateplant(1,2)= 0;    %% State start time
Info_twostateplant(1,3)= 0;    %% State end time, this value is unknown
Info_twostateplant(1,4)= 6;    %% capacity in this state 
Total_Dur = 0;
iter =1;
                              
                              
while Total_Dur < mission_time 
  
  state_no = Info_twostateplant(iter,1);  
  
  if  twostateplant_rates(state_no,1) >0 
    a1 = log(rand(1));
    b1 =  twostateplant_rates(state_no,1);  %% Transition rate 
    Dur1 = - a1 /b1;%Ö»ÓÐ2×´Ì¬£¬¿¼ÂÇËð»µ
  else    
    Dur1 = 10^10;      
  end    
    
  if  twostateplant_rates(state_no,2) >0 
    a2 = log(rand(1));
    b2 =  twostateplant_rates(state_no,2);  %% Transition rate 
    Dur2 = - a2 /b2;%ÐÞ¸´
  else    
    Dur2 = 10^10;      
  end 
  
  if Dur1 < Dur2  %% It will go to the next state  
    
    Info_twostateplant(iter,3) = Info_twostateplant(iter,2) + Dur1;  
    Total_Dur = Info_twostateplant(iter,3);   
   
    if Total_Dur < mission_time   
      iter = iter +1;
      Info_twostateplant(iter,1)= Info_twostateplant((iter-1),1)+1;    %% State number, go to next state
      Info_twostateplant(iter,2)= Info_twostateplant((iter-1),3);    %% State start time
      Info_twostateplant(iter,3)= 0;    %% State end time, this value is unknown
      Info_twostateplant(iter,4)= twostateplant_state(Info_twostateplant(iter,1));    %% available capacity  
    end
   
   
  else            %% It will go to the previous state
    
    Info_twostateplant(iter,3) = Info_twostateplant(iter,2) + Dur2;  
    Total_Dur = Info_twostateplant(iter,3);   
   
    if Total_Dur < mission_time   
      iter = iter +1;
      Info_twostateplant(iter,1)= Info_twostateplant((iter-1),1)-1;    %% State number, go to previous state
      Info_twostateplant(iter,2)= Info_twostateplant((iter-1),3);    %% State start time
      Info_twostateplant(iter,3)= 0;    %% State end time, this value is unknown
      Info_twostateplant(iter,4)= twostateplant_state(Info_twostateplant(iter,1));    %% available capacity  
    end  
             
  end

Info_twostateplant(end,3) = mission_time;
end 









