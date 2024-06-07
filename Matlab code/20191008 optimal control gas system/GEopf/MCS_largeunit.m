function [Info_largeunit] = MCS_largeunit(mission_time)

%   Detailed explanation goes here
%input:mission_time,
%matrix:largeunitInfo(unit_no,unit_capacity_state,Lamda,

Unit_No= 4;  %% the number of multistate unit 
unit_capacity_state  =  [0,247,482,575]; 

Lamda = [-0.0933, 0.0800, 0.0133, 0;
          0.0294,-0.3823, 0.3235, 0.0294;
               0, 0.0288, -0.3846,0.3558;
           0.0002, 0.0001, 0.0007,-0.0010];
unitstateno = 4; 


%Suppose it starts from the best state: state 4
Info_largeunit(1,1)= 4;    %% State number
Info_largeunit(1,2)= 0;    %% State start time
Info_largeunit(1,3)= 0;    %% State end time
Info_largeunit(1,4)= 575;  %% State capacity
Total_Dur = 0;
iter =1;

while Total_Dur < mission_time 
    
    for i = 1:4 
       
      a = log(rand(1));
      state_no = Info_largeunit(iter,1);
      b =  Lamda(state_no,i);  %% Transition rate 
      Dur(i) = - a /b;
       if Dur(i) <0   
           Dur(i) = 10^10; 
       end
    end
    
   [stateDur, Nextstate] = min(Dur);
   Info_largeunit(iter,3) = Info_largeunit(iter,2) + stateDur;  
   Total_Dur = Info_largeunit(iter,3);
   
   if Total_Dur < mission_time 
    iter = iter +1;
    Info_largeunit(iter,1)= Nextstate;
    Info_largeunit(iter,2)= Info_largeunit((iter-1),3);    %% State start time
    Info_largeunit(iter,3)= 0;    %% State end time
    Info_largeunit(iter,4)= unit_capacity_state(Nextstate);    %% State capacity
   end
   
end

a =1; 









