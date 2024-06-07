function rts = Case24ReliabillityDatav3()
%%---------  Reliability Data for the IEEE-RTS-------- %%
%% generation unit data
%1BUS %2Forced Outage Rate %3MTTF(hour) %4MTTR(hour) %5Scheduled
%Maintenance (wks/year) % 6 number of states
rts.gen=[
450 50
450 50
1960 40
1960 40
450 50  %5
450 50 
1960 40 
1960 40 
1200 50 
1200 50  %10
1200 50 
950 50 
950 50 
950 50 
999999 0.0000001  %15
2940 60 
2940 60 
2940 60 
2940 60 
2940 60  %20
960 40 
960 40 
1100 150 
1100 150 
1980 20 %25
1980 20 
1980 20 
1980 20 
1980 20 
1980 20  %30
960 40 
960 40 
1150 100 
    ];

%% gas source failure data
% suppose it is similiar with gen
%等分的划成multistate.第三列是小气源的个数（一个小气源2Mm3/day左右比较好）
rts.Gsou = [
1960 40
1960 40
1960 40
1960 40
1960 40
1200 50
1200 50
1200 50
1200 50
2940 60 
2940 60 
9999999999999999999999999 0.00000000000000000000000000001 % regarded as never broken
960 40 
1980 20 ];
rts.Gsou(:,1) = rts.Gsou(:,1)/2;% double the failure rate


end