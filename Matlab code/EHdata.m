function [summerLoad1,EHpara] = EHdata()
%% data input
% winterLoad = [  58.5	70.2	81.9	105.3	117	140.4	198.9	245.7	257.4	245.7	163.8	152.1	140.4	152.1	163.8	187.2	198.9	234	256.5	222.3	187.2	140.4	105.3	58.5
%                 192.2	207.6	239.3	253.2	383.4	253.2	336.8	398.5	504	523.2	492.7	520.6	613.6	540	582.2	444	360	300.5	286.4	270.6	242.4	227.9	163.2	208.2
%                 23.2	24.1	23.8	24.4	25.2	27.5	40.5	45	49	54.4	67.5	62.5	70	65	55	43.5	35.5	26.5	24.1	23.8	24.3	25.2	21	21.5];
% summerLoad = [  58.5	70.2	81.9	105.3	117	140.4	198.9	245.7	257.4	245.7	163.8	152.1	140.4	152.1	163.8	187.2	198.9	234	256.5	222.3	187.2	140.4	105.3	58.5
summerLoad = [58.5000000000000,70.2000000000000,81.9000000000000,105.300000000000,117,140.400000000000,198.900000000000,245.700000000000,257.400000000000,245.700000000000,183.800000000000,177.100000000000,170.400000000000,177.100000000000,183.800000000000,197.200000000000,203.900000000000,234,256.500000000000,222.300000000000,187.200000000000,140.400000000000,105.300000000000,58.5000000000000;
                35.5	22.9	15.9	15.9	25.5	63.7	57.3	25.5	79.6	101.9	108.3	76.4	89.2	89.2	106.4	75.2	67.1	77.4	91.4	126.3	119.2	91.2	56.1	35.1
                212.4	222	202.7	173.7	135.1	193.1	154.4	270.3	463.3	559.8	637.1	675.7	666	700.6	748.1	637.1	579.2	540.5	386.1	328.2	308.9	289.6	270.3	251];
% summer load normalization
% summerLoad(3,:) = summerLoad(3,:) / max(summerLoad(3,:)) * 380;

[HA,EA,HB,EB,HC,EC,HD,ED] = deal(0,250,110,210,90,50,0,100);
[HO2MAX,HO2MIN,HO3MAX,HO3MIN,CO3MAX,CO3MIN,CO4MAX,CO4MIN] = deal(300,0,600,0,600,0,300,0);
% [HO2MAX,HO2MIN,HO3MAX,HO3MIN,CO3MAX,CO3MIN,CO4MAX,CO4MIN] = deal(250,0,450,0,450,0,300,0);
% [a,b,c,d,e,f] = deal(0.00216,0.90625,0.00188,0.26250,0.00188,16.56);
[eta1e,eta1h,eta2,COP3hSummer,COP3hWinter,COP3c,COP4] = deal(0.3,0.4,0.8,4,3,3,0.7);
% nOMEGA = 13;
% load curtailment data
[alpha,beta] = deal(0.2,0.1);
[CDF] = deal(22600);% 22.6$/kWh
% % scale the load curve
% loadCurve = summerLoad / max(summerLoad(1,:));
%% price
% electricityPriceAll.time = [1  6   11  15  19  22  24];
% electricityPriceAll.price = [40  60  80    60  80  60  40];
% gasPrice = 48;
%%
% para.HA = HA;para.EA = EA;para.HB = HB; para.EB = EB; para.HC = HC; para.EC = EC; para.HD = HD; para.ED = ED;
% para.HO2MAX = HO2MAX; para.HO2MIN = HO2MIN; 
% para.HO3MAX = HO3MAX; para.HO3MIN = HO3MIN; para.CO3MAX = CO3MAX; 
% para.CO3MIN = CO3MIN; para.CO4MAX = CO4MAX; para.CO4MIN = CO4MIN; 
% % para.a = a; para.b = b; para.c = c; para.d = d; para.e = e; para.f = f; 
% para.eta1e = eta1e; para.eta1h = eta1h; para.eta2 = eta2; para.eta3 = eta3; para.COP4hSummer = COP4hSummer; 
% para.COP4hWinter = COP4hWinter; para.COP4c = COP4c; para.COP5 = COP5;
% para.alpha = alpha; para.beta = beta;
% para.CDF = CDF;
EHpara.capacity = [HA,EA,HB,EB,HC,EC,HD,ED,HO2MAX,HO2MIN,HO3MAX,HO3MIN,CO3MAX,CO3MIN,CO4MAX,CO4MIN];
EHpara.efficiency = [eta1e,eta1h,eta2,COP3hSummer,COP3hWinter,COP3c,COP4];
EHpara.schedule = [alpha, beta, CDF];
%% swap electircity and cooling load pattern
Emean = mean(summerLoad(1,:));Cmean = mean(summerLoad(3,:));
summerLoad1 = summerLoad;
summerLoad1(1,:) = summerLoad(3,:)/Cmean*Emean;
summerLoad1(3,:) = summerLoad(1,:)/Emean*Cmean;
end

