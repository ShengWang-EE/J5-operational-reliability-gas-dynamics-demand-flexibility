function [nodalEHpara] = scaleEH(mpc,loadCurve,EHpara)
% [HA,EA,HB,EB,HC,EC,HD,ED,HO2MAX,HO2MIN,HO3MAX,HO3MIN,HO4MAX,HO4MIN,CO4MAX,CO4MIN,CO5MAX,CO5MIN,a,b,c,d,e,f,eta1e,eta1h,eta2,eta3,COP4hSummer,COP4hWinter,COP4c,COP5,alpha,beta,CDF] = unpackEHpara(EHpara);
electricityMax = max(loadCurve(1,:));
% EHboundary = [HA,EA,HB,EB,HC,EC,HD,ED,HO2MAX,HO2MIN,HO3MAX,HO3MIN,HO4MAX,HO4MIN,CO4MAX,CO4MIN,CO5MAX,CO5MIN];
nEH = size(mpc.EHlocation,1);
nodalEHpara = cell(nEH,1);
% convert into MW
% a = a*1000; c = c*1000; e = e*1000; f = f/1000;
for i = 1:nEH
    iEH = mpc.EHlocation(i,2);%the electricity bus of EH
    ratio = mpc.bus(iEH,3) / electricityMax * 0.25; % %认为EH最大电力终端负荷占节点总电力负荷的1/4
    nodalEHpara{i}.load = loadCurve * ratio;
    nodalEHpara{i}.capacity = EHpara.capacity * ratio;
    nodalEHpara{i}.efficiency = EHpara.efficiency;
    nodalEHpara{i}.schedule = EHpara.schedule;
end

end