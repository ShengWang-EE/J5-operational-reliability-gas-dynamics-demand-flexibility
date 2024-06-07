function [result] = opfWithoutConstraint(mpc)
% only curtailments are concerned (not caring about the generators, or gas
% sources)
sumOfGen = sum(mpc.gen(1:mpc.originalGenNumber,9));
sumOfGas = sum(mpc.Gsou(:,4));
sumOfElectricityLoad = sum(mpc.bus(:,3));
sumOfGasLoad = sum(mpc.Gbus(:,3));

electricityLC = zeros(size(mpc.gen,1)-mpc.originalGenNumber,1);
gasLC = zeros(size(mpc.Gbus(:,3)));
UBelectricityLC = mpc.gen(mpc.originalGenNumber+1:end,9);
UBgasLC = mpc.Gbus(:,3);
% no upboundaries for this uncervergecible opf
if sumOfGen < sumOfElectricityLoad
    % electricity curtailment is needed
    electricityLC = (sumOfElectricityLoad - sumOfGen) / sum(UBelectricityLC) * UBelectricityLC;
end
if sumOfGas < sumOfGasLoad
    gasLC = (sumOfGasLoad - sumOfGas) / UBgasLC * UBgasLC;
end

result = mpc;
result.gen(mpc.originalGenNumber+1:end,2) = electricityLC;
result.Gbus(:,10) = gasLC;
% result.Gsou(:,5) = result.Gsou(:,2);
end