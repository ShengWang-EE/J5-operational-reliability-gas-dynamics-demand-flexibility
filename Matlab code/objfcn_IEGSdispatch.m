function [totalCost,sumGenAndLCeCost,sumGasPurchasingCost,sumGasCurtailmentCost] = objfcn_IEGSdispatch(Pg,PGs,LCg,mpc,NK)

% totalCost = electricity generation cost + gas purchasing cost + CDFe + CDFg
% 
CDFg = calculateGasCDF(mpc)*1;
Pg = mpc.baseMVA*Pg;

    genAndLCeCost = sum(Pg * mpc.gencost(:,6));
%     sum(totcost_yalmip(mpc.gencost, Pg(:,:)')); % Pg includes GFU, but the gencost are 0.
    gasPurchasingCost = sum(PGs * mpc.Gcost);
    gasCurtailmentCost = sum(LCg) * CDFg;

sumGenAndLCeCost = sum(genAndLCeCost)/4;
sumGasPurchasingCost = sum(gasPurchasingCost)/4;%?15min
sumGasCurtailmentCost = sum(gasCurtailmentCost)/4;
totalCost = sumGenAndLCeCost + sumGasPurchasingCost + sumGasCurtailmentCost;


end



