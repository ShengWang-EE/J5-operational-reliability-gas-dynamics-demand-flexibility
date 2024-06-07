function totalCost = objfcn_SP(lc,CDF,NK)
curtailmentCost = 0; 
for l = 1:3 % 3 energies
    for ift = 1:NK % from time
        curtailmentCost = curtailmentCost + lc(ift,l) * CDF;
    end
end
totalCost = curtailmentCost; % sum all
end