function [J,dJ,d2J] = objfcn_LaCMS_mips(x,mpc,possibleResults, probability,nx,NK,baseMVA)
CDFg = calculateGasCDF(mpc);
[Pg,Prs,PGs,LCg] = deal(x{:});
[nPg,nPrs,nPGs,nLCg] = deal(size(Pg,1),size(Prs,1),size(PGs,1),size(LCg,1));
nvar = nPg+nPrs+nPGs+nLCg;
Nx = nvar/NK;
[Pg,Prs0, PGs, LCg] = deal(reshape(Pg,[],NK),reshape(Prs,[],NK),reshape(PGs,[],NK),reshape(LCg,[],NK));%ÁÐÊÇk
% obtain the nodal pressure
nGb = size(mpc.Gbus,1);
nGl = size(mpc.Gline,1);
PrsNodal = zeros(nGb,1);
lastPrs = mat2cell(Prs0(:,NK),nx+1,1);
for i = 1:nGl
    PrsNodal(mpc.Gline(i,1)) = lastPrs{i}(1);
    PrsNodal(mpc.Gline(i,2)) = lastPrs{i}(end);
end

[J_dynamic,J_endpoint,~,~,~] = LaCMScost(Pg,LCg,PGs,PrsNodal,mpc,possibleResults,probability,NK);

J = J_dynamic + J_endpoint;

%% dJ
[dJ] = deal(zeros(Nx,1));
iPg = 1:nPg;
iPrs = mat2cell(nPg+nPrs-sum(nx+1)+1:nPg+nPrs,1,nx+1);
iPGs = 1+nPg+nPrs:nPg+nPrs+nPGs;
iLCg = 1+nPg+nPrs+nPGs:nPg+nPrs+nPGs+nLCg;


for k = 1:NK
    for i = 1:nPg/NK
        dJ(iPg((k-1)*nPg/NK+i)) = baseMVA * polycost(mpc.gencost(i,:), Pg(i,k), 1);
    end
    for i = 1:nLCg/NK
        dJ(iLCg((k-1)*nLCg/NK+i)) = CDFg;
    end
    for i = 1:nPGs/NK
        dJ(iPGs((k-1)*nPGs/NK+i)) = mpc.Gcost(i);
    end
end
d2J = zeros(nvar);
%     
end