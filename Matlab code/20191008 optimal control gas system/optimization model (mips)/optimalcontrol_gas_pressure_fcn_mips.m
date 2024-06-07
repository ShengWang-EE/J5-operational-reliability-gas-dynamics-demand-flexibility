function [g,dg] = optimalcontrol_gas_pressure_fcn_mips(x,pp,mpc,NK,nx)
nGb = size(mpc.Gbus,1); nGl = size(mpc.Gline,1);
[Prs] = deal(x{:});
[nPrs] = deal(size(Prs,1));
cumnvar = cumsum([nPrs]);
nvar = nPrs;
[Prs0] = deal(reshape(Prs,[],NK));%ÁÐÊÇk
g = zeros(NK*(2*nGl-nGb),1); dg = zeros(NK*(2*nGl-nGb),nvar);[g0,dg0] = deal(cell(NK,1));
for k = 1:NK
    Prs = Prs0(:,k);
%     [g((k-1)*(2*nGl-nGb)+1:k*(2*nGl-nGb)),  dg((k-1)*sum(nx+1)+1:k*sum(nx+1),(k-1)*(2*nGl-nGb)+1:k*(2*nGl-nGb))] ...
%         = opf_gas_pressure_fcn(Prs,pp,mpc);
    [g0{k},dg0{k}] = opf_gas_pressure_fcn(Prs,pp,mpc,nx);
    g((k-1)*(2*nGl-nGb)+1:k*(2*nGl-nGb)) = g0{k};   
end
for i = 1:nvar
    itype = min(find(cumnvar>=i));
    k = mod(i,NK);
    if k==0
        k=NK;
    end
    dg((k-1)*(2*nGl-nGb)+1:k*(2*nGl-nGb),1) = dg0{k}(:,itype);
end

end
function [f_equalP,df_equalP] = opf_gas_pressure_fcn(Prs,pp,mpc,nx)
nPrs = size(Prs,1);
nGl = size(mpc.Gline,1);
nGb = size(mpc.Gbus,1);
iPrs = mat2cell(1:nPrs,1,nx+1);
Prs1 = zeros(nGl,2);

Prs1(:,1) = Prs(pp(:,1));
Prs1(:,2) = Prs(pp(:,2));

counter1 = 1; 
f_equalP = zeros(2*nGl-nGb,1);
df_equalP = zeros(nPrs,2*nGl-nGb);
for i = 1:nGb
    % equal P
    connectLineAsFromBus = find(mpc.Gline(:,1)==i);
    connectLineAsToBus = find(mpc.Gline(:,2)==i);
    connectLine = [connectLineAsFromBus; connectLineAsToBus];

    if size(connectLine,1) > 1
        for j = 2:size(connectLineAsFromBus,1) %take the first P in from line as a base
            f_equalP(counter1) = Prs1(connectLineAsFromBus(1),1) - Prs1(connectLineAsFromBus(j),1);
            df_equalP(iPrs{connectLineAsFromBus(1)}(1),counter1) = 1;
            df_equalP(iPrs{connectLineAsFromBus(j)}(1),counter1) = -1;
            counter1 = counter1 + 1;
        end
        for j = 1:size(connectLineAsToBus,1)
            f_equalP(counter1) = Prs1(connectLineAsFromBus(1),1) - Prs1(connectLineAsToBus(j),2);
            df_equalP(iPrs{connectLineAsFromBus(1)}(1),counter1) = 1;
            df_equalP(iPrs{connectLineAsToBus(j)}(end),counter1) = -1;
            counter1 = counter1 + 1;
        end  
    end
end
df_equalP = df_equalP';
end