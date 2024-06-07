function [g,dg] = optimalcontrol_gas_pressure_fcn(p0,pp,mpc,NK)
for k = 1:NK
    p = p0(k,:)';
    [g(k,:)] = opf_gas_pressure_fcn(p,pp,mpc)';
%     [~,dg{k}] = opf_gas_pressure_fcn(p,pp,mpc);
end
end
function [f_equalP,df_equalP] = opf_gas_pressure_fcn(p,pp,mpc)
nGl = size(mpc.Gline,1);
nGb = size(mpc.Gbus,1);
% iPrs = mat2cell(vv.i1.Prs:vv.iN.Prs,1,nx+1);
Prs = zeros(nGl,2);
for i = 1:nGl
    Prs(:,1) = p(pp(:,1));
    Prs(:,2) = p(pp(:,2));
end
counter1 = 1; f_equalP = zeros(2*nGl-nGb,1);
% df_equalP = zeros(vv.N.all,2*nGl-nGb);
for i = 1:nGb
    % equal P
    connectLineAsFromBus = find(mpc.Gline(:,1)==i);
    connectLineAsToBus = find(mpc.Gline(:,2)==i);
    connectLine = [connectLineAsFromBus; connectLineAsToBus];

    if size(connectLine,1) > 1
        for j = 2:size(connectLineAsFromBus,1) %take the first P in from line as a base
            f_equalP(counter1) = Prs(connectLineAsFromBus(1),1) - Prs(connectLineAsFromBus(j),1);
%             df_equalP(iPrs{connectLineAsFromBus(1)}(1),counter1) = 1;
%             df_equalP(iPrs{connectLineAsFromBus(j)}(1),counter1) = -1;
            counter1 = counter1 + 1;
        end
        for j = 1:size(connectLineAsToBus,1)
            f_equalP(counter1) = Prs(connectLineAsFromBus(1),1) - Prs(connectLineAsToBus(j),2);
%             df_equalP(iPrs{connectLineAsFromBus(1)}(1),counter1) = 1;
%             df_equalP(iPrs{connectLineAsToBus(j)}(end),counter1) = -1;
            counter1 = counter1 + 1;
        end  
    end
end
end