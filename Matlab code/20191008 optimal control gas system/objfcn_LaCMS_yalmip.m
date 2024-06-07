function [J] = objfcn_LaCMS_yalmip(Pg,PrsRaw,PGs,LCg,mpc,vv,possibleResults, probability,nx,NK,baseMVA)
% obtain the nodal pressure
nGb = size(mpc.Gbus,1);
nGl = size(mpc.Gline,1);

p = mat2cell(PrsRaw(NK,:),1,nx+1);
for i = 1:nGl
    Prs(mpc.Gline(i,1)) = p{i}(1);
    Prs(mpc.Gline(i,2)) = p{i}(end);
end

[J_dynamic,J_endpoint,~,~,~] = LaCMScost(Pg,LCg,PGs,Prs,mpc,possibleResults,probability,NK);

J = J_dynamic + J_endpoint;
J = J_dynamic;
J = sum(sum(LCg));
% J = -sum(sum(Pg)); % 经测试，不是目标函数的问题
%% dJ
% [dJ,dJdPg,dJdLCg,dJdPGs] = deal(zeros(nvar,1));
% iPg = zeros(NK,vv.N.Pg);iLCg = zeros(NK,vv.N.LCg);iPGs = zeros(NK,vv.N.PGs);
% 
% for k = 1:NK
%     iPg(k,:) = Nx*(k-1) + vv.i1.Pg:vv.iN.Pg;
%     iLCg(k,:) = Nx*(k-1) + vv.i1.LCg:vv.iN.LCg;
%     iPGs(k,:) = Nx*(k-1) + vv.i1.PGs:vv.iN.PGs;
% end
% for k = 1:NK
%     for i = 1:vv.N.Pg
%         dJdPg(iPg(k,i)) = baseMVA * polycost(mpc.gencost(i, :), Pg(k,i), 1);
%     end
%     for i = 1:vv.N.LCg
%         dJdLCg(iLCg(k,i)) = CDFg;
%     end
%     for i = 1:vv.N.PGs
%         dJdPGs(iPGs(k,i)) = mpc.Gcost(i);
%     end
% end
% dJ = dJ + dJdPg + dJdLCg + dJdPGs;
%     
end