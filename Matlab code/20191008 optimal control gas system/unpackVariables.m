function [Pg,Qg,Va,Vm,Prs,Gf,PGs,LCg] = unpackVariables(x,vv,NK)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
Nx = size(x,1) / NK;
Pg = zeros(NK,vv.N.Pg);Qg = zeros(NK,vv.N.Qg);Va = zeros(NK,vv.N.Va);Vm = zeros(NK,vv.N.Vm);
Prs = zeros(NK,vv.N.Prs);Gf = zeros(NK,vv.N.Gf);PGs = zeros(NK,vv.N.PGs);LCg = zeros(NK,vv.N.LCg);
for k = 1:NK
    xk = x((k-1)*Nx+1:k*Nx);
    
    Pg(k,:) = xk(vv.i1.Pg:vv.iN.Pg)';Qg(k,:) = xk(vv.i1.Qg:vv.iN.Qg)';
    Va(k,:) = xk(vv.i1.Va:vv.iN.Va)';Vm(k,:) = xk(vv.i1.Vm:vv.iN.Vm)';
    Prs(k,:) = xk(vv.i1.Prs:vv.iN.Prs)';Gf(k,:) = xk(vv.i1.Gf:vv.iN.Gf)';
    PGs(k,:) = xk(vv.i1.PGs:vv.iN.PGs)';LCg(k,:) = xk(vv.i1.LCg:vv.iN.LCg)';
end
end

