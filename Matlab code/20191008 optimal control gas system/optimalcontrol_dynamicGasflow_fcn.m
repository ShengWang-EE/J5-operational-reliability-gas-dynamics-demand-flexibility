function [f,df] = optimalcontrol_dynamicGasflow_fcn(p0,q0,mpc,gtd,initialCondition,para,dx,nx,dt,NK)
% nvar = vv.N.all * NK; 
[C_q,R,T,Z_T,Z,B,~,~,rhon] = unpackPara2(para);
PrsInitial = cell2mat(initialCondition.P')';
GfInitial =  cell2mat(initialCondition.Q')';
nGl = size(mpc.Gline,1);
Prs = mat2cell([PrsInitial';p0*1e5],NK+1,nx+1);
gasflow = mat2cell([GfInitial';q0],NK+1,nx+1);
% get index
[f_continuity, f_motion] = deal(cell(nGl,1));
[iContinuity,iMotion,iPrs,iGf] = cell(nGl,1);
for k = 1:NK
    for i = 1:nGl
        for m = 1:nx(i)
            iContinuity{i}(k,m) = sum(nx)*(k-1) + sum(nx(1:i-1)) + m;
            iMotion{i}(k,m) = sum(nx)*NK + sum(nx)*(k-1) + sum(nx(1:i-1)) + m;
            iPrs{i}(k,m) = vv.N.all*(k-1) + vv.i1.Prs-1 + sum(nx(1:i-1)) + m;
            iGf{i}(k,m) = vv.N.all*(k-1) + vv.i1.Gf-1 + sum(nx(1:i-1)) + m;
        end
    end
end

for i = 1:nGl
    p = Prs{i}; q = gasflow{i};
%     [f_continuity{i}, f_motion{i}] = deal(zeros(nx(i),NK));
%     [df_continuity{i}, df_motion{i}] = deal(zeros(nvar,nx(i)*NK));
    % parameter
    r = gtd.Gline(i,3); A = pi*r^2; D = 2*r;
    F = gtd.Gline(i,5);
    dl = dx(i); 
    for m = 2:nx(i)+1
        for k = 2:NK+1
            P = p(k-1,m-1); P_X = p(k-1,m); P_T = p(k,m-1); P_XT = p(k,m);
            Q = q(k-1,m-1); Q_X = q(k-1,m); Q_T = q(k,m-1); Q_XT = q(k,m);
            
            P_AV = 2/3 * (P + P_X - P*P_X/(P+P_X));
            P_TAV = 2/3 * (P_T + P_XT - P_T*P_XT/(P_T+P_XT));
            
            f_continuity{i}(m-1,k-1) = C_q*R*T/A * dt/dl * (Q_XT - Q_T) + P_TAV/Z_T - P_AV/Z;
            f_motion{i}(m-1,k-1) = (P_XT^2-P_T^2)/dl * (0.5) +  ...
                (Q_T+Q_XT)*abs(Q_T+Q_XT)*C_q^2*B^2/(2*F^2*D*A^2) + C_q^2*B^2/A^2 * (Q_XT^2-Q_T^2)/dl + ...
                C_q*P/A * (Q_T+Q_XT-Q-Q_X)/(2*dt);
%             %            ---------- calculate jacobi ----------
%             df_continuity_dP = -2/(3*Z) * (1 - P_X^2/(P+P_X)^2);
%             df_continuity_dP_X = -2/(3*Z) * (1 - P^2/(P+P_X)^2);
%             df_continuity_dP_T = 2/(3*Z_T) * (1 - P_XT^2/(P_T+P_XT)^2);
%             df_continuity_dP_XT = 2/(3*Z_T) * (1 - P_T^2/(P_T+P_XT)^2);
%             df_continuity_dQ = 0;
%             df_continuity_dQ_X = 0;
%             df_continuity_dQ_T = -(C_q*R*T)/A * dt/dl;
%             df_continuity_dQ_XT = (C_q*R*T)/A * dt/dl;
% 
%             df_motion_dP = C_q*(Q_T+Q_XT-Q-Q_X)/(2*A*dt);
%             df_motion_dP_X = 0;
%             df_motion_dP_T = -P_T/dl;
%             df_motion_dP_XT = P_XT/dl;
%             df_motion_dQ = - C_q*P/(2*A*dt);
%             df_motion_dQ_X = - C_q*P/(2*A*dt);
%             df_motion_dQ_T = sign(Q_T+Q_XT) * (Q_T+Q_XT)*C_q^2*B^2 / (F^2*D*A^2) + ...
%                 C_q^2*B^2/A^2 * (-2*Q_T/dl) + C_q*P/A * 1/(2*dt);
%             df_motion_dQ_XT = sign(Q_T+Q_XT) * (Q_T+Q_XT)*C_q^2*B^2 / (F^2*D*A^2) + ...
%                 C_q^2*B^2/A^2 * (2*Q_XT/dl) + C_q*P/A * 1/(2*dt);   
%             
%             df_continuity{i}(iPrs(k-1,m-1)) = df_continuity_dP; df_continuity{i}(iPrs(k-1,m)) = df_continuity_dP_X;
%             df_continuity{i}(iPrs(k,m-1)) = df_continuity_dP_T; df_continuity{i}(iPrs(k-1,m-1)) = df_continuity_dP_XT;
%             df_continuity{i}(iGf(k-1,m-1)) = df_continuity_dQ; df_continuity{i}(iGf(k-1,m)) = df_continuity_dQ_X;
%             df_continuity{i}(iGf(k,m-1)) = df_continuity_dQ_T; df_continuity{i}(iGf(k-1,m-1)) = df_continuity_dQ_XT;
%             df_motion{i}(iPrs(k-1,m-1)) = df_motion_dP; df_motion{i}(iPrs(k-1,m)) = df_motion_dP_X;
%             df_motion{i}(iPrs(k,m-1)) = df_motion_dP_T; df_motion{i}(iPrs(k-1,m-1)) = df_motion_dP_XT;
%             df_motion{i}(iGf(k-1,m-1)) = df_motion_dQ; df_motion{i}(iGf(k-1,m)) = df_motion_dQ_X;
%             df_motion{i}(iGf(k,m-1)) = df_motion_dQ_T; df_motion{i}(iGf(k-1,m-1)) = df_motion_dQ_XT;
        end
    end
end
f = [cell2mat(f_continuity); cell2mat(f_motion)]';
end
