function f = optimalcontrol_dynamicGasflow_fcn2(p0,q0,mpc,gtd,initialCondition,para,dx,nx,dt,NK)

[C_q,R,T,Z_T,Z,B,~,~,rhon] = unpackPara2(para);
PrsInitial = cell2mat(initialCondition.P')';
GfInitial =  cell2mat(initialCondition.Q')';
nGl = size(mpc.Gline,1);
Prs = mat2cell([PrsInitial';p0*1e5],NK+1,nx+1);
gasflow = mat2cell([GfInitial';q0],NK+1,nx+1);

[f_continuity, f_motion] = deal(cell(nGl,1));

% for i = 1:nGl
%     p = Prs{i}; q = C_q * gasflow{i};
%     [f_continuity{i}, f_motion{i}] = deal(zeros(nx(i),NK));
%     % parameter
%     r = gtd.Gline(i,3); A = pi*r^2; D = 2*r;
%     F = gtd.Gline(i,5);
%     dl = dx(i); 
%     for m = 2:nx(i)+1
%         for k = 2:NK+1
%             f_continuity{i}(m-1,k-1) = 2*rhon*dt*B^2/(A) * (q(k,m)-q(k,m-1)) + ...
%                 dl * (p(k,m)-p(k-1,m)+p(k,m-1)-p(k-1,m-1));
%             f_motion{i}(m-1,k-1) = dt*(p(k,m)-p(k,m-1)) + ...
%                 dl*rhon/(2*A) * (q(k,m)-q(k-1,m)+q(k,m-1)-q(k-1,m-1)£© + ...
%                 (dt*dl*rhon^2*B^2)/(F^2*D*A^2) * ( q(k,m-1)*abs(q(k,m-1))/p(k,m-1) + q(k,m)*abs(q(k,m))/p(k,m)) ;
%         end
%     end
% end
for i = 1:nGl
    p = Prs{i}; q = gasflow{i};
    [f_continuity{i}, f_motion{i}] = deal(zeros(nx(i),NK));
    % parameter
    r = gtd.Gline(i,3); A = pi*r^2; D = 2*r;
    F = gtd.Gline(i,5);
    dl = dx(i); 
    for m = 2:nx(i)+1
        for k = 2:NK+1
            P = p(k-1,m-1);
            P_X = p(k-1,m);
            P_T = p(k,m-1);
            P_XT = p(k,m);
            Q = q(k-1,m-1);
            Q_X = q(k-1,m);
            Q_T = q(k,m-1);
            Q_XT = q(k,m);
            
            P_AV = 2/3 * (P + P_X - P*P_X/(P+P_X));
            P_TAV = 2/3 * (P_T + P_XT - P_T*P_XT/(P_T+P_XT));
            f_continuity{i}(m-1,k-1) = C_q*R*T/A * dt/dl * (Q_XT - Q_T) + P_TAV/Z_T - P_AV/Z;
            f_motion{i}(m-1,k-1) = (P_XT^2-P_T^2)/dl * (0.5) +  ...
                (Q_T+Q_XT)*abs(Q_T+Q_XT)*C_q^2*B^2/(2*F^2*D*A^2) + C_q^2*B^2/A^2 * (Q_XT^2-Q_T^2)/dl + ...
                C_q*P/A * (Q_T+Q_XT-Q-Q_X)/(2*dt);
        end
    end
end
f = [cell2mat(f_continuity); cell2mat(f_motion)]';
%%
    % ---------- calculate jacobi ----------
%     %evaluate the value
%     df_continuity_dP_T = 2/(3*Z_T) * (1 - P_XT^2/(P_T+P_XT)^2);
%     df_continuity_dP_XT = 2/(3*Z_T) * (1 - P_T^2/(P_T+P_XT)^2);
%     df_continuity_dQ_T = -(C_q*R*T)/A * dt/dx;
%     df_continuity_dQ_XT = (C_q*R*T)/A * dt/dx;
%     
%     df_motion_dP_T = -P_T/dx;
%     df_motion_dP_XT = P_XT/dx;
%     df_motion_dQ_T = sign(Q_T+Q_XT) * (Q_T+Q_XT)*C_q^2*B^2 / (F^2*D*A^2) + ...
%         C_q^2*B^2/A^2 * (-2*Q_T/dx) + C_q*P/A * 1/(2*dt);
%     df_motion_dQ_XT = sign(Q_T+Q_XT) * (Q_T+Q_XT)*C_q^2*B^2 / (F^2*D*A^2) + ...
%         C_q^2*B^2/A^2 * (2*Q_XT/dx) + C_q*P/A * 1/(2*dt);   
%     % substitute into
%     dfcontinuitydP(i,i) = df_continuity_dP_T;
%     dfcontinuitydP(i+1,i) = df_continuity_dP_XT;
%     dfcontinuitydQ(i,i) = df_continuity_dQ_T;
%     dfcontinuitydQ(i+1,i) = df_continuity_dQ_XT;
%     
%     dfmotiondP(i,i) = df_motion_dP_T;
%     dfmotiondP(i+1,i) = df_motion_dP_XT;
%     dfmotiondQ(i,i) = df_motion_dQ_T;
%     dfmotiondQ(i+1,i) = df_motion_dQ_XT;
% jacobi = [dfcontinuitydP, dfmotiondP;
%           dfcontinuitydQ, dfmotiondQ];
end
