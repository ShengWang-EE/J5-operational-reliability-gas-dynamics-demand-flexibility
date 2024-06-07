function [f] = pipeSectionForOneTimeStep2(allP,allQ,allP_T,allQ_T,diameter,F,nxi,dx,dt,para)
[C_q,R,T,Z_T,Z,B,g,alpha,~] = unpackPara2(para);
D = diameter;
A = (D/2)^2 * pi;
nx = nxi;
%% solve euqations using Matlab built in solver
f_continuity = zeros(nx,1);
f_motion = zeros(nx,1);
% jacobi = zeros(2*(nx+1),2*nx);
% dfcontinuitydP = zeros(nx+1,nx);
% dfcontinuitydQ = zeros(nx+1,nx);
% dfmotiondP = zeros(nx+1,nx);
% dfmotiondQ = zeros(nx+1,nx);
for i = 1:nx
    % variables known:(P,P_X,Q,Q_X) and two of (P_T,P_XT,Q_T,Q_XT)
    % variables unknow: other two of (P_T,P_XT,Q_T,Q_XT)
    % give values
    P = allP(i); Q = allQ(i); 
    P_X = allP(i+1); Q_X = allQ(i+1);
    P_T = allP_T(i); Q_T = allQ_T(i);
    P_XT = allP_T(i+1); Q_XT = allQ_T(i+1);
    % 
    P_AV = 2/3 * (P + P_X - P*P_X/(P+P_X));
    P_TAV = 2/3 * (P_T + P_XT - P_T*P_XT/(P_T+P_XT));
    % ----------- calculate f ----------
    f_continuity(i) = C_q*R*T/A * dt/dx * (Q_XT - Q_T) + P_TAV/Z_T - P_AV/Z;
%     f_motion(i) = (P_XT^2-P_T^2)/dx * (0.5 - C_q^2*B^2*(Q_T+Q_XT)^2/(8*A^2*P_TAV^2)) + P_TAV^2*g*sin(alpha)/B^2 + ...
%         (Q_T+Q_XT)*abs(Q_T+Q_XT)*C_q^2*B^2/(2*F^2*D*A^2) + C_q^2*B^2/A^2 * (Q_XT^2-Q_T^2)/dx + ...
%         C_q*P/A * (Q_T+Q_XT-Q-Q_X)/(2*dt);
    f_motion(i) = (P_XT^2-P_T^2)/dx * (0.5) + P_TAV^2*g*sin(alpha)/B^2 + ...
        (Q_T+Q_XT)*abs(Q_T+Q_XT)*C_q^2*B^2/(2*F^2*D*A^2) + C_q^2*B^2/A^2 * (Q_XT^2-Q_T^2)/dx + ...
        C_q*P/A * (Q_T+Q_XT-Q-Q_X)/(2*dt);
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
end
% test
f_motion = f_motion / 1e5;
f = [f_continuity; f_motion;];
% jacobi = [dfcontinuitydP, dfmotiondP;
%           dfcontinuitydQ, dfmotiondQ];
end
