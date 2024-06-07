function [f] = pipeSectionForOneTimeStep(x,boundaryCondition,dt,para)
    [C_q,R,T,A,Z_T,Z,B,g,alpha,F,D,~,dx,nx] = unpackPara(para);
    %% unpack x considering boundaries
    switch boundaryCondition.flag
        case 1 % set the Q of beginning of the pipe, and the P of end of the pipe                
             P_Tall = [x(1:nx-1), boundaryCondition.P_T]; 
%              P_Tall = [boundaryCondition.P_T, x(1:nx-1)]; 
             Q_Tall = [boundaryCondition.Q_T, x((nx):(2*(nx-1)))];
             % nx = number of pipe cell + 1
    end
    %% solve euqations using Matlab built in solver
    f_continuity = zeros(nx-1,1);
    f_motion = zeros(nx-1,1);
    for i = 1:nx-1
        % variables known:(P,P_X,Q,Q_X) and two of (P_T,P_XT,Q_T,Q_XT)
        % variables unknow: other two of (P_T,P_XT,Q_T,Q_XT)
        P = boundaryCondition.P(i); Q = boundaryCondition.Q(i); 
        P_X = boundaryCondition.P(i+1); Q_X = boundaryCondition.Q(i+1);
        P_T = P_Tall(i); Q_T = Q_Tall(i);
        P_XT = P_Tall(i+1); Q_XT = Q_Tall(i+1);
        % 
        P_AV = 2/3 * (P + P_X - P*P_X/(P+P_X));
        P_TAV = 2/3 * (P_T + P_XT - P_T*P_XT/(P_T+P_XT));
        %
        f_continuity(i) = C_q*R*T/A * dt/dx * (Q_XT - Q_T) + P_TAV/Z_T - P_AV/Z;
        f_motion(i) = (P_XT^2-P_T^2)/dx * (0.5 - C_q^2*B^2*(Q_T+Q_XT)^2/(8*A^2*P_TAV^2)) + P_TAV^2*g*sin(alpha)/B^2 + ...
            (Q_T+Q_XT)*abs(Q_T+Q_XT)*C_q^2*B^2/(2*F^2*D*A^2) + C_q^2*B^2/A^2 * (Q_XT^2-Q_T^2)/dx + ...
            C_q*P/A * (Q_T+Q_XT-Q-Q_X)/(2*dt);
%         f_motion(i) = (P_XT^2-P_T^2)/dx * (0.5) + P_TAV^2*g*sin(alpha)/B^2 + ...
%             (Q_T+Q_XT)*abs(Q_T+Q_XT)*C_q^2*B^2/(2*F^2*D*A^2) + C_q^2*B^2/A^2 * (Q_XT^2-Q_T^2)/dx + ...
%             C_q*P/A * (Q_T+Q_XT-Q-Q_X)/(2*dt);
    end

    f = [f_continuity; f_motion;];
end
