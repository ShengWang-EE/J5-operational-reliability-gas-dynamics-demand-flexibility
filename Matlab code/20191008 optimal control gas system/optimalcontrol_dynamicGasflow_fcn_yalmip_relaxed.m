function [f] = optimalcontrol_dynamicGasflow_fcn_yalmip_relaxed(p0,q0,mpc,gtd,initialCondition,para,vv,dx,nx,dt,NK)

[C_q,R,T,Z_T,Z,B,~,~,rhon] = unpackPara2(para);
PrsInitial = cell2mat(initialCondition.P')';
GfInitial =  cell2mat(initialCondition.Q')';
nGl = size(mpc.Gline,1);
Prs = mat2cell([PrsInitial';p0*1e6],NK+1,nx+1);
gasflow = mat2cell([GfInitial';q0],NK+1,nx+1);
% get index

[iContinuity,iMotion,iPrs,iGf] = deal(cell(nGl,1));
for k = 1:NK
    for i = 1:nGl
        for m = 1:nx(i)
            iContinuity{i}(k,m) = sum(nx)*(k-1) + sum(nx(1:i-1)) + m;
            iMotion{i}(k,m) = sum(nx)*NK + sum(nx)*(k-1) + sum(nx(1:i-1)) + m;
        end
    end
end
f_continuityAndMontion = []; counter = 1;
for i = 1:nGl
    p = Prs{i}; q = gasflow{i};
    % parameter
    r = gtd.Gline(i,3); A = pi*r^2; D = 2*r;
    F = gtd.Gline(i,5);
    dl = dx(i); 
    for m = 2:nx(i)+1
        for k = 2:NK+1
            P = p(k-1,m-1); P_X = p(k-1,m); P_T = p(k,m-1); P_XT = p(k,m);
            Q = q(k-1,m-1); Q_X = q(k-1,m); Q_T = q(k,m-1); Q_XT = q(k,m);
            
            P_AV = (1/2*P + 1/2*P_X);
            P_TAV = (1/2*P_T + 1/2*P_XT);
            
            % relax the motion equation: (P_XT^2-P_T^2)/dl * (0.5) + (Q_T+Q_XT)*signGf*(Q_T+Q_XT)*C_q^2*B^2/(2*F^2*D*A^2) 
            % + C_q^2*B^2/A^2 * (Q_XT^2-Q_T^2)/dl + C_q*P/A * (Q_T+Q_XT-Q-Q_X)/(2*dt);
            k1 = 1/dl * (0.5); k2 = *C_q^2*B^2/(2*F^2*D*A^2); k3 = C_q^2*B^2/A^2 /dl; k4 = C_q/A/(2*dt);
            signGf = sign(mean(initialCondition.Q{i}));
            f_continuity = C_q*R*T/A * dt/dl * (Q_XT - Q_T) + P_TAV/Z_T - P_AV/Z;
            if signGf >= 0
                f_motion1 = norm([sqrt(k1)*P_XT, sqrt(k2)*(Q_T+Q_XT), sqrt(k3)*Q_XT, sqrt(k4/2)*(P+Q_T), ...
                    sqrt(k4/2)*(P+Q_XT), sqrt(k4)*P, sqrt(k4/2)*Q, sqrt(k4/2)*Q_X],2) - z(counter);
                f_motion2 = norm([sqrt(k1)*P_T, sqrt(k3)*Q_T, sqrt(k4)*P, sqrt(k4/2)*Q_T, ...
                    sqrt(k4/2)*Q_XT, sqrt(k4/2)*(P+Q), sqrt(k4/2)*(P+Q_X)],2) - z(counter);
            else
                f_motion1 = norm([sqrt(k1)*P_XT, sqrt(k3)*Q_XT, sqrt(k4/2)*(P+Q_T), ...
                    sqrt(k4/2)*(P+Q_XT), sqrt(k4)*P, sqrt(k4/2)*Q, sqrt(k4/2)*Q_X],2) - z(counter);
                f_motion2 = norm([sqrt(k1)*P_T, sqrt(k2)*(Q_T+Q_XT), sqrt(k3)*Q_T, sqrt(k4)*P, sqrt(k4/2)*Q_T, ...
                    sqrt(k4/2)*Q_XT, sqrt(k4/2)*(P+Q), sqrt(k4/2)*(P+Q_X)],2) - z(counter);
            end
            f_continuityAndMontion = [f_continuityAndMontion;f_continuity;f_motion1; f_motion2];
            counter = counter + 1;
        end
    end
end
% f = [cell2mat(f_continuity); cell2mat(f_motion)]';
f = f_continuityAndMontion'/1e10;
end
