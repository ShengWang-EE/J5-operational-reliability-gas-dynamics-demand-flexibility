function initialCondition = setInitialConditionForGasSystem(result,gtd,nx,dx)
%the initial condition including P and Q of each pipeline cell
%   Detailed explanation goes here
%% initialization
% grab dimenssions
nGl = size(result.Gline,1);
initialCondition.P = cell(nGl,1); initialCondition.Q = cell(nGl,1);
[para] = initializeParameters2();
[~,~,~,~,~,B,~,~,rhon] = unpackPara2(para);
%%
for i = 1:nGl
    P = zeros(1,nx(i)+1);Q = zeros(1,nx(i)+1);
    Q(:) = result.Gline(i,6) * 1000000 / 86400; % convert to m3/s;
    fb = result.Gline(i,1); tb = result.Gline(i,2);
    P(1) = result.Gbus(fb,7) * 100000; %convert to Pa
    P(nx(i)+1) = result.Gbus(tb,7) * 100000;
    
    D = gtd.Gline(i,3);
    A = (D/2)^2 * pi;
    F = gtd.Gline(i,5);
    dp2dx = -4 * B^2 * rhon^2 * Q(1) * abs(Q(1)) /( F^2 * A^2 * D );
    P(2:nx(i)+1) = sqrt(P(1)^2 + dp2dx * (1:nx(i)) * dx(i));
    % 这里使用暂态的参数算的，和稳态的结果应该是一样的，这就说明暂态和稳态的参数统一起来了
    initialCondition.P{i} = P;% Pa
    initialCondition.Q{i} = Q*86400/1000000; %Mm3/day
end
% 4*rhon^2*B^2*4000*[6.02400000000000]^2/(D*A^2*864^2)

end

