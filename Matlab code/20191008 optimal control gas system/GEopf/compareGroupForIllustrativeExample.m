%% compare group
clear
t = 0:1:168;
lamda = 1/960; mu = 1/40;
nGs = 3;
gasLoad = 5; electricityLoad = 200; GFUefficiency = 200; 
pSuccess = mu/(lamda+mu) + lamda/(lamda+mu) * exp(-(lamda+mu)*t);
pFailure = lamda/(lamda+mu) * (1 - exp(-(lamda+mu)*t));
p(:,1) = pSuccess .^ 3;
p(:,2) = pSuccess.^2 .* pFailure * 3;
p(:,3) = pSuccess .* pFailure.^2 .* 3;
p(:,4) = pFailure.^3;
expectedGasProduction = p * [6,4,2,0]';
gasLOLP_steady = 1-p(:,1);
gasEENS_steady = (5-4) * p(:,2) + (5-2) * p(:,3) + (5-0) * p(:,4);
electricityLOLP_steady = 1-p(:,1);
electricityEENS_steady = 200 * (p(:,2)+p(:,3)+p(:,4));

% plot(t,y)