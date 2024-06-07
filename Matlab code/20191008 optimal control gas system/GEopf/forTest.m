%% for test
for i = 1:nGl
	a(i) = transientResult.Q{i}(end,1);
end
%%
gasReliability{1} = load('transient4.mat',gasReliability);
electricityReliability{1} = load('transient4.mat',electricityReliability);
gasReliability{2} = load('transient5.mat',gasReliability);
electricityReliability{2} = load('transient5.mat',electricityReliability);
gasReliability{3} = load('transient6.mat',gasReliability);
electricityReliability{3} = load('transient6.mat',electricityReliability);

sumGasReliability.system = zeros(size(gasReliability.system{1}));
sumGasReliability.nodal = zeros(size(gasReliability.nodal{1}));

for i = 1:max(size(gasReliability))
    sumGasReliability.system = sumGasReliability.system + gasReliability{i}.system;
    sumElectricityReliability.system = sumElectricityReliability.system + gasReliability{i}.system;
end
averageGasReliability.system = sumGasReliability.system / max(size(gasReliability));
averageElectricityReliability.system = sumElectricityReliability.system / max(size(gasReliability));