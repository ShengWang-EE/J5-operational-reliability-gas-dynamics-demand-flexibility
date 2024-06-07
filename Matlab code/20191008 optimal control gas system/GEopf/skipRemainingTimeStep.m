function [YES] = skipRemainingTimeStep(transientResult,nGl,timeStepCounter)
alpha = 0.0001;
difference0 = ones(nGl,2);
differenceP = ones(nGl,1);differenceQ = ones(nGl,1);
for m = 1:nGl
    P = transientResult.P{m}(timeStepCounter-2:timeStepCounter,:);
    Q = transientResult.Q{m}(timeStepCounter-2:timeStepCounter,:);
    differenceP(m) = max((P(1,:) - P(2,:)) ./ P(1,:));
    differenceQ(m) = max((Q(1,:) - Q(2,:)) ./ Q(1,:));
    difference0(m,1) = max(differenceP(m),differenceQ(m));
    differenceP(m) = max((P(2,:) - P(3,:)) ./ P(2,:));
    differenceQ(m) = max((Q(2,:) - Q(3,:)) ./ Q(2,:));
    difference0(m,2) = max(differenceP(m),differenceQ(m));
end
difference = max(difference0);
% if difference(2)/difference(1) < alpha
if difference(2) < alpha
    YES = 1;
else
    YES = 0;
end
end