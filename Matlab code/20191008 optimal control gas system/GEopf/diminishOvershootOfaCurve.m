function newCurve = diminishOvershootOfaCurve(curve)
% NK = max(size(curve));
% threshold = 0.01;
% nStall = 10;
% newCurve = curve;
% steadyValue = curve(1);
% k=0;
% while k <= NK-nStall
%     k = k + 1;
%     if abs((curve(k) - steadyValue) / steadyValue) > threshold
%         stallValue = curve(k:k+nStall-1);
%         steadyValue2 = stallValue(end);
% 
%         stallValue(stallValue>max(steadyValue,steadyValue2)) = max(steadyValue,steadyValue2);
%         stallValue(stallValue<min(steadyValue,steadyValue2)) = min(steadyValue,steadyValue2);
%         newCurve(k:k+nStall-1) = stallValue;
%         
%         steadyValue = steadyValue2;
%         k = k+nStall-1;
%     end         
% end
%% method 2
NK = max(size(curve));
amplitudeThreshold = 0.01;
durationThreshold = NK/30;
pointer = 1;
i = 1;
while pointer < NK
    allValueInfo(i,3) = curve(pointer);
    allValueInfo(i,1) = pointer;
    allValueInfo(i,2) = max(find(abs(curve(pointer:NK)-allValueInfo(i,3)) / (allValueInfo(i,3)+0.000000001) < amplitudeThreshold))+pointer-1;

    pointer = allValueInfo(i,2)+1;
    i = i+1;
end
allValueInfo(:,4) = allValueInfo(:,2) - allValueInfo(:,1);

transientIndex = find(allValueInfo(:,4) < durationThreshold);
steadyIndex = find(allValueInfo(:,4) >= durationThreshold);
for j = 1:size(transientIndex,1)
    i = transientIndex(j);
    steadyValue1Index = steadyIndex(max(find(steadyIndex<i)));
    steadyValue2Index = steadyIndex(min(find(steadyIndex>i)));
    if isempty(steadyValue1Index)
        steadyValue1Index = 1;
    end
    if isempty(steadyValue2Index)
        steadyValue2Index = size(allValueInfo,1);
    end 
    steadyValue1 = allValueInfo(steadyValue1Index,3);
    steadyValue2 = allValueInfo(steadyValue2Index,3);
    if allValueInfo(i,3) > max(steadyValue1,steadyValue2)
         allValueInfo(i,3) = max(steadyValue1,steadyValue2);
    end
    if allValueInfo(i,3) < min(steadyValue1,steadyValue2)
         allValueInfo(i,3) = min(steadyValue1,steadyValue2);
    end
end
%
newCurve = curve;
for i = 1:size(allValueInfo,1)
    newCurve(allValueInfo(i,1):allValueInfo(i,2)) = allValueInfo(i,3);
end
    
end
