function new_transientResult = diminishOvershoot(transientResult,nGl)
new_transientResult = transientResult;
for m = 1:nGl
    for j = 1:size(transientResult.P{m},2)
        curve = transientResult.P{m}(:,j);
        newCurve = diminishOvershootOfaCurve(curve);
        new_transientResult.P{m}(:,j) = newCurve;
        %
        curve = transientResult.Q{m}(:,j);
        newCurve = diminishOvershootOfaCurve(curve);
        new_transientResult.Q{m}(:,j) = newCurve;
    end
end
end
%%
