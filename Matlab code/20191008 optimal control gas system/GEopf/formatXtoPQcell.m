function [stuctureName] = formatXtoPQcell(x,nGl,nx)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
stuctureName.P = cell(nGl,1); stuctureName.Q = cell(nGl,1);
nP = sum(nx+1); nQ = sum(nx+1);
PforAllPipeline = x(1:nP); QforAllPipeline = x(nP+1:nP+nQ); % first sort out P and Q
cumnx = cumsum(nx+1);
for i = 1:nGl
    stuctureName.P{i} = zeros(1,nx(i)+1);stuctureName.Q{i} = zeros(1,nx(i)+1);
    if i==1
        stuctureName.P{i}(1:nx(i)+1) = PforAllPipeline(1:nx(i)+1);stuctureName.Q{i}(1:nx(i)+1) = QforAllPipeline(1:nx(i)+1);
    else
        stuctureName.P{i}(1:nx(i)+1) = PforAllPipeline(cumnx(i-1)+1:cumnx(i)); stuctureName.Q{i}(1:nx(i)+1) = QforAllPipeline(cumnx(i-1)+1:cumnx(i));
    end    
end
end

