function [dx, nx] = decidePipelineCell(gtd)
%DECIDEPIPELINECELL Summary of this function goes here
%   Detailed explanation goes here
dx0 = 10000; % the aproximate dx = 2 km
nx = floor(gtd.Gline(:,4)/dx0) + 1;
dx = gtd.Gline(:,4) ./ nx;
end

