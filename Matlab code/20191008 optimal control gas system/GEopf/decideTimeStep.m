function [dt, nt] = decideTimeStep(interruptionTime)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
controlMethod = 1;
if controlMethod == 1
    dt0 = 600; %default 30 min
    interruptionTime = interruptionTime * 3600;
    nt = floor(interruptionTime/dt0) + 1;
    dt = interruptionTime ./ nt;
end
% if controlMethod == 2
%     % if is not gas source failure (will not trigger transient)
%     if InfoComponents == 
%         dt = interruptionTime;
%         nt = 1;
%     else
%         timeUB = 1800; timeLB = 300; % uper and lower boundaries for time steps
%         % select the pressure change for last two times
%         % if the change decrease, the time step should be longer
%         failedGSindex = 

