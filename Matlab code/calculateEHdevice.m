function [ei,gi,CHP,GB,EHP,AB] = calculateEHdevice(ei, gi, eee, ee3, e13, e1e, gg1, gg2, h1h, h14, h24, h2h, c3c, c4c)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
ei = ei;
gi = gi;

CHP.HO = h1h+h14;
CHP.EO = e13+e1e;
GB = h24+h2h;
EHP = c3c;
AB = c4c;
end

