function [C_q,R,T,Z_T,Z,B,g,alpha,para,pn,rhon] = unpackPara3(para)
%UNPACKPARA2 Summary of this function goes here
%   Detailed explanation goes here
%% new
C_q = para.rhon*1000000/86400;
R = para.R;
T = para.T;
Z_T = para.Z; % not sure
Z = para.Z;
B = para.e;
g = 9.98;
alpha = 0;
pn = para.pn;
rhon = para.rhon;
end

