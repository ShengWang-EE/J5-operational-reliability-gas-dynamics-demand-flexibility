function [C_q,R,T,A,Z_T,Z,B,g,alpha,F,D,rhon,dx,nx] = unpackPara(para)
%UNPACKPARA Summary of this function goes here
%   Detailed explanation goes here
C_q = para.rhon*1000000/86400 ;
R = para.R;
T = para.T;
A = para.A;
Z_T = para.Z; % not sure
Z = para.Z;
B = para.e;
g = para.g;
alpha = para.alpha;
F = para.F;
D = para.D;
rhon = para.rhon;

dx = para.dx;
nx = para.nx;
end

