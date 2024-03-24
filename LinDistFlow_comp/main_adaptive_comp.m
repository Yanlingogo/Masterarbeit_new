clc;
clear;
close all;

mpc = loadcase("case33_org.m");
%% vertexes
[vertexes, dispatch] = 
%% Jacobian at z_0
data_basepoint = jacobian_mpc(mpc);
%% compensation
for i = size()
