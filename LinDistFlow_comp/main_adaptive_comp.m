clc;
clear;

mpc = loadcase("case33_org.m");
%% edit the parameters of case

%% non-linear solver
% result of AC model
AC_points = sample_AC(mpc);

%% Adaptive compensation
% vertexes
[vertexes, dispatch] = vertex_ipopt(mpc);
% Jacobian at z_0
data_basepoint = jacobian_mpc(mpc);
% compensation
Resoultion  = 30;
line_points = cell(size(vertexes,1),1);
numlines = size(vertexes,1);
for i = 1: numlines
    power_dispatch = [dispatch{i} dispatch{mod(i,numlines)+1}];
    endpoints = [vertexes(i,:)' vertexes(mod(i,numlines)+1,:)'];
    line_points{i} = compensation_line_cor(data_basepoint, endpoints, power_dispatch, Resoultion);
end

%% visulization
% AC model
lightzs = [181, 180, 214]/255;
violett = [158, 49, 252]/255;
h1 = fill(AC_points(:,1), AC_points(:,2), lightzs);
set(h1, 'facealpha', 0.9, 'EdgeColor', violett,'LineWidth',2);
hold all
% LinDistFlow 
lightorange = [250, 188, 113]/255;
orange = [255, 99, 0]/255;
h2 = fill(vertexes(:,1), vertexes(:,2), lightorange, 'FaceAlpha',0.7);
set(h2, 'facealpha', 0.6, 'EdgeColor', orange,'LineWidth', 2);
% compensation result
visualization_lines(line_points)