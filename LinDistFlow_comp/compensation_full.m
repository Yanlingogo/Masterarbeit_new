clc
clear;


mpc = loadcase("case33_org.m");
%%  solver
% result of AC model
AC_points = sample_AC(mpc);
% result of LinDistFlow model
LDF_points = sample_LinDistFlow(mpc);
% result of Compensation 
[Comp_points, filtered_points] = compensation_correction(LDF_points,mpc);
%% visulization 
% LinDistFlow
figure;
lightorange = [250, 188, 113]/255;
orange = [255, 99, 0]/255;
h1 = fill(LDF_points(:,1), LDF_points(:,2), lightorange, 'FaceAlpha',0.7);
set(h1, 'facealpha', 0.6, 'EdgeColor', orange,'LineWidth', 2);
hold all
% AC model
lightzs = [181, 180, 214]/255;
violett = [158, 49, 252]/255;
h2 = fill(AC_points(:,1), AC_points(:,2), lightzs);
set(h2, 'facealpha', 0.9, 'EdgeColor', violett,'LineWidth',2);
% Compensation
scatter(Comp_points(:,1),Comp_points(:,2),5,'filled');
scatter(filtered_points(:,1),filtered_points(:,2),5,'filled');