clc
clear;


mpc = loadcase("case33_org.m");
%%  solver
% result of AC model
AC_points = sample_AC(mpc);
% result of LinDistFlow model
LDF_points = sample_LinDistFlow(mpc);
% result of Compensation 
[Comp_points,filtered_points] = compensation(LDF_points,mpc);
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

% shp = alphaShape(filtered_points(:,1),filtered_points(:,2));
% shp.Alpha = 0.1;
% % 绘制 Alpha Shape
% % plot(shp)
% 
% 
% % 获取边界点
% [bndPoints, bnd] = boundaryFacets(shp);
% bnd(end+1,:) = bnd(1,:);
% plot(bnd(:,1), bnd(:,2), 'k-', 'LineWidth', 2) % 绘制为黑色线条
% 
% % bndPoints 给出了边界上的点坐标
% disp('边界点坐标:')
% disp(bndPoints);
% 
% % bnd 给出了连接边界点的顺序，可以用来绘制多边形或进一步处理
% disp('边界连接顺序:')
% disp(bnd);