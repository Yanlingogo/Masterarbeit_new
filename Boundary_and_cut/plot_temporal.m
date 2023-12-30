load("intersections.mat")
figure;

% 设置颜色地图，以便每个多边形有不同的颜色
colors = jet(length(intersections)); % 'jet' 是 MATLAB 的一个颜色地图

% 循环通过每个 cell 绘制多边形
for i = 1:length(intersections)
    % 获取当前 cell 的坐标点
    xy = (intersections{i})';
    centroid = mean(xy, 1);
    angles = atan2(xy(:,2) - centroid(2), xy(:,1) - centroid(1));
    [~, order] = sort(angles);
    sortedPoints = xy(order, :);
    % 假设 xy 是一个 n x 2 的矩阵，其中第一列是 x 坐标，第二列是 y 坐标
    x = sortedPoints(:, 1);
    y = sortedPoints(:, 2);
    
    % 在三维空间中添加多边形。Z 坐标由 i 控制，以将它们堆叠起来
    z = ones(size(x)) * i;
    patch(z, x, y, colors(i, :), 'EdgeColor', 'none');
    
    % 可以选择添加一些透明度
    alpha(0.5);

    hold on; % 保持当前图形，以便在上面绘制
    plot3(z, x, y, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', colors(i, :));

end

ylim([-0.2 0.7]);
zlim([-1.2 1.2]);

set(gca, 'YDir','reverse')
view(3);

xlabel('T(time period)');
ylabel('P(Mvar)');
zlabel('Q(Mvar)');
title('Feasible region with temporal coupling');

% 优化图形显示
grid on;
%axis tight;
rotate3d on; % 允许使用鼠标旋转视图