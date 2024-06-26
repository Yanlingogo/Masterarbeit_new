% 生成随机的三维点
num_points = 30; % 你可以改变点的数量
points = rand(num_points, 3) + [zeros(num_points, 1), zeros(num_points, 1), 0.3*ones(num_points,1)]; % 确保z值至少为2，保持在x-y平面上方

% 计算点集的凸包
[K, v] = convhull(points);

% 获取所有顶点的x和y坐标，计算其在x-y平面的投影的凸包
xyPoints = points(:, 1:2);
[K2, v2] = convhull(xyPoints);

% 创建一个新的图形窗口
fig1=figure; box on; grid on; hold all; set(fig1, 'Position', [100, 100, 650, 550])


% 定义颜色
deepBlue = [0, 0, 0.5]; % 深蓝色轮廓线
lightBlue = [0.7, 0.7, 1]; % 浅蓝色填充
lightRed = [1, 0.5, 0.5]; % 浅红色轮廓线
veryLightRed = [1, 0.7, 0.7];% 更浅的红色填充
lightGreen = [117, 170, 107]/255; % 浅绿色填充
Green = [0.1, 0.4, 0.1];

% 第一个子图：三维凸包及其投影（包括虚线）

trisurf(K, points(:,1), points(:,2), points(:,3), 'FaceColor', lightBlue, 'EdgeColor', deepBlue);
axis equal;
hold on;
fill(xyPoints(K2,1), xyPoints(K2,2), lightGreen, 'EdgeColor', Green);
for i = 1:length(K2)
    idx = K2(i); % 凸包顶点的索引
    plot3([points(idx, 1), points(idx, 1)], [points(idx, 2), points(idx, 2)], [points(idx, 3), 0], 'k--','LineWidth', 1.5);
end
xlabel('X-axis', 'FontName', 'Times New Roman','FontSize', 14);
ylabel('Y-axis', 'FontName', 'Times New Roman','FontSize', 14);
zlabel('Z-axis', 'FontName', 'Times New Roman','FontSize', 14);
% title('3D Convex Hull with Projection', 'FontName', 'Times New Roman','FontSize', 14);
grid on;
view(140, 30);

% 设置坐标轴的刻度间隔为0.2，但不超过坐标轴的最大值
set(gca, 'XTick', 0:0.2:1, 'YTick', 0:0.2:1, 'ZTick', 0:0.2:2);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);

% 调整坐标轴范围以设置与凸包的距离
axisPadding = 0.1; % 设置坐标轴的padding值
xlim([min(points(:,1))-axisPadding, max(points(:,1))+axisPadding]);
ylim([min(points(:,2))-axisPadding, max(points(:,2))+axisPadding]);
zlim([0, max(points(:,3))+axisPadding]);

% 第二个子图：只有三维凸包，没有虚线
fig2=figure; box on; grid on; hold all; set(fig2, 'Position', [100, 100, 650, 550])
trisurf(K, points(:,1), points(:,2), points(:,3), 'FaceColor', lightGreen, 'EdgeColor', Green);
axis equal;
xlabel('X-axis', 'FontName', 'Times New Roman','FontSize', 14);
ylabel('Y-axis', 'FontName', 'Times New Roman','FontSize', 14);
zlabel('Z-axis', 'FontName', 'Times New Roman','FontSize', 14);
% title('Top view of a 3D convex hull', 'FontName', 'Times New Roman','FontSize', 14);
grid on;
view(180, 90); % 设置为俯视图

% 设置坐标轴的刻度间隔为0.2，但不超过坐标轴的最大值
set(gca, 'XTick', 0:0.2:1, 'YTick', 0:0.2:1, 'ZTick', 0:0.2:1);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);

% 同样调整坐标轴范围
xlim([min(points(:,1))-axisPadding, max(points(:,1))+axisPadding]);
ylim([min(points(:,2))-axisPadding, max(points(:,2))+axisPadding]);
zlim([0, max(points(:,3))+axisPadding]);

% 调整视角以更好地观察
% exportgraphics(fig, 'output.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none');