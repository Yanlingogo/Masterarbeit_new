function stochastic_test(random_data,robust_data)
% 示例：生成并绘制100组随机二维坐标的凸包
numGroups = 100; % 组数
N = 30; % 每组数据的点数

figure;
hold on; % 保持当前图形，以便叠加新图形
grid on;

% 使用循环生成和绘制每组数据的凸包
for i = 1:numGroups
    % 生成每组数据的随机二维坐标
    data = rand(N, 2) * 10; % 生成位于0到10之间的随机坐标
    
    % 计算凸包
    K = convhull(data(:,1), data(:,2));
    
    % 使用fill函数填充凸包区域，设置颜色和透明度
    fill(data(K,1), data(K,2), [1 0.6 0], 'FaceAlpha', 0.03,'EdgeColor', 'none'); % 浅橙色，透明度设置为0.1
end

% 设置图形属性
axis equal; % 保持横纵轴比例一致
xlim([0, 10]); % 设置X轴的范围
ylim([0, 10]); % 设置Y轴的范围
xlabel('X Axis');
ylabel('Y Axis');
title('Overlay of 100 Convex Hulls with Color Intensity Based on Overlap');

hold off; % 解除hold on
