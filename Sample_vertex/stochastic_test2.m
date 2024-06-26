load("stochastic_25.mat");
numDatasets = size(stochastic_AF,1);
% 初始化一个网格，用于统计重叠次数
xRange = linspace(-1.8, -0.4, 1000); % X范围和精度
yRange = linspace(1, 2.8, 1000); % Y范围和精度
[X, Y] = meshgrid(xRange, yRange);
overlapCount = zeros(size(X)); % 重叠计数

% 生成并绘制每个凸包，更新重叠次数
for i = 1:200
    % 生成随机凸包
    data = stochastic_AF{i};
    K = convhull(data(:,1), data(:,2));
    
    % 使用inpolygon判断网格中的点是否在凸包内
    in = inpolygon(X, Y, data(K,1), data(K,2));
    
    % 更新重叠计数
    overlapCount = overlapCount + in;
end

% 绘制重叠计数图
% figure;
% imagesc(xRange, yRange, overlapCount);
% axis xy; % 确保Y轴方向正确
% colormap jet; % 使用jet颜色图
% colorbar; % 显示颜色条
% clim([1 numDatasets]); % 调整颜色轴的范围
% title('Overlap Count of Convex Hulls');

% 初始化绘图
fig=figure; box on; grid on; hold all; set(fig, 'Position', [100, 100, 650, 550])
imagesc(xRange, yRange, overlapCount/max(max(overlapCount)) * 100);
axis xy; % 确保Y轴方向正确

% Apply custom color mapping
blueToYellow = [linspace(0, 1, 256)' linspace(0, 1, 256)' linspace(1, 0, 256)'];
colormap("parula"); 
hColorbar = colorbar; % Show color bar
clim([0 100]); % Set the color axis range to 0-100
% Custom color bar scale labels, add percent signs
hColorbar.Ticks = [0 25 50 75 100]; 
hColorbar.TickLabels = {'0%', '25%', '50%', '75%', '100%'}; 


% Special marking of points with more than 95 overlaps
hold on; 
line_red = [242, 0, 26]/255;
[~, h] = contour(X, Y, overlapCount, [200, 200], 'LineColor', 'k', 'LineWidth', 2); 


load("vert22.mat")
Points_exact = vert;
lightbl = [66, 148, 249]/255;
h2 = fill(Points_exact(:,1), Points_exact(:,2), lightbl);
set(h2, 'facealpha', 0.6, 'EdgeColor', 'none');

load("vert23.mat")
Points_exact = vert;
lightbl = [79, 108, 221]/255;
h3 = fill(Points_exact(:,1), Points_exact(:,2), lightbl);
set(h3, 'facealpha', 0.6, 'EdgeColor', 'none');

xlim([-2, 2.5]);
ylim([-3, 3]);
set(gca, 'FontSize', 14,'FontName', 'Times New Roman');

x_label = xlabel('$p^{\mathrm{pcc}}$/(p.u.)'); % 修正了大括号
set(x_label, 'Interpreter', 'latex', 'FontSize', 20, 'FontName', 'Times New Roman');

y_label = ylabel('$q^{\mathrm{pcc}}$/(p.u.)'); % 修正了大括号
set(y_label, 'Interpreter', 'latex', 'FontSize', 20, 'FontName', 'Times New Roman');
