load('Tv.mat');
T = Tv.T;
v = Tv.v;
[valid_intersections] = intersection(T,v);
x = linspace(-2, 2, 400);
y = linspace(-2.2, 2.2, 400);
[X, Y] = meshgrid(x, y);

% 初始化绘图
fig=figure; box on; grid on; hold all; set(fig, 'Position', [100, 100, 650, 550])


% 设置坐标轴范围和网格线间隔
xlim([-2, 2]);
ylim([-2.2, 2.2]);
xticks(-2:0.5:2);
yticks(-2:0.5:2);
set(gca, 'FontSize', 14,'FontName', 'Times New Roman');
grid on;

% 绘制每个不等式定义的线
lightorange = [250, 188, 113]/255;
orange = [252, 128, 2]/255;
for i = 1:size(T,1)
    if ismember(v(i), [1.5,1.6,2])
        continue;
    elseif T(i,1) == 0 
        line(xlim, [v(i)/T(i,2) v(i)/T(i,2)], 'Color', orange,'LineWidth', 1.5);
    elseif T(i,2) == 0 
        line([v(i)/T(i,1) v(i)/T(i,1)], ylim, 'Color', orange,'LineWidth', 1.5);
    else
        plot(x, (v(i) - T(i,1)*x)/T(i,2), 'Color',orange,'LineWidth', 1.5);
    end
end
deepgreen = [114, 156, 97]/255;
deepyellow = [226, 201, 93]/255;
yline(2, 'Color', deepgreen,'LineWidth',2.5)
yline(-2, 'Color', deepgreen,'LineWidth',2.5,'LineStyle','--')
xline(1.6, 'Color', deepyellow,'LineWidth',2.5)
xline(-1.5, 'Color', deepyellow,'LineWidth',2.5,'LineStyle','--')

fill(valid_intersections(:,1), valid_intersections(:,2), lightorange, 'FaceAlpha',0.7, 'EdgeColor', 'none');

%% results considering uncertainty
load("AF_uncertainty.mat")
Points_exact = sortedPoints;
lightbl = [66, 148, 249]/255;
h1 = fill(Points_exact(:,1), Points_exact(:,2), lightbl);
set(h1, 'facealpha', 0.6, 'EdgeColor', 'none');