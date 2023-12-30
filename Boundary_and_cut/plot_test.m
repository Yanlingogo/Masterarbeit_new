% 设置 x 轴和 y 轴的范围
x = linspace(-0.2, 0.7, 400);
y = linspace(-1.2, 1.2, 400);
[X, Y] = meshgrid(x, y);

% 初始化绘图
fig=figure; box on; hold all; set(fig, 'Position', [100, 100, 850, 650]);

% 设置坐标轴范围和网格线间隔
xlim([-0.2, 0.7]);
ylim([-1.2, 1.2]);
xticks(-0.2:0.1:0.7);
yticks(-1.2:0.5:1.2);
grid on;

% 绘制每个不等式定义的线
for i = 1:size(T,1)
    if T(i,1) == 0
        line(xlim, [v(i)/T(i,2) v(i)/T(i,2)], 'Color', 'r');
    elseif T(i,2) == 0
        line([v(i)/T(i,1) v(i)/T(i,1)], ylim, 'Color', 'r');
    else
        plot(x, (v(i) - T(i,1)*x)/T(i,2), 'r');
    end
end

% 检查每个点是否满足所有不等式
inside = all((A * [X(:), Y(:)]' <= b)');

% 绘制可行区域
fill(X(inside), Y(inside), 'g', 'FaceAlpha', 0.3);

% 设置坐标轴标签和标题
xlabel('x');
ylabel('y');
title('Linear Inequality Constraints');
hold off;