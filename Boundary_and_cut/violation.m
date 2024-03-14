fig=figure; box on; grid on; hold all; set(fig, 'Position', [100, 100, 650, 550])
load("violation.mat")
x = 1:1:size(violation,2);
plot(x, violation, '-o', ... % '-o'表示线和圆点标记
    'LineWidth', 2, ... % 设置线宽
    'Color', [68, 114, 196]/255, ... % 设置线颜色
    'MarkerSize', 8, ... % 设置标记大小
    'MarkerEdgeColor', 'none', ... % 设置标记边缘颜色
    'MarkerFaceColor', [68, 114, 196]/255); % 设置标记填充颜色
xlim([1 12]);
ylim([0 1]);
set(gca, 'FontSize', 14,'FontName', 'Times New Roman');

x_label = xlabel('Iteration number');
y_label = ylabel('Violation of the cut region');
set(x_label, 'Interpreter', 'latex', 'FontSize', 20, 'FontName', 'Times New Roman');
set(y_label, 'Interpreter', 'latex', 'FontSize', 20, 'FontName', 'Times New Roman');
exportgraphics(gca, 'violation.pdf', 'ContentType', 'vector');