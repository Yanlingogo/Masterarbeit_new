% visualizeZonotope.m
% 引入 CORA 工具箱
addpath(genpath('path_to_CORA')); % 替换 'path_to_CORA' 为 CORA 工具箱的实际路径

% 定义 zonotope
center = [1; 2];
generators = [1, 0; 0, 1];
Z = zonotope([center, generators]);

% 可视化 zonotope
figure;
plot(Z, [1, 2], 'Filled', true, 'FaceColor', [0.5, 0.5, 0.5], 'EdgeColor', 'b');
grid on;
xlabel('x');
ylabel('y');
title('Zonotope Visualization');
