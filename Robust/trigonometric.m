% 定义theta_min为-60度，转换为弧度
theta_min = -180 * (pi / 180);
theta_max = 180 * (pi / 180);
% 定义theta的范围，例如从-60度到0度，转换为弧度
theta = linspace(theta_min, theta_max, 1000);

% 计算特定函数的值
f_theta = theta + ((sin(theta_min) - theta_min) / theta_min^2) .* theta.^2;
f_theta2 = theta + ((sin(theta_max) - theta_max) / theta_max^2) .* theta.^2;

% 计算原始的sin函数值
sin_theta = sin(theta);

% 绘制特定函数的图像
plot(theta, f_theta, 'b', 'LineWidth', 2)
hold on % 保持当前图像，以便在上面添加新的图像
plot(theta, f_theta2, 'g', 'LineWidth', 2)
hold on 
% 在同一图像中绘制原始的sin函数
plot(theta, sin_theta, 'r--', 'LineWidth', 2)

% 添加图例
legend('Modified Function', 'sin(\theta)')

% 添加坐标轴标签和标题
xlabel('\theta (radians)')
ylabel('Function Value')
title('Comparison of Modified Function and sin(\theta)')

% 开启网格
grid on

% 取消保持状态
hold off
