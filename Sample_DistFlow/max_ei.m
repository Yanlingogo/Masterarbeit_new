function [E] = max_ei(points)
% 使用convhull生成凸包
K = convhull(points);

% 绘制凸包
plot(points(:,1), points(:,2), 'o');
hold on;
plot(points(K,1), points(K,2), 'r-');

% 计算凸包
K = convhull(points(:,1), points(:,2));

% 构造线性不等式
A = [];
b = [];

for i = 1:length(K)-1
    point1 = points(K(i), :);
    point2 = points(K(i+1), :);
    
    % 使用两点形式的线的方程来计算 a 和 b
    normal = [point2(2) - point1(2), -(point2(1) - point1(1))];
    A = [A; normal];
    b = [b; normal * point1'];
end

disp('A = ')
disp(A)
disp('b = ')
disp(b)

yalmip('clear');

% 假设您已经有了 m x 1 的 A 和 b
m = length(A);

% 定义优化变量
B = sdpvar(2, 2);
d = sdpvar(2, 1);

% 构建约束
constraints = [];
constraints = [constraints, B >= eye(2) * 1e-4];
for i = 1:m
     constraints = [constraints, norm(B * A(i, :)', 2) + A(i, :) * d <= b(i)];
end

% 定义目标函数
objective = -logdet(B);

% 解决问题
options = sdpsettings('solver', 'sdpt3'); % 使用 SDPT3 求解器，您可以选择其他的
sol = optimize(constraints, objective, options);

% 提取结果
B_optimal = value(B);
d_optimal = value(d);
% 定义 u 为单位圆上的点
theta = linspace(0, 2*pi, 1000);
u = [cos(theta); sin(theta)];

% 应用仿射变换
E = B_optimal*u + d_optimal;

end
% 画图
% plot(E(1,:), E(2,:), 'b');
% axis equal;
% title('Transformed Ellipse: Bu + d');