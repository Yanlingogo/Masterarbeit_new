% 定义 x1 和 x2 的范围
x1 = linspace(-2, 2, 400);
x2 = linspace(-2, 2, 400);

% 生成网格
[X1, X2] = meshgrid(x1, x2);

% 计算函数值
F = X1.^2 + X1.*X2 + X2.^2;

% 绘制等高线，其中 F = 1
contour(X1, X2, F, [1 1])
title('Contour of x1^2 + x1*x2 + x2^2 = 1')
xlabel('x1')
ylabel('x2')
