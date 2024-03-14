% 定义二次规划问题的数据
Q = [2, 1;    % 二次目标函数系数矩阵
     1, 2];
c = [-1; -1]; % 线性目标函数系数

% 线性约束条件
A = [1, 1];  % 约束系数矩阵
b = [1];    % 约束边界

% 构造 MOSEK 格式的问题
prob.c = c;
prob.qosubi = [1; 1; 2];
prob.qosubj = [1; 2; 2];
prob.qoval = [Q(1,1); Q(1,2); Q(2,2)]; % 对于对称矩阵，只需提供上三角部分
prob.a = A;
prob.blc = -inf(size(b));
prob.buc = b;
prob.blx = [];  % 变量的下界（默认为空，表示没有下界）
prob.bux = [];  % 变量的上界（默认为空，表示没有上界）

% 使用 MOSEK 来求解
[r, res] = mosekopt('minimize', prob);

% 提取解决方案
if r == 0
    x = res.sol.itr.xx;
    fprintf('Optimal objective value: %f\n', c'*x + 0.5*x'*Q*x);
    disp('Optimal solution:');
    disp(x);
else
    fprintf('Solver reported an error. Error code: %d\n', r);
end
