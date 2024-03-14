clear;
clc;
% 定义二次规划问题的数据
Q = [2, 1;    % 二次目标函数系数矩阵
     1, 2];
c = [-1; -1]; % 线性目标函数系数

% 线性约束条件
A = [1, 1];  % 约束系数矩阵
b = [1];    % 约束边界

% 二次约束条件
Qc = [1, 0;   % 二次约束系数矩阵
      0, 1];
d = [1; 1];  % 线性部分系数
e = 2;       % 二次约束边界

% 构造 Gurobi 模型
model.Q = sparse(0.5 * Q);
model.obj = c;
model.A = sparse(A);
model.rhs = b;
model.sense = '<';

% 添加二次约束
model.quadcon(1).Qc = sparse(Qc);
model.quadcon(1).q = d;
model.quadcon(1).rhs = e;

% 使用 Gurobi 来求解
result = gurobi(model);

% 显示结果
if strcmp(result.status, 'OPTIMAL')
    fprintf('Optimal objective value: %f\n', result.objval);
    disp('Optimal solution:');
    disp(result.x);
else
    fprintf('Solver reported status: %s\n', result.status);
end
