time_points = 0:5:20; % 每隔5小时的时刻

% 对于每个时间点，我们有一个功率（P-Q）的可行域，这里我们用随机数据来模拟
for k = 1:length(time_points)
    % 生成P-Q的可行域，这里使用随机数据来模拟
    P = linspace(-50, 50, 10); % 有功功率的范围
    Q = rand(10, 1) * 20; % 无功功率的范围，随机数
    
    % 创建一个网格，这样每个P值都会对应一个Q值数组
    [P, Q] = meshgrid(P, Q);
    
    % 时刻对应的Z轴坐标值
    T = time_points(k) * ones(size(P));
    
    % 绘制三维表面图，每个表面图对应一个时刻
    surf(T, P, Q)
    
    % 保持图形，以便在同一图中绘制多个图层
    hold on
end