function [A_extended, distances] = Cal_dist(A)
    
    % 在A的末尾添加A的第一行
    A_extended = [A; A(1,:)];
    
    % 计算每对相邻坐标之间的差分
    diffA = diff(A_extended, 1, 1);  % 对行进行操作
    
    % 计算欧几里得距离
    distances = sqrt(sum(diffA.^2, 2));
end

