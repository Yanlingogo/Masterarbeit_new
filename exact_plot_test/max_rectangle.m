function [minC, maxC, minR, maxR, maxArea] = max_rectangle(matrix)
    [m, n] = size(matrix);
    left = zeros(m, n);

    % 计算每个元素左侧连续1的数量
    for i = 1:m
        for j = 1:n
            if matrix(i, j) == 1
                if j == 1
                    left(i, j) = 1;
                else
                    left(i, j) = left(i, j - 1) + 1;
                end
            end
        end
    end

    % 初始化最大矩形的位置和面积
    minC = -1; maxC = -1; minR = -1; maxR = -1; maxArea = 0;

    % 处理每一列
    for j = 1:n
        up = -ones(m, 1);
        down = ones(m, 1) * (m + 1);

        % 计算 up
        stack = java.util.LinkedList();
        for i = 1:m
            while ~stack.isEmpty() && left(stack.peek() + 1, j) >= left(i, j)
                stack.pop();
            end
            if ~stack.isEmpty()
                up(i) = stack.peek() + 1;
            end
            stack.push(i - 1);
        end

        % 计算 down
        stack.clear();
        for i = m:-1:1
            while ~stack.isEmpty() && left(stack.peek() + 1, j) >= left(i, j)
                stack.pop();
            end
            if ~stack.isEmpty()
                down(i) = stack.peek() + 1;
            end
            stack.push(i - 1);
        end

        % 计算最大矩形
        for i = 1:m
            height = down(i) - up(i) - 1;
            area = height * left(i, j);
            if area > maxArea
                maxArea = area;
                minC = up(i) + 1;
                maxC = down(i) - 1;
                minR = j - left(i, j) + 1;
                maxR = j;
            end
        end
    end

    % 调整返回值的顺序以匹配 Groovy 代码
    [minR, maxR, minC, maxC, maxArea] = deal(minR, maxR, minC, maxC, maxArea);
end

