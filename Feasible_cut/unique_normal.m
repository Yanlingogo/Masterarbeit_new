function filteredMatrixB = unique_normal(matrixA, matrixB, tolerance)
    % 初始化一个逻辑向量，用于标记matrixB中需要保留的行
    toKeep = true(size(matrixB, 1), 1);
    
    for i = 1:size(matrixB, 1)
        % 对于matrixB中的每一行
        currentRowB = matrixB(i, :);
        for j = 1:size(matrixA, 1)
            % 计算与matrixA中每一行的差的绝对值
            diff = abs(currentRowB - matrixA(j, :));
            % 如果差的最大值在所有维度上都小于阈值，则认为行向量接近相等
            if all(diff <= tolerance, 2)
                toKeep(i) = false; % 标记为不保留
                break; % 已找到接近相等的行，无需进一步比较
            end
        end
    end
    
    % 返回去除接近重复行之后的matrixB
    filteredMatrixB = matrixB(toKeep, :);
end


