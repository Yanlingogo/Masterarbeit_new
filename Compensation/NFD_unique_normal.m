% function filteredMatrixB = NFD_unique_normal(matrixA, matrixB, tolerance)
%     % 初始化一个逻辑向量，用于标记matrixB中需要保留的行
%     toKeep = true(size(matrixB, 1), 1);
% 
%     for i = 1:size(matrixB, 1)
%         % 对于matrixB中的每一行
%         currentRowB = matrixB(i, :);
%         for j = 1:size(matrixA, 1)
%             % 计算与matrixA中每一行的差的绝对值
%             diff = abs(currentRowB - matrixA(j, :));
%             % 如果差的最大值在所有维度上都小于阈值，则认为行向量接近相等
%             if all(diff <= tolerance, 2)
%                 toKeep(i) = false; % 标记为不保留
%                 break; % 已找到接近相等的行，无需进一步比较
%             end
%         end
%     end
% 
%     % 返回去除接近重复行之后的matrixB
%     filteredMatrixB = matrixB(toKeep, :);
% end

function filteredMatrixB = NFD_unique_normal(matrixA, matrixB, tolerance)
    % 初始化一个逻辑向量，用于标记matrixB中需要保留的行
    toKeep = true(size(matrixB, 1), 1);
    
    for i = 1:size(matrixB, 1)
        % 对于matrixB中的每一行
        currentRowB = matrixB(i, :);
        for j = 1:size(matrixA, 1)
            % 提取matrixA中的对应行
            currentRowA = matrixA(j, :);
            % 计算两行的比例因子，仅在非零分量上计算
            if currentRowA(1) ~= 0 && currentRowB(1) ~= 0
                k = currentRowB(1) / currentRowA(1);
            elseif currentRowA(2) ~= 0 && currentRowB(2) ~= 0
                k = currentRowB(2) / currentRowA(2);
            else
                continue; % 如果对应的分量均为0，则跳过
            end
            
            % 计算比例因子后的理论值
            scaledRowA = currentRowA * k;
            % 计算与理论值的差的绝对值
            diff = abs(currentRowB - scaledRowA);
            % 如果差的最大值在所有维度上都小于阈值，则认为两向量方向相同
            if all(diff <= tolerance, 2)
                toKeep(i) = false; % 标记为不保留
                break; % 已找到方向相同的向量，无需进一步比较
            end
        end
    end
    
    % 返回去除方向相同行之后的matrixB
    filteredMatrixB = matrixB(toKeep, :);
end
