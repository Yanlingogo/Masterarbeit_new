function [uniqueCoords] = unique_coor(coords, tolerance)
    % 初始化一个逻辑向量，用于标记需要保留的行
    toKeep = true(size(coords, 1), 1);
    
    for i = 1:size(coords, 1)
        if toKeep(i)
            % 对于每个标记为保留的行，计算与其余行的差的绝对值
            diff = abs(coords - coords(i, :));
            % 计算差的最大值是否在所有维度上都小于阈值
            isClose = all(diff <= tolerance, 2);
            % 标记当前行之后的“接近”行为不保留（因为它们视为重复）
            toKeep(i+1:end) = toKeep(i+1:end) & ~isClose(i+1:end);
        end
    end
    
    % 返回去除接近重复行之后的坐标数组
    uniqueCoords = coords(toKeep, :);
end
