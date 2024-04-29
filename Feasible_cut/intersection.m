function [sortedPoints] = intersection(A,b)
    [A_x,~] = size(A);
    intersections = [];

    for i = 1:A_x-1
        for j = i+1:A_x
            A_sub = [A(i, :); A(j, :)];
            B = [b(i); b(j)]; 
            if rank(A_sub) == 2 % 检查A是否满秩，以确保有唯一解
                intersection = A_sub \ B;
                intersections = [intersections,intersection];
            else 
                disp(['no intersection:',num2str(i),',',num2str(j)])
            end
        end
    end

    % Check if the intersection point is within the constraints
    valid = true(1, size(intersections, 2));

    for i = 1:size(intersections, 2)
        if any(A * intersections(:, i) > b+1e-6)
            valid(i) = false;
        end
    end
    valid_intersections = intersections(:, valid);
    [valid_intersections, ~, ~] = unique(valid_intersections', 'rows', 'stable');

    % order the points by angles
    centroid = mean(valid_intersections, 1);
    angles = atan2(valid_intersections(:,2) - centroid(2), valid_intersections(:,1) - centroid(1));
    [~, order] = sort(angles);
    sortedPoints = valid_intersections(order, :);
    % 假设 xy 是一个 n x 2 的矩阵，其中第一列是 x 坐标，第二列是 y 坐标

end

