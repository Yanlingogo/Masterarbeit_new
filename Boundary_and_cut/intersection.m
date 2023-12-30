function [valid_intersections] = intersection(A,b)
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
        if any(A * intersections(:, i) > b)
            valid(i) = false;
        end
    end
    valid_intersections = intersections(:, valid);
    [valid_intersections, ~, ~] = unique(valid_intersections', 'rows', 'stable');
    valid_intersections = valid_intersections';
end

