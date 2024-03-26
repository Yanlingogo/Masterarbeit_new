function [points] = vector_segment(A,resolution)
    [row, ~] = size(A);
    points = zeros(row, resolution);
    for i = 1:row
        points(i,:) = linspace(A(i,1),A(i,2),resolution);
    end
end

