function [sortedPoints] = Orderpoints(Points)
    [Points, ~, ~] = unique(Points, 'rows', 'stable');
    centroid = mean(Points, 1);
    angles = atan2(Points(:,2) - centroid(2), Points(:,1) - centroid(1));
    [~, order] = sort(angles);
    sortedPoints = Points(order, :);
end

