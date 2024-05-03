function normals = NFD_Norm_vector(points,vert_c)
    % add the first point to the end to form a closed area
    if points(1,:) ~= points(end,:)
        points = [points; points(1,:)]; 
    end
    
    % intialization
    normals = zeros(size(points,1)-1, 2);
    
    for i = 1:size(points,1)-1
        % direction of current line
        d = points(i+1, :) - points(i, :);
        
        % Calculate the outer normal vector (rotate 90 degrees counterclockwise)
        n = [d(2), -d(1)];
        
        % standardization
        n = n / abs((n*(points(i,:)-vert_c)'));
        
        normals(i, :) = n;
    end
end

