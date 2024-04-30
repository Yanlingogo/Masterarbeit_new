function [lost_flexibility,risk_area,Area_data,Area_ref] = compensation_error(data,data_ref)
    
    % create polyshape 
    poly = polyshape(data(:,1), data(:,2));
    poly_ref = polyshape(data_ref(:,1), data_ref(:,2)); % 基准多边形
    
    % calculate the union of polys
    intersectPoly = intersect(poly, poly_ref);
    
    % poly - union 
    mismatch_1 = subtract(poly, intersectPoly);
    mismatch_2 = subtract(poly_ref, intersectPoly);
    
    % calculate the area and error factor
    Area_data = area(poly);
    Area_ref = area(poly_ref);
    lost_flexibility = area(mismatch_2)/Area_ref;
    risk_area = area(mismatch_1)/Area_ref;

   
    fprintf('Lost flexibility: %.2f%%, risk area: %.2f%%\n', lost_flexibility*100, risk_area * 100);

end

