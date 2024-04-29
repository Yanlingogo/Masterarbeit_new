function [power] = func_speed_power(speed,vi,vr,v0,Pr)
    [x,y] = size(speed);
    power = zeros(x,y);
    for i = 1:x
        for j = 1:y
            if speed(i,j) < vi || speed(i,j)> v0
                power(i,j) = 1e-6; % ensure the covariance are positive definite
            elseif speed(i,j) >= vi && speed(i,j) <= vr
                power(i,j) = Pr * ((speed(i,j)^2-vi^2)/(vr^2-vi^2));
            else
                power(i,j) = Pr;
            end
        end
    end
end

