function [Points] = Inner_appr(mpc) % equality and inequality constraints
    i = 1;
    Points = zeros(8,2);
    for c1 = -1:1
        for c2 = -1:1
            if c1== 0 && c2 ==0
                continue;
            else
                ffun = c1*z1+c2*z2;

                Points(i,:) = [z1,z2];
                i = i+1;
            end
        end
    end
        
end
runopf()
