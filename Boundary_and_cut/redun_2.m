function [A,b] = redun_2(A,b)
% given a constrained set X={x| Ax<=b}
% remove redundant constraints and return  X={x| Ax<=b}

% if remove constraint Cx<=d, 
% the optimal value of   max Cx  s.t. x in X does not increase
% then, this constraint is redundent
    for i = size(A,1):-1:1
        tempA = A;  tempA(i,:) = [];
        tempb = b;  tempb(i,:) = [];

        options = optimoptions('linprog', 'Display', 'none');
        [~,temp]  = linprog( -A(i,:), A, b, [], [], [], [],options);
        [~,temp1] = linprog( -A(i,:), tempA, tempb,[], [], [], [],options);  % if it is redundant, temp1 does not increase

        if -temp1 <= -temp
            A = tempA;
            b = tempb;
        end
    end
end
