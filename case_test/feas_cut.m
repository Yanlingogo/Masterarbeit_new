function [D,d] = feas_cut(A, B, b, D, d, M)
    f = 1;
    violation = [];
    while ~(f<=1e-5)
        yalmip('clear')
        
        Numb = numel(b);
        Numd = numel(d);
        
        z = sdpvar(Numb,1);
        v = sdpvar(Numd,1);
    
        alpha = binvar(Numd,1);
        [~, D_y] = size(D);
        x = sdpvar(D_y,1);
        
        Cons = A'*z + D'*v == 0;
        Cons = [Cons, B' * z == 0];
        Cons = [Cons, z >= -1, z <=0];
        Cons = [Cons, d - D*x >= 0, d - D*x <= M*alpha];
        Cons = [Cons, v >= 0, v <= M*(1-alpha)];
        
        Obj = -(z' * b + d' * v);
        options = sdpsettings('solver', 'gurobi','verbose',0);
        optimize(Cons, Obj, options);

        f = -value(Obj);
        disp(['Exceeding the constraints: f = ', num2str(f)]);
        violation  = [violation,f];
        value_z = value(z);
        % value_v = value(v);
        % value_alpha = value(alpha);

        if f ~= 0
            if any(-value_z'*A ~= 0)
                D = [D;-value_z'*A];
                d = [d;-value_z'*b];
            else
                error('no feasible region')
            end
        end
    end
    % [D, d] = E_UCI(D, d);
 
end

