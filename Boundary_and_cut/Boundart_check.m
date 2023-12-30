function [Q, h_s] = Boundart_check(b0, beta, gamma, z_b)
    Numb0 = numel(b0);

    h = sdpvar(Numb0, 1);
    
    Cons = [beta' * h == 0];
    Cons = [Cons, h>= -1, h <= 0];

    Obj = -(h' * (b0 - gamma * z_b));
    
    %options = sdpsettings('solver', 'gurobi','gurobi.FeasibilityTol', 1e-8);
    options = sdpsettings('verbose', 0, 'solver', 'gurobi');
    sol = optimize(Cons, Obj, options);

    h_s = value(h);
    Q = value(-Obj);
end

