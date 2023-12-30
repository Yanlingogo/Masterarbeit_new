function [Q, z_s] = Boundary_search(b0, T, v, gamma, beta, M)
    Numb0 = numel(b0);
    Numv = numel(v);

    lambda = sdpvar(Numb0,1);
    mu = sdpvar(Numv,1);
    delta = binvar(Numv,1);
    [~, T_y] = size(T);
    z = sdpvar(T_y,1);
        
    Cons = [T' * mu + gamma' * lambda == 0];
    Cons = [Cons, beta' * lambda == 0];
    Cons = [Cons, lambda >= -1, lambda <=0];
    Cons = [Cons, v - T*z >= 0, v - T*z <= M*delta];
    Cons = [Cons, mu >= 0, mu <= M*(1-delta)];
    
    Obj = -(lambda' * b0 + mu' * v);

    options = sdpsettings('solver', 'gurobi','verbose', 0);
    %options = sdpsettings('verbose', 0, 'solver', 'gurobi', 'gurobi.DisplayInterval', 0);
    sol = optimize(Cons, Obj, options);

    Q = -value(Obj);
    z_s = value(z);

    value_lambda = value(lambda);
    value_mu = value(mu);
    value_delta = value(delta);
end

