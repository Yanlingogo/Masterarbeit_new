function [A_u, b_u] = E_UCI(A, b)
    Numb0 = numel(b);
    beta = binvar(Numb0, 1);
    [~, A_y] = size(A);
    x = sdpvar(A_y, 1);
    s = sdpvar(Numb0, 1);

    flag = 0;
    S = 1e3;
    idx = [];

    Cons = [A*x <= b];
    Cons = [Cons, A*x + s >= b];
    Cons = [Cons, s >= 0];
    Cons = [Cons, beta - s/S >= 0];

    while flag == 0
        Cons = [Cons, beta(idx) == 1];
        Obj = sum(beta);
        %options = sdpsettings('verbose', 0, 'solver', 'gurobi', 'gurobi.DisplayInterval', 0);
        options = sdpsettings('verbose', 0,'solver','gurobi');
        sol = optimize(Cons, Obj, options);

        beta_value = value(beta);
        idx = vertcat(idx, find(beta_value <= 1e-5));
        idx = sort(idx);

        if sum(beta_value) == Numb0
            flag = 1;
        else
            flag = 0;
        end
    end
  A_u = A(idx,:);
  b_u = b(idx,:);
  disp(['The original number of constraints: ', num2str(Numb0),', Number of compressed model constraints: ',num2str(numel(idx))]);
end

