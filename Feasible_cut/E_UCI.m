function [A_u, b_u] = E_UCI(A, b)
    yalmip('clear')

    Numb0 = numel(b);
    beta = binvar(Numb0, 1);
    [~, A_y] = size(A);
    x = sdpvar(A_y, 1);
    s = sdpvar(Numb0, 1);

    Cons0 = A*x<= b;
    Obj0 = sum(x);
    options = sdpsettings('verbose', 0,'solver','gurobi');
    sol0 = optimize(Cons0, Obj0, options);

    if sol0.problem == 0
        disp('Optimization Succeeded!');

    else
        disp(['Solver reported: ', yalmiperror(sol0.problem)]);
        error('Optimization Failed: The problem has no solution or other errors occurred.');
    end


    flag = 0;
    S = 1e2;
    idx = [];

    Cons = A*x <= b;
    Cons = [Cons, A*x + s == b];
    Cons = [Cons, s >= 0];
    Cons = [Cons, beta - s/S >= 0];

    while flag == 0
        Cons = [Cons, beta(idx) == 1];
        Obj = sum(beta);
        %options = sdpsettings('verbose', 0, 'solver', 'gurobi', 'gurobi.DisplayInterval', 0);
        sol = optimize(Cons, Obj, options);

        beta_value = value(beta);
        idx = vertcat(idx, find(beta_value <= 1e-6));
        idx = sort(idx);

        if sum(beta_value) >= (Numb0 - 1e-6)
            flag = 1;
        else
            flag = 0;
        end
    end
  A_u = A(idx,:);
  b_u = b(idx,:);
  disp(['The original number of constraints: ', num2str(Numb0),', Number of compressed model constraints: ',num2str(numel(idx))]);
end

