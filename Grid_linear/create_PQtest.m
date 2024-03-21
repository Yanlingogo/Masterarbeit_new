function [PQ_test] = create_PQtest(P,Q,Pg0,Qg0,id)

    Test_P = P(id) - Pg0(id);
    Test_Q = Q(id) - Qg0(id);
    PQ_test = vertcat(Test_P,Test_Q);

end

