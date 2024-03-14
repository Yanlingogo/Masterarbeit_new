e_st = zeros(Nbus,1);
e_st(1) = 1;

from_bus   = mpc.branch(:, F_BUS);                           
to_bus     = mpc.branch(:, T_BUS);
Cf         = sparse(1:Nbranch,from_bus,ones(Nbranch,1),Nbranch,Nbus);
Ct         = sparse(1:Nbranch,to_bus,ones(Nbranch,1),Nbranch,Nbus);

C = Cf - Ct;

R   = diag(mpc.branch(:,BR_R));
X   = diag(mpc.branch(:,BR_X));

A = [e_st' zeros(1,2*Nbus);
     C -2*R -2*X zeros(2,1);
     zeros(Nbus) -C' zeros(Nbus,Nbranch) e_st zeros(Nbus,1);
     zeros(Nbus,Nbus+Nbranch) -C' zeros(Nbus,1) e_st];
B = [zeros(Nbus,2*(Ngen-1));
     Cg_ns zeros(Nbus,Ngen-1);
     zeros(Nbus,Ngen-1), Cg_ns];
b = [1;zeros(Nbranch,1);Pd;Qd];

%% x = A^(-1)*(b - B*u)
% u = ([p^pcc; q^pcc] - [sum(Pd);sum(Qd)])/(Ngen-1);
% x = A^(-1)*(b - B*u);
% U = x(1:Nbus);
% Pij = x(Nbus+1:Nbus+Nbranch);
% Qij = x(Nbus+Nbranch+1:Nbus+2*Nbranch);
% l = (Pij.^2 + Qij.^2)./U;
% P_loss = sum(R*l);
% Q_loss = sum(Q*l);



