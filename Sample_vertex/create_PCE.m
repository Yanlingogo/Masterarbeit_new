function [PCE] = create_PCE(Pst,Pnd,Ust,Und, Keps,Pgmax,Pgmin,Umax,Umin)
    un1 = Pst(2:end) + Keps*sqrt(Pnd.^2) - Pgmax(2:end);
    un2 = Pst(2:end) - Keps*sqrt(Pnd.^2) - Pgmax(2:end);
    un3 = -Pst(2:end) + Keps*sqrt(Pnd.^2) + Pgmin(2:end);
    un4 = -Pst(2:end) - Keps*sqrt(Pnd.^2) + Pgmin(2:end);

    un5 = Ust(2:end) + Keps*sqrt(Und.^2) - Umax(2:end);
    un6 = Ust(2:end) - Keps*sqrt(Und.^2) - Umax(2:end);
    un7 = -Ust(2:end) + Keps*sqrt(Und.^2) + Umin(2:end);
    un8 = -Ust(2:end) - Keps*sqrt(Und.^2) + Umin(2:end);

    PCE = vertcat(un1,un2,un3,un4,un5,un6,un7,un8);
    % PCE = vertcat(un1,un2,un3,un4);
end

