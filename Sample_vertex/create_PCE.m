function [PCE] = create_PCE(Pst,Pnd,Qst,Qnd,Ust,Und,Keps,Pgmax,Pgmin,Qgmax,Qgmin,Umax,Umin)
    un1 = Pst(2:end) + Keps*sqrt(Pnd.^2) - Pgmax(2:end);
    un2 = Pst(2:end) - Keps*sqrt(Pnd.^2) - Pgmax(2:end);
    un3 = -Pst(2:end) + Keps*sqrt(Pnd.^2) + Pgmin(2:end);
    un4 = -Pst(2:end) - Keps*sqrt(Pnd.^2) + Pgmin(2:end);

    un5 = Qst(2:end) + Keps*sqrt(Qnd.^2) - Qgmax(2:end);
    un6 = Qst(2:end) - Keps*sqrt(Qnd.^2) - Qgmax(2:end);
    un7 = -Qst(2:end) + Keps*sqrt(Qnd.^2) + Qgmin(2:end);
    un8 = -Qst(2:end) - Keps*sqrt(Qnd.^2) + Qgmin(2:end);

    un9 = Ust(2:end) + Keps*sqrt(Und.^2) - Umax(2:end);
    un10 = Ust(2:end) - Keps*sqrt(Und.^2) - Umax(2:end);
    un11 = -Ust(2:end) + Keps*sqrt(Und.^2) + Umin(2:end);
    un12 = -Ust(2:end) - Keps*sqrt(Und.^2) + Umin(2:end);

    PCE = vertcat(un1,un2,un3,un4,un5,un6,un7,un8,un9,un10,un11,un12);
    % PCE = vertcat(un1,un2,un3,un4);
end

