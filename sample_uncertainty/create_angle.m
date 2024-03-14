function [cons_ang] = create_angle(vang,C,Phimin,Phimax)
    ineq1 = Phimin - C*vang;
    ineq2 = C*vang - Phimax;
    cons_ang = vertcat(ineq1,ineq2);
end

