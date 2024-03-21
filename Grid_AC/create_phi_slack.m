function [phi_slack] = create_phi_slack(vang,slack_phi,C,Phimax,Phimin)
    phi1 = C*vang + slack_phi - Phimax;
    phi2 = -C*vang + slack_phi + Phimin;

    phi_slack = vertcat(phi1,phi2);

    %2Nbranch
end

