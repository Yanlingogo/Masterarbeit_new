function [cossin] = create_ineq_cossin(Dcuu,Dcul,Dclu,Dcll,Dsuu,Dsul,Dslu,Dsll,Dphiu,Dphil,Phi0)
% <=0
    c1 = -Dcuu-sin(Phi0).*Dphiu+1/2*(Dphiu).^2;
    c2 = -Dcul-sin(Phi0).*Dphil+1/2*(Dphil).^2;
    c3 =  Dclu+sin(Phi0).*Dphiu+1/2*(Dphiu).^2;
    c4 =  Dcll+sin(Phi0).*Dphil+1/2*(Dphil).^2;

    s1 = -Dsuu+cos(Phi0).*Dphiu+1/2*(Dphiu).^2;
    s2 = -Dsul+cos(Phi0).*Dphil+1/2*(Dphil).^2;
    s3 =  Dslu-cos(Phi0).*Dphiu+1/2*(Dphiu).^2;
    s4 =  Dsll-cos(Phi0).*Dphil+1/2*(Dphil).^2;
    
    cossin = vertcat(c1,c2,c3,c4,s1,s2,s3,s4);
    %8*Nbranch
end

