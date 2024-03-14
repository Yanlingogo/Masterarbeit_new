function [Psivv] = create_Psi(Psi_u_vvcos,Psi_l_vvcos,Psi_u_vvsin,Psi_l_vvsin,Psi_u_vv,Psi_l_vv, Dvvu_uu,Dvvl_uu,Dvvl_ll, Dcuu,Dcul,Dclu,Dcll,Dsuu,Dsll, v_u,v_l,Dvu,Dvl,Psi_vvcos0,Psi_vvsin0,Psi_vv0,Phi0,v0,idx_fr,idx_to)
    Psi1  = -Psi_u_vvcos + Psi_vvcos0 + cos(Phi0).*Dvvu_uu + v0(idx_fr).*v0(idx_to).*Dcuu + 1/4*(Dvvu_uu+Dcuu).^2;
    Psi2  = -Psi_u_vvcos + Psi_vvcos0 + cos(Phi0).*Dvvu_uu + v0(idx_fr).*v0(idx_to).*Dcul + 1/4*(Dvvu_uu+Dcul).^2;
    Psi3  =  Psi_l_vvcos - Psi_vvcos0 - cos(Phi0).*Dvvl_ll - v0(idx_fr).*v0(idx_to).*Dclu + 1/4*(Dvvl_ll-Dclu).^2;
    Psi4  =  Psi_l_vvcos - Psi_vvcos0 - cos(Phi0).*Dvvl_ll - v0(idx_fr).*v0(idx_to).*Dcll + 1/4*(Dvvl_ll-Dcll).^2;
    
    Psi5  = -Psi_u_vvsin + Psi_vvsin0 + sin(Phi0).*Dvvu_uu + v0(idx_fr).*v0(idx_to).*Dsuu + 1/4*(Dvvu_uu+Dsuu).^2;
    Psi6  = -Psi_u_vvsin + Psi_vvsin0 + sin(Phi0).*Dvvl_uu + v0(idx_fr).*v0(idx_to).*Dsuu + 1/4*(Dvvl_uu+Dsuu).^2;
    Psi7  =  Psi_l_vvsin - Psi_vvsin0 - sin(Phi0).*Dvvl_ll - v0(idx_fr).*v0(idx_to).*Dsll + 1/4*(Dvvl_ll-Dsll).^2;
    Psi8  =  Psi_l_vvsin - Psi_vvsin0 - sin(Phi0).*Dvvu_uu - v0(idx_fr).*v0(idx_to).*Dsll + 1/4*(Dvvu_uu-Dsll).^2;
    
    Psi9  = -Psi_u_vv + v_u.^2; 
    Psi10 = -Psi_u_vv + v_l.^2;

    Psi11 = Psi_l_vv - Psi_vv0 - 2*v0.*Dvu;
    Psi12 = Psi_l_vv - Psi_vv0 - 2*v0.*Dvl;

    Psivv = vertcat(Psi1,Psi2,Psi3,Psi4,Psi5,Psi6,Psi7,Psi8,Psi9,Psi10,Psi11,Psi12);
    %8*Nbranch+4*Nbus
end

