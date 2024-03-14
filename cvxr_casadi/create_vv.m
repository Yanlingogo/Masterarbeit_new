function [vv] = create_vv(Dvvu_uu,Dvvu_ll,Dvvu_ul,Dvvu_lu,Dvvl_uu,Dvvl_ll,Dvvl_ul,Dvvl_lu, Dvu,Dvl,v0,idx_fr,idx_to)
% checked 
    vv1 = -Dvvu_uu+Dvu(idx_fr).*v0(idx_to)+v0(idx_fr).*Dvu(idx_to)+1/4*(Dvu(idx_fr)+Dvu(idx_to)).^2;
    vv2 = -Dvvu_ll+Dvl(idx_fr).*v0(idx_to)+v0(idx_fr).*Dvl(idx_to)+1/4*(Dvl(idx_fr)+Dvl(idx_to)).^2;
    vv3 = -Dvvu_ul+Dvu(idx_fr).*v0(idx_to)+v0(idx_fr).*Dvl(idx_to)+1/4*(Dvu(idx_fr)+Dvl(idx_to)).^2;
    vv4 = -Dvvu_lu+Dvl(idx_fr).*v0(idx_to)+v0(idx_fr).*Dvu(idx_to)+1/4*(Dvl(idx_fr)+Dvu(idx_to)).^2;
    
    vv5 =  Dvvl_uu-Dvu(idx_fr).*v0(idx_to)-v0(idx_fr).*Dvu(idx_to)+1/4*(Dvu(idx_fr)-Dvu(idx_to)).^2;
    vv6 =  Dvvl_ll-Dvl(idx_fr).*v0(idx_to)-v0(idx_fr).*Dvl(idx_to)+1/4*(Dvl(idx_fr)-Dvl(idx_to)).^2;
    vv7 =  Dvvl_ul-Dvu(idx_fr).*v0(idx_to)-v0(idx_fr).*Dvl(idx_to)+1/4*(Dvu(idx_fr)-Dvl(idx_to)).^2;
    vv8 =  Dvvl_lu-Dvl(idx_fr).*v0(idx_to)-v0(idx_fr).*Dvu(idx_to)+1/4*(Dvl(idx_fr)-Dvu(idx_to)).^2;

    vv= vertcat(vv1,vv2,vv3,vv4,vv5,vv6,vv7,vv8);
    %8*Nbranch
end


