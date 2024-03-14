function obj_p = create_obj_p(vmag, vang, Pg, Qg, Gbus, Bbus, Pd, Qd, Cg,id_slack)%,...
    %% power flow equation
   
   length_vang = length(vang);
   rep_vang = repmat(vang,1,length_vang);
   
   cos_vang = cos(rep_vang - rep_vang');
   sin_vang = sin(rep_vang - rep_vang');

   P = vmag .* (Gbus .* cos_vang + Bbus .* sin_vang)* vmag;
   Q = vmag .* (Gbus .* sin_vang - Bbus .* cos_vang)* vmag;
   
    obj_p = P(id_slack);
end