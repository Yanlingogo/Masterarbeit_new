function pf_eq = create_local_power_flow_equation_pol(vmag, vang, Pg, Qg, Gbus, Bbus, Pd, Qd, Cg)%,...
    %% power flow equation
   
   length_vang = length(vang);
   rep_vang = repmat(vang,1,length_vang);
   
   cos_vang = cos(rep_vang - rep_vang');
   sin_vang = sin(rep_vang - rep_vang');

   P = vmag .* (Gbus .* cos_vang + Bbus .* sin_vang)* vmag + Pd - Cg*Pg;
   Q = vmag .* (Gbus .* sin_vang - Bbus .* cos_vang)* vmag + Qd - Cg*Qg;
   
    pf_eq = vertcat(P,Q);
end