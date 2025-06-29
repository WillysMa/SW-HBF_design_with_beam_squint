function [S_og,S_pg,S_pgr]=dimension_reduction(params_pga,H_all,Pb,Ns,N_rf,InI,noise_power,obj_fun_pg,delta)
Set_fsb=[0 1];

[Nr,Nt,Nc]=size(H_all);

power_avr=Pb/Ns;
F_ini=rand(Nt,N_rf);
F_ini_int = mapping(F_ini,Ns);
coef=power_avr/noise_power;

for c=1:Nc
    Hc=H_all(:,:,c);
    HBF_tx.H_eff{c}=Hc'*InI.proj_comb{c}*Hc;
end


[SE_pg,F_pg,~]=PGA_new(params_pga,F_ini_int,HBF_tx,coef,obj_fun_pg);
diff_pg=setdiff(F_pg,Set_fsb);
S_pg=length(diff_pg);

F_sub= rounding(diff_pg,delta);
diff_refined=setdiff(F_sub,Set_fsb);
S_pgr=length(diff_refined);

S_og=Nt*N_rf;

ccc=1;
end

