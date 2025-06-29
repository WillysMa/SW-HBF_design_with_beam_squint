function [SE_Apga_TS,SE_Rdn_TS,SE_rs]=Multi_InI(params_pga,params_ts,params_TbRs,Pb,Ns,N_rf,...
    Hcc,noise_power,episilon,delta,obj_fun,obj_fun_pg,Flag_Identity_Ini,Proj_Comb,ii_max)


parfor ii=1:ii_max
 [SW_Apga_TS(ii),~]=HBF_APGA_TS(params_pga,params_ts,Pb,Ns,N_rf,Hcc,noise_power,episilon,delta,obj_fun,obj_fun_pg,Flag_Identity_Ini,Proj_Comb);
 [SW_Rdn_TS(ii),~]=HBF_TS(params_ts,Pb,Ns,N_rf,Hcc,noise_power,episilon,obj_fun,Flag_Identity_Ini,Proj_Comb);      
 [SW_rs(ii),~]=HBF_APGA_TbRs(params_pga,params_TbRs,Pb,Ns,N_rf,Hcc,noise_power,episilon,delta,obj_fun,obj_fun_pg,Flag_Identity_Ini,Proj_Comb); 

end


SE_Apga_TS=max(SW_Apga_TS);
SE_Rdn_TS=max(SW_Rdn_TS);
SE_rs=max(SW_rs);