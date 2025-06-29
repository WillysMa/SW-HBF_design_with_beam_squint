function [SE_pgats,SE_es,SE_rdn]=TwoStepHBF_design_MISO(params_ts,params_pga,Hall,Wgtu,Nrf,Pb,noise_power,iter_max,tol)
[Nt,Nu,Nc]=size(Hall);
delta=0.1;
[Frf_pga,Fbb_all,~,~]=FP_relaxed_HBF(params_pga,Hall,Wgtu,Nrf,Pb,noise_power,iter_max,tol);

Frf_int = mapping(Frf_pga,Nrf);%random
SE_rdn=rate_calcu_MISO(Hall,Frf_int,Fbb_all,noise_power,Wgtu);

[Frf_pgr,Vec_pgr,Index_bf,Npart1]=TS_pga_params(Frf_pga,Nu,delta);

lenVar=length(Index_bf);
space_size=2^lenVar;
sBest=0;
vBest=0;
% for ii=space_size-1:-1:0
% %     ii
%    sp=dec2bin(ii,lenVar);
%    sp_binary=sp-'0';
%    Vec_sp=Vec_pgr;
%    Vec_sp(Index_bf)=sp_binary;
%    X_rshp=reshape(Vec_sp,[Nt,Nrf]);
%    Rate=rate_calcu_MISO(Hall,X_rshp,Fbb_all,noise_power,Wgtu);
% 
%    if Rate>=vBest
%       sBest=X_rshp;
%       vBest=Rate;
%    end
%    
%    ccc=1;
% end
SE_es=vBest;


[Obj_opt,Frf_opt,iter]=TS_pga_MU_new(params_ts,Frf_pgr,Nu,Vec_pgr,Index_bf,Npart1,Hall,Fbb_all,noise_power,Wgtu);
SE_pgats=Obj_opt;
% Frf=Frf_opt;

end