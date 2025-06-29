function [W_rf_ini,Vec_pgr,Index_bf,Npart1]=TS_pga_params(W_pg,Ns,delta)
[N,N_rf]=size(W_pg);
W_pgr = rounding(W_pg,delta);
rank_pgr=rank(W_pgr,1e-4);
if rank_pgr<Ns
        W_pgr=W_pg;        
end
    
% W_pgr
% W_pg
Vec_pgr=vec(W_pgr);
W_rf_ini = mapping(W_pgr,Ns);
    
V_rf_ini=vec(W_rf_ini);
Diff_bf=abs(Vec_pgr-V_rf_ini);
[Index_bf,~]=find(vec(Diff_bf)>0);% Index_bf denotes the index in Vec_pgr

Determined_ele=setdiff((1:N*N_rf),Index_bf);
Vec_pgr_determined=Vec_pgr(Determined_ele);
Npart1=sum(Vec_pgr_determined);% number of total active antennas so far

ccc=1;
end