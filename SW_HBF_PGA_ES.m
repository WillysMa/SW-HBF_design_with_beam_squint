function [SE_ts,W_rf]=SW_HBF_PGA_ES(params_pga,Pb,Ns,N_rf,H_all,noise_power,episilon,delta,obj_fun,obj_fun_pg)

[Nr,Nt,Nc]=size(H_all);
Np=10;
power_avr=Pb/Ns;
F_ini=rand(Nt,N_rf);
F_ini_int = mapping(F_ini,Ns);
coef=power_avr/noise_power;
Ht=0;
for c=1:Nc
    Hc=H_all(:,:,c);
    HBF_tx.H_eff{c}=Hc'*Hc;
    Rc=rank(Hc,1e-4);
    [~,~,V]=svd(Hc);
    Vc=V(:,1:Rc);
    Ht=Ht+Vc*Vc';
end
Ht=Ht/Nc; 
[Vt,Dt]=eig(Ht);
Rt = Vt*sqrt(Dt)*Vt'; 

H_avr=zeros(Nt,Nt);
for c=1:Nc
    H_avr=H_avr+H_all(:,:,c)'*H_all(:,:,c);
end
H_avr=1/Nc*H_avr;
Vtx=block_power_method(H_avr,N_rf,Np);
Ax_tx=eye(Nr)-Vtx*Vtx';
Coef_tx=Ax_tx'*Ax_tx;

Mix_tx=Ht-Coef_tx;



func_name=func2str(obj_fun);
index=str2double(func_name(end));
switch index
    case 0
        Srrgt=HBF_tx;
    case 1
        Srrgt=H_avr;
    case 2
        Srrgt=Coef_tx;
    case 3
        Srrgt=Rt;
    case 4       
        Srrgt=Mix_tx;
end

[SE_pg,F_pg,~]=PGA_new(params_pga,F_ini_int,HBF_tx,coef,obj_fun_pg);
if rank(F_pg,1e-4)<Ns
    F_pg=F_ini_int;
end
[F_rf_ini,Vec_pgr,Index_bf,Npart1]=TS_pga_params(F_pg,Ns,delta);
if length(Index_bf)>=1
    F_rf=PGA_Ehaustive_Search(Nt,N_rf,Ns,Index_bf,Vec_pgr,Npart1,Srrgt,coef,obj_fun) ;

else
   F_rf= F_pg;
end
% F_rf
Fp=F_rf'*F_rf;
[Uq,Sq]=eig(Fp);
Froot=Uq*sqrt(inv(Sq))*Uq';% (F'F)^(-0.5)
Hr=0;
F_avr=zeros(Nr,Nr);
for c=1:Nc
    Hc=H_all(:,:,c);

    Qc=Hc*F_rf*Froot;
    [~,Sc,Vc]=svd(Qc);
    if N_rf>1
        rank_apx=min(Ns,trace(Sc(1:N_rf,:).^2>episilon));
    else
        rank_apx=Ns;
    end
    HBF_info.approx_rank{c}=rank_apx;
    baseband_bf=Froot*Vc(:,1:rank_apx)*sqrt(power_avr);  
    precoder=F_rf*baseband_bf;
    HBF_info.transmitter_bf{c}=precoder;   
    HBF_info.H_eff{c}=Hc*precoder*(Hc*precoder)';
    
    Rc=rank(HBF_info.H_eff{c},1e-4);
    [~,~,V]=svd(HBF_info.H_eff{c});
    Vc=V(:,1:Rc);
    Hr=Hr+Vc*Vc';
    
    F_avr=F_avr+HBF_info.H_eff{c};   
end
F_avr=1/Nc*F_avr;
Hr=Hr/Nc;
[Vr,Dr]=eig(Hr);
Rr = Vr*sqrt(Dr)*Vr'; 

HBF_info.H_eff_avr=F_avr;

Vrx=block_power_method(F_avr,N_rf,Np);
Ax_rx=eye(Nr)-Vrx*Vrx';
Coef_rx=Ax_rx'*Ax_rx;

Mix_rx=Hr-Coef_rx;


switch index
    case 0
        Srrgt=HBF_info;
    case 1
        Srrgt=F_avr;
    case 2
        Srrgt=Coef_rx;
    case 3
        Srrgt=Rr;
    case 4
        Srrgt=Mix_rx;    
end


W_ini=rand(Nr,N_rf);
W_ini_int = mapping(W_ini,Ns);
coef=1/noise_power;


[SE_pg,W_pg,~]=PGA_new(params_pga,W_ini_int,HBF_info,coef,obj_fun_pg);
if rank(W_pg,1e-4)<Ns
    W_pg=W_ini_int;
end
[W_rf_ini,Vec_pgr,Index_bf,Npart1]=TS_pga_params(W_pg,Ns,delta);
if length(Index_bf)>=1
   W_rf=PGA_Ehaustive_Search(Nr,N_rf,Ns,Index_bf,Vec_pgr,Npart1,Srrgt,coef,obj_fun) ;
   ccc=1;

else
    W_rf=W_pg;
end


SE_ts=compute_rate_rx(HBF_info,W_rf,coef);




ccc=1;
end

%% defined function
function [W_rf_ini,Vec_pgr,Index_bf,Npart1]=TS_pga_params(W_pg,Ns,delta)
[N,N_rf]=size(W_pg);
W_pgr = rounding(W_pg,delta);
rank_pgr=rank(W_pgr,1e-4);
if rank_pgr<Ns
        W_pgr=W_pg;        
end
    
% W_pgr= W_pg;
Vec_pgr=vec(W_pgr);
W_rf_ini = mapping(W_pgr,Ns);
    
V_rf_ini=vec(W_rf_ini);
Diff_bf=abs(Vec_pgr-V_rf_ini);
[Index_bf,~]=find(vec(Diff_bf)>0);% Index_bf denotes the index in Vec_pgr

Determined_ele=setdiff((1:N*N_rf),Index_bf);
Vec_pgr_determined=Vec_pgr(Determined_ele);
Npart1=sum(Vec_pgr_determined);
end