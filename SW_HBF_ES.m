function [SE_ts,W_rf,TC]=SW_HBF_ES(Pb,Ns,N_rf,H_all,noise_power,episilon,obj_fun)
tStart = tic;
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
%     Rc=rank(Hc,1e-4);
%     [~,~,V]=svd(Hc);
%     Vc=V(:,1:Rc);
%     Ht=Ht+Vc*Vc';
end
% Ht=Ht/Nc; 
% [Vt,Dt]=eig(Ht);
% Rt = Vt*sqrt(Dt)*Vt'; 
% 
% H_avr=zeros(Nt,Nt);
% for c=1:Nc
%     H_avr=H_avr+H_all(:,:,c)'*H_all(:,:,c);
% end
% H_avr=1/Nc*H_avr;
% Vtx=block_power_method(H_avr,N_rf,Np);
% Ax_tx=eye(Nr)-Vtx*Vtx';
% Coef_tx=Ax_tx'*Ax_tx;
% 
% Mix_tx=Ht-Coef_tx;
% 
% 
% 
% func_name=func2str(obj_fun);
% index=str2double(func_name(end));
% switch index
%     case 0
%         Srrgt=HBF_tx;
%     case 1
%         Srrgt=H_avr;
%     case 2
%         Srrgt=Coef_tx;
%     case 3
%         Srrgt=Rt;
%     case 4       
%         Srrgt=Mix_tx;
% end
F_rf=Exhaustive_search(Nt,N_rf,Ns,HBF_tx,coef,obj_fun);
% [SE_tx,F_rf,iter_tx]=Tabu_search(params_ts,F_ini_int,Ns,Srrgt,coef,obj_fun);

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
    
%     Rc=rank(HBF_info.H_eff{c},1e-4);
%     [~,~,V]=svd(HBF_info.H_eff{c});
%     Vc=V(:,1:Rc);
%     Hr=Hr+Vc*Vc';
%     
%     F_avr=F_avr+HBF_info.H_eff{c};   
end
% F_avr=1/Nc*F_avr;
% Hr=Hr/Nc;
% [Vr,Dr]=eig(Hr);
% Rr = Vr*sqrt(Dr)*Vr'; 
% 
% HBF_info.H_eff_avr=F_avr;
% 
% Vrx=block_power_method(F_avr,N_rf,Np);
% Ax_rx=eye(Nr)-Vrx*Vrx';
% Coef_rx=Ax_rx'*Ax_rx;
% 
% Mix_rx=Hr-Coef_rx;
% 
% 
% switch index
%     case 0
%         Srrgt=HBF_info;
%     case 1
%         Srrgt=F_avr;
%     case 2
%         Srrgt=Coef_rx;
%     case 3
%         Srrgt=Rr;
%     case 4
%         Srrgt=Mix_rx;    
% end



coef=1/noise_power;

W_rf=Exhaustive_search(Nr,N_rf,Ns,HBF_info,coef,obj_fun);

% vBest=0;
% sBest=0;
% N=Nr*N_rf;
% space_size=2^N;
% Nlb=N_rf;
% Nub=Nr+(N_rf-1)*(Nr-1);
% for i=1:space_size
%     symi=dec2bin(i,N);
%     vi=symi-'0';
%     Ni=sum(vi);
%     if Ni>=Nlb && Ni<=Nub
%        Wi=reshape(vi,[Nr,N_rf]);
%        if rank(Wi,1e-4)>=Ns
%            SEi=obj_fun(Srrgt,Wi,coef);
%            if SEi>=vBest
%                vBest=SEi;
%                sBest=Wi;
%            end
%        end    
%     end        
%     ccc=1;
% end


% W_ini=rand(Nr,N_rf);
% W_ini_int = mapping(W_ini,Ns);
% [~,W_rf,~]=Tabu_search(params_ts,W_ini_int,Ns,Srrgt,coef,obj_fun);
% W_rf=sBest;
SE_ts=compute_rate(HBF_info,W_rf,coef);
TC=toc(tStart);


ccc=1;
end