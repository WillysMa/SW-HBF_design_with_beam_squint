function [SE_ts,W_rf,TC]=HBF_TS(params_ts,Pb,Ns,N_rf,H_all,noise_power,episilon,obj_fun,Flag,InI)
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
    if Flag
        HBF_tx.H_eff{c}=Hc'*Hc;

    else
        HBF_tx.H_eff{c}=Hc'*InI(:,:,c)*Hc;
    end
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
[SE_tx,F_rf,data_ts]=Tabu_search(params_ts,F_ini_int,Ns,HBF_tx,coef,obj_fun);
% save('TS_data.mat','data_ts');
Fp=F_rf'*F_rf;
[Uq,Sq]=eig(Fp);
Froot=Uq*sqrt(inv(Sq))*Uq';% (F'F)^(-0.5)
% Hr=0;
% F_avr=zeros(Nr,Nr);
for c=1:Nc
    Hc=H_all(:,:,c);
    if Flag
        Qc=Hc*F_rf*Froot;
    else
        Qc=InI(:,:,c)*Hc*F_rf*Froot;
    end

    [~,Sc,Vc]=svd(Qc);
%     if N_rf>1
%         rank_apx=min(Ns,trace(Sc(1:N_rf,:).^2>episilon));
%     else
%         rank_apx=Ns;
%     end
    rank_apx=Ns;
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

W_ini=rand(Nr,N_rf);
W_ini_int = mapping(W_ini,Ns);
coef=1/noise_power;
[~,W_rf,~]=Tabu_search(params_ts,W_ini_int,Ns,HBF_info,coef,obj_fun);

SE_ts=compute_rate(HBF_info,W_rf,coef);


TC=toc(tStart);

% SE_all=zeros(Nc,N_rf);
% for i=1:N_rf
%     wb=W_rf(:,i);
%     for c=1:Nc
%         SE_all(c,i)=real(log2(det(eye(i)+coef*wb'*HBF_info.H_eff{c}*wb/(wb'*wb))));
%     end
% end
% altstyles = {'g-'; 'r-'; 'b-'; 'k-'; 'mo-.'; 'g*--'; 'r+--'};
% figure
% for i=1:N_rf
%     plot(SE_all(:,i),altstyles{i})
%     hold on
% end
% h_all=flipud(get(gca,'Children'));
% legend([h_all(1),h_all(2),h_all(3),h_all(4)], '$N_{RF} $ 1','$N_{RF} $ 2','$N_{RF} $ 3','$N_{RF} $ 4','Interpreter','latex')
% xlabel('Carrier');
% ylabel('SE');
% grid on
ccc=1;
end