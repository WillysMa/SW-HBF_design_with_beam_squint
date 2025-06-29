function tmpSE=PSHBF_MISO_WeiYu(Hall,Wgtu,Nrf,Pb,noise_power,PhaseSet,iter_max,tol)
[Nt,Nu,Nc]=size(Hall);
acc=1e-4;
Fbb_all=zeros(Nrf,Nu,Nc);
%% design analog precoder
Havr=0;
for nc=1:Nc
    Hc=Hall(:,:,nc);
    Havr=Havr+Hc*Hc';
end
Havr=Havr/Nc;
coef_t=Pb/(Nu*noise_power);
Frf=analog_bf_compute(coef_t,Nt,Nrf,Havr);
Frf=angle_qt(Frf,PhaseSet);

%% design digital precoder
for nc=1:Nc
    Hc=Hall(:,:,nc);
    Heff=Frf'*Hc*Hc'*Frf;
    [U,D,V]=svd(Heff);
    Fbb=V(:,1:Nu);    
    Fbb_all(:,:,nc)=sqrt(Pb)*Fbb/norm(Frf*Fbb,'fro');
end


Uall=zeros(Nu,Nc);
Wall=zeros(Nu,Nc);
Dall=zeros(Nu,Nc);
Yall=zeros(Nt,Nt,Nc);


SE=[];
tmpSE=1;
for iter=1:iter_max
SE_old=tmpSE;

for nc=1:Nc
    tmpY=0;
    for ii=1:Nu          
        Interference_power=norm(Hall(:,ii,nc)'*Frf*Fbb_all(:,:,nc),2)^2+noise_power;
        target=Hall(:,ii,nc)'*Frf*Fbb_all(:,ii,nc);
        target_power=abs(target)^2;
        
        Uall(ii,nc)=target/Interference_power;

        Wall(ii,nc)=1+target_power/(Interference_power-target_power);

        tmpY=tmpY+Wgtu(ii)*Wall(ii,nc)*abs(Uall(ii,nc))^2*Hall(:,ii,nc)*Hall(:,ii,nc)';
        Dall(ii,nc)=Wgtu(ii)*Wall(ii,nc)*Uall(ii,nc);
    end
    Yall(:,:,nc)=tmpY;

end

for nc=1:Nc
    
    J=Frf'*Yall(:,:,nc)*Frf;
    tmpF=Frf'*Hall(:,:,nc)*diag(Dall(:,nc));
    Fnc=pinv(J)*tmpF;
    if  norm(Frf*Fnc,'fro')^2<= Pb
          Fnc=pinv(J)*tmpF;
          mu=0;
    else
    
        mu_ub=10*norm(tmpF,'fro')/sqrt(Pb);
        mu_lb=0;
        while abs(mu_ub-mu_lb)/abs(mu_ub)>=acc
            mu=(mu_ub+mu_lb)/2;
            Fnc=pinv(J+mu*Frf'*Frf)*tmpF;
            if norm(Frf*Fnc,'fro')^2<= Pb
                mu_ub=mu;
            else
                mu_lb=mu;
            end
    
        end
    end
%     power=norm(Frf*Fnc,'fro')^2;
Fbb_all(:,:,nc)=Fnc;

end

temp_x1=0;
for nc=1:Nc
    for ii=1:Nu
        Interference_power=norm(Hall(:,ii,nc)'*Frf*Fbb_all(:,:,nc),2)^2+noise_power;
        target=Hall(:,ii,nc)'*Frf*Fbb_all(:,ii,nc);
        target_power=abs(target)^2;
        
        SINR=target_power/(Interference_power-target_power);

        temp_x1=temp_x1+Wgtu(ii)*log2(1+SINR);
    end
end
tmpSE=temp_x1/Nc;

SE=[SE tmpSE];


err_ratlv=abs(tmpSE-SE_old)/abs(SE_old);
if err_ratlv<=tol
    break
end

ccc=1;
end
% tmpSE

% figure
% plot(SE)
% xlabel('Iterations')
% ylabel('Sum rate (bits/s/Hz)')
% grid on
% box on
% title('PS-WMMSE')
ccc=1;
end