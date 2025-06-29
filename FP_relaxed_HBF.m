function [Frf,Fbb_all,SE_list,iter_data]=FP_relaxed_HBF(params_pga,Hall,Wgtu,Nrf,Pb,noise_power,iter_max,tol)
% iter_max=1e+3;
acc=1e-4;
[Nt,Nu,Nc]=size(Hall);
Frf=rand([Nt,Nrf]);
Fbb_all=zeros(Nrf,Nu,Nc);
for nc=1:Nc
    Hc=Hall(:,:,nc);
    Heff=Frf'*Hc*Hc'*Frf;
    [U,D,V]=svd(Heff);
    Fbb=V(:,1:Nu);    
    Fbb_all(:,:,nc)=sqrt(Pb)*Fbb/norm(Frf*Fbb,'fro');
end

Rall=zeros(Nu,Nc);
Qall=zeros(Nu,Nc);
Zall=zeros(Nrf,Nrf,Nc);
Yall=zeros(Nt,Nt,Nc);
Dall=zeros(Nu,Nc);

Cache=[];SE_list=[];
obj_iter=1;tmpSE=1;
for iter=1:iter_max
obj_old=obj_iter;
SE_old=tmpSE;
Xi=0;
tmpZY=0;
for nc=1:Nc
    tmpZ=0;
    tmpY=0;
    for ii=1:Nu          
        Interference_power=norm(Hall(:,ii,nc)'*Frf*Fbb_all(:,:,nc),2)^2+noise_power;
        target=Hall(:,ii,nc)'*Frf*Fbb_all(:,ii,nc);
        target_power=abs(target)^2;
        
        Rall(ii,nc)=target_power/(Interference_power-target_power);

        temp_coef=sqrt(Wgtu(ii)*(Rall(ii,nc)+1));
        Qall(ii,nc)=temp_coef*target/Interference_power;

        Xi=Xi+temp_coef*Qall(ii,nc)*Hall(:,ii,nc)*Fbb_all(:,ii,nc)';
        tmpZ=tmpZ+Fbb_all(:,ii,nc)*Fbb_all(:,ii,nc)';

        tmpY=tmpY+abs(Qall(ii,nc))^2*Hall(:,ii,nc)*Hall(:,ii,nc)';
        Dall(ii,nc)=Qall(ii,nc)*temp_coef;
    end
    tmpZY=tmpZY+kron(tmpZ.',tmpY);
    Zall(:,:,nc)=tmpZ;
    Yall(:,:,nc)=tmpY;

end

[Frf,iter_data]=PGA_MU(params_pga,Frf,Xi,Yall,Zall);

% cvx_begin quiet
%     variable Frf(Nt,Nrf) %nonnegative
%     expression term1
%     for nc=1:Nc
%         tempS=Yall(:,:,nc)^(0.5)*Frf*Zall(:,:,nc)^(0.5);
%         term1=term1+sum_square_abs(vec(tempS));
%     end    
%     minimize( term1-2*real(trace(Xi*Frf')) )
%     for rr=1:Nt
%         for col=1:Nrf
%             Frf(rr,col)<=1;
%         end
%     end
% 
%     cvx_end
% Frf_cvx=Frf;

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
temp_x2=0;
temp_x3=0;
for nc=1:Nc
    for ii=1:Nu
        Interference_power=norm(Hall(:,ii,nc)'*Frf*Fbb_all(:,:,nc),2)^2+noise_power;
        target=Hall(:,ii,nc)'*Frf*Fbb_all(:,ii,nc);
        target_power=abs(target)^2;
        
        SINR=target_power/(Interference_power-target_power);

        temp_x1=temp_x1+Wgtu(ii)*log2(1+SINR);
        temp_x2=temp_x2+Wgtu(ii)*SINR;
        

        temp_x3=temp_x3+Wgtu(ii)*(SINR+1)*target_power/Interference_power;
    end
end
tmpSE=temp_x1/Nc;
obj_iter=tmpSE-(temp_x2-temp_x3)/(Nc*log(2));

SE_list=[SE_list tmpSE];
Cache=[Cache obj_iter];

err_ratlv=abs(tmpSE-SE_old)/abs(SE_old);
if err_ratlv<=tol
    break
end

ccc=1;
end
% figure
% plot(Cache)
% xlabel('Iterations')
% ylabel('Objective')
% grid on
% box on

% figure
% plot(SE_list)
% xlabel('Iterations')
% ylabel('Sum rate (bits/s/Hz)')
% grid on
% box on
% title('FP')
ccc=1;
end