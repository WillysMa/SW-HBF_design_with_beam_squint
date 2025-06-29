function SE_iter=WMMSE_MISO(Hall,Wgtu,Pb,noise_power,iter_max,tol)
[Nt,Nu,Nc]=size(Hall);
Fd_all=zeros(Nt,Nu,Nc);
for nc=1:Nc
    Hc=Hall(:,:,nc);
    [Uc,Dc,Vc]=svd(Hc);
    Fd_all(:,:,nc)=sqrt(Pb/Nu)*Uc(:,1:Nu);

end
Uall=zeros(Nu,Nc);
Wall=zeros(Nu,Nc);
Dall=zeros(Nu,Nc);

% Wgtu=ones(1,Nu);
SE_iter=1;SE=[];
% iter_max=1e+3;
acc=1e-4;
for iter=1:iter_max
SE_old=SE_iter;

tmpSE=0;
for nc=1:Nc
    tmpY=0;
    for ii=1:Nu          
        Interference_power=norm(Hall(:,ii,nc)'*Fd_all(:,:,nc),2)^2+noise_power;
        target=Hall(:,ii,nc)'*Fd_all(:,ii,nc);
        target_power=abs(target)^2;
        
        Uall(ii,nc)=target/Interference_power;

        Wall(ii,nc)=1+target_power/(Interference_power-target_power);


        Dall(ii,nc)=Wgtu(ii)*Wall(ii,nc)*Uall(ii,nc);

        tmpY=tmpY+Wgtu(ii)*Wall(ii,nc)*abs(Uall(ii,nc))^2*Hall(:,ii,nc)*Hall(:,ii,nc)';

        tmpSE=tmpSE+Wgtu(ii)*log2(1+Wall(ii,nc));
    end%->ii=1:Nu
%     Yall(:,:,nc)=tmpY;

    J=tmpY;
    tmpF=Hall(:,:,nc)*diag(Dall(:,nc));
    Fnc=pinv(J)*tmpF;
    if  norm(Fnc,'fro')^2<= Pb
          Fnc=pinv(J)*tmpF;
          mu=0;
    else
    
        mu_ub=10*norm(tmpF,'fro')/sqrt(Pb);
        mu_lb=0;
        while abs(mu_ub-mu_lb)/abs(mu_ub)>=acc
            mu=(mu_ub+mu_lb)/2;
            Fnc=pinv(J+mu*eye(Nt))*tmpF;
            if norm(Fnc,'fro')^2<= Pb
                mu_ub=mu;
            else
                mu_lb=mu;
            end
    
        end
    end
    power=norm(Fnc,'fro')^2;
    Fd_all(:,:,nc)=Fnc;
    
end%->nc=1:Nc

   
SE_iter=tmpSE/Nc;
SE=[SE SE_iter];

err_ratlv=abs(SE_iter-SE_old)/abs(SE_old);
if err_ratlv<=tol
    break
end

end
% figure
% plot(SE)
% xlabel('Iterations')
% ylabel('Sum rate (bits/s/Hz)')
% grid on
% box on
% title('WMMSE')
cccc=1;
end