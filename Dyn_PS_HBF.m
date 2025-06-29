function [SE,TC]=Dyn_PS_HBF(Hcc,N_rf,Ns,Pb,noise_power,PhaseSet,iter_max,tol_obj)
tStart = tic;
[Nr,Nt,Nc]=size(Hcc);
rhot=100/Nt;rhor=100/Nr;acc=1e-4;
% iter_max=50;tol_obj=1e-6;
%% initialization

Havr=0;
for nc=1:Nc
    Havr=Havr+Hcc(:,:,nc);
end
Havr=Havr/Nc;

[U,D,V]=svd(Havr);

Ntset=int32(Nt/N_rf);Nrset=int32(Nr/N_rf);
Frf_sub=zeros(Ntset,N_rf);
Wrf_sub=zeros(Nrset,N_rf);
for n=1:N_rf
    Frf_sub(:,n)=V((n-1)*Ntset+1:n*Ntset,n);
    Wrf_sub(:,n)=U((n-1)*Nrset+1:n*Nrset,n);
end
Frf_sub_qtz=angle_qt(Frf_sub,PhaseSet);
Wrf_sub_qtz=angle_qt(Wrf_sub,PhaseSet);

Frf=zeros(Nt,N_rf);
Wrf=zeros(Nr,N_rf);
for n=1:N_rf
    Frf((n-1)*Ntset+1:n*Ntset,n)=Frf_sub_qtz(:,n);
    Wrf((n-1)*Nrset+1:n*Nrset,n)=Wrf_sub_qtz(:,n);
end


Xall=zeros(Nt,Ns,Nc);
Yall=zeros(Nr,Ns,Nc);
Fbb_all=zeros(N_rf,Ns,Nc);
Wbb_all=zeros(N_rf,Ns,Nc);
for nc=1:Nc
    Heff=Wrf'*Hcc(:,:,nc)*Frf;
    [Ue,De,Ve]=svd(Heff);
    Fbb_all(:,:,nc)=sqrt(Pb)*Ve(:,1:Ns)/norm(Frf*Ve(:,1:Ns),'fro');
    Wbb_all(:,:,nc)=Ue(:,1:Ns);
    
    Xall(:,:,nc)=Frf*Fbb_all(:,:,nc);
    Yall(:,:,nc)=Wrf*Wbb_all(:,:,nc); 

end

%% main loop

Mall=zeros(Ns,Ns,Nc);

SE_iter=0;Obj_iter=0;
SE_cache=[];
Obj_cache=[];
for ii=1:iter_max
    SE_old=SE_iter;
    Obj_old=Obj_iter;
    X_set=[];
    Fbb_set=[];
    for nc=1:Nc
        Hc=Hcc(:,:,nc);
        Xc=Xall(:,:,nc); 
        Yc=Yall(:,:,nc);

        tmpA=Hc*Xc*(Hc*Xc)'+noise_power*eye(Nr);
        Mc=inv(Yc'*tmpA*Yc-2*real(Yc'*Hc*Xc)+eye(Ns));
        Mall(:,:,nc)=Mc;

        tmpX=Hc'*Yc*Mc+1/(2*rhot)*Frf*Fbb_all(:,:,nc);
        J=Hc'*Yc*Mc*Yc'*Hc+1/(2*rhot)*eye(Nt);
        Xc=inv(J)*tmpX;
        if  norm(Xc,'fro')^2<= Pb
              Xc=inv(J)*tmpX;
        else
        
            mu_ub=norm(tmpX,'fro')/sqrt(Pb);
            mu_lb=0;
            while abs(mu_ub-mu_lb)/abs(mu_ub)>=acc
                mu=(mu_ub+mu_lb)/2;
                Xc=inv(J+mu*eye(Nt))*tmpX;
                if norm(Xc,'fro')^2<= Pb
                    mu_ub=mu;
                else
                    mu_lb=mu;
                end
        
            end
        end
        Xall(:,:,nc)=Xc;
%         ccc=1;
        Fbb_all(:,:,Nc)=pinv(Frf)*Xc;

        X_set=[X_set Xc];     
        Fbb_set=[Fbb_set Fbb_all(:,:,Nc)];
    end


    Frf=find_BF(Fbb_set,X_set,PhaseSet,Nt,N_rf);
    
    Y_set=[];
    Wbb_set=[];
    for nc=1:Nc
        Hc=Hcc(:,:,nc);
        Xc=Xall(:,:,nc);
        Mc=Mall(:,:,nc);
        Ac=Hc*Xc*(Hc*Xc)'+noise_power*eye(Nr);

        tem_med=Hc*Xc*Mc+ 1/(2*rhor)*Wrf*Wbb_all(:,:,nc);
        Coef_mtx=kron(Mc.',Ac)+1/(2*rhor)*eye(Ns*Nr);
        yvec=inv(Coef_mtx)*vec(tem_med);

        for ns=1:Ns
            Yall(:,ns,nc)=yvec((ns-1)*Nr+1:ns*Nr);
        end

        Wbb_all(:,:,nc)=pinv(Wrf)*Yall(:,:,nc);

        Y_set=[Y_set Yall(:,:,nc)];
        Wbb_set=[Wbb_set  Wbb_all(:,:,nc)];

    end
    
    Wrf=find_BF(Wbb_set,Y_set,PhaseSet,Nr,N_rf);

    SE_iter=achievable_rate(Hcc,Frf,Fbb_all,Wrf,Wbb_all,noise_power,Ns);
    Obj_iter=Obj_cal(Hcc,Frf,Fbb_all,Wrf,Wbb_all,Mall,Xall,Yall,rhot,rhor);
    
    SE_cache=[SE_cache SE_iter];
    Obj_cache=[Obj_cache Obj_iter];
    err=abs(Obj_iter-Obj_old)/abs(Obj_old);
    if err<=tol_obj
        break
    end
    ccc=1;
end

for nc=1:Nc
    Fbb_all(:,:,nc)=sqrt(Pb)*Fbb_all(:,:,nc)/norm(Frf*Fbb_all(:,:,nc),'fro');
end
SE=achievable_rate(Hcc,Frf,Fbb_all,Wrf,Wbb_all,noise_power,Ns);

% figure
% subplot(121)
% plot(SE_cache)
% xlabel('Iteration')
% ylabel('SE')
% box on
% grid on
% 
% subplot(122)
% plot(Obj_cache)
% xlabel('Iteration')
% ylabel('Objective value')
% box on
% grid on
TC=toc(tStart);
ccc=1;
end
%% defined functions
function BF=find_BF(Fbb_set,X_set,PhaseSet,Nt,N_rf)
%     Delta=2*pi/2^b;
    Varpi=zeros(Nt,N_rf);
%     angle_all=zeros(Nt,N_rf);
    for m=1:N_rf
        for n=1:Nt        
            Varpi(n,m)=Fbb_set(m,:)*X_set(n,:)';
%             angle_tmp=-angle(Varpi(n,m));
%             if angle_tmp>=0
%                 angle_all(n,m)=round(angle_tmp/Delta)*Delta+2*pi;
% 
%             else
%                 angle_all(n,m)=round(angle_tmp/Delta)*Delta;
% 
%             end
           
        end
    end
    
    Varpi_tmp=exp(-angle(Varpi)*1i);
    Varpi_qtz=angle_qt(Varpi_tmp,PhaseSet);

    Ind_cache=zeros(Nt,N_rf);
    for n=1:Nt
        tmp_cache=zeros(1,N_rf);
        for m=1:N_rf
%             tmp_cache(m)=abs(Varpi(n,m))*cos(angle(Varpi(n,m))+ angle_all(n,m));
            tmp_cache(m)=abs(Varpi(n,m))*cos(angle(Varpi(n,m))+ angle(Varpi_qtz(n,m)));
        end
        [v,id]=max(tmp_cache);
        Ind_cache(n,id)=1;
    end
    BF=Varpi_qtz.*Ind_cache;
%     BF=exp(1i*angle_all).*Ind_cache;

end

function Obj=Obj_cal(Hall,Frf,Fbb_all,Wrf,Wbb_all,Mall,Xall,Yall,rhot,rhor)
[Nr,Nt,Nc]=size(Hall);
obj=0;
for nc=1:Nc
    Hc=Hall(:,:,nc);
    Xc=Xall(:,:,nc);
    Mc=Mall(:,:,nc);
    Yc=Yall(:,:,nc);

    obj=obj+real(log(det(Mc)))-1/(2*rhot)*norm(Xc-Frf*Fbb_all(:,:,nc),'fro')^2-...
        1/(2*rhor)*norm(Yc-Wrf*Wbb_all(:,:,nc),'fro')^2;
end
Obj=obj/Nc;
end

function SE=achievable_rate(Hall,Frf,Fbb_all,Wrf,Wbb_all,noise_power,Ns)
[Nr,Nt,Nc]=size(Hall);
rate=0;
for nc=1:Nc
    Hc=Hall(:,:,nc);
    Wc=Wrf*Wbb_all(:,:,nc);
    Fc=Frf*Fbb_all(:,:,nc);
    tmp=Hc*Fc;
    rate=rate+real(log2(det( eye(Ns)+ 1/noise_power*pinv(Wc)*tmp*tmp'*Wc  )));
end
SE=rate/Nc;
end