function [SE,SE0,TC]=HBF_FTTD(Hcc,DBF,N_rf,Nfdt,Nfdr,noise_power,AtnDis,B, ofdm_paras)

tStart = tic;
fc=ofdm_paras.fc;% carrier frequency
ofdm_paras.B=B;% bandwidth
[Nr,Nt,Nc]=size(Hcc);
stopping_criteria.iter_max=1e+3;
stopping_criteria.iter_acc=1e-4;
         
Delay_ub=(Nt-1)*AtnDis/fc;
%         Delay_ub=(Nt-2)*AtnDis/fc;

BF_opt=DBF.Precoder_all;
[~,Ns,~]=size(BF_opt);
[Slct_mtrx_tx,TTD_mtrx_tx,BF_digital_tx,Array_effciency_tx]=FTTD_Alg(BF_opt,N_rf,Nfdt,Delay_ub,ofdm_paras,stopping_criteria);
for nc=1:Nc
    precoder_opt=BF_opt(:,:,nc);
    BF_digital_tx(:,:,nc)=norm(precoder_opt,'fro')/norm(Slct_mtrx_tx*TTD_mtrx_tx(:,:,nc)*BF_digital_tx(:,:,nc),'fro')*BF_digital_tx(:,:,nc);     
end

Combiner_opt=zeros(Nr,Ns,Nc);
for nc=1:Nc
    Hc=Hcc(:,:,nc)*Slct_mtrx_tx*TTD_mtrx_tx(:,:,nc)*BF_digital_tx(:,:,nc);
    [U,S,V]=svd(Hc);
    Combiner_opt(:,:,nc)=U(:,1:Ns);
end

Delay_ub=(Nr-1)*AtnDis/fc;
[Slct_mtrx_rx,TTD_mtrx_rx,BF_digital_rx,Array_effciency_rx]=FTTD_Alg(Combiner_opt,N_rf,Nfdr,Delay_ub,ofdm_paras,stopping_criteria);

rate=zeros(1,Nc);
rate0=zeros(1,Nc);
for nc=1:Nc
    Fk=Slct_mtrx_tx*TTD_mtrx_tx(:,:,nc)*BF_digital_tx(:,:,nc);
    Wk=Slct_mtrx_rx*TTD_mtrx_rx(:,:,nc)*BF_digital_rx(:,:,nc);
    rate(nc)=real(log2(det( eye(Ns)+1/noise_power* pinv(Wk)*Hcc(:,:,nc)*Fk*Fk'*Hcc(:,:,nc)'*Wk  )));
     rate0(nc)=real(log2(det( eye(Nr)+1/noise_power*Hcc(:,:,nc)*Fk*Fk'*Hcc(:,:,nc)' )));

end
SE=mean(rate);
SE0=mean(rate0);
TC=toc(tStart);
ccc=1;
end
%% defined function
function [Slct_mtrx,TTD_mtrx,BF_digital,Array_effciency]=FTTD_Alg(BF_opt,N_rf,Q,Delay_ub,ofdm_paras,stopping_criteria)
[N,Ns,Nc]=size(BF_opt);
fc=ofdm_paras.fc;% carrier frequency
B=ofdm_paras.B;
iter_max=stopping_criteria.iter_max;
iter_acc=stopping_criteria.iter_acc;

Slct_mtrx=zeros(N,N_rf*Q);
index_vec=randi(N_rf*Q,[1,N]);
for nt=1:N
    Slct_mtrx(nt,index_vec(nt))=1;
end
Delay_all=Delay_ub/(Q-1)*(0:Q-1)';
TTD_mtrx=zeros(N_rf*Q,N_rf,Nc);
BF_digital=zeros(N_rf,Ns,Nc);
for nc=1:Nc
    f=fc+(nc-(Nc+1)/2)*B/Nc;
    for nn=1:N_rf
        TTD_mtrx((nn-1)*Q+1:nn*Q,nn,nc)=exp(-1i*2*pi*f*Delay_all);
    end
    precoder_opt=BF_opt(:,:,nc);
    tmp=precoder_opt'*Slct_mtrx*TTD_mtrx(:,:,nc);
    [Ue,Se,Ve]=svd(tmp);
    BF_digital(:,:,nc)=Ve(:,1:Ns)*Ue';

%     BF_digital(:,:,nc)=pinv(Slct_mtrx*TTD_mtrx(:,:,nc))*precoder_opt;
end

% iter_max=1e+3;
% iter_acc=1e-6;
Err_cache=[];
SE_cache=[];
for iter=1:iter_max
                               
    for nt=1:N
         Sr=zeros(1,N_rf*Q);
         Cache_vec=0;
        for nc=1:Nc
            precoder_opt=BF_opt(:,:,nc);
            Cache_vec=Cache_vec-2*real(TTD_mtrx(:,:,nc)*BF_digital(:,:,nc)*precoder_opt(nt,:)')...
                        +real(diag(TTD_mtrx(:,:,nc)*BF_digital(:,:,nc)*BF_digital(:,:,nc)'*TTD_mtrx(:,:,nc)'));
        end
        [value,index]=min(Cache_vec);
        Sr(index)=1;
        Slct_mtrx(nt,:)=Sr;     
    end
    err_sum=0;
    SE_iter=0;
    for nc=1:Nc
        precoder_opt=BF_opt(:,:,nc);

        tmp=precoder_opt'*Slct_mtrx*TTD_mtrx(:,:,nc);
        [Ue,Se,Ve]=svd(tmp);
        BF_digital(:,:,nc)=Ve(:,1:Ns)*Ue';

%         BF_digital(:,:,nc)=pinv(Slct_mtrx*TTD_mtrx(:,:,nc))*precoder_opt;

         
%         Heff=Hcc(:,:,nc)* Slct_mtrx*TTD_mtrx(:,:,nc)*BF_digital(:,:,nc);
%         SE_iter=SE_iter+real(log2(det( eye(Ns)+ 1/noise_power*pinv(combiner_opt)*Heff*Heff'*combiner_opt  )));

        err_sum=err_sum+norm(precoder_opt-Slct_mtrx*TTD_mtrx(:,:,nc)*BF_digital(:,:,nc),'fro')^2;
    end
%     SE_cache=[SE_cache SE_iter/Nc];
    Err_cache=[Err_cache err_sum];
    if iter>=2
        diff=abs(Err_cache(end)-Err_cache(end-1));
        if (diff/abs(Err_cache(end))<=iter_acc) || (iter>=iter_max)
            break
        end

    end

end%->iter=1:iter_max


AntnNum_perTTD=sum(Slct_mtrx);
NumAntn_eff=length(find(AntnNum_perTTD>0));
Array_effciency=NumAntn_eff/N;

% figure
% plot(Err_cache)
% xlabel('iter')
% ylabel('Error')
% grid on


ccc=1;
end