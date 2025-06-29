function SE=compute_rate(HBF_info,BeaMtx,coef)
[~,Nc]=size(HBF_info.H_eff);
[~,N_rf]=size(BeaMtx);

% tic
% rate2=zeros(1,Nc);
% for k=1:Nc
%     rate2(k)=real(log2(det(  eye(N_rf)+ coef* pinv(BeaMtx)*HBF_info.H_eff{k}*BeaMtx  )));
% end
% SE=mean(rate2);
% toc

% tic
rate=zeros(1,Nc);
[U,~]=qr(BeaMtx,0);
for k=1:Nc
    Surrgate=U'*HBF_info.H_eff{k}*U;
    Eigen_all=eig(Surrgate);
    rate(k)=real(sum(log2(  1+ coef* Eigen_all  )));
end

SE=mean(rate);

% if abs(SE-SE_cmp)>=1e-4
%     disp('Warning')
% end
% toc
ccc=1;
end


% function eig_value=block_power_eig(A,Ns)
% [N,~]=size(A);
% V=rand(N,Ns);
% acc=1e-2;
% iter_max=1e+8;
% for iter=1:iter_max
%     B=A*V;[Q,R]=qr(B,0);
%     V=Q(:,1:Ns);Lambda=R(1:Ns,:);
%     err=norm(A*V-V*Lambda,'fro')^2/norm(A*V,'fro')^2;
%     if err<=acc
%         break
%     end
% end
% eig_value=diag(Lambda);
% end