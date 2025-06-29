function [SE_new,W_new,iter_data]=PGA_new(params_pga,W,Srrt,coef,obj_fun)
tol=params_pga.acc;
iter_max=params_pga.iter_max;
alpha=0.1;beta=0.3;
[~,Nc]=size(Srrt.H_eff);
[N,N_rf]=size(W);
fw=obj_fun(Srrt,W,coef);
fw_all=[];
for iter=1:iter_max

    fw_old=fw;

    accu=zeros(N,N_rf);
    for k=1:Nc
    Ak=eye(N) + coef*Srrt.H_eff{k};
    accu=accu+2*Ak*W*inv(W'*Ak*W);
    end
    grad=accu/Nc-2*W*inv(W'*W);
    grad_unit=grad/norm(grad,'fro');

    step_size=1;counter=0;
    while 1 && counter<=30
        counter=counter+1;
        DeltaX=grad_unit;
        x_new=W+step_size*DeltaX;
        obj_new=compute_rate(Srrt,x_new,coef);
        obj_benchmark=fw+alpha*step_size*norm(grad_unit,"fro")^2;

        if obj_new>obj_benchmark
            break
        else
            step_size=beta*step_size;
        end
    end%backtracking line search


    W=W+step_size*real(grad_unit);% gradient ascend
    W=project(W);% projection
    fw=compute_rate(Srrt,W,coef);

    fw_all=[fw_all,fw];
%     err=fw-fw_old;
    err=abs(fw-fw_old)/abs(fw);
    if err<=tol
        break
    end
end
W_new=W;
SE_new=compute_rate(Srrt,W,coef);
iter_data=fw_all;
% figure
% plot(fw_all)
% ylabel('SE (bps/Hz)')
% xlabel('iteration')
% grid on
% box on
% xlim([1,length(fw_all)])
% title('PGA Algorithm')
cccc=1;
end
%% defined functions
function W_proj=project(W)
[r,c]=size(W);
W_proj=zeros(r,c);
for i=1:r
    for j=1:c
        x=W(i,j);
        if x<=0
            W_proj(i,j)=0;
        elseif x>=1
            W_proj(i,j)=1;
        else
            W_proj(i,j)=x;
        end
        
    end
end

end


