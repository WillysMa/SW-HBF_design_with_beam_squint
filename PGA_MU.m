function [W_new,iter_data]=PGA_MU(params_pga,W,Xi,Yall,Zall)
acc=params_pga.acc;
iter_max=params_pga.iter_max;
c=0.1;alpha=0.1;beta=0.1;
[N,~,Nc]=size(Yall);
[N,N_rf]=size(W);
fw=1;
fw_old=0;iter=0;
fw_all=[];err=1;
for iter=1:iter_max 
    fw_old=fw;

    temp_sum=0;
    for nc=1:Nc
    temp_sum=temp_sum+Yall(:,:,nc)*W*Zall(:,:,nc);
    end
    grad=2*Xi-2*temp_sum;
    grad_norm=grad/norm(grad,"fro");
    step_size=1;
    counter=0;
    while 1 && counter<=100
        counter=counter+1;
        DeltaX=grad_norm;
        x_new=W+step_size*DeltaX;
        obj_new=Calcu_FrfMU(x_new,Xi,Yall,Zall);
        obj_benchmark=fw+alpha*step_size*trace(grad_norm'*DeltaX);

        if obj_new>obj_benchmark
            break
        else
            step_size=beta*step_size;
        end
    end%backtracking line search

    W=W+step_size*real(grad_norm);% gradient ascend
    W=project(W);% projection
    fw=Calcu_FrfMU(W,Xi,Yall,Zall);

    fw_all=[fw_all,fw];
    
    err=abs(fw-fw_old)/abs(fw_old);
    if err<=acc
       break 
    end
%     err=abs(fw-fw_old)/abs(fw)
end
W_new=W;

iter_data=fw_all;
% figure
% plot(fw_all)
% ylabel('Objective')
% xlabel('iteration')
% grid on
% box on
% xlim([1,length(fw_all)])
% title('PGA-MU')
cccc=1;
end
%% defined functions
function W_proj=project(W)
[r,c]=size(W);
Max_v=1;
W_proj=zeros(r,c);
for i=1:r
    for j=1:c
        x=real(W(i,j));
        if x<=0
            W_proj(i,j)=0;
        elseif x>=Max_v
            W_proj(i,j)=Max_v;
        else
            W_proj(i,j)=x;
        end
        
    end
end

end

function obj_value=Calcu_FrfMU(W,Xi,Yall,Zall)
[N,~,Nc]=size(Yall);
 tmp=0;
 for nc=1:Nc
    tmp=tmp+trace(W'*Yall(:,:,nc)*W*Zall(:,:,nc));
 end
obj_value=2*real(trace(Xi*W'))-real(tmp);
end
