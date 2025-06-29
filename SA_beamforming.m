function [SE,BF]=SA_beamforming(Surrgate,F_ini_int,coef,Ns,params)
T=params.T;
tol=params.tol;
kappa=params.kappa;
Obj_list=[];

Obj_old=compute_rate(Surrgate,F_ini_int,coef);
Fsmp=F_ini_int;
while T>=tol

[Neighbor_set,num_neighbor]=neighbor_generate(Fsmp,Ns);
smp_id=randi(num_neighbor);
F_new=Neighbor_set.W{smp_id};
Obj_new=compute_rate(Surrgate,F_new,coef);
if Obj_new>Obj_old
    Fsmp=F_new;
else
    prob=exp(-(Obj_old-Obj_new)/T);
    if prob>=rand
        Fsmp=F_new;
    end

end
Obj_old=Obj_new;
T=kappa*T;

Obj_list=[Obj_list Obj_new];
ccc=1;
end
SE=Obj_new;
BF=F_new;
% figure
% plot(Obj_list)
% xlabel('Iteration')
% ylabel('Objective')
% grid on 
% box on
ccc=1;
end

%% defined functions
function [Neighbor_set,num_neighbor]=neighbor_generate(W_int_ini,Ns)
[N,N_rf]=size(W_int_ini);
tol=1e-4;
W_vec=vec(W_int_ini);
neighbor_size=N*N_rf;
num_neighbor=0;
Neighbor_set.W=[]; 
Neighbor_set.SE=[];
    for s=1:neighbor_size
        Wx=W_vec;
        Wx(s)=~W_vec(s);
        Nsum=sum(Wx);
        if Nsum>=N_rf && Nsum<=(N-1)*N_rf+1% narrow down search space
            Wx=reshape(Wx,[N,N_rf]);
            if rank(Wx,tol)>=Ns
                num_neighbor=num_neighbor+1;
                Neighbor_set.W{num_neighbor}=Wx;
            end
        else
            continue
        end
    end

end