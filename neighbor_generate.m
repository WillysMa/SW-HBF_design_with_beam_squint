function Neighbor_set=neighbor_generate(W_int_ini,Ns,X,coef,obj_fun)
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
                Neighbor_set.SE{num_neighbor}=obj_fun(X,Wx,coef);
            end
        else
            continue
        end
    end

end