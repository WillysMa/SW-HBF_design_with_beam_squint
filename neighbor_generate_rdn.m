function Neighbor_set=neighbor_generate_rdn(W_rf_ini,Vec_pgr,Index_bf,Npart1,N,Ns,N_rf,X,coef,obj_fun,NumNeighbor)
tol=1e-6;
V_rf_ini=vec(W_rf_ini);
BF_int_ini=V_rf_ini(Index_bf);
W_vec=BF_int_ini;
neighbor_size=length(W_vec);
num_neighbor=0;
Neighbor_set.W=[]; 
Neighbor_set.SE=[];
counter=0;  
List_visited=[];
    while counter<NumNeighbor % randomly generate NumNeighbor candidate as the neighbor set.       
        while 1
            s=randi(neighbor_size);
            if ~ismember(s,List_visited)
                break
            end
        end
        Vec_sp=Vec_pgr;
        Wx=W_vec;
        Wx(s)=~W_vec(s);
        Nsum=sum(Wx)+Npart1;
        if Nsum>=max(N*N_rf/5,N_rf) && Nsum<=(N-1)*N_rf+1% narrow down search space
            Vec_sp(Index_bf)=Wx;
            Wx_rshp=reshape(Vec_sp,[N,N_rf]);
            NumActiveAntennas_perRFchain=sum(Wx_rshp);
            if min(sum(NumActiveAntennas_perRFchain))>=N/5 && rank(Wx_rshp,tol)>=Ns
                
                num_neighbor=num_neighbor+1;
                Neighbor_set.W{num_neighbor}=Wx_rshp;
                Neighbor_set.SE{num_neighbor}=obj_fun(X,Wx_rshp,coef);

                counter=counter+1; 
               
            end
        else
            continue
        end

        List_visited=[List_visited s];
        ccc=1;
    end

ccc=1;
end