function BF=PGA_Ehaustive_Search(N,N_rf,Ns,Index_bf,Vec_pgr,Npart1,Srrgt,coef,obj_fun)

vBest=0;
sBest=0;
Nbit=length(Index_bf);
space_size=2^Nbit-1;

    for s=0:space_size
        s_str=dec2bin(s,Nbit);
        s_value=s_str-'0';
        Vec_sp=Vec_pgr;
        
        Nsum=sum(s_value)+Npart1;
        if Nsum>=N_rf && Nsum<=(N-1)*N_rf+1% narrow down search space
            Vec_sp(Index_bf)=s_value;
            Wx_rshp=reshape(Vec_sp,[N,N_rf]);
            if rank(Wx_rshp,1e-4)>=Ns
                SEs=obj_fun(Srrgt,Wx_rshp,coef);
                if SEs>=vBest
                    sBest=Wx_rshp;
                    vBest=SEs;
                end
                
            end
        else
            continue
        end
        ccc=1;
    end
    BF=sBest;
end