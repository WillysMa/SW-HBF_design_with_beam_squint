function BF=Exhaustive_search(N,N_rf,Ns,Srrgt,coef,obj_fun)
vBest=0;
sBest=0;
Nd=N*N_rf;
space_size=2^Nd-1;
Nlb=N_rf;
Nub=N+(N_rf-1)*(N-1);
for i=1:space_size
    symi=dec2bin(i,Nd);
    vi=symi-'0';
    Ni=sum(vi);
    if Ni>=Nlb && Ni<=Nub
       Wi=reshape(vi,[N,N_rf]);
       if rank(Wi,1e-4)>=Ns
           SEi=obj_fun(Srrgt,Wi,coef);
           if SEi>=vBest
               vBest=SEi;
               sBest=Wi;
           end
       end    
    end        
    ccc=1;
end
BF=sBest;

end