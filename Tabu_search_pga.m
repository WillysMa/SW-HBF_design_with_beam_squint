function [SE_SW_tx,sBest,SE_best_all]=Tabu_search_pga(params_ts,BF_int_ini,Ns,Vec_pgr,Index_bf,Npart1,Surrgate,coef,obj_fun,NumNeighbor)

local_iter_max=params_ts.local_iter_max;%TS stop criteria
iter_max=params_ts.iter_max;
Len_ts=params_ts.Len_tl;% length of tabu list
%% only short memory, one-elemnet difference neighbor

[N,N_rf]=size(BF_int_ini);
SE_ini=obj_fun(Surrgate,BF_int_ini,coef);

sBest=BF_int_ini;
vBest=SE_ini;
Best_candidate_s=sBest;

tabu.s{1}=sBest;
tabu.v{1}=vBest;

iter=0;record=0;SE_best_all=[];
while iter<=iter_max && record<=local_iter_max
    iter=iter+1;
    vBest_old=vBest;

    if NumNeighbor>=1e+3
        Neighbors=neighbor_generate_pga(Best_candidate_s,Vec_pgr,Index_bf,Npart1,N,Ns,N_rf,Surrgate,coef,obj_fun);

    else
        Neighbors=neighbor_generate_rdn(Best_candidate_s,Vec_pgr,Index_bf,Npart1,N,Ns,N_rf,Surrgate,coef,obj_fun,NumNeighbor);

    end
[~,Nb]=size(Neighbors.W);
[~,num_s]=size(tabu.s);

if Nb>=1
    Best_candidate_s=Neighbors.W{1};
    Best_candidate_v=Neighbors.SE{1};
    for i=1:Nb
        Candidate_s=Neighbors.W{i};
        Candidate_v=Neighbors.SE{i};
            clen=min(num_s,Len_ts);
            diff=zeros(1,clen);
            for j=1:clen
                diff(j)=norm( Candidate_s - tabu.s{j},'fro');
            end
            diff_med=~(diff>0); 
            if sum(diff_med)==0 % the generated point didn't appeared in Tabu list.
                if Candidate_v>=Best_candidate_v
                    Best_candidate_s=Candidate_s;
                    Best_candidate_v=Candidate_v;
                end
            end
    end

    if Best_candidate_v>=vBest
        sBest=Best_candidate_s;
        vBest=Best_candidate_v;  
    end

    if iter>Len_ts
        pos=mod(iter-1,Len_ts);
        tabu.s{pos+1}=Best_candidate_s;
        tabu.v{pos+1}=Best_candidate_v;
    else
        tabu.s{iter}=Best_candidate_s;
        tabu.v{iter}=Best_candidate_v; 
    end
end
    if vBest-vBest_old==0
        record=record+1;
    else
        record=0;
    end
    SE_best_all=[SE_best_all,vBest];
end

SE_SW_tx=obj_fun(Surrgate,sBest,coef);
% figure
% plot(SE_best_all)
% ylabel('SE-ts')
% xlabel('iterations')
% grid on
ccc=1;
end



