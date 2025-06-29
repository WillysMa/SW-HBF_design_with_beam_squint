function [SE_SW_tx,F_rf,iter_data]=Tabu_search(params_ts,BF_int_ini,Ns,Surrate,coef,obj_fun)

local_iter_max=params_ts.local_iter_max;%TS stop criteria
iter_max=params_ts.iter_max;
Len_ts=params_ts.Len_tl;% length of tabu list
%% only short memory, one-elemnet difference neighbor

SE_ini=obj_fun(Surrate,BF_int_ini,coef);

sBest=BF_int_ini;
vBest=SE_ini;
Best_candidate_s=sBest;

tabu.s{1}=sBest;
tabu.v{1}=vBest;

iter=0;record=0;SE_best_all=[];
while iter<=iter_max && record<=local_iter_max
    iter=iter+1;
    vBest_old=vBest;
Neighbors=neighbor_generate(Best_candidate_s,Ns,Surrate,coef,obj_fun);
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

F_rf=sBest;
SE_SW_tx=obj_fun(Surrate,F_rf,coef);
iter_data=SE_best_all;
% figure
% plot(SE_best_all)
% ylabel('Objective value')
% xlabel('iteration')
% grid on
% box on
% xlim([1 length(SE_best_all)])
% title('TS Algorithm')
ccc=1;
end



