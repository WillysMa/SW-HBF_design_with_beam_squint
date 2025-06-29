function [vBest,sBest,iter_data]=Tabu_RS(params_TbRs,W_pg,Ns,Srrgt,coef,obj_fun)
[N,N_rf]=size(W_pg);
sp_max=params_TbRs.iter_max;
Ls=params_TbRs.Len_tl;
local_iter_max=params_TbRs.local_iter_max;
iter_max=5;

vBest=0;
sBest=0;
tabu.s{1}=sBest;
tabu.v{1}=vBest;
id_sp=1;
num_tb=1;
record=0;
SE_all=[];
    for sp=1:sp_max
        vBest_old=vBest;
        W_rf = mapping(W_pg,Ns);
        Nall=sum(sum(W_rf));
        if Nall>=N_rf && Nall<=(N-1)*N_rf+1
            if sp==1
                SE_sp=obj_fun(Srrgt,W_rf,coef);
                Candidate.Wrf{sp}=W_rf;
                Candidate.SE(sp)=SE_sp;
                
                sBest=W_rf;
                vBest=SE_sp;
                
                    tabu.s{num_tb}=W_rf;
                    tabu.v{num_tb}=SE_sp; 
                    num_tb=num_tb+1;
            else
                flag=1;iter_id=0;
                while flag
                   
                    W_rf = mapping(W_pg,Ns);
                    [~,num_s]=size(tabu.s);
                    clen=min(num_s,Ls);
                    diff=zeros(1,clen);
                    for j=1:clen
                        diff(j)=norm( W_rf - tabu.s{j},'fro');
                    end
                    diff_ind=diff>0;
                    
                    if isempty(find(diff_ind==0))% we have found a new feasible point
                        
                        if num_tb>Ls
                            pos=mod(num_tb-1,Ls);
                            tabu.s{pos+1}=Best_candidate_s;
                            tabu.v{pos+1}=Best_candidate_v;
                        else
                            tabu.s{num_tb}=W_rf;
                            tabu.v{num_tb}=SE_sp; 
                        end
                        num_tb=num_tb+1;
                
                        break                       
                    else
                        iter_id=iter_id+1;
%                         iter_id=0; 
                    end
                    if iter_id>=iter_max
                        flag=0;% it is a feasible point
                    end
                end
                Candidate.Wrf{sp}=W_rf;
                SE_sp=obj_fun(Srrgt,W_rf,coef);
                Candidate.SE(sp)=SE_sp;
                
                if SE_sp>=vBest
                    sBest=W_rf;
                    vBest=SE_sp; 
                    id_sp=id_sp+1;
                 end
    

            end
            
        else
            continue
        end
        
        if vBest-vBest_old==0
            record=record+1;
        else
            record=0;
        end
        
        if record>=local_iter_max
            break
        end
    SE_all=[SE_all vBest];
    end
iter_data=SE_all;
% figure
% plot(SE_all)
% ylabel('Objective')
% xlabel('Iteration')
% grid on
% box on
% title('APGA-TbRs Algorithm')

ccc=1;
end