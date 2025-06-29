function Wq=angle_qt(W,phase_Set)

Phase_threshold=(phase_Set(1:end-1)+phase_Set(2:end))/2;
[N,N_rf]=size(W);
Wq=zeros(N,N_rf);
    for row_id=1:N
        for col_id=1:N_rf


            angle_opt=angle(W(row_id,col_id));

            for t=1:length(Phase_threshold)
                if angle_opt<=Phase_threshold(t)
                    angle_q=phase_Set(t);
                    break
                else
                    angle_q=phase_Set(end);
                end
            end

%             [~,id]=find(phase_Set>=angle_opt);
%             if length(id)>=1
%                 if id(1)>=2
%                     intval_low=phase_Set(id(1)-1);
%                     intval_up=phase_Set(id(1));
%                     med_val=(intval_low+intval_up)/2;
%                     if  angle_opt>=med_val
%                         angle_q=intval_up;
%                     else
%                         angle_q=intval_low;
%                     end
%                 else
%                     angle_q=phase_Set(1);
%                 end
% 
%             else
%                 angle_q=phase_Set(end);% angle is larger than all values in angle_dictionary
%             end

            Wq(row_id,col_id)=abs(W(row_id,col_id))*exp(1j*angle_q);%1/sqrt(N)*
        end
    end
end