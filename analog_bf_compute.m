function analog_bf=analog_bf_compute(coef,N,N_rf,F_eff)
% coef=1/noise_power^2/Nr;
W_rf=ones(N,N_rf);
for j=1:N_rf
    sub_W_rf=[W_rf(:,1:j-1),W_rf(:,j+1:end)];
    C=eye(N_rf-1)+coef*sub_W_rf'*F_eff*sub_W_rf;
    G=coef*F_eff-coef^2*F_eff*sub_W_rf*inv(C)*sub_W_rf'*F_eff;
    for i=1:N
        eta_ij=G(i,:)*W_rf(:,j);
        if abs(eta_ij)<=1e-4
            W_rf(i,j)=1;
        else
            W_rf(i,j)=eta_ij/abs(eta_ij);
        end      
    end
end
analog_bf=W_rf;
end