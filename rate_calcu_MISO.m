function Rate=rate_calcu_MISO(Hall,Frf,Fbb_all,noise_power,Wgtu)
[Nt,Nu,Nc]=size(Hall);
temp=0;
for nc=1:Nc
    for ii=1:Nu
        Interference_power=norm(Hall(:,ii,nc)'*Frf*Fbb_all(:,:,nc),2)^2+noise_power;
        target=Hall(:,ii,nc)'*Frf*Fbb_all(:,ii,nc);
        target_power=abs(target)^2;
        
        SINR=target_power/(Interference_power-target_power);

        temp=temp+Wgtu(ii)*log2(1+SINR);        
    end
end
Rate=temp/Nc;
end